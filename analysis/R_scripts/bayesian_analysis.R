# Author: Kevin See
# Purpose: Calculate adult overshoots at Priest Rapids dam, using Bayesian model
# Created: 5/8/20
# Last Modified: 9/3/20
# Notes: 

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(readxl)
library(magrittr)
library(msm)
library(ggrepel)
library(rjags)
library(postpack)


theme_set(theme_bw())

#-----------------------------------------------------------------
# read in data
obs_ovrsht_org = read_excel('data/overshoot estimates.xlsx',
                            range = 'A3:C10',
                            col_names = c('Year',
                                          'juv_tags',
                                          'pom_overshoots')) %>%
  mutate(origin = 'Wild',
         location = 'downstream_PRA') %>%
  bind_rows(read_excel('data/overshoot estimates.xlsx',
                       range = 'A19:B26',
                       col_names = c('Year',
                                     'juv_tags')) %>%
              mutate(origin = 'Wild',
                     location = 'at_PRA')) %>%
  bind_rows(read_excel('data/overshoot estimates.xlsx',
                       range = 'A31:C38',
                       col_names = c('Year',
                                     'juv_tags',
                                     'pom_overshoots')) %>%
              mutate(origin = 'Hatchery',
                     location = 'downstream_PRA')) %>%
  bind_rows(read_excel('data/overshoot estimates.xlsx',
                       range = 'A47:B54',
                       col_names = c('Year',
                                     'juv_tags')) %>%
              mutate(origin = 'Hatchery',
                     location = 'at_PRA')) %>%
  select(Year,
         Origin = origin,
         location,
         everything()) %>%
  arrange(location, Origin, Year)

# data on individual overshoots
obs_ovrsht = read_excel("data/known overshoot fallback locations.xlsx") %>%
  janitor::clean_names() %>%
  rename(notes = x5,
         site = downstream_obs_site) %>%
  tidyr::fill(total) %>%
  mutate(year = year + 1)

# estimated overshoots
est_ovrsht = obs_ovrsht %>%
  group_by(Year = year, site) %>%
  summarise(obs_ovrst_tags = n_distinct(pit_tag)) %>%
  ungroup() %>%
  # get estimates of downstream detection prob
  left_join(as.list(2011:2018) %>%
              rlang::set_names() %>%
              map_df(.id = 'Year',
                     .f = function(x) {
                       read_excel(paste0('data/DABOM_results/PRA_Steelhead_', x[1], '_20190610.xlsx'),
                                  "Detection") %>%
                         janitor::clean_names()
                     }) %>%
              mutate_at(vars(Year),
                        list(as.numeric)) %>%
              filter(!grepl("A0$", node)) %>%
              mutate(site = str_remove(node, "B0$")) %>%
              select(Year, site, 
                     det_est = estimate, 
                     det_se = se)) %>%
  # drop sites with no detection probability: they are all not the head of branches
  filter(!is.na(det_est)) %>%
  mutate(beta_alpha = ((1 - det_est) / det_se^2 - 1 / det_est) * det_est^2,
         beta_alpha = if_else(beta_alpha < 0, 0.01, beta_alpha),
         beta_beta = beta_alpha * (1 / det_est - 1)) %>%
  mutate(beta_alpha = if_else(det_est == 1, 1, beta_alpha),
         beta_beta = if_else(det_est == 1, 0, beta_beta)) %>%
  mutate(est_ovrst_tags = obs_ovrst_tags / det_est,
         # est_ovrst_tags = floor(est_ovrst_tags),
         est_ovrst_tags = round(est_ovrst_tags),
         rho = 1 / (beta_alpha + beta_beta + 1),
         var_ovrst_tags = est_ovrst_tags * det_est * (1 - det_est) * (1 + (est_ovrst_tags - 1)*rho),
         se_ovrst_tags = sqrt(var_ovrst_tags)) %>%
  mutate(Origin = 'Wild') %>%
  select(Year, Origin, everything())

# average detection at sites
as.list(2011:2018) %>%
  rlang::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           read_excel(paste0('data/DABOM_results/PRA_Steelhead_', x[1], '_20190610.xlsx'),
                      "Detection") %>%
             janitor::clean_names()
         }) %>%
  mutate_at(vars(Year),
            list(as.numeric)) %>%
  filter(grepl('PRO', node) |
           grepl('ICH', node) |
           grepl('PRV', node) |
           grepl('TMF', node) |
           grepl('^JD1', node)) %>%
  filter(!grepl('A0$', node)) %>%
  mutate(cv = se / estimate) %>%
  group_by(node) %>%
  summarise(mean_tags = mean(n_tags),
            mean = mean(estimate),
            sd_est = sd(estimate),
            mean_se = mean(se),
            mean_cv = mean(cv, na.rm = T))

# get estimates of downstream escapement
dwnstrm_est = as.list(2011:2018) %>%
  rlang::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           read_excel(paste0('data/DABOM_results/PRA_Steelhead_', x[1], '_20190610.xlsx'),
                      1) %>%
             filter(Population == 'BelowPriest')
         }) %>%
  mutate_at(vars(Year),
            list(as.numeric)) %>%
  mutate(Origin = recode(Origin,
                         'Natural' = 'Wild'))

# join escapement estimates and estimated overshoot tags, using all branches and all overshoot tags
# Andrew's preferred method
mod_data = est_ovrsht %>%
  group_by(Year, Origin) %>%
  summarise(tags_obs = sum(obs_ovrst_tags),
            det_prob = mean(det_est),
            det_prob_wgt = weighted.mean(det_est, w = est_ovrst_tags),
            tags_est = sum(est_ovrst_tags),
            tags_se = sqrt(sum(se_ovrst_tags^2))) %>%
  ungroup() %>%
  full_join(dwnstrm_est %>%
              filter(Origin == "Wild") %>%
              select(-Population) %>%
              rename(escp_est = Estimate,
                     escp_se = SE,
                     escp_lwr = lowerCI,
                     escp_upr = upperCI))

#-------------------------------------------------
# turn it into data for JAGS
jags_data = mod_data %>%
  select(Year, Origin, tags_est:escp_se) %>%
  mutate(tags_prec = 1 / tags_se^2,
         escp_prec = 1 / escp_se^2,
         escp_est_log = log(escp_est)) %>%
  select(-ends_with("se"),
         -Origin) %>%
  select(-escp_est,
         -escp_prec,
         -Year) %>%
  as.list()

# add data for predictions
# known overshoot tags observed at PRA, their origin, and POM estimates of downstream escapement
jags_data = c(jags_data,
              obs_ovrsht_org %>%
                filter(location == "at_PRA") %>%
                select(Year, Origin, ovrst_tags = juv_tags) %>%
                left_join(dwnstrm_est %>%
                            select(Year, Origin, dwnstrm_escp = Estimate, dwnstrm_se = SE)) %>%
                # filter(Origin == 'Wild') %>%
                # mutate(ovrst_org = as.numeric(as.factor(Origin))) %>%
                select(-Year,
                       -Origin) %>%
                as.list())

# # fix some of the SE to close to 0
# jags_data$tags_prec = rep(1000, length(jags_data$tags_prec))
# jags_data$dwnstrm_se = rep(0.001, length(jags_data$dwnstrm_se))

#----------------------------------------
# function to write JAGS model
jags_simple_model = function() {
  "  # PRIORS
  for(i in 1:2) {
    beta[i] ~ dt(0, 0.01, 1)
  }
  sigma ~ dt(0, 0.01, 1)T(0,)
  tau <- pow(sigma, -2)
  
  # MODEL
  for(i in 1:length(tags_est)) {
    mu[i] <- beta[1] + beta[2] * log(tags_est[i])
    escp_est_log[i] ~ dnorm(mu[i], tau)
  }
  
  for(i in 1:length(ovrst_tags)) {
    pred_mu_log[i] <- beta[1] + beta[2] * log(ovrst_tags[i])
    pred_ovrshts_log[i] ~ dnorm(pred_mu_log[i], tau)T(log(dwnstrm_escp[i]),)
    pred_ovrshts[i] <- round(exp(pred_ovrshts_log[i]))
    
    phi[i] <- dwnstrm_escp[i] / pred_ovrshts[i]
  }"
}


#----------------------------------------
jags_model = function() {
  "  # PRIORS
  for(i in 1:2) {
    beta[i] ~ dt(0, 0.01, 1)
  }
  sigma ~ dt(0, 0.01, 1)T(0,)
  tau <- pow(sigma, -2)
  
  # MODEL
  for(i in 1:length(tags_est)) {
    n_tags_org[i] ~ dnorm(tags_est[i], tags_prec[i])
    n_tags[i] <- round(n_tags_org[i])
    # n_escp_log[i] ~ dlnorm(escp_est[i], escp_prec[i])
    
    mu[i] <- beta[1] + beta[2] * log(n_tags[i])
    # mu[i] <- beta[1] + beta[2] * n_tags[i]
    
    # assuming downstream escapement estimates are known
    escp_est_log[i] ~ dnorm(mu[i], tau)
  }
  
  for(i in 1:length(ovrst_tags)) {
    # deal with uncertainty in downstream escapement estimates
    est_dwnstrm_org[i] ~ dnorm(dwnstrm_escp[i], 1 / (dwnstrm_se[i]^2))
    est_dwnstrm[i] <- round(est_dwnstrm_org[i])
    
    # predict the number of overshoot fish at Priest
    pred_mu_log[i] <- beta[1] + beta[2] * log(ovrst_tags[i])
    # pred_mu_log[i] <- beta[1] + beta[2] * ovrst_tags[i]
    pred_ovrshts_log[i] ~ dnorm(pred_mu_log[i], tau)T(log(est_dwnstrm[i]),log(1e4))
    pred_ovrshts[i] <- round(exp(pred_ovrshts_log[i]))

    # estimate survival of overshoots
    phi[i] <- est_dwnstrm[i] / pred_ovrshts[i]
  }"
}

#----------------------------------------
jags_hier_model = function() {
"  # PRIORS
  for(i in 1:2) {
    beta[i] ~ dt(0, 0.01, 1)
  }
  sigma ~ dt(0, 0.01, 1)T(0,)
  tau <- pow(sigma, -2)
  
  # hyperparameters for survival
  for(i in 1:2) {
    u[i] ~ dbeta(1,1)
    v[i] ~ dgamma(1, 0.05)
  }
  
  for(i in 1:length(ovrst_tags)) {
    # phi[i] ~ dbeta(1,1)
    phi[i] ~ dbeta(u[ovrst_org[i]] * v[ovrst_org[i]], (1 - u[ovrst_org[i]]) * v[ovrst_org[i]])T(0.001, 0.999)
  }
  
  # MODEL
  for(i in 1:length(tags_est)) {
    n_tags[i] ~ dnorm(tags_est[i], tags_prec[i])
    # n_escp_log[i] ~ dlnorm(escp_est[i], escp_prec[i])
    
    mu[i] <- beta[1] + beta[2] * log(n_tags[i])
    escp_est_log[i] ~ dnorm(mu[i], tau)
  }
  
  for(i in 1:length(ovrst_tags)) {
    pred_mu_log[i] <- beta[1] + beta[2] * log(ovrst_tags[i])
    pred_ovrshts_log[i] ~ dnorm(pred_mu_log[i], tau)T(log(dwnstrm_escp[i]),)
    pred_ovrshts[i] <- round(exp(pred_ovrshts_log[i]))

    dwnstrm_escp[i] ~ dbin(phi[i], pred_ovrshts[i])
    
  }"
}



# write model to a text file
jags_file = "model.txt"
write_model(jags_model, jags_file)

# which paramters to track?
jags_params = c("beta",
                "sigma",
                "mu",
                "n_tags",
                "phi",
                "est_dwnstrm",
                "pred_ovrshts")

# using rjags package
set.seed(3)
jags = jags.model(jags_file,
                  data = jags_data,
                  inits = list(pred_ovrshts_log = log(jags_data$dwnstrm_escp + 10)),
                  n.chains = 4,
                  n.adapt = 1000)

# burnin
update(jags, n.iter = 50000)
# posterior sampling
post = coda.samples(jags,
                    jags_params,
                    n.iter = 50000,
                    thin = 50)

# convert posteriors into long tibble
post_df = as.matrix(post,
                    chains = T,
                    iters = T) %>%
  as_tibble() %>%
  pivot_longer(cols = c(-CHAIN, -ITER),
               names_to = 'param',
               values_to = 'value') %>%
  mutate(param_fam = str_split(param, "\\[", simplify = T)[,1],
         param_num = str_extract(param, "[:digit:]+"),
         param_num = as.numeric(param_num))

param_summ = post_summ(post,
                       jags_params,
                       Rhat = T,
                       ess = T) %>%
  t() %>%
  as_tibble(rownames = "param") %>%
  mutate(cv = sd / mean)

param_summ

param_summ %>%
  mutate(param_fam = str_remove(param, "\\[[:digit:]+\\]")) %>%
  group_by(param_fam) %>%
  summarise(n_params = n_distinct(param),
            big_Rhat = n_distinct(param[Rhat > 1.05]),
            perc_prob = big_Rhat / n_params)

param_summ %>%
  filter(grepl('beta', param) |
           grepl('sigma', param))

# try to calculate posterior of R2
r2 = post_df %>%
  filter(param_fam == "mu") %>%
  left_join(mod_data %>%
              mutate(obs = log(escp_est)) %>%
              select(Year, Origin,
                     obs) %>%
              mutate(param_num = 1:n())) %>%
  group_by(CHAIN, ITER) %>%
  nest() %>%
  ungroup() %>%
  # slice(1:1000) %>%
  mutate(r2 = map_dbl(data,
                      .f = function(x) {
                        x %>%
                          select(value, obs) %>%
                          corrr::correlate(quiet = T) %>%
                          corrr::stretch(na.rm = T,
                                         remove.dups = T) %>%
                          mutate(r2 = r^2) %>%
                          pull(r2)
                      }))
as.mcmc(r2$r2) %>%
  summary()


r2 %>%
  summarise_at(vars(r2),
               list(mean = mean,
                    sd = sd,
                    `50%` = median))

quantile(r2$r2, probs = c(0.025, 0.975))


cor_df = param_summ %>%
  filter(grepl('mu', param)) %>%
  mutate(obs = jags_data$escp_est_log) %>%
  select(mean, obs) %>%
  corrr::correlate()

cor_df %>%
  corrr::stretch()

cor_mat[1,2]^2

# param_summ %>%
#   filter(grepl("^u", param) |
#            grepl('^v', param))

param_summ %>%
  filter(grepl('phi', param)) %>%
  bind_cols(obs_ovrsht_org %>%
              filter(location == "at_PRA") %>%
              select(Year, Origin)) %>%
  group_by(Origin) %>%
  summarise_at(vars(mean, sd, cv),
               list(mean))
param_summ %>%
  filter(grepl('n_tags', param))

dodge_width = 0.5
param_summ %>%
  filter(grepl('phi', param)) %>%
  bind_cols(obs_ovrsht_org %>%
              filter(location == "at_PRA") %>%
              # filter(Origin == "Wild") %>%
              select(Year, Origin, ovrst_tags = juv_tags) %>%
              left_join(dwnstrm_est %>%
                          select(Year, Origin, 
                                 dwnstrm_escp = Estimate, 
                                 dwnstrm_se = SE))) %>%
  mutate_at(vars(Year, Origin),
            list(as.factor)) %>%
  filter(Origin == "Wild") %>%
  ggplot(aes(x = Year,
             y = mean,
             color = Origin)) + 
  geom_errorbar(aes(ymin = `2.5%`,
                    ymax = `97.5%`),
                position = position_dodge(width = dodge_width),
                width = 0.1) +
  geom_point(size = 4,
             # aes(size = 1/Rhat),
             position = position_dodge(width = dodge_width)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = 'Year',
       y = 'Survival to Natal Trib',
       title = 'Overshoot Survival')


dodge_width = 0.5
param_summ %>%
  filter(grepl('pred_ovrshts', param)) %>%
  bind_cols(obs_ovrsht_org %>%
              filter(location == "at_PRA") %>%
              # filter(Origin == "Wild") %>%
              select(Year, Origin, ovrst_tags = juv_tags) %>%
              left_join(dwnstrm_est %>%
                          select(Year, Origin, 
                                 dwnstrm_escp = Estimate, 
                                 dwnstrm_se = SE))) %>%
  mutate_at(vars(Year, Origin),
            list(as.factor)) %>%
  ggplot(aes(x = Year,
             y = mean,
             color = Origin)) + 
  geom_errorbar(aes(ymin = `2.5%`,
                    ymax = `97.5%`),
                position = position_dodge(width = dodge_width),
                width = 0.1) +
  geom_point(#aes(size = 1/Rhat),
             size = 5,
             position = position_dodge(width = dodge_width)) +
  geom_point(aes(y = dwnstrm_escp),
             shape = 1,
             size = 3,
             position = position_dodge(width = dodge_width)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = 'Year',
       y = 'Overshoots at PRA',
       title = 'Overshoot Abundance')


param_summ %>%
  filter(grepl('pred_ovrshts', param)) %>%
  bind_cols(obs_ovrsht_org %>%
              filter(location == "at_PRA") %>%
              select(Year, Origin)) %>%
  bind_rows(param_summ %>%
              filter(grepl('est_dwnstrm', param)) %>%
              bind_cols(obs_ovrsht_org %>%
                          filter(location == "at_PRA") %>%
                          select(Year, Origin))) %>%
  mutate(param = str_split(param, "\\[", simplify = T)[,1]) %>%
  pivot_longer(cols = mean:cv,
               names_to = "stat",
               values_to = "value") %>%
  pivot_wider(names_from = "param",
              values_from = "value") %>%
  filter(stat %in% c("mean",
                     "50%",
                     "2.5%",
                     "97.5")) %>%
  mutate(phi = est_dwnstrm / pred_ovrshts)



post_df %>%
  filter(param_fam == 'phi') %>%
  # filter(param_fam == 'est_dwnstrm') %>%
  # filter(param_fam == 'pred_ovrshts') %>%
  left_join(dwnstrm_est %>%
              select(Year, Origin) %>% 
              arrange(Origin, Year) %>%
              mutate_all(as.factor) %>%
              mutate(param_num = 1:n())) %>%
  ggplot(aes(x = Year,
             y = value,
             fill = Origin)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  # coord_cartesian(ylim = c(0, 1.5e4)) +
  stat_summary(fun = mean,
               position = position_dodge(width = 0.9)) +
  scale_fill_brewer(palette = "Set1")


# re-create Table 3 in manuscript
param_summ %>%
  filter(grepl('pred_ovrshts', param)) %>%
  bind_cols(obs_ovrsht_org %>%
              filter(location == "at_PRA") %>%
              select(Year, Origin, ovrst_tags = juv_tags) %>%
              left_join(dwnstrm_est %>%
                          select(Year, Origin, dwnstrm_escp = Estimate, dwnstrm_se = SE))) %>%
  select(Year, Origin, 
         ovst_est = mean,
         ovst_lci = `2.5%`,
         ovst_uci = `97.5%`) %>%
  left_join(param_summ %>%
              filter(grepl('phi', param)) %>%
              bind_cols(obs_ovrsht_org %>%
                          filter(location == "at_PRA") %>%
                          select(Year, Origin, ovrst_tags = juv_tags) %>%
                          left_join(dwnstrm_est %>%
                                      select(Year, Origin, dwnstrm_escp = Estimate, dwnstrm_se = SE))) %>%
              select(Year, Origin, 
                     phi_est = mean,
                     phi_lci = `2.5%`,
                     phi_uci = `97.5%`)) %>%
  filter(Origin == "Wild") %>%
  write_csv("outgoing/Wild_Overshoots_Survival.csv")


plot_df = post_df %>%
  filter(param_fam == "mu") %>%
  left_join(mod_data %>%
              mutate(obs = log(escp_est)) %>%
              select(Year, Origin,
                     tags_est,
                     obs) %>%
              mutate(param_num = 1:n())) %>%
  arrange(CHAIN, ITER,
          tags_est) %>%
  filter(ITER == 51050,
         CHAIN == 1) %>%
  ggplot(aes(x = tags_est,
             y = value)) +
  geom_line(color = 'lightgray',
            alpha = 0.5) +
  geom_point()

mod_data %>%
  ggplot(aes(x = log(tags_est),
             y = log(escp_est))) +
  geom_abline(data = post_df %>%
                filter(param_fam == "beta") %>%
                select(CHAIN:value) %>%
                pivot_wider(names_from = "param",
                            values_from = "value"),
              aes(intercept = `beta[1]`,
                  slope = `beta[2]`),
              color = 'lightgray',
              alpha = 0.1) +
  geom_point() +
  geom_abline(data = param_summ %>%
                filter(grepl('beta', param)),
              aes(intercept = mean[param == "beta[1]"],
                  slope = mean[param == "beta[2]"]),
              color = 'blue',
              lwd = 2) +
  geom_smooth(method = lm,
              formula = y ~ x,
              color = 'red',
              se = F,
              fullrange = T) +
  labs(x = 'Adults Tagged as Juveniles\nthat Returned Downstream',
       y = 'POM Sucessful Overshoots')

#--------------------------------------------------
# diagnostic plots
#--------------------------------------------------
library(ggmcmc)


my_fam = "beta"
my_fam = 'sigma'
my_fam = "^u"
my_fam = "^v"
my_fam = "phi"
my_fam = "pred_ovrshts"

my_ggs = ggs(post,
             family = my_fam)

ggs_traceplot(my_ggs) +
  facet_wrap(~ Parameter)
ggs_density(my_ggs) +
  facet_wrap(~ Parameter,
             scales = 'free')
ggs_Rhat(my_ggs)
ggs_geweke(my_ggs)

ggs_diagnostics(ggs(post))
ggs_diagnostics(my_ggs) %>%
  filter(Diagnostic == "Rhat")

# black lines are priors

# ggs(post,
#     family = "^u") %>%
#   ggs_density() +
#   facet_wrap(~ Parameter) +
#   stat_function(fun = dbeta,
#                 args = list(shape1 = 1,
#                             shape2 = 1),
#                 color = 'black',
#                 linewidth = 2)
# 
# ggs(post,
#     family = "^v") %>%
#   ggs_density() +
#   facet_wrap(~ Parameter) +
#   stat_function(fun = dgamma,
#                 args = list(shape = 1,
#                             scale = 20),
#                 color = 'black',
#                 linewidth = 2)

ggs(post,
    family = "beta") %>%
  ggs_density() +
  facet_wrap(~ Parameter) +
  stat_function(fun = dt,
                args = list(df = 1),
                color = 'black',
                linewidth = 2)

ggs(post,
    family = "sigma") %>%
  ggs_density() +
  facet_wrap(~ Parameter) +
  stat_function(fun = dt,
                args = list(df = 1),
                color = 'black',
                linewidth = 2)


ggs(post,
    family = "phi") %>%
  mutate(param_num = str_extract(Parameter, "[:digit:]+"),
         param_num = as.integer(param_num)) %>%
  left_join(obs_ovrsht_org %>%
              filter(location == "at_PRA") %>%
              select(Year, Origin) %>%
              mutate(param_num = 1:n())) %>%
  ggs_density() +
  # facet_wrap(~ Parameter) +
  facet_wrap(~ Origin + Year) +
  stat_function(fun = dunif,
                args = list(min = 0,
                            max = 1),
                color = 'black',
                linewidth = 2)

#--------------------------------------------------
# compare counts at Priest and estimates across UC pops
#--------------------------------------------------
uc_pops = as.list(2011:2018) %>%
  rlang::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           read_excel(paste0('data/DABOM_results/PRA_Steelhead_', x[1], '_20190610.xlsx'),
                      "Population Escapement") %>%
             janitor::clean_names()
         }) %>%
  filter(population %in% c("Wenatchee",
                           "Entiat",
                           "Methow",
                           "Okanogan")) %>%
  group_by(Year, Origin = origin) %>%
  summarise_at(vars(uc_pops = estimate),
               list(sum))


pra_total = as.list(2011:2018) %>%
  rlang::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           read_excel(paste0('data/DABOM_results/PRA_Steelhead_', x[1], '_20190610.xlsx'),
                      "All Escapement") %>%
             janitor::clean_names()
         }) %>%
  filter(param %in% c("past_RIA",
                      'dwnStrm',
                      "PRA_bb")) %>%
  group_by(Year, Origin = origin) %>%
  summarise_at(vars(PRA = estimate),
               list(sum))

pra_total %>%
  full_join(uc_pops) %>%
  mutate_at(vars(Year),
            list(as.numeric)) %>%
  mutate(Origin = recode(Origin,
                         "Natural" = "Wild")) %>%
  full_join(param_summ %>%
              filter(grepl('pred_ovrshts', param)) %>%
              bind_cols(obs_ovrsht_org %>%
                          filter(location == "at_PRA") %>%
                          select(Year, Origin)) %>%
              select(Year, Origin,
                     pred_ovrst = mean,
                     pred_median = `50%`)) %>%
  mutate(uc_plus_ovrst = uc_pops + pred_ovrst,
         leftover = PRA - uc_plus_ovrst,
         rel_diff = leftover / PRA) %>%
  arrange(Origin, Year) %>%
  # group_by(Origin) %>%
  # summarise_at(vars(leftover, rel_diff),
  #              list(mean))
  ggplot(aes(x = Year)) +
  geom_line(aes(y = PRA,
                color = "At Priest")) +
  geom_point(aes(y = PRA,
                 color = "At Priest")) +
  geom_line(aes(y = uc_plus_ovrst,
                color = "UC Pops plus Est. Overshoots")) +
  geom_point(aes(y = uc_plus_ovrst,
                 color = "UC Pops plus Est. Overshoots")) +
  scale_color_brewer(palette = "Set1",
                     name = "Source") +
  theme(legend.position = "bottom") +
  facet_wrap(~ Origin,
             scales = 'fixed') +
  labs(y = "Total Steelhead")


pra_total %>%
  mutate_at(vars(Year),
            list(as.numeric)) %>%
  mutate(Origin = recode(Origin,
                         "Natural" = "Wild")) %>%
  full_join(param_summ %>%
              filter(grepl('pred_ovrshts', param)) %>%
              bind_cols(obs_ovrsht_org %>%
                          filter(location == "at_PRA") %>%
                          select(Year, Origin)) %>%
              select(Year, Origin,
                     pred_ovrst = mean,
                     pred_median = `50%`)) %>%
  ungroup() %>%
  filter(Origin == "Wild") %>%
  mutate_at(vars(starts_with("pred")),
            list(~ . / PRA)) %>%
  summarise_at(vars(starts_with("pred")),
               list(mean = mean,
                    sd = sd,
                    min = min,
                    max = max))


param_summ %>%
  filter(grepl('phi', param)) %>%
  bind_cols(obs_ovrsht_org %>%
              filter(location == "at_PRA") %>%
              select(Year, Origin)) %>%
  select(Year, Origin,
         pred_mean = mean,
         pred_median = `50%`,
         pred_cv = cv) %>%
  filter(Origin == "Wild") %>%
  summarise_at(vars(pred_cv),
               list(mean = mean,
                    sd = sd,
                    min = min,
                    max = max))

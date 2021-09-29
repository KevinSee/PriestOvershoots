# Author: Kevin See
# Purpose: Calculate adult overshoots at Priest Rapids dam
# Created: 9/19/19
# Last Modified: 10/4/19
# Notes: 

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(readxl)
library(magrittr)
library(msm)
library(ggrepel)

theme_set(theme_bw())

#-----------------------------------------------------------------
# read in data
data = read_excel('data/overshoot estimates.xlsx',
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
                     location = 'at_PRA'))

data

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
            list(as.numeric))
dwnstrm_est %>%
  filter(Origin == 'Natural') %>%
  # mutate(CV = SE / Estimate)
  summarise(mean_var = mean(SE^2),
            mean_se = mean(SE),
            mean_cv = mean(SE / Estimate))


data %<>%
  select(-pom_overshoots) %>%
  left_join(dwnstrm_est %>%
              mutate(origin = recode(Origin,
                                     'Natural' = 'Wild')) %>%
              select(-Population, -Origin) %>%
              rename(pom_overshoots = Estimate,
                     pom_se = SE,
                     lwr95 = lowerCI,
                     upr95 = upperCI) %>%
              mutate(location = 'downstream_PRA')) %>%
  select(origin, location, Year, juv_tags, pom_overshoots, pom_se, lwr95, upr95)


p = data %>%
  filter(location == 'downstream_PRA') %>%
  ggplot(aes(x = juv_tags,
             y = pom_overshoots)) +
  geom_point() +
  geom_smooth(method = lm,
              formula = y ~ x - 1,
              se = F,
              fullrange = T) +
  geom_smooth(method = lm,
              formula = y ~ x,
              color = 'red',
              se = F,
              fullrange = T) +
  facet_wrap(~ origin,
             scales = 'free') +
  labs(x = 'Adults Tagged as Juveniles',
       y = 'POM Overshoots')
p

ggsave('figures/Regressions.pdf',
       p,
       width = 7,
       height = 7)



# get estimates of downstream escapement to specific arrays
site_df = PITcleanr::writePRDNodeNetwork()
site_df %>% 
  filter(Step2 == 'BelowPriest') %>%
  filter(Step4 == '')

# use PRV instead of breaking out by HST and MDR
dwnstrm_sites = as.list(2011:2018) %>%
  rlang::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           read_excel(paste0('/Users/seek/Documents/GitProjects/MyProjects/DABOM_PriestRapids_Sthd_old/outgoing/PITcleanr/PRA_Steelhead_', x[1], '_20190610.xlsx'),
                      'All Escapement') %>%
             filter(param %in% c('BelowJD1',
                                 'past_RSH',
                                 'past_PRH',
                                 'past_JD1',
                                 'past_TMF',
                                 'past_ICH',
                                 'past_PRO',
                                 'past_PRV'))
         }) %>%
  mutate_at(vars(Year),
            list(as.numeric))

dwnstrm_sites %>%
  write_csv('outgoing/DownstreamEstimates.csv')

#--------------------------------------------------------------
# focus on wild fish, using data from downstream of PRA
mod_df = data %>%
  filter(origin == 'Wild',
         location == 'downstream_PRA') %>%
  group_by(origin) %>%
  nest() %>%
  mutate(lin_mod = map(.x = data,
                       .f = function(x) {
                         lm(pom_overshoots ~ -1 + juv_tags,
                         data = x)
                       }),
         lin_mod_int = map(.x = data,
                           .f = function(x) {
                             lm(pom_overshoots ~ juv_tags,
                                data = x)
                           }),
         summ = map(.x = lin_mod,
                    .f = summary)) %>%
  mutate(intercept = map_dbl(.x = lin_mod,
                             .f = function(x) coef(x)[match('(Intercept)', names(coef(x)))]),
         slope = map_dbl(.x = lin_mod,
                         .f = function(x) coef(x)[match('juv_tags', names(coef(x)))]),
         slope_se = map_dbl(.x = summ,
                            .f = function(x) coefficients(x)[1,2]),
         R2 = map_dbl(.x = summ,
                      .f = 'r.squared'),
         adj_R2 = map_dbl(.x = summ,
                          .f = 'adj.r.squared'))

mod_df

wild_p = mod_df %>%
  ungroup() %>%
  select(data) %>%
  unnest(cols = c(data)) %>%
  ggplot(aes(x = juv_tags,
             y = pom_overshoots)) +
  geom_point() +
  geom_smooth(method = lm,
              formula = y ~ x - 1,
              se = T,
              fullrange = T) +
  geom_text_repel(aes(label = Year)) +
  geom_text(x = 5,
            y = 1700,
            vjust = 0, 
            hjust = 0,
            parse = T,
            label = paste("r^2 ==", round(mod_df$R2, 3))) +
  labs(x = 'Successful Overshoot Adults Tagged as Juveniles\nDetected at Priest Rapids Dam',
       y = 'Total Successful Overshoots (from POM)')

ggsave('figures/wild_regressions.pdf',
       wild_p,
       width = 7,
       height = 7)

# using only the model for wild downstream of PRA
pred_df = data %>%
  select(origin, location, Year, juv_tags, pom_overshoots, lwr95, upr95)

# add some additional uncertainty from the estimates of overshoots
# this is the residual variance from the model
var_org = summary(mod_df$lin_mod[[1]])$sigma^2

# this is the average variance in our estimates of overshoots that make it downstream
# var_ovst = 0
var_ovst = dwnstrm_est %>%
  # filter(Origin == 'Natural') %>%
  summarise(mean_var = mean(SE^2)) %>%
  pull(mean_var)


pred_int = predict(mod_df$lin_mod[[1]],
                   newdata = pred_df,
                   # interval = 'prediction',
                   interval = 'confidence',
                   level = 0.95,
                   se.fit = T,
                   scale = sqrt(var_org + var_ovst),
                   df = summary(mod_df$lin_mod[[1]])$df[2])

all_preds = pred_df %>%
  bind_cols(pred_int$fit %>%
              as_tibble() %>%
              mutate(se_overshoot = pred_int$se.fit) %>%
              select(fit, se_overshoot, lwr, upr)) %>%
  rename(pred_overshoot = fit) %>%
  mutate_at(vars(pred_overshoot, se_overshoot, lwr, upr),
            list(round)) %>%
  mutate_at(vars(lwr),
            list(~ if_else(. < 0, 0, .)))


all_preds %>%
  # filter(location == 'downstream_PRA')
  filter(location == 'at_PRA')

all_preds %>%
  filter(location == 'downstream_PRA') %>%
  mutate(inCI = if_else(lwr <= pom_overshoots & upr >= pom_overshoots, T, F))


pred_cone = tibble(juv_tags = seq(0, max(data$juv_tags))) %>%
  bind_cols(predict(mod_df$lin_mod[[1]],
                    newdata = .,
                    interval = 'prediction',
                    # interval = 'confidence',
                    level = 0.95,
                    se.fit = F,
                    scale = sqrt(var_org + var_ovst),
                    df = summary(mod_df$lin_mod[[1]])$df[2]) %>%
              as_tibble()) %>%
  mutate_at(vars(fit, lwr),
            list(~ if_else(. < 0, 0, .)))



all_preds %<>%
  mutate(location = recode(location,
                           'at_PRA' = 'At Priest',
                           'downstream_PRA' = 'Downstream'))

all_preds %>%
  # filter(origin == 'Wild') %>%
  ggplot(aes(x = juv_tags)) +
  # geom_smooth(method = lm,
  #             aes(y = pom_overshoots),
  #             formula = y ~ x - 1,
  #             se = T,
  #             color = 'green',
  #             fill = 'green',
  #             fullrange = T) +
  geom_ribbon(data = pred_cone,
              aes(ymin = lwr,
                  ymax = upr),
              color = 'lightgray',
              alpha = 0.5) +
  geom_abline(slope = mod_df$slope,
              intercept = 0,
              color = 'blue') +
  geom_errorbar(aes(ymin = lwr95,
                    ymax = upr95,
                    color = 'POM')) +
  geom_point(aes(y = pom_overshoots,
                 color = 'POM',
                 shape = 'POM')) +
  geom_point(aes(y = pred_overshoot,
                 shape = 'regression',
                 color = 'regression')) +
  geom_text_repel(aes(y = pred_overshoot,
                      label = Year)) +
  geom_text_repel(data = dwnstrm_est %>%
                    mutate(location = 'Downstream',
                           origin = recode(Origin,
                                           'Natural' = 'Wild')) %>%
                    # filter(origin == 'Wild') %>%
                    select(Year, location, origin, Estimate:upperCI) %>%
                    left_join(all_preds %>%
                                select(Year, location, origin, juv_tags)) %>%
                    filter(Year %in% c(2013, 2015, 2017)),
                  aes(y = Estimate,
                      label = Year,
                      color = 'POM')) +
  geom_errorbar(aes(ymin = lwr,
                    ymax = upr)) +
  scale_color_manual(values = c('POM' = 'red',
                                'regression' = 'black')) +
  scale_shape_manual(values = c('POM' = 1,
                                'regression' = 19)) +
  facet_grid(origin ~ location,
             scales = 'free') +
  labs(x = 'Overshoot Adults Tagged as Juveniles\nDetected at Priest Rapids Dam',
       y = 'Estimated Overshoots at Priest',
       color = 'Source',
       shape = 'Source')

# save a csv
all_preds %>%
  rename(pred_lwr95 = lwr,
         pred_upr95 = upr) %>%
  select(-se_overshoot) %>%
  write_csv('outgoing/Overshoot_Estimates.csv')

#-----------------
# now estimate conversion (survival) rate from Priest to downstream spawning locations (grouped)

conv_df = all_preds %>%
  filter(location == 'At Priest') %>%
  select(Year, location, origin, est = pred_overshoot) %>%
  bind_rows(dwnstrm_est %>%
              mutate(location = 'Downstream',
                     origin = recode(Origin,
                                     'Natural' = 'Wild')) %>%
              select(Year, location, origin,
                     est = Estimate)) %>%
  spread(location, est) %>%
  rename(PRA = `At Priest`,
         Dwn = Downstream) %>%
  left_join(all_preds %>%
              filter(location == 'At Priest') %>%
              select(Year, location, origin, se = se_overshoot) %>%
              bind_rows(dwnstrm_est %>%
                          mutate(location = 'Downstream',
                                 origin = recode(Origin,
                                                 'Natural' = 'Wild')) %>%
                          select(Year, location, origin,
                                 se = SE)) %>%
              spread(location, se) %>%
              rename(PRA_se = `At Priest`,
                     Dwn_se = `Downstream`)) %>%
  mutate(Conv_rate = Dwn / PRA) %>%
  rowwise() %>%
  mutate(Conv_rate_se = deltamethod(~ x1 / x2,
                                    mean = c(Dwn, PRA),
                                    cov = diag(c(Dwn_se, PRA_se)^2))) %>%
  ungroup() %>%
  mutate(lwr = Conv_rate + qnorm(0.025) * Conv_rate_se,
         upr = Conv_rate + qnorm(0.975) * Conv_rate_se) #%>%
# mutate(lwr = if_else(lwr < 0, 0, lwr),
#        upr = if_else(upr > 1, 1, upr),
#        Conv_rate = if_else(Conv_rate > 1, 1, Conv_rate))

conv_df %>%
  filter(Conv_rate > 1)

# save a csv
conv_df %>%
  write_csv('outgoing/Conversion_Rate.csv')

conv_df %>%
  mutate(Conv_rate_cv = Conv_rate_se / Conv_rate) %>%
  group_by(origin) %>%
  summarise_at(vars(Conv_rate, Conv_rate_cv),
               list(mean = mean,
                    sd = sd))

#--------------------------------------------------------------
# Look at detection estimates
# calculate joint detection probability for downstream sites
#--------------------------------------------------------------
# These averages go in Table 2

dwn_det = as.list(2011:2017) %>%
  rlang::set_names() %>%
  map_df(.id = "Year",
         .f = function(x) {
           read_excel(paste0('data/DABOM_results/PRA_Steelhead_', x[1], '_20190610.xlsx'),
                      "Detection")
         }) %>%
  mutate_at(vars(Year),
            list(as.numeric)) %>%
  filter(grepl('BelowJD', Node) |
           grepl('JD1', Node) |
           grepl('TMF', Node) |
           grepl('PRV', Node) |
           grepl('HST', Node) |
           grepl('MDR', Node) |
           grepl('ICH', Node) |
           grepl('PRO', Node) |
           grepl('RSH', Node) |
           grepl('PRH', Node)) %>%
  mutate(Site = str_remove(Node, 'B0$'),
         Site = str_remove(Site, 'A0$'))

dwn_tags = as.list(2011:2017) %>%
  rlang::set_names() %>%
  map_df(.id = "Year",
         .f = function(x) {
           load(paste0('data/DABOMready/UC_Steelhead_', x, '.rda'))
           proc_list$ProcCapHist %>%
             filter(AutoProcStatus) %>%
             filter(grepl('BelowJD', Node) |
                      grepl('JD1', Node) |
                      grepl('TMF', Node) |
                      grepl('PRV', Node) |
                      grepl('HST', Node) |
                      grepl('MDR', Node) |
                      grepl('ICH', Node) |
                      grepl('PRO', Node) |
                      grepl('RSH', Node) |
                      grepl('PRH', Node)) %>%
             select(TagID, Origin, SiteID, Node) %>%
             distinct()
         }) %>%
  mutate_at(vars(Year),
            list(as.numeric)) %>%
  mutate(Site = str_remove(Node, 'B0$'),
         Site = str_remove(Site, 'A0$'))

dwn_tags %>%
  group_by(Site, Year) %>%
  summarise(n_tags = n_distinct(TagID)) %>%
  ungroup() %>%
  spread(Year, n_tags,
         fill = 0)


joint_det = dwn_det %>%
  arrange(Year, Site, desc(Node)) %>%
  group_by(Year, Site) %>%
  mutate(node_num = paste("node", 1:n(), sep = "_")) %>%
  select(Year, Site, node_num, Estimate) %>%
  spread(node_num, Estimate,
         fill = 0) %>%
  left_join(dwn_det %>%
              arrange(Year, Site, desc(Node)) %>%
              group_by(Year, Site) %>%
              mutate(node_num = paste("node", 1:n(), "se", sep = "_")) %>%
              select(Year, Site, node_num, SE) %>%
              spread(node_num, SE,
                     fill = 0)) %>%
  rowwise() %>%
  mutate(joint_det = 1 - ( (1 - node_1) * (1 - node_2) ),
         joint_det_se = msm::deltamethod(~ 1 - ((1 - x1) * (1 - x2)),
                                         mean = c(node_1, node_2),
                                         cov = diag(c(node_1_se, node_2_se)^2))) %>%
  ungroup() %>%
  left_join(dwn_tags %>%
              group_by(Site, Year) %>%
              summarise(n_tags = n_distinct(TagID)) %>%
              ungroup()) %>%
  mutate(n_tags = if_else(is.na(n_tags),
                          as.integer(0),
                          n_tags)) %>%
  mutate(valid_est = if_else(n_tags > 0 &
                       (node_1 > 0 & node_2 > 0),
                       T, F)) %>%
  select(Year, Site, n_tags, valid_est, starts_with("joint"))




joint_det %>%
  filter(valid_est) %>%
  filter(Site %in% c('BelowJD1',
                     'JD1',
                     'TMF',
                     'PRV',
                     'ICH',
                     'PRO',
                     'RSH',
                     'PRH')) %>%
  group_by(Year) %>%
  summarise(yr_tags = sum(n_tags),
            mean_det = mean(joint_det),
            wgt_det = weighted.mean(joint_det, w = n_tags)) %>%
  ungroup() %>%
  summarise_at(vars(mean_det, wgt_det),
               list(mean))

joint_det %>%
  filter(Site == 'PRV') %>%
  filter(n_tags > 0) %>%
  group_by(Site) %>%
  summarise_at(vars(n_tags, starts_with("joint")),
               list(mean))

dwn_det %>%  
  filter(Site == 'PRV') %>%
  filter(n_tags > 0) %>%
  group_by(Site) %>%
  summarise_at(vars(n_tags, Estimate, SE),
               list(mean))

dwn_det %>%  
  filter(Site == 'TMF') %>%
  filter(n_tags > 0) %>%
  group_by(Site) %>%
  summarise_at(vars(n_tags, Estimate, SE),
               list(mean))

joint_det %>%
  mutate(est_tags = n_tags / joint_det) %>%
  filter(Site %in% c('BelowJD1',
                     'JD1',
                     'TMF',
                     'PRV',
                     'ICH',
                     'PRO',
                     'RSH',
                     'PRH')) %>%
  filter(!is.na(est_tags)) %>%
  group_by(Year) %>%
  summarise_at(vars(n_tags, est_tags),
               list(sum)) %>%
  mutate(diff = est_tags - n_tags,
         rel_diff = diff / n_tags)

joint_det %>%
  select(Year, Site, starts_with("joint")) %>%
  left_join(data %>%
              filter(location == 'downstream_PRA') %>%
              select(origin, Year, juv_tags))

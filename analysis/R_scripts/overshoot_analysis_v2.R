# Author: Kevin See
# Purpose: Calculate adult overshoots at Priest Rapids dam
# Created: 9/1/20
# Last Modified: 9/8/20
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
#-----------------------------------------------------------------
# original data of known overshoots
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
         everything())

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

# # get branch specific estimates of downstream escapement
# dwnstrm_loc_est = as.list(2011:2018) %>%
#   rlang::set_names() %>%
#   map_df(.id = 'Year',
#          .f = function(x) {
#            read_excel(paste0('data/DABOM_results/PRA_Steelhead_', x[1], '_20190610.xlsx'),
#                       2) %>%
#              filter(param %in% c("BelowJD1",
#                                  'past_RSH',
#                                  'past_PRH',
#                                  'past_JD1',
#                                  'past_TMF',
#                                  'past_ICH',
#                                  'past_PRO',
#                                  'past_PRV')) %>%
#              rename(site = param) %>%
#              mutate(site = str_remove(site, "^past_"))
#          }) %>%
#   mutate_at(vars(Year),
#             list(as.numeric)) %>%
#   mutate(Origin = recode(Origin,
#                          'Natural' = 'Wild')) %>%
#   rename(pom_est = Estimate,
#          pom_se = SE,
#          pom_ci_low = lowerCI,
#          pom_ci_high = upperCI)


# # 
# est_ovrsht %>%
#   full_join(dwnstrm_loc_est %>%
#               filter(Origin == 'Wild')) %>%
#   filter(pom_est == 0 | is.na(obs_ovrst_tags)) %>%
#   group_by(Year) %>%
#   summarise_at(vars(obs_ovrst_tags,
#                     pom_est),
#                list(sum),
#                na.rm = T)


# missing overshoot tags?
est_ovrsht %>%
  group_by(Year, Origin) %>%
  summarise(n_obs_tags = sum(obs_ovrst_tags),
            n_est_tags = sum(est_ovrst_tags),
            se_est_tags = sqrt(sum(se_ovrst_tags^2))) %>%
  ungroup() %>%
  full_join(obs_ovrsht_org %>%
              filter(Origin == 'Wild',
                     location == "downstream_PRA") %>%
              select(Year, Origin,
                     tot_obs_tags = juv_tags)) %>%
  mutate(n_miss_tags = tot_obs_tags - n_obs_tags)

# look at which fish ended up being dropped, because of their downstream detection location
# these match up with the "n_miss_tags" in the data.frame above
read_excel("data/known overshoot fallback locations.xlsx") %>%
  janitor::clean_names() %>%
  rename(notes = x5,
         site = downstream_obs_site) %>%
  tidyr::fill(total) %>%
  mutate(year = year + 1) %>%
  filter(!is.na(notes)) %>%
  arrange(year, site)
# only issue appears to be in 2015. Original data says there were 32 downstream tags, updated data only has 31


# # join escapement estimates and estimated overshoot tags, but only for places that have both estimates in a given year
# mod_data = est_ovrsht %>%
#   inner_join(dwnstrm_loc_est %>%
#               filter(Origin == "Wild")) %>%
#   filter(pom_est > 0) %>%
#   group_by(Year, Origin) %>%
#   summarise(tags_obs = sum(obs_ovrst_tags),
#             det_prob = mean(det_est),
#             det_prob_wgt = weighted.mean(det_est, w = est_ovrst_tags),
#             tags_est = sum(est_ovrst_tags),
#             tags_se = sqrt(sum(se_ovrst_tags^2)),
#             escp_est = sum(pom_est),
#             escp_se = sqrt(sum(pom_se^2))) %>%
#   ungroup()

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
  # full_join(dwnstrm_loc_est %>%
  #              filter(Origin == "Wild") %>%
  #              group_by(Year, Origin) %>%
  #              summarise(escp_est = sum(pom_est),
  #                        escp_se = sqrt(sum(pom_se^2))))
  



# quick plot
p = mod_data %>%
  ggplot(aes(x = tags_est,
             y = escp_est)) +
  geom_errorbar(aes(ymin = escp_est + escp_se * qnorm(0.025),
                    ymax = escp_est + escp_se * qnorm(0.975)),
                width = 0) +
  geom_errorbarh(aes(xmin = tags_est + tags_se * qnorm(0.025),
                     xmax = tags_est + tags_se * qnorm(0.975)),
                 height = 0) +
  geom_point() +
  geom_smooth(method = lm,
              formula = y ~ x - 1,
              se = F,
              fullrange = T) +
  # geom_smooth(method = lm,
  #             formula = y ~ x,
  #             color = 'red',
  #             se = F,
  #             fullrange = T) +
  scale_x_continuous(trans = "log",
                     breaks = scales::pretty_breaks()) +
  scale_y_continuous(trans = "log",
                     breaks = scales::pretty_breaks()) +
  labs(x = 'Adults Tagged as Juveniles\nthat Returned Downstream',
       y = 'POM Sucessful Overshoots')
p

# on a regular scale, showing log-log curve
p2 = mod_data %>%
  ggplot(aes(x = tags_est,
             y = escp_est)) +
  geom_errorbar(aes(ymin = escp_est + escp_se * qnorm(0.025),
                    ymax = escp_est + escp_se * qnorm(0.975)),
                width = 0) +
  geom_errorbarh(aes(xmin = tags_est + tags_se * qnorm(0.025),
                     xmax = tags_est + tags_se * qnorm(0.975)),
                 height = 0) +
  geom_point() +
  geom_smooth(method = glm,
              formula = y ~ log(x),
              method.args = list(family = gaussian(link = 'log')),
              se = F,
              fullrange = T) +
  geom_smooth(method = lm,
              formula = y ~ -1 + x,
              color = "red",
              se = F) +
  labs(x = 'Adults Tagged as Juveniles\nthat Returned Downstream',
       y = 'POM Sucessful Overshoots')

p2

#--------------------------------------------------------------
# focus on wild fish, using data from downstream of PRA


# # fit a linear model
# mod_df = mod_data %>%
#   group_by(Origin) %>%
#   nest() %>%
#   mutate(lin_mod = map(.x = data,
#                        .f = function(x) {
#                          lm(escp_est ~ -1 + tags_est,
#                             data = x)
#                        }),
#          lin_mod_int = map(.x = data,
#                            .f = function(x) {
#                              lm(escp_est ~ tags_est,
#                                 data = x)
#                            }),
#          summ = map(.x = lin_mod,
#                     .f = summary)) %>%
#   mutate(intercept = map_dbl(.x = lin_mod_int,
#                              .f = function(x) coef(x)[match('(Intercept)', names(coef(x)))]),
#          slope_int = map_dbl(.x = lin_mod_int,
#                              .f = function(x) coef(x)[match('tags_est', names(coef(x)))]),
#          slope = map_dbl(.x = lin_mod,
#                          .f = function(x) coef(x)[match('tags_est', names(coef(x)))]),
#          slope_se = map_dbl(.x = summ,
#                             .f = function(x) coefficients(x)[1,2]),
#          R2 = map_dbl(.x = summ,
#                       .f = 'r.squared'),
#          adj_R2 = map_dbl(.x = summ,
#                           .f = 'adj.r.squared'))
# 
# 
# mod_log_df = mod_data %>%
#   mutate(escp_est_log = log(escp_est),
#          tags_est_log = log(tags_est)) %>%
#   group_by(Origin) %>%
#   nest() %>%
#   mutate(lin_mod = map(.x = data,
#                        .f = function(x) {
#                          lm(escp_est_log ~ tags_est_log,
#                             data = x)
#                        }),
#          glm_mod = map(.x = data,
#                        .f = function(x) {
#                          glm(escp_est ~ log(tags_est),
#                              data = x,
#                              family = gaussian(link = "log"))
#                        }),
#          lin_summ = map(.x = lin_mod,
#                         .f = summary),
#          glm_summ = map(.x = glm_mod,
#                         .f = summary)) %>%
#   mutate(coefs = map(.x = lin_mod,
#                      .f = coefficients),
#          int = map_dbl(coefs,
#                        .f = function(x) x[1]),
#          slope = map_dbl(coefs,
#                          .f = function(x) x[2]),
#          slope_se = map_dbl(.x = lin_summ,
#                             .f = function(x) coefficients(x)[1,2]),
#          R2 = map_dbl(.x = lin_summ,
#                       .f = 'r.squared'),
#          adj_R2 = map_dbl(.x = lin_summ,
#                           .f = 'adj.r.squared'))


mod_df = mod_data %>%
  group_by(Origin) %>%
  nest() %>%
  ungroup() %>%
  mutate(lin_mod = map(.x = data,
                       .f = function(x) {
                         lm(escp_est ~ tags_est,
                            data = x)
                       }),
         log_mod = map(.x = data,
                        .f = function(x) {
                          lm(log(escp_est) ~ log(tags_est),
                             data = x)
                        })) %>%
  mutate(glm_mod = map(.x = data,
                       .f = function(x) {
                         glm(escp_est ~ log(tags_est),
                            data = x,
                            family = gaussian(link = "log"))
                       })) %>%
  pivot_longer(cols = ends_with('mod'),
               names_to = "type",
               values_to = "model") %>%
  mutate(summ = map(model,
                    summary),
         coefs = map(model,
                     coef)) %>%
  mutate(pred_data = list(obs_ovrsht_org %>%
                            filter(location == "at_PRA") %>%
                            select(Year, Origin, tags_est = juv_tags) %>%
                            left_join(dwnstrm_est %>%
                                        select(Year, Origin, dwnstrm_escp = Estimate, dwnstrm_se = SE)) %>%
                            arrange(Origin, Year))) %>%
  # mutate(pred_data = list(obs_ovrsht_org %>%
  #                           # assume 100% detection at Priest
  #                           rename(tags_est = juv_tags))) %>%
  mutate(preds = map2(model,
                      pred_data,
                      .f = function(x, y) {
                        p = predict(x,
                                newdata = y,
                                type = "response",
                                se.fit = T)
                        tibble(pred = p$fit,
                               pred_se = p$se.fit)
                      }))

mod_df %>%
  mutate(coefs = map(coefs, enframe)) %>%
  select(type, coefs) %>%
  unnest(cols = coefs) %>%
  pivot_wider(names_from = "type",
              values_from = "value")


mod_df %>%
  select(type, pred_data, preds) %>%
  unnest(cols = c(pred_data, preds)) %>%
  select(-pred_se) %>%
  mutate(pred = if_else(type == 'log_mod',
                        exp(pred),
                        pred)) %>%
  # filter(pred < dwnstrm_escp) %>%
  pivot_wider(names_from = "type",
              values_from = "pred") %>%
  arrange(Origin, Year) %>%
  rename(obs = dwnstrm_escp) %>%
  select(-dwnstrm_se) %>%
  pivot_longer(cols = -(Year:tags_est),
               names_to = c("source"),
               values_to = "overshoots") %>%
  ggplot(aes(x = tags_est,
             y = overshoots,
             color = source)) +
  geom_point() +
  geom_smooth(se = F) +
  facet_wrap(~ Origin)
  
mod_df %>%
  select(type, pred_data, preds) %>%
  unnest(cols = c(pred_data, preds)) %>%
  select(-pred_se) %>%
  mutate(pred = if_else(type == 'log_mod',
                        exp(pred),
                        pred)) %>%
  pivot_wider(names_from = "type",
              values_from = "pred") %>%
  mutate(lin_surv = dwnstrm_escp / lin_mod,
         log_surv = dwnstrm_escp / log_mod,
         glm_surv = dwnstrm_escp / glm_mod)
  

library(ggfortify)
autoplot(mod_df$model[[1]], which = 1:6)
autoplot(mod_df$model[[2]], which = 1:6)
autoplot(mod_df$model[[3]], which = 1:6)

mod_df$coefs[[1]]
mod_df$coefs[[2]]
mod_df$coefs[[3]]

exp(mod_df$coefs[[1]][1])

# curve(2*x^4.55, xlim = c(0, 35))

curve(exp(mod_df$coefs[[2]][1])*x^mod_df$coefs[[2]][2], 
      xlim = c(0, 35),
      col = "blue")

curve(exp(mod_df$coefs[[3]][1])*x^mod_df$coefs[[3]][2], 
      add = T,
      col = "red")

#--------------------------------------------------------------
# predict number of overshoots at Priest
#--------------------------------------------------------------
# get detections at Priest
pred_df = obs_ovrsht_org %>%
  # assume 100% detection at Priest
  rename(tags_est = juv_tags) %>%
  mutate(tags_est_log = log(tags_est),
         escp_est_log = log(pom_overshoots))

# add some additional uncertainty from the estimates of overshoots
# this is the residual variance from the model
var_org = summary(mod_log_df$lin_mod[[1]])$sigma^2

# this is the average variance in our estimates of overshoots that make it downstream
# var_ovst = 0
var_ovst = dwnstrm_est %>%
  # filter(Origin == 'Natural') %>%
  mutate(escp_est_log = log(Estimate),
         escp_var_log = log(SE^2)) %>%
  summarise(mean_var = mean(escp_var_log)) %>%
  pull(mean_var)


pred_int = predict(mod_log_df$lin_mod[[1]],
                   newdata = pred_df,
                   # interval = 'prediction',
                   interval = 'confidence',
                   level = 0.95,
                   se.fit = T,
                   scale = sqrt(var_org + var_ovst),
                   df = summary(mod_log_df$lin_mod[[1]])$df[2])

all_preds = pred_df %>%
  bind_cols(pred_int$fit %>%
              as_tibble() %>%
              mutate(se_overshoot = pred_int$se.fit) %>%
              select(fit, se_overshoot, 
                     lwr, upr)) %>%
  rename(pred_overshoot = fit)

all_preds %>%
  # filter(Origin == "Wild",
  #        location == "downstream_PRA") %>%
  ggplot(aes(x = tags_est_log)) +
  # geom_ribbon(data = pred_cone,
  #             aes(ymin = lwr,
  #                 ymax = upr),
  #             color = 'lightgray',
  #             alpha = 0.5) +
  geom_point(aes(y = escp_est_log,
                 color = 'POM',
                 shape = 'POM')) +
  geom_point(aes(y = pred_overshoot,
                 shape = 'regression',
                 color = 'regression')) +
  geom_text_repel(aes(y = pred_overshoot,
                      label = Year)) +
  scale_color_manual(values = c('POM' = 'red',
                                'regression' = 'black')) +
  scale_shape_manual(values = c('POM' = 1,
                                'regression' = 19)) +
  facet_grid(Origin ~ location,
             scales = 'free') +
  scale_x_continuous(trans = "exp") +
  scale_y_continuous(trans = "exp") +
  labs(x = 'Overshoot Adults Tagged as Juveniles\nDetected at Priest Rapids Dam',
       y = 'Estimated Overshoots at Priest',
       color = 'Source',
       shape = 'Source')


pred_cone = tibble(tags_est = seq(1, max(obs_ovrsht_org$juv_tags)),
                   tags_est_log = log(tags_est)) %>%
  bind_cols(predict(mod_log_df$lin_mod[[1]],
                    newdata = .,
                    interval = 'prediction',
                    # interval = 'confidence',
                    level = 0.95,
                    se.fit = F,
                    scale = sqrt(var_org + var_ovst),
                    df = summary(mod_log_df$lin_mod[[1]])$df[2]) %>%
              as_tibble()) %>%
  mutate_at(vars(fit, lwr, upr),
            list(exp)) %>%
  mutate_at(vars(fit, lwr),
            list(~ if_else(. < 0, 0, .)))

all_preds %>%
  # filter(Origin == 'Wild') %>%
  ggplot(aes(x = tags_est)) +
  geom_ribbon(data = pred_cone,
              aes(ymin = lwr,
                  ymax = upr),
              color = 'lightgray',
              alpha = 0.5) +
  geom_abline(slope = mod_log_df$slope,
              intercept = 0,
              color = 'blue') +
  # geom_errorbar(aes(ymin = lwr95,
  #                   ymax = upr95,
  #                   color = 'POM')) +
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
                           Origin = recode(Origin,
                                           'Natural' = 'Wild')) %>%
                    # filter(Origin == 'Wild') %>%
                    select(Year, location, Origin, Estimate:upperCI) %>%
                    left_join(all_preds %>%
                                select(Year, location, Origin, tags_est)) %>%
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
  facet_grid(Origin ~ location,
             scales = 'free') +
  labs(x = 'Overshoot Adults Tagged as Juveniles\nDetected at Priest Rapids Dam',
       y = 'Estimated Overshoots at Priest',
       color = 'Source',
       shape = 'Source')


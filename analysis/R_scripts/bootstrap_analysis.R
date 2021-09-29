# Author: Kevin See
# Purpose: Calculate adult overshoots at Priest Rapids dam
# Created: 9/2/20
# Last Modified: 9/2/20
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

#-----------------------------------------------------------------
# set up some bootstrapped datasets
#-----------------------------------------------------------------
# data to make predictions on
pred_df = obs_ovrsht_org %>%
  # assume 100% detection at Priest
  rename(tags_est = juv_tags) %>%
  filter(location == "at_PRA") %>%
  select(-pom_overshoots)

# %>%
#   left_join(dwnstrm_est %>%
#               select(Year, Origin,
#                      dwnstm_escp = Estimate,
#                      dwnstm_se = SE))

n_boot = 10
set.seed(5)

dwnstrm_boot = dwnstrm_est %>%
  nest(data = everything()) %>%
  crossing(iter = 1:n_boot) %>%
  select(iter, data) %>%
  mutate(escp_df = map(data,
                       .f = function(x) {
                         x %>%
                           rowwise() %>%
                           mutate(escp_est = rnorm(1, Estimate, SE)) %>%
                           select(Year, Origin, escp_est)
                       })) %>%
  select(iter, escp_df) %>%
  unnest(cols = escp_df)

boot_df = est_ovrsht %>%
  select(Year:obs_ovrst_tags, beta_alpha, beta_beta) %>%
  nest(data = everything()) %>%
  crossing(iter = 1:n_boot) %>%
  select(iter, data) %>%
  mutate(tags_df = map(data,
                       .f = function(x) {
                         x %>%
                           # rowwise() %>%
                           mutate(det_p = rbeta(1, beta_alpha, beta_beta),
                                  n_tags = obs_ovrst_tags / det_p,
                                  n_tags = round(n_tags)) %>%
                           group_by(Year) %>%
                           summarise(tags_est = sum(n_tags),
                                     .groups = "drop")
                       })) %>%
  select(iter, tags_df) %>%
  unnest(cols = tags_df) %>%
  left_join(dwnstrm_boot %>%
              filter(Origin == "Wild")) %>%
  group_by(iter) %>%
  nest() %>%
  left_join(pred_df %>%
              full_join(dwnstrm_boot) %>%
              group_by(iter) %>%
              nest() %>%
              rename(pred_data = data)) %>%
  mutate(log_mod = map(data,
                       .f = function(x) {
                         lm(log(escp_est) ~ log(tags_est),
                            data = x)
                       })) %>%
  # make predictions of total overshoots at Priest
  mutate(preds = map2(log_mod,
                      pred_data,
                      .f = function(x, y) {
                        p = predict(x,
                                    newdata = y,
                                    type = "response",
                                    se.fit = F)
                        y %>%
                          mutate(pred = p) %>%
                          mutate(pred_ovrst = exp(pred)) %>%
                          mutate(ovrst_surv = escp_est / pred_ovrst)
                      })) %>%
  # fit a GLM with log link
  mutate(glm_mod = map(data,
                       .f = function(x) {
                         glm(escp_est ~ log(tags_est),
                             data = x,
                             family = gaussian(link = "log"))
                       })) %>%
  # make predictions of total overshoots at Priest
  mutate(preds_glm = map2(glm_mod,
                      pred_data,
                      .f = function(x, y) {
                        p = predict(x,
                                    newdata = y,
                                    type = "response",
                                    se.fit = F)
                        y %>%
                          mutate(pred_ovrst = p) %>%
                          mutate(ovrst_surv = escp_est / pred_ovrst)
                      })) %>%
  ungroup()
  
  
  boot_df %>%
    select(iter, preds = preds_glm) %>%
    # select(iter, preds) %>%
    unnest(cols = preds) %>%
    # filter(ovrst_surv <= 1) %>%
    # group_by(Year, Origin) %>%
    # summarise(n_iter = n_distinct(iter)) %>%
    # arrange(n_iter)
    group_by(Year, Origin, tags_est) %>%
    summarise_at(vars(escp_est, pred_ovrst, ovrst_surv),
                 list(mean)) %>%
  filter(Origin == "Wild")
  summarise_at(vars(ovrst_surv),
               list(min = min,
                    mean = mean,
                    median = median,
                    max = max)) %>%
  arrange(Origin, Year)

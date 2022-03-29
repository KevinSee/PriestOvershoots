# Author: Kevin See
# Purpose: Logistic regression on overshoot success based on number of dams crossed
# Created: 9/15/20
# Last Modified: 3/28/2022
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(magrittr)
library(ggpubr)
library(here)
library(readxl)
library(boot)
library(ResourceSelection) # for goodness of fit tests

#-----------------------------------------------------------------
# here's the data, provided by Andrew Murdoch.
# All the numbers with dams > 0 are from fish tagged at PRD that return to a downstream trib
# for dam == 0, Andrew used known Yakima fish that crossed McNary.
# because slightly different data, we decided to only use PRD fish for the model

# dam_df = tibble(dams = c(5:3, 1, 0),
#                 fail = c(38,
#                          9,
#                          14,
#                          16,
#                          14),
#                 success = c(11,
#                             14,
#                             31,
#                             113,
#                             262)) %>%
#   mutate(n_tot = success + fail,
#          suc_perc = success / n_tot) %>%
#   arrange(dams) %>%
#   filter(dams > 0)

# data by year
dam_yr = read_excel(here("analysis/data/raw_data",
                         "dan's request by year.xlsx"),
                    range = "E2:L6") %>%
  mutate(dams = c(5:3, 1)) %>%
  pivot_longer(cols = -dams,
               names_to = "year",
               values_to = "fail") %>%
  left_join(read_excel(here("analysis/data/raw_data",
                            "dan's request by year.xlsx"),
                       range = "E10:L13",
                       col_names = as.character(2010:2017)) %>%
              mutate(dams = c(5:3, 1)) %>%
              pivot_longer(cols = -dams,
                           names_to = "year",
                           values_to = "success")) %>%
  mutate(across(year,
                as_factor)) %>%
  mutate(n_tot = success + fail,
         suc_perc = success / n_tot)

# expand to one row per individual
dam_individ <- dam_yr %>%
  select(dams:success) %>%
  pivot_longer(cols = c(fail, success),
               names_to = "outcome",
               values_to = "freq") %>%
  uncount(freq) %>%
  mutate(outcome_fct = recode(outcome,
                              "fail" = 0,
                              "success" = 1),
         across(outcome_fct,
                as.numeric))

# pool years
dam_df <- dam_yr %>%
  group_by(dams) %>%
  summarize(across(c(success, fail),
                   sum)) %>%
  mutate(n_tot = success + fail,
         suc_perc = success / n_tot) %>%
  arrange(dams)

# fit a binomial logistic model
# pooling years
binom_mod = glm(suc_perc ~ dams,
                data = dam_df,
                weights = dam_df$n_tot,
                family = "binomial")
null_mod = glm(suc_perc ~ 1,
               data = dam_df,
               weights = dam_df$n_tot,
               family = "binomial")
# summarizing within years
binom_mod = glm(suc_perc ~ dams,
                data = dam_yr,
                family = "binomial",
                weights = dam_yr$n_tot)
null_mod = glm(suc_perc ~ 1,
               data = dam_yr,
               family = "binomial",
               weights = dam_yr$n_tot)
# using individual outcomes
binom_mod = glm(outcome_fct ~ dams,
                data = dam_individ,
                family = "binomial")
null_mod = glm(outcome_fct ~ 1,
                data = dam_individ,
                family = "binomial")




summary(binom_mod)
hoslem.test(binom_mod$y, fitted(binom_mod), g = 4)
pchisq(deviance(binom_mod), binom_mod$df.null, lower = F)

# looking at goodness of fit
anova(binom_mod,
      null_mod,
      test = "Chisq")

lmtest::lrtest(binom_mod,
               null_mod)
pscl::pR2(binom_mod)
MuMIn::r.squaredGLMM(binom_mod)

exp(coef(binom_mod)[2])
# odds of successful return get multiplied by 44% for every dam you cross

# predicted success if dams == 0
boot::inv.logit(coef(binom_mod)[1])
# matches observed 0.949 from Yakima fish crossing McNary

# looking at odds vs. prob of success
tibble(dams = seq(0, 5)) %>%
  mutate(odds = predict(binom_mod,
                        newdata = .,
                        type = "link"),
         odds = exp(odds)) %>%
  mutate(suc_perc = predict(binom_mod,
                            newdata = .,
                            type = "response"))

# make predictions of probability of success, with confidence intervals
pred_df <- tibble(dams = seq(0, 5, by = 0.1)) %>%
  bind_cols(predict(binom_mod,
                    newdata = .,
                    type = "link",
                    se = T) %>%
              magrittr::extract(1:2) %>%
              map_df(.id = 'pred',
                     .f = identity) %>%
              pivot_longer(cols = -pred,
                           names_to = "dams",
                           values_to = "log_odds") %>%
              pivot_wider(names_from = "pred",
                          values_from = "log_odds") %>%
              select(-dams)) %>%
  mutate(lowCI = fit + se.fit * qnorm(0.025),
         uppCI = fit + se.fit * qnorm(0.975)) %>%
  mutate_at(vars(pred_suc = fit, lowCI, uppCI),
            list(boot::inv.logit)) %>%
  select(dams,
         pred_suc,
         ends_with('CI'))

pred_df %>%
  filter(dams %in% c(0:5)) %>%
  mutate_at(vars(pred_suc:uppCI),
            list(round),
            digits = 3) %>%
  write_csv("outgoing/Dams_LogisticPredictions.csv")

# plot it
logis_p = dam_df %>%
  ggplot(aes(x = dams,
             y = suc_perc)) +
  geom_ribbon(data = pred_df,
              aes(ymin = lowCI,
                  ymax = uppCI,
                  y = pred_suc),
              color = NA,
              fill = 'gray90') +
  geom_line(data = pred_df,
            aes(y = pred_suc),
            color = "black",
            linetype = 2) +
  # geom_smooth(data = pred_df,
  #             method = loess,
  #             aes(y = pred_suc),
  #             color = "black",
  #             linetype = 2) +
  geom_point(size = 4) +
  # geom_point(aes(size = n_tot)) +
  theme_pubr(base_family = "serif",
             base_size = 10) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Number of Columbia River Dams",
       y = "Estimated Overshoot Fallback Probability")

ggsave(here("analysis/figures/Dams_Logistic.jpeg"),
       logis_p,
       width = 5,
       height = 5)

#-----------------------------------------------------------------
# using a mixed-effects model with random effect for year
library(lme4)
library(broom.mixed)
library(ggeffects)

# random intercept
m1 = glmer(cbind(success, fail) ~ dams + (1 | year),
           data = dam_yr,
           family = "binomial")

# random intercept and slope
m2 = glmer(cbind(success, fail) ~ dams + (1 + dams | year),
           data = dam_yr,
           family = "binomial")

# m2 = glmer(outcome_fct ~ dams + (1 + dams | year),
#            data = dam_individ,
#            family = "binomial")

# random slope
m3 = glmer(cbind(success, fail) ~ 1 + (dams - 1 | year),
           data = dam_yr,
           family = "binomial")

lapply(list(m1, m2, m3), isSingular)

anova(m1, m2, m3)

exp(fixef(m2)[2])
# odds of successful return get multiplied by 47% for every dam you cross

# predicted success if dams == 0
boot::inv.logit(fixef(m2)[1])
# matches observed 0.949 from Yakima fish crossing McNary

# pseudo-R2
MuMIn::r.squaredGLMM(m2)
# Hoslem GOF
hoslem.test(m2@resp$y, fitted(m2), g = 10)

tidy(m2,
     conf.int = T) %>%
  filter(term == "(Intercept)") %>%
  mutate(across(c(estimate, conf.low, conf.high),
                inv.logit))

tibble(year = rownames(ranef(m1)$year),
       b0 = as_vector(ranef(m1)$year + fixef(m1)[1])) %>%
  mutate(succ_dam0 = boot::inv.logit(b0))

ranef(m2) %>%
  as_tibble() %>%
  left_join(fixef(m2) %>%
              enframe(name = "term",
                      value = "fixed")) %>%
  mutate(beta = fixed + condval) %>%
  select(term, grp, beta) %>%
  pivot_wider(names_from = term,
              values_from = beta) %>%
  mutate(across(2,
                boot::inv.logit),
         across(dams,
                exp))

tibble(dams = seq(0, 5)) %>%
  mutate(odds = predict(m2,
                        newdata = .,
                        re.form = NA,
                        type = "link"),
         odds = exp(odds)) %>%
  mutate(suc_perc = predict(m2,
                            newdata = .,
                            re.form = NA,
                            type = "response"))

ctab <- confint(m2,
              parm = "beta_") %>%
  as_tibble(rownames = "param") %>%
  mutate(est = fixef(m2)) %>%
  relocate(est,
           .after = 1)
ctab %>%
  mutate(across(-param,
                exp))

tidy(m1,
     conf.int = T,
     exponentiate = F,
     effects = "fixed")

# make predictions of probability of success, with confidence intervals
pred_df <- ggpredict(m2,
                     terms = c("dams [0:5 by = 0.1]"),
                     type = "fixed",
                     ci.lvl = 0.95) %>%
  as_tibble() %>%
  rename(dams = x)

# plot it
logis_p = pred_df %>%
  ggplot(aes(x = dams)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              color = NA,
              fill = 'gray80') +
  geom_line(aes(y = predicted),
            color = "black",
            linetype = 2) +
  geom_point(data = dam_df,
             aes(y = suc_perc,
                 shape = "All Years",
                 size = "All Years")) +
  geom_point(data = dam_yr,
             aes(y = suc_perc,
                 shape = "Individ. Years",
                 size = "Individ. Years")) +
  scale_shape_manual(values = c("All Years" = 19,
                                "Individ. Years" = 1),
                     name = "Data Source") +
  scale_size_manual(values = c("All Years" = 4,
                               "Individ. Years" = 2),
                    name = "Data Source") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_pubr(base_family = "serif",
             base_size = 10) +
  theme(legend.position = c(0.2, 0.2)) +
  labs(x = "Number of Columbia River Dams",
       y = "Estimated Successful Overshoot Fallback Probability")

logis_p

ggsave(here("analysis/figures/Dams_Logistic_RE_glmer.jpeg"),
       logis_p,
       width = 5,
       height = 5)

# include random effects
pred_df_re <- ggpredict(m2,
                        terms = c("dams [0:5 by = 0.1]"),
                        type = "random") %>%
  as_tibble() %>%
  rename(dams = x)


logis_p2 = pred_df %>%
  ggplot(aes(x = dams)) +
  geom_ribbon(data = pred_df_re,
              aes(ymin = conf.low,
                  ymax = conf.high),
              color = NA,
              fill = 'gray90') +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              color = NA,
              fill = 'gray70') +
  geom_line(aes(y = predicted),
            color = "black",
            linetype = 2) +
  geom_point(data = dam_df,
             aes(y = suc_perc,
                 shape = "All Years",
                 size = "All Years")) +
  geom_point(data = dam_yr,
             aes(y = suc_perc,
                 shape = "Individ. Years",
                 size = "Individ. Years")) +
  scale_shape_manual(values = c("All Years" = 19,
                                "Individ. Years" = 1),
                     name = "Data Source") +
  scale_size_manual(values = c("All Years" = 4,
                               "Individ. Years" = 2),
                    name = "Data Source") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_pubr(base_family = "serif",
             base_size = 10) +
  theme(legend.position = c(0.2, 0.2)) +
  labs(x = "Number of Columbia River Dams",
       y = "Estimated Successful Overshoot Fallback Probability")

logis_p2

#-----------------------------------------------------------------
# Bayesian approach
library(brms)
library(tidybayes)

# no random effects
# fit1 <- brm(success | trials(success + fail) ~ dams,
#             data = dam_yr,
#             family = binomial)

fit1 <- brm(outcome_fct ~ dams,
            data = dam_individ,
            family = bernoulli)


# random intercept and slope
# fit2 <- brm(success | trials(success + fail) ~ dams + (1 + dams | year),
#             data = dam_yr,
#             family = binomial,
#             control = list(adapt_delta = 0.9))

fit2 <- brm(outcome_fct ~ dams + (1 + dams | year),
            data = dam_individ,
            family = bernoulli,
            control = list(adapt_delta = 0.9))


# random intercept
# fit3 <- brm(success | trials(success + fail) ~ dams + (1 | year),
#             data = dam_yr,
#             family = binomial,
#             control = list(adapt_delta = 0.9))

fit3 <- brm(outcome_fct ~ dams + (1 | year),
            data = dam_individ,
            family = bernoulli,
            control = list(adapt_delta = 0.9))

# random slope
# fit4 <- brm(success | trials(success + fail) ~ dams + (dams | year),
#             data = dam_yr,
#             family = binomial,
#             control = list(adapt_delta = 0.9))

fit4 <- brm(outcome_fct ~ dams + (dams | year),
            data = dam_individ,
            family = bernoulli,
            control = list(adapt_delta = 0.9))


fit1 <- add_criterion(fit1,
                      criterion = c("bayes_R2",
                                    "loo_R2"))
mean(fit1$criteria$bayes_R2)
mean(fit1$criteria$loo_R2)


fit2 <- add_criterion(fit2,
                       criterion = c("bayes_R2",
                                     "loo_R2"))
mean(fit2$criteria$bayes_R2)
mean(fit2$criteria$loo_R2)

# list(re_no = fit1,
#      re_si = fit2,
#      re_i = fit3,
#      re_s = fit4) %>%
#   map(.id = "model",
#       .f = loo)
# loo(fit1)
# loo(fit1)

# pull out all fixed effects
fix_est <- map(list(re_no = fit1,
                    re_si = fit2,
                    re_i = fit3,
                    re_s = fit4),
               fixef) %>%
  map_df(.id = "model",
         .f = as_tibble,
         rownames = "param")

fix_est %>%
  filter(param == "Intercept") %>%
  select(-Est.Error) %>%
  mutate(across(Estimate:Q97.5,
                boot::inv.logit))

fix_est %>%
  filter(param == "dams") %>%
  select(-Est.Error) %>%
  mutate(across(Estimate:Q97.5,
                exp))

posterior_predict(fit2,
        newdata = tibble(dams = seq(0, 5, by = 0.2)),
        re_formula = NA)

stanplot(fit2,
         type = "trace")

mod_coef <- fix_est %>%
  filter(model == "re_si")

logit_p <- fit2 %>%
  spread_draws(b_Intercept, b_dams) %>%
  mutate(dams = list(seq(0, 5, by = 0.1))) %>%
  unnest(dams) %>%
  mutate(pred_suc = inv.logit(b_Intercept + b_dams * dams)) %>%
  group_by(dams) %>%
  summarize(pred = mean(pred_suc),
            lci = quantile(pred_suc, 0.025),
            uci = quantile(pred_suc, 0.975),
            .groups = "drop") %>%
  ggplot(aes(x = dams)) +
  geom_ribbon(aes(ymin = lci,
                  ymax = uci),
              color = NA,
              fill = 'gray80') +
  geom_line(aes(y = pred),
            color = "black",
            linetype = 2) +
  geom_point(data = dam_df,
             aes(y = suc_perc,
                 shape = "All Years",
                 size = "All Years")) +
  geom_point(data = dam_yr,
             aes(y = suc_perc,
                 shape = "Individ. Years",
                 size = "Individ. Years")) +
  scale_shape_manual(values = c("All Years" = 19,
                                "Individ. Years" = 1),
                     name = "Data Source") +
  scale_size_manual(values = c("All Years" = 4,
                               "Individ. Years" = 2),
                    name = "Data Source") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_pubr(base_family = "serif",
             base_size = 10) +
  theme(legend.position = c(0.2, 0.2)) +
  labs(x = "Number of Columbia River Dams",
       y = "Estimated Successful Overshoot Fallback Probability")

logit_p

ggsave(here("analysis/figures/Dams_Logistic_RE_brms.jpeg"),
       logit_p,
       width = 5,
       height = 5)

# make predictions of probability of success, with confidence intervals
pred_df <- ggpredict(fit2,
                     terms = c("dams [0:5 by = 0.1]"),
                     type = "fixed",
                     ci.lvl = 0.95) %>%
  as_tibble() %>%
  rename(dams = x)

pred_df %>%
  ggplot(aes(x = dams)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              color = NA,
              fill = 'gray80') +
  geom_line(aes(y = predicted),
            color = "black",
            linetype = 2) +
  geom_point(data = dam_df,
             aes(y = suc_perc,
                 shape = "All Years",
                 size = "All Years")) +
  geom_point(data = dam_yr,
             aes(y = suc_perc,
                 shape = "Individ. Years",
                 size = "Individ. Years")) +
  scale_shape_manual(values = c("All Years" = 19,
                                "Individ. Years" = 1),
                     name = "Data Source") +
  scale_size_manual(values = c("All Years" = 4,
                               "Individ. Years" = 2),
                    name = "Data Source") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_pubr(base_family = "serif",
             base_size = 10) +
  theme(legend.position = c(0.2, 0.2)) +
  labs(x = "Number of Columbia River Dams",
       y = "Estimated Successful Overshoot Fallback Probability")

# Author: Kevin See
# Purpose: Logistic regression on overshoot success based on number of dams crossed
# Created: 9/15/20
# Last Modified: 9/16/20
# Notes: 

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(ResourceSelection) # for goodness of fit tests

#-----------------------------------------------------------------
# here's the data, provided by Andrew Murdoch. 
# All the numbers with dams > 0 are from fish tagged at PRD that return to a downstream trib
# for dam == 0, Andrew used known Yakima fish that crossed McNary.
# because slightly different data, we decided to only use PRD fish for the model
dam_df = tibble(dams = c(5:3, 1, 0),
                fail = c(38,
                         9,
                         14,
                         16,
                         14),
                success = c(11,
                            14,
                            31,
                            113,
                            262)) %>%
  mutate(n_tot = success + fail,
         suc_perc = success / n_tot) %>%
  arrange(dams) %>%
  filter(dams > 0)

# fit a binomial logistic model
binom_mod = glm(suc_perc ~ dams,
                data = dam_df,
                weights = dam_df$n_tot,
                family = "binomial")

null_mod = glm(suc_perc ~ 1,
               data = dam_df,
               weights = dam_df$n_tot,
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
pred_df = predict(binom_mod,
                  newdata = tibble(dams = 0:5),
                  type = "link",
                  se = T) %>%
  extract(1:2) %>%
  map_df(.id = 'pred',
         .f = identity) %>%
  pivot_longer(cols = -pred,
               names_to = "dams",
               values_to = "log_odds") %>%
  mutate(dams = as.integer(dams) - 1) %>%
  pivot_wider(names_from = "pred",
              values_from = "log_odds") %>%
  mutate(lowCI = fit + se.fit * qnorm(0.025),
         uppCI = fit + se.fit * qnorm(0.975)) %>%
  mutate_at(vars(pred_suc = fit, lowCI, uppCI),
            list(boot::inv.logit)) %>%
  select(dams,
         pred_suc, 
         ends_with('CI'))

pred_df %>%
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
  geom_smooth(data = pred_df,
              method = loess,
              aes(y = pred_suc),
              color = "black",
              linetype = 2) +
  geom_point(size = 4) +
  theme_pubr(base_family = "Times New Roman",
             base_size = 12) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Number of Columbia River Dams",
       y = "Overshoot Return Rate")

ggsave("outgoing/Dams_Logistic.jpeg",
       logis_p,
       width = 5,
       height = 5)

# Author: Kevin See
# Purpose: compare timing of overshoot tags to all tags at Prosser Dam
# Created: 1/6/2022
# Last Modified: 1/10/2022
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(magrittr)
library(here)
library(ggpubr)
library(readxl)
library(lubridate)
library(STADEM)
# library(ResourceSelection) # for goodness of fit tests

#-----------------------------------------------------------------
# here's the data, provided by Andrew Murdoch.
pro_tags <- read_excel(here("analysis/data/raw_data/",
                            "PIT tags at Prosser.xlsx"),
                       1) %>%
  mutate(year = year(Date),
         month = month(Date))

# can also get dam counts at Prosser
pro_cnts = tibble(spawn_year = 2011:2018) %>%
  mutate(win_cnts = map(spawn_year,
                        .f = quietly(function(yr) {
                          queryWindowCnts(dam = "PRO",
                                          spp_code = "fsw",
                                          start_date = paste0(yr-1, "0701"),
                                          end_date = paste0(yr, "0630"))
                        })),
         win_cnts = map(win_cnts,
                        "result")) %>%
  unnest(win_cnts) %>%
  mutate(month = month(Date))

pro_julian_cnts <- pro_cnts %>%
  mutate(Date = ymd(paste(if_else(month >= 7,
                                   2019,
                                   2020),
                           month(Date),
                           day(Date)))) %>%
  group_by(Date,
           month) %>%
  summarize(across(Wild_Steelhead,
                   sum),
            .groups = "drop") %>%
  uncount(Wild_Steelhead) %>%
  mutate(source = "Dam Counts") %>%
  mutate(day = difftime(Date,
                        min(Date),
                        units = "days")) %>%
  mutate(across(day,
                as.numeric))

# convert tag counts to individual rows
ind_tags = pro_tags %>%
  select(-year) %>%
  pivot_longer(cols = c(Yakima, Priest),
               names_to = "source",
               values_to = "n_tags") %>%
  uncount(n_tags) %>%
  mutate(day = difftime(Date, min(pro_julian_cnts$Date),
                        units = "days"),
         day = as.numeric(day))

# compare CDFs of arrival timing
ind_tags %>%
  bind_rows(pro_julian_cnts) %>%
  mutate(source = fct_recode(source,
                             "Yakima Tags" = "Yakima",
                             "Priest Tags" = "Priest")) %>%
  # filter(source != "Priest Tags") %>%
  mutate(across(Date,
                as.Date)) %>%
  ggplot(aes(x = Date,
             fill = source,
             color = source)) +
  scale_color_brewer(palette = "Set1",
                     name = "Data Source") +
  scale_fill_brewer(palette = "Set1",
                    name = "Data Source") +
  stat_ecdf() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_date(date_breaks = "2 months",
               date_labels =  "%b") +
  labs(y = "Cumulative Probability")

# density plot
ind_tags %>%
  bind_rows(pro_julian_cnts) %>%
  mutate(source = fct_recode(source,
                             "Yakima Tags" = "Yakima",
                             "Priest Tags" = "Priest")) %>%
  # filter(source != "Priest Tags") %>%
  mutate(across(Date,
                as.Date)) %>%
  mutate(across(Date,
                as.Date)) %>%
  ggplot(aes(x = Date,
             fill = source,
             color = source)) +
  scale_color_brewer(palette = "Set1",
                     name = "Data Source") +
  scale_fill_brewer(palette = "Set1",
                    name = "Data Source") +
  geom_density(alpha = 0.2) +
  scale_x_date(date_breaks = "2 months",
               date_labels =  "%b") +
  theme_bw()

# pull out vectors of days by source
yak_tags <- ind_tags %>%
  filter(source == "Yakima")
prd_tags <- ind_tags %>%
  filter(source == "Priest")

# Kolmogorov-Smirnov test
# using day
ks.test(yak_tags$day,
        prd_tags$day)
# using month
ks.test(yak_tags$month,
        prd_tags$month)

# looking at dam counts and tags
ks.test(pro_julian_cnts$day,
        prd_tags$day)

# some subsets, trying to find a combination that looks the same
ks.test(pro_julian_cnts$day[pro_julian_cnts$Date <= ymd("20200430")],
        yak_tags$day[yak_tags$Date <= ymd("20200430")])

ks.test(pro_julian_cnts$day[pro_julian_cnts$Date >= ymd("20191101") &
                              pro_julian_cnts$Date <= ymd("20200430")],
        yak_tags$day[yak_tags$Date >= ymd("20191101") &
                       yak_tags$Date <= ymd("20200430")])

ks.test(pro_julian_cnts$month[pro_julian_cnts$month >= 10 |
                                pro_julian_cnts$month <= 4],
        yak_tags$month[yak_tags$month >= 10 |
                         yak_tags$month <= 4])

#-----------------------------------------------------------------
#-----------------------------------------------------------------
# recreate figure 3
# get temperature data
temp_df = read_excel(here('analysis/data/raw_data/temps for kevin.xlsx'),
                     range = "A1:D367") %>%
  mutate(Month = month(Date,
                       label = T),
         month = month(Date),
         Month = fct_reorder(Month,
                             .x = Date),
         month_num = as.numeric(Month)) %>%
  pivot_longer(cols = c(Kiona,
                        PRD),
               names_to = 'source',
               values_to = 'temp') %>%
  group_by(Month,
           month,
           month_num,
           source) %>%
  summarize(across(temp,
                   mean,
                   na.rm = T),
            .groups = "drop")

temp_df %>%
  ggplot(aes(x = month_num,
             y = temp,
             linetype = source)) +
  # geom_point() +
  geom_line() +
  theme_pubr(base_family = "Times New Roman",
             base_size = 12) +
  scale_x_continuous(name = "Date",
                     breaks = 1:nlevels(temp_df$Month),
                     labels = levels(temp_df$Month))

# divide the temperature by this coefficient to put on same plot
coeff <- 70

fig_3 <- pro_tags %>%
  pivot_longer(cols = c(Yakima, Priest),
               names_to = "source",
               values_to = "n_tags") %>%
  group_by(year,
           month,
           source) %>%
  summarize(across(n_tags,
                   sum),
            .groups = "drop") %>%
  group_by(source) %>%
  mutate(tag_prop = n_tags / sum(n_tags)) %>%
  ungroup() %>%
  left_join(temp_df %>%
              pivot_wider(names_from = source,
                          values_from = temp),
            by = c("month")) %>%
  mutate(Date = ymd(paste(year, month, "01"))) %>%
  mutate(source = as_factor(source),
         source = fct_relevel(source,
                              "Priest",
                              after = 0)) %>%
  ggplot(aes(x = Date,
             y = tag_prop,
             fill = source)) +
  scale_fill_grey(name = "Tag\nSource") +
  geom_col(position = "dodge",
           color = "black") +
  geom_line(aes(y = Kiona / coeff,
                linetype = "Yakima"),
            lwd = 1) +
  geom_line(aes(y = PRD / coeff,
                linetype = "Columbia"),
            lwd = 1) +
  scale_x_date(breaks = "1 month",
               date_labels = "%B") +
  scale_y_continuous(name = "Mean proportion at Prosser Dam",
                     # limits = c(0,0.45),
                     sec.axis = sec_axis(name = "Water Temperature (C)",
                                         trans = ~ . * coeff,
                                         breaks = seq(0, 25, length.out = 6))) +
  theme_pubr(base_family = "serif",
             base_size = 10) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2,
                             title.position = "left"),
         linetype = guide_legend(nrow = 2,
                                 title.position = "left")) +
  theme(legend.position = c(0.75, 0.8)) +
  # theme(legend.background = element_rect(fill = "white", color = "black")) +
  labs(linetype = "Temperature\nSource")


fig_3

ggsave(here("analysis/figures/Prosser_Timing.jpeg"),
       fig_3,
       width = 5,
       height = 5)


# what proportion of each group of tags arrived by Jan 1?
pro_tags %>%
  pivot_longer(cols = c(Yakima, Priest),
               names_to = "source",
               values_to = "n_tags") %>%
  group_by(year,
           month,
           source) %>%
  summarize(across(n_tags,
                   sum),
            .groups = "drop") %>%
  group_by(source) %>%
  mutate(tag_prop = n_tags / sum(n_tags)) %>%
  mutate(cum_prop = cumsum(tag_prop)) %>%
  ungroup() %>%
  filter(month %in% c(12,1))


#-------------------------------------------------------------
# time-series analysis
library(forecast)

yak_ts = ts(pro_tags$Yakima)
prd_ts = ts(pro_tags$Priest)

plot(yak_ts)

yak_arima = auto.arima(yak_ts)
prd_arima = auto.arima(prd_ts)

yak_arima = Arima(yak_ts,
                  order = c(1,1,1))
prd_arima = Arima(prd_ts,
                  order = c(1,1,1))


yak_arima
prd_arima

anova(yak_arima,
      prd_arima)

summary(yak_arima)
coef(yak_arima)

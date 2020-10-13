# Chart 1:

## Load seroprevalence data 
sero <- read.csv("source-data/sero-surveys.csv")

# Notes:
# We collected the sero-surveys from medical journals, government websites, or, in a few instances, press reports (which were triangulated). All were conducted by medical researchers and/or government agencies. We only included serosurveys taken in general populations (i.e. not only health-care workers). 
# Because sero-surveys are conducted over an interval of time, one has to make a decision about how to match it to case tallies. We decided to use the midpoint of the period in which samples were collected.
# Because sero-surveys measure the prevalence of antibodies rather than active infections, they must be brought backwards in time to correspond to percentages of those with a past or current likely active infection. We selected an interval of two weeks. 
# As we discuss in the briefing, sero-surveys may not reflect the number of people who have been infected with covid-19. In some world regions, the circulation of viruses similar to covid-19 may trigger antibody-tests, even if a person has not been infected. Tests may suffer from poor sensitivity, or specificity (where available, we used results which corrected for this). Test subjects may not be fully random, and methods to make them closer to random or as-if random (through e.g. re-weighting by demographics), which we used were available, may not always be successful. Finally, the extent to which covid-19 produces antibodies was at the time of writing not fully understood. It may be that some or many, especially those who have suffered only mild infections, do not produce antibodies, or do so to an extent that is picked up by antibody tests.   

# NB: As we make clear in the chart footnote and briefing text, this model extrapolates from data in places where serosurveys are available to places they are not. This means they must be taken with more than a grain of salt. Some countries may not fit the pattern. There are fewer available serosurveys from countries with lower average incomes, and fewer serosurveys from recent dates. This means that all these estimates are very uncertain, especially for times and places where serosurveys are less available. It may be that patterns found among these serosurveys and cases and deaths do not extrapolate well. It is however crucial to get better estimates of the scale of the pandemic than those provided by official case counts alone, which as of September 15 suggested 30m cases worldwide. This presents one effort to do so with the data limitations we have. To improve our sense of the scale of the pandemic, more serosurveys conducted according to scientific standards, especially in developing countries, are crucial.

## Variables: 
# sero_ir : seroprevalence of covid-19 antibodies as indicated by serosurvey.
# sample.size : Number of samples in serosurvey
# case_ir : cumulative reported cases of covid-19 / population
# death_r : cumulative reported covid-19 deaths / population
# log_gdp_ppp : Log(GDP per capita, adjusted for purchasing power).

# Notes: 
# GDP per capita is at the national level, case_ir and death_r are for the geographical units in which the serosurvey was conducted (national, province/state, region, or metropolitan level). In a few very rare instances, a sero-survey was conducted across a large geographical region within a country which could plausibly similar to the country average, e.g. Niger state in Nigeria, where no case or death data was available. In these instances, we used data from the national level instead. Smaller surveys (e.g. cities) which could not be matched to time-series data on cases and deaths were not used.

## Model seroprevalence through case and death data:
lm_model <- glm(sero_ir ~ (log(case_ir) + log(death_r + 1e-7)) * log_gdp_ppp,
                data = sero, family = binomial, weights = sample.size)

# Notes:
# The model assumes that seroprevalence is some multiple of case rates and death rates and their interaction, scaled by a proxy for testing capacity (log_gdp_per_capita). Furthermore, it assumes that as diagnosed deaths and cases rise as a share of population, they come to better reflect underlying cases and deaths (hence taking the natural logs of both). We weight observations by the log of sample size, so that larger serosurveys, which have less uncertainty, are more influential.
# If available for all countries, we would have preferred to use tests conducted and their positive rate in these models. However, while this data is available for many countries (Our World in Data provides an excellent resource), they are often incomparable and more importantly, not available for all countries. Furthermore, even where they are available, they are usually not available throughout the pandemic, especially in the crucial early months. 
# Diagnostics of this model revealed Afghanistan's national survey to be especially high-leverage. We therefore contacted the WHO to confirm the estimate provided of seroprevalence in that survey. The WHO Press Office confirmed it was the correct preliminary estimate.
# Many other modelling strategies are possible, some of which we tried (e.g. beta regression and a random forest), their estimates were similar and landed within the interval of our predictions. A linear model was thought to be the most straightforward, interpretable, and conservative choice.

## Predict infection rate from offical cases, deaths, scaled by testing capacity (gdp per capita)
prediction_df <- read.csv("source-data/cases_deaths_gdp.csv")
prediction_df$date <- as.Date(prediction_df$date)
prediction_df$pred_ir <- predict(lm_model, newdata = prediction_df, type = 'response')
prediction_df$pred_ir_low <- qbinom(0.05, 1000, prediction_df$pred_ir)/1000
prediction_df$pred_ir_high <- qbinom(0.95, 1000, prediction_df$pred_ir)/1000

# We know cumulative infection rates have lower bound of zero and are strictly increasing, and use this to improve our predictions. We first use a 10-day average to avoid day-to-day variation driving changes.
library(runner)
make_average <- function(x, n){
  temp <- x
  for(i in 1:length(x)){
    x[i] <- mean(temp[max(c(1, i-n)):min(c(length(x), i+n))], na.rm = T)
  }
  x
}
prediction_df$pred_ir[prediction_df$pred_ir < 0] <- 0 
prediction_df$pred_ir_low[prediction_df$pred_ir_low < 0] <- 0
prediction_df$pred_ir_high[prediction_df$pred_ir_high < 0] <- 0
prediction_df <- prediction_df[order(prediction_df$date), ]
prediction_df$pred_ir <- ave(prediction_df$pred_ir, prediction_df$iso3c, FUN = function(x) max_run(make_average(x, n = 5)))
prediction_df$pred_ir_low <- ave(prediction_df$pred_ir_low, prediction_df$iso3c, FUN = function(x) max_run(make_average(x, n = 5)))
prediction_df$pred_ir_high <- ave(prediction_df$pred_ir_high, prediction_df$iso3c, FUN = function(x) max_run(make_average(x, n = 5)))

# From predicted sero-prevalence rates, we extract implied case counts
prediction_df$pred_cases <- prediction_df$pred_ir*prediction_df$population
prediction_df$pred_cases_low <- prediction_df$pred_ir_low*prediction_df$population
prediction_df$pred_cases_high <- prediction_df$pred_ir_high*prediction_df$population

# The following five lines creates world totals:
prediction_df$world_cases <- ave(prediction_df$cases, prediction_df$date, FUN = sum)
prediction_df$pred_world_cases <- ave(prediction_df$pred_cases, prediction_df$date, FUN = sum)
prediction_df$pred_world_cases_low <- ave(prediction_df$pred_cases_low, prediction_df$date, FUN = sum)
prediction_df$pred_world_cases_high <- ave(prediction_df$pred_cases_high, prediction_df$date, FUN = sum)

# This re-produces our inset plot:
library(ggplot2)
ggplot(prediction_df[!duplicated(prediction_df$date), ], 
       aes(x=as.Date(date), ymin = 0))+geom_ribbon(aes(ymax=pred_world_cases, fill = "Probably Infected, World"))+
  theme_minimal()+geom_ribbon(aes(ymin = 0, ymax = world_cases, fill = "Reported Cases, World"))+
  geom_line(aes(y=pred_world_cases_low), col = "white")+
  geom_line(aes(y=pred_world_cases_high), col = "black")+
  scale_y_continuous(labels = scales::comma)+
  theme(legend.title = element_blank(), legend.position = "bottom")+xlab("")+ylab("")

# To generate the large plot, we first-differences and by continent and a select few large countries:

# These lines define our groups (colors in the large plot)
prediction_df$continents_plus <- prediction_df$continent
prediction_df$continents_plus[prediction_df$country == "United States"] <- "United States"
prediction_df$continents_plus[prediction_df$country == "China"] <- "China"
prediction_df$continents_plus[prediction_df$country == "India"] <- "India"
prediction_df$continents_plus[prediction_df$country == "Brazil"] <- "Brazil"

# This function takes first differences by country and sums countries together by day and group: 
big_chart_data <- function(prediction_df, grouping = "continent_plus"){
  
  prediction_df$region <- prediction_df[, grouping]
  
  regions <- prediction_df
  
  regions$region_cases <- ave(regions$cases, paste0(regions$date, "_", regions$region), FUN = sum)
  regions$pred_region_cases <- ave(regions$pred_cases, paste0(regions$date, "_", regions$region), FUN = sum)
  regions$pred_region_cases_low <- ave(regions$pred_cases_low, paste0(regions$date, "_", regions$region), FUN = sum)
  regions$pred_region_cases_high <- ave(regions$pred_cases_high, paste0(regions$date, "_", regions$region), FUN = sum)
  
  new_cases_fun <- function(x) {
    x <- x - c(0, x)[1:length(x)]
    x <- make_average(x, n = 10)
    x }
  regions <- regions[!duplicated(paste0(regions$date, "_", regions$region)), ]
  
  regions$pred_region_cases <- ave(regions$pred_region_cases, regions$region, FUN = function(x) make_average(x, n = 10))
  regions$pred_region_cases_low <- ave(regions$pred_region_cases_low, regions$region, FUN = function(x) make_average(x, n = 10))
  regions$pred_region_cases_high <- ave(regions$pred_region_cases_high, regions$region, FUN = function(x) make_average(x, n = 10))
  
  regions$new_cases <- ave(regions$region_cases, regions$region, FUN = new_cases_fun)
  regions$new_pred_cases <- ave(regions$pred_region_cases, regions$region, FUN = new_cases_fun)
  regions$new_pred_cases_high <- ave(regions$pred_region_cases_high, regions$region, FUN = new_cases_fun)
  regions$new_pred_cases_low <- ave(regions$pred_region_cases_low, regions$region, FUN = new_cases_fun)
  return(regions)}
continents_plus <- big_chart_data(prediction_df, grouping = "continents_plus") # This runs the above function

# This reproduces our first large plot:
ggplot(continents_plus, 
       aes(x=date, ymin = 0))+geom_area(aes(y=-new_pred_cases, fill = continents_plus), col = "white")+
  theme_minimal()+geom_ribbon(aes(ymin = 0, ymax = -new_cases, fill = "Reported Cases"))+
  #  geom_line(aes(y=new_pred_cases_high), col = "white")+
  # geom_line(aes(y=new_pred_cases_low), col = "black")+
  scale_y_continuous(labels = scales::comma)+ylab("")+
  theme(legend.title = element_blank(), legend.position = "bottom")+xlab("")

# Notes: As we write and show through our confidence interval above, newer dates are less certain, as sero-surveys are only released some time after data is collected. The labels and order in which groups appear and other visual details were tweaked by our graphical designers. 
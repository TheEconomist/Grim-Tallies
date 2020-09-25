# Chart 2, Panel A:

# This chart is a simple thermometer plot of the following data:
thermometer_chart <- read.csv("source-data/serosurveys_selected.csv")

# Note: This also features a few surveys which we could not consistently match to case and death counts, and thus did not feature in the main model data or inform this model (see Chart 1). The multiplier noted is the average ratio of sero-survey prevalence percent to cumulative cases / population for each country. 
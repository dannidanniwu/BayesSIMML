#got data after running interim.R
library("collapse")
library("dplyr")
library("data.table")
library("tidyverse")

##sim16
for ( i in 1:7){
load(paste0("./scenarios",i,".rda"))
results.aggregated2[[i]] <- res
}

save(results.aggregated2, file = "C:/Users/Danni/OneDrive - NYU Langone Health/BayesSIMML/results/all_24_scenarios.rda")

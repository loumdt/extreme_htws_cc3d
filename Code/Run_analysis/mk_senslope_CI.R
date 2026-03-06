#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(MannKendallTrends)
library(jsonlite)
library(readr)
# Load correct path from argument
file_path <- args[1]
# Load data to compute trend on frequency
data_frequency <- read_csv(paste(file_path,"/frequency.csv", sep=""), col_types = cols(Year = col_datetime(format = "%Y")))
# Load data to compute other trends
df_htws_step2 <- read_csv(paste(file_path,"/df_htws_step3.csv", sep="" ),
    col_types = cols(...1 = col_skip(), Year = col_datetime(format = "%Y")))
# Select subset of columns
cols_to_exclude <- names(df_htws_step2) %in% c('Start Date', 'End Date', 'EM-DAT heatwaves', 'EM-DAT Total Deaths', 'HWMId_pop_mean', 'HWMId_mean', 'Temp_sum', 'Temp_pop', 'Total Exposed_population')
sub_df <- df_htws_step2[!cols_to_exclude]
# Create output data.frame
result_mk <- data.frame(slope = numeric(), UCL = numeric(), LCL = numeric(), P = numeric(), ss = numeric())
# Compute Mann-Kendall trend test on frequency
res <- MK.tempAggr(data = data_frequency, resolution = 0.9)
result_mk['Frequency',] <- res
# Compute other Mann-Kendall trend tests
for (i in 2:length((colnames(sub_df))))
{
res <- MK.tempAggr(data=sub_df[,c(1,i)], resolution = 0.9)
index_name <- colnames(sub_df[,i]) # Get index name and use it as output row name
result_mk[index_name,] <- res
}
# Write result to json file
write_json(result_mk,paste(file_path,"/mk_result.json",sep=""))

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
df_htws_step3 <- read_csv(paste(file_path,"/df_htws_step3.csv", sep="" ),
    col_types = cols(...1 = col_skip(), Year = col_datetime(format = "%Y")))
# Select subset of columns
cols_to_exclude <- names(df_htws_step3) %in% c('Start Date', 'End Date', 'EM-DAT heatwaves', 'EM-DAT Total Deaths', 'HWMId_pop_mean', 'HWMId_mean', 'Temp_sum', 'Temp_pop', 'Total Exposed_population')
sub_df <- df_htws_step3[!cols_to_exclude]

sub_df$HWMId_pop <- sub_df$HWMId_pop/1e6 # Divide by 1e6 to avoid memory issues
sub_df$HWMId_sum <- sub_df$HWMId_sum/1e5 # Divide by 1e3 to avoid memory issues
sub_df$'Spatial extent' <- sub_df$'Spatial extent'/1e3 # Divide by 1e3 to avoid memory issues
# Create output data.frame
result_mk <- data.frame(slope = numeric(), UCL = numeric(), LCL = numeric(), P = numeric(), ss = numeric())
# Compute Mann-Kendall trend test on frequency
res <- MK.tempAggr(data = data_frequency, resolution = 0.9)
result_mk['Frequency',] <- res
# Compute other Mann-Kendall trend tests
last_col <- length((colnames(sub_df))) - 1
for (i in 2:last_col)
{
print(colnames(sub_df[,i]))
res <- MK.tempAggr(data=sub_df[,c(1,i)], resolution = 0.9)
index_name <- colnames(sub_df[,i]) # Get index name and use it as output row name
result_mk[index_name,] <- res
print(index_name)
}

i <- i+1
print(colnames(sub_df[,i]))
print('One')
res <- MK.tempAggr(data=sub_df[,c(1,i)], resolution = 0.01)
print('two')
index_name <- colnames(sub_df[,i]) # Get index name and use it as output row name
print('three')
result_mk[index_name,] <- res
print(index_name)
print("\n")


result_mk['HWMId_pop','slope'] <- result_mk['HWMId_pop','slope']*1e6 # Set back to orginal unit
result_mk['HWMId_pop','UCL'] <- result_mk['HWMId_pop','UCL']*1e6 # Set back to orginal unit
result_mk['HWMId_pop','LCL'] <- result_mk['HWMId_pop','LCL']*1e6 # Set back to orginal unit

result_mk['HWMId_sum','slope'] <- result_mk['HWMId_sum','slope']*1e5 # Set back to orginal unit
result_mk['HWMId_sum','UCL'] <- result_mk['HWMId_sum','UCL']*1e5 # Set back to orginal unit
result_mk['HWMId_sum','LCL'] <- result_mk['HWMId_sum','LCL']*1e5 # Set back to orginal unit

result_mk['Spatial extent','slope'] <- result_mk['Spatial extent','slope']*1e3 # Set back to orginal unit
result_mk['Spatial extent','UCL'] <- result_mk['Spatial extent','UCL']*1e3 # Set back to orginal unit
result_mk['Spatial extent','LCL'] <- result_mk['Spatial extent','LCL']*1e3 # Set back to orginal unit
# Write result to json file
write_json(result_mk,paste(file_path,"/mk_result.json",sep=""))

# RUN SCRUB for decontaminating microbiome data

# Set directory
setwd('C:\\Users\\mspg\\Documents\\Projects\\ADRC')

# Load libraries
libs <- c('tidyverse','SCRuB','here')

# Install packages if not already installed
install.packages(setdiff(libs, rownames(installed.packages())))
# Load libraries and dependencies
sapply(libs,require,character.only=TRUE)

# Run tutorial
set.seed(1)
library(SCRuB)

folders <- c('194890','194939')

for (i in 1:2){

# Load microbiome data
microbiomePath <- here('ADRC','study_15448_031424-023713','BIOM',folders[i],'free.txt')
data <- read.delim(microbiomePath,check.names = FALSE,row.names=1) %>% as.matrix() %>% t()

# Load metadata
metadataPath <- here('ADRC','study_15448_031424-023713','mapping_files',str_c(folders[i],'_mapping_file.txt'))
metadata <- read.delim(metadataPath,check.names = FALSE,row.names=1)

# Sort samples
data <- data[order(rownames(data)), ]
metadata <- metadata[order(rownames(metadata)), ]

# Add binary variable indicating negative controls and variable for control and feces
metadata <- metadata %>% 
  mutate(is_control = str_detect(metadata$run_prefix,'BLANK'),
         sample_type = if_else(str_detect(metadata$run_prefix,'BLANK'),'control blank','feces')) %>%
  select(is_control,sample_type,sample_well)

# Run scrub
scr_out <- SCRuB(data,metadata)

# Investigate decontaminated samples
decontaminated_samples <- scr_out$decontaminated_samples

decontaminated_samples <- as_tibble(decontaminated_samples,rownames = "ID")

# # Visualise the level of contamination
# (1 - scr_out$p) %>%
#   boxplot(ylab='Level of contamination')
# 
# # Check the level of well leakage
# boxplot( 1 - scr_out$inner_iterations$`control blank`$alpha[, ncol(scr_out$inner_iterations$`control blank`$alpha)],
#          ylab='Level of well leakage in process controls')
# 
# # Check the estimate abundance of the contamination community
# scr_out$inner_iterations$`control blank`$gamma %>% 
#   plot(ylab='Relative abundance', xlab='Different OTUs in contamination source')

# Save results
write.table(decontaminated_samples,file=here('ADRC',str_c('reads_',folders[i],'_decontaminated.txt')),sep = '\t',row.names = FALSE)
}
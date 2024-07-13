library(Seurat)
library(readr)

#Set working directory
setwd("/parallel_scratch/mp01950/raw_data/")

# Define the organ directories
organs <- c("Liver", "Colon", "Lung","Thyroid","Breast")

# Initialize an empty list to store the loaded data
data_list <- list()

#Loading .rds files into environment
# Loop over each organ directory
for (organ in organs) {
  # Define the path to the saved_RDS directory
  saved_rds_dir <- file.path(organ, "saved_RDS")
  
  # Get the list of all .rds files with the pattern *_md.rds
  rds_files <- list.files(path = saved_rds_dir, pattern = "_mmd\\.rds$", full.names = TRUE)
  
  # Loop over each file and load the data
  for (rds_file in rds_files) {
    # Load the .rds file
    data <- readRDS(rds_file)
    
    # Extract the name of the file without the extension for naming the list entry
    file_name <- tools::file_path_sans_ext(basename(rds_file))
    
    # Store the data in the list
    data_list[[file_name]] <- data
  }
}

for (i in seq_along(data_list)) {
  
  #Normalizaion
  data_list[[i]] <- NormalizeData(data_list[[i]])
  
  #Find HVGs
  data_list[[i]] <- FindVariableFeatures(data_list[[i]],selection.method = "vst", nfeatures = 2000)
  
}

#Find integration anchors
anchors <- FindIntegrationAnchors(object.list = data_list)

#Save anchors
saveRDS(anchors, file = "anchors.rds")

#Integration for batch effect
integrated_data <- IntegrateData(anchorset = anchors)

#Save anchors
saveRDS(integrated_data, file = "Integrated_dataset.rds")

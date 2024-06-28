#This is for GSE154826, metadata was provided for batches instead of single cells
#All single cells in a batch will have the same metadata

# Load necessary libraries
library(Seurat)
library(dplyr)

#Define the organ/tissue to build the data path
study_subject <- "lung/"
#modify the path according to your data storage location
raw_data_path <- "/parallel_scratch/mp01950/raw_data/"

subject_data_path <- paste(raw_data_path, study_subject, sep = "")
setwd(subject_data_path)
directory_path <- file.path(subject_data_path, "saved_RDS")
metadata_path <- file.path(subject_data_path, "saved_RDS", "metadata", "GSE154826_meta.csv.gz")

# List all .rds files in the directory
all_files <- list.files(directory_path, pattern = "\\.rds$", full.names = TRUE)

# Filter files matching the pattern "GSE154826_batch_XXX_filtered"
pattern <- "GSE154826_batch_\\d+_filtered\\.rds"
filtered_files <- grep(pattern, all_files, value = TRUE)

# Initialize a list to store the Seurat objects
seurat_objects <- list()

# Read each filtered .rds file and store it in the list
for (file in filtered_files) {
  seurat_obj <- readRDS(file)
  
  # Extract the base name of the file without the directory path
  base_name <- basename(file)
  
  # Store the Seurat object in the list with the base name as the key
  seurat_objects[[base_name]] <- seurat_obj
}

# Load the metadata CSV file
metadata <- read.csv(metadata_path)

# Verify the column names
print(colnames(metadata))

# Iterate through each Seurat object and add metadata
for (name in names(seurat_objects)) {
  # Extract the batch ID from the object name
  batch_id <- as.numeric(gsub("GSE154826_batch_|_filtered.*", "", name))
  
  # Find the matching row in the metadata
  meta_row <- metadata %>% filter(amp_batch_ID == batch_id)
  
  # Check if the metadata row is found
  if (nrow(meta_row) == 1) {
    # Remove the amp_batch_ID column
    meta_row <- meta_row %>% select(-amp_batch_ID)
    
    # Convert the metadata row to a list
    meta_list <- as.list(meta_row)
    
    # Get the Seurat object
    seurat_obj <- seurat_objects[[name]]
    
    # Add metadata to the Seurat object
    for (meta_name in names(meta_list)) {
      seurat_obj[[meta_name]] <- meta_list[[meta_name]]
    }
    
    # Store the updated Seurat object back in the list
    seurat_objects[[name]] <- seurat_obj
  } else {
    message(paste("No matching metadata found for", name))
  }
}

# List all .rds files in the directory
all_files <- list.files(directory_path, pattern = "\\.rds$", full.names = TRUE)

# Filter files matching the pattern "GSE154826_batch_XXX_filtered"
pattern <- "GSE154826_batch_\\d+_filtered\\.rds"
filtered_files <- grep(pattern, all_files, value = TRUE)

# Merge all Seurat objects into one
batch_ids <- gsub("GSE154826_batch_(\\d+)_filtered\\.rds", "\\1", names(seurat_objects))
merged_seurat <- merge(x = seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = batch_ids)

# Save the merged Seurat object as .rds file
merged_file_path <- file.path(subject_data_path, "saved_RDS", "GSE154826_filtered_md.rds")
saveRDS(merged_seurat, file = merged_file_path)

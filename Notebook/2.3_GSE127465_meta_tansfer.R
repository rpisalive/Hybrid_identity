#This script is to rename all _filtered.rds files for GSE127465 and transfer metadata.
#Reason for individual workflow: the matching of the metadata relies on both "Library" and "Barcode"
#Original count data in corresponding directories are required to run this.
#Seuat objects with only 1 cell will (and must) be removed due to the design logic of Seurat Objects.

# Load necessary libraries
library(Seurat)

#Define the organ/tissue to build the data path
study_subject <- "lung/"
#modify the path according to your data storage location
raw_data_path <- "/parallel_scratch/mp01950/raw_data/"

subject_data_path <- paste(raw_data_path, study_subject, sep = "")
setwd(subject_data_path)
# Get the list of GSM directories
gsm_dirs <- list.dirs(".", recursive = FALSE, full.names = TRUE)
gsm_dirs <- gsm_dirs[grepl("GSM3635\\d{3}$", gsm_dirs)]

# Define the path to the saved_RDS directory
saved_rds_dir <- file.path(".", "saved_RDS")

# Iterate over each GSM directory
for (gsm_dir in gsm_dirs) {
  # Extract the base GSM name (e.g., "GSM3635296")
  base_gsm_name <- basename(gsm_dir)
  
  # Find the .tsv.gz file in the current GSM directory
  tsv_file <- list.files(gsm_dir, pattern = "_raw_counts.tsv.gz$", full.names = TRUE)
  
  if (length(tsv_file) == 1) {
    # Extract the new name from the .tsv.gz file name
    new_name <- gsub("_raw_counts.tsv.gz$", "", basename(tsv_file))
    
    # Define the current .rds file path in the saved_RDS directory
    current_rds_path <- file.path(saved_rds_dir, paste0(base_gsm_name, ".rds"))
    
    # Define the new .rds file path
    new_rds_path <- file.path(saved_rds_dir, paste0(new_name, ".rds"))
    
    # Rename the .rds file
    if (file.exists(current_rds_path)) {
      file.rename(current_rds_path, new_rds_path)
      cat("Renamed:", current_rds_path, "to", new_rds_path, "\n")
    } else {
      cat("File does not exist:", current_rds_path, "\n")
    }
  } else {
    cat("No unique .tsv.gz file found in directory:", gsm_dir, "\n")
  }
}

# Get a list of all .rds files in the saved_RDS directory
rds_files <- list.files(saved_rds_dir, pattern = "^GSM.*_filtered\\.rds$", full.names = TRUE)

# Create an empty list to store the loaded Seurat objects
seurat_objects <- list()

# Load each .rds file and store it in the seurat_objects list
for (rds_file in rds_files) {
  seurat_object <- readRDS(rds_file)
  object_name <- gsub("\\.rds$", "", basename(rds_file))
  seurat_objects[[object_name]] <- seurat_object
}


# Set the path to the metadata file
metadata_file_path <- file.path("saved_RDS", "metadata", "GSE127465_meta.tsv.gz")

# Read the metadata file into a dataframe
metadata_df <- read.delim(metadata_file_path, sep = "\t", header = TRUE)

# Iterate over each row in the metadata dataframe
for (i in 1:nrow(metadata_df)) {
  # Extract the Library and Barcode values
  library_value <- metadata_df$Library[i]
  barcode_value <- metadata_df$Barcode[i]
  
  # Find the matching Seurat object based on the Library value (YYY)
  matching_objects <- grep(paste0("_", library_value, "$"), names(seurat_objects), value = TRUE)
  
  # If there is a matching Seurat object
  if (length(matching_objects) == 1) {
    matching_object_name <- matching_objects[1]
    matching_seurat_object <- seurat_objects[[matching_object_name]]
    
    # Check if the Barcode exists in the rownames of the matching Seurat object's @meta.data
    if (barcode_value %in% rownames(matching_seurat_object@meta.data)) {
      # Transfer all other data from the metadata dataframe row to the Seurat object's @meta.data
      row_index <- which(rownames(matching_seurat_object@meta.data) == barcode_value)
      matching_seurat_object@meta.data[row_index, names(metadata_df)[!names(metadata_df) %in% c("Library", "Barcode")]] <- metadata_df[i, !names(metadata_df) %in% c("Library", "Barcode")]
      
      # Update the Seurat object in the list
      seurat_objects[[matching_object_name]] <- matching_seurat_object
    }
  }
}

rm(metadata_df)
rm(seurat_object)

# Filter the list to remove Seurat objects with only one cell in the RNA assay
seurat_objects <- lapply(seurat_objects, function(seurat_obj) {
  if (ncol(seurat_obj) > 1) {
    return(seurat_obj)
  } else {
    return(NULL)
  }
})

# Remove NULL elements from the list
seurat_objects <- seurat_objects[!sapply(seurat_objects, is.null)]

# Loop through each Seurat object and save it with the modified name
for (i in seq_along(seurat_objects)) {
  # Extract the original name of the Seurat object
  original_name <- names(seurat_objects)[i]
  
  # Create the new file name by appending "_md" before ".rds"
  new_name <- paste0("saved_RDS/",original_name, "_md.rds")
  
  # Save the Seurat object with the new name
  saveRDS(seurat_objects[[i]], file = new_name)
}

# Load necessary libraries
library(Seurat)
library(readr)
library(data.table)

#Define the organ/tissue to build the data path
study_subject <- "liver/"
#modify the path according to your data storage location
raw_data_path <- "/parallel_scratch/mp01950/raw_data/"

subject_data_path <- paste(raw_data_path, study_subject, sep = "")
setwd(subject_data_path)

# Define the paths to the directories
saved_rds_dir <- "saved_RDS"
metadata_dir <- "saved_RDS/metadata"

# Get a list of all .rds files in the saved_RDS directory
rds_files <- list.files(path = saved_rds_dir, pattern = "_filtered\\.rds$", full.names = TRUE)

# Initialize lists to store the Seurat objects and their corresponding metadata
seurat_objects <- list()
metadata_objects <- list()

# Function to load metadata file based on its extension
load_metadata <- function(file) {
  if (grepl("\\.txt\\.gz$", file)) {
    return(fread(file))
  } else if (grepl("\\.csv\\.gz$", file)) {
    return(fread(file))
  }
}

# Loop through each .rds file and load the Seurat object
for (file in rds_files) {
  # Extract the name of the file without extension
  object_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the Seurat object from the .rds file
  seurat_object <- readRDS(file)
  
  # Assign the Seurat object to a variable in the environment
  assign(object_name, seurat_object)
  
  # Store the Seurat object in the list
  seurat_objects[[object_name]] <- seurat_object
  
  # Extract the prefix to match metadata files
  prefix <- sub("_filtered$", "", object_name)
  
  # Find metadata files corresponding to the prefix
  metadata_files <- list.files(path = metadata_dir, pattern = paste0("^", prefix, ".*\\.gz$"), full.names = TRUE)
  
  # Load each metadata file
  for (meta_file in metadata_files) {
    meta_data <- load_metadata(meta_file)
    metadata_objects[[basename(meta_file)]] <- meta_data
    metadata_objects[[basename(meta_file)]] <- as.data.frame(metadata_objects[[basename(meta_file)]])
  }
}
# Define object names to remove
objects_to_remove <- c("meta_data", "seurat_object", "seurat_objects")

# Remove objects from the global environment
rm(list = objects_to_remove)

#Replace all "-" in the first column in metadata objects with ".".
#If the the cell_ID in your metadata is consistent with the rownames of your meta.data of your seurat objects, do not run this part.
# Loop through each data frame in the list
for (i in seq_along(metadata_objects)) {
  # Replace "-" with "." in the first column of each data frame
  metadata_objects[[i]][, 1] <- gsub("-", ".", metadata_objects[[i]][, 1])
}

#Metadata transfer
# Iterate over each metadata data frame
for (meta_name in names(metadata_objects)) {
  # Extract the prefix from the metadata frame name
  meta_prefix <- gsub("_meta.*", "", meta_name)
  
  # Find the corresponding Seurat object in the environment
  seurat_name <- paste0(meta_prefix, "_filtered")
  
  if (seurat_name %in% ls(pattern = seurat_name)) {
    # Retrieve Seurat object
    seurat_obj <- get(seurat_name)
    
    # Retrieve cell names from Seurat object
    seurat_cell_names <- rownames(seurat_obj@meta.data)
    
    # Retrieve metadata data frame
    metadata_df <- metadata_objects[[meta_name]]
    
    # Check for matching cell names
    matching_cells <- metadata_df[,1] %in% seurat_cell_names
    
    # Iterate over matching cell names
    for (cell in metadata_df[,1][matching_cells]) {
      # Find corresponding row index in metadata data frame
      row_index <- which(metadata_df[,1] == cell)
      
      # Add metadata to Seurat object's @meta.data if cell name exists
      if (cell %in% seurat_cell_names) {
        seurat_obj@meta.data[cell, names(metadata_df)[-1]] <- metadata_df[row_index, names(metadata_df)[-1]]
      }
    }
    
    # Update Seurat object in the environment
    assign(seurat_name, seurat_obj, envir = .GlobalEnv)
    
    cat("Added metadata to", seurat_name, "\n")
  } else {
    cat("Seurat object", seurat_name, "not found in the environment. Skipping...\n")
  }
}

# Define object names to remove
objects_to_remove <- c("metadata_df", "metadata_objects", "seurat_obj")

# Remove objects from the global environment
rm(list = objects_to_remove)

# List all objects in the environment
all_objects <- ls()
  
# Filter objects that are of class Seurat
seurat_objects <- all_objects[sapply(all_objects, function(x) inherits(get(x), "Seurat"))]
  
# Iterate over each Seurat object
for (seurat_name in seurat_objects) {
  # Retrieve the Seurat object
  seurat_obj <- get(seurat_name)
  
  # Construct the filename
  file_name <- paste0(seurat_name, "_md.rds")
    
  # Construct the full file path
  file_path <- file.path("saved_RDS/", file_name)
    
  # Save the Seurat object as an .rds file
  saveRDS(seurat_obj, file = file_path)
    
  # Print a message indicating the object has been saved
  cat("Saved", seurat_name, "as", file_path, "\n")
}

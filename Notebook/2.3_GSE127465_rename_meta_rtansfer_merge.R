#This script is to rename all .rds files for GSE127465, transfer metadata and merging
#Reason: The filtering pipeline would yield seurat objects with only 1 cell,
#which cannot be merged with the standard merging pipeline.
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
rds_files <- list.files(saved_rds_dir, pattern = "\\.rds$", full.names = TRUE)

# Create an empty list to store the loaded Seurat objects
seurat_objects <- list()

# Load each .rds file and store it in the seurat_objects list
for (rds_file in rds_files) {
  seurat_object <- readRDS(rds_file)
  object_name <- gsub("\\.rds$", "", basename(rds_file))
  seurat_objects[[object_name]] <- seurat_object
}

# Print the names of all loaded Seurat objects
cat("Loaded Seurat objects:\n")
print(names(seurat_objects))

# Remove any Seurat objects with names starting with "GSE"
seurat_objects <- seurat_objects[!grepl("^GSE", names(seurat_objects))]

# Remove Seurat objects with names containing the "filtered" string
seurat_objects <- seurat_objects[!grepl("filtered", names(seurat_objects))]

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
rm(matching_seurat_object)
rm(metadata_df)
rm(seurat_object)

#Merging
# Extract the first Seurat object
GSM3635278_human_p1t1 <- seurat_objects[[1]]
# Remove the first Seurat object from the list
seurat_objects <- seurat_objects[-1]
batch_ids <- c("GSM3635278_human_p1t1", names(seurat_objects))
merged_obj <- merge(GSM3635278_human_p1t1, y = seurat_objects, add.cell.ids = batch_ids)

rm(GSM3635278_human_p1t1)
rm(seurat_objects)

# Function to perform quality control on a Seurat object
filter_seurat_object <- function(seurat_obj, seurat_version) {
  # Filter cells with less than 500 genes
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 500)
  
  # Calculate the percentage of mitochondrial reads
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Filter cells with more than 20% mitochondrial reads
  seurat_obj <- subset(seurat_obj, subset = percent.mt <= 20)
  
  # Define features to check
  features_to_check <- c("CD3D", "CD3E", "CD3G")
  
  # Check if the Seurat object has the features CD3D, CD3E, and CD3G
  if (!all(features_to_check %in% rownames(seurat_obj))) {
    cat("Seurat object does not have either CD3D, CD3E, or CD3G feature.\n")
    return(NULL)
  }
  
  # Fetch data for CD3D, CD3E, and CD3G
  gene_data <- FetchData(seurat_obj, vars = features_to_check)
  
  # Identify cells that meet the filtering criteria
  cells_to_keep <- rownames(gene_data)[gene_data$CD3D > 0 & gene_data$CD3E > 0 & gene_data$CD3G > 0]
  
  # Check if cells_to_keep has length 0
  if (length(cells_to_keep) == 0) {
    cat("No cells pass the CD3D, CD3E, and CD3G filter.\n")
    return(NULL)
  }
  
  # Subset the Seurat object to keep only the selected cells
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
  
  return(seurat_obj)
}

# Get Seurat version
seurat_version <- as.numeric(substr(packageVersion("Seurat"), 1, 1))

# Perform quality control
filtered_seurat_obj <- filter_seurat_object(merged_obj, seurat_version)

# Save the Seurat object as an .rds file
saveRDS(filtered_seurat_obj, file = "saved_RDS/GSE127465_filtered_md.rds")


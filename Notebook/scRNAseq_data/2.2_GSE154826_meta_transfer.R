#This is for GSE154826, metadata were provided in three .csv files:
#1) At sample level - named as "GSE154826_meta.csv.gz"
#2) At single cell level, sample_ID is the linkage to 1 - named as "GSE154826_cell_metadata.csv.gz"
#3) At cluster level, cluster no. is the linkage to 2 - named as "GSE154826_annots_list.csv.gz"

# Load necessary libraries
library(Seurat)
library(dplyr)
library(stringr)
library(tibble)

#Define the organ/tissue to build the data path
study_subject <- "Lung/"
#modify the path according to your data storage location
raw_data_path <- "/parallel_scratch/mp01950/raw_data/"

subject_data_path <- paste(raw_data_path, study_subject, sep = "")
setwd(subject_data_path)
directory_path <- file.path(subject_data_path, "saved_RDS")
metadata_path <- file.path(subject_data_path, "saved_RDS", "metadata", "GSE154826_meta.csv.gz")

# List all .rds files in the directory
all_files <- list.files(directory_path, pattern = "\\.rds$", full.names = TRUE)

# Filter files matching the pattern "GSE154826_amp_batch_XXX_filtered"
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

rm(seurat_obj)

# Load the CSV files
meta_data <- read.csv(gzfile("saved_RDS/metadata/GSE154826_meta.csv.gz"))
cell_metadata <- read.csv(gzfile("saved_RDS/metadata/GSE154826_cell_metadata.csv.gz"))
annots_list <- read.csv(gzfile("saved_RDS/metadata/GSE154826_annots_list.csv.gz"))

# Convert amp_batch_ID to character in meta_data
meta_data$amp_batch_ID <- as.character(meta_data$amp_batch_ID)

# Iterate over each Seurat object
for (i in seq_along(seurat_objects)) {
  # Extract the Seurat object and its name
  seurat_obj <- seurat_objects[[i]]
  seurat_name <- names(seurat_objects)[i]
  
  #Create a column for cell_ID
  seurat_obj@meta.data$cell_ID <- rownames(seurat_obj@meta.data)
  
  # Extract the XXX number from the Seurat object name
  amp_batch_ID <- str_extract(seurat_name, "(?<=GSE154826_batch_)\\d+(?=_filtered)")
  
  #Add amp_batch_ID & sample_ID columns into meta.data & list of column names to initialize with NA
  columns_to_initialize <- c(
    "sample_ID", "cluster", "lineage", "sub_lineage", "norm_group", "lig_rec_group",
    "patient_ID", "old_lib_name", "HTO", "tissue", "disease", "Species",
    "library_chemistry", "prime", "vdj_kit", "prep"
  )
  
  # Initialize amp_batch_ID
  seurat_obj@meta.data$amp_batch_ID <- amp_batch_ID
  
  # Loop through each column name and initialize it with NA
  for (col in columns_to_initialize) {
    seurat_obj@meta.data[[col]] <- NA
  }
  
  # Find rows in meta_data with the same amp_batch_ID
  meta_rows <- meta_data %>% filter(amp_batch_ID == !!amp_batch_ID)
  
  for (j in seq_len(nrow(meta_rows))) {
    sample_ID <- meta_rows$sample_ID[j]
    
    # Check if the sample_ID exists in cell_metadata
    if (sample_ID %in% cell_metadata$sample_ID) {
      # Get the cell_IDs from cell_metadata for the sample_ID
      cell_IDs <- cell_metadata %>% filter(sample_ID == !!sample_ID)
      
      # Extract the part of cell_ID after "_" and check if it matches rownames in Seurat object
      cell_IDs_to_check <- str_extract(cell_IDs$cell_ID, "(?<=_).*")
      matching_cell_IDs <- cell_IDs_to_check[cell_IDs_to_check %in% seurat_obj@meta.data$cell_ID]
      
      if (length(matching_cell_IDs) > 0) {
        # Create a temporary data frame for the matching cells with amp_batch_ID and sample_ID
        temp_meta <- data.frame(
          cell_ID = matching_cell_IDs,
          amp_batch_ID = amp_batch_ID,
          sample_ID = sample_ID,
          stringsAsFactors = FALSE
        )
        
        # Merge data frames based on cell_ID to add sample_ID to the matching cells in Seurat object's metadata
        seurat_obj@meta.data <- left_join(seurat_obj@meta.data, temp_meta, by = "cell_ID", suffix = c("", ".b")) %>%
          mutate(sample_ID = coalesce(sample_ID, sample_ID.b)) %>%
          select(-sample_ID.b, -matches("\\.b$"))
        
        # Add the meta_data information to Seurat object metadata
        meta_info <- meta_rows %>% filter(sample_ID == !!sample_ID)
        
        # Add the meta_data information to Seurat object metadata
        seurat_obj@meta.data <- left_join(seurat_obj@meta.data, meta_info, by = "sample_ID", suffix = c("", ".b")) %>%
          mutate(patient_ID = coalesce(patient_ID, patient_ID.b),
                 old_lib_name = coalesce(old_lib_name, old_lib_name.b),
                 HTO = coalesce(HTO, HTO.b),
                 tissue = coalesce(tissue, tissue.b),
                 disease = coalesce(disease, disease.b),
                 Species = coalesce(Species, Species.b),
                 library_chemistry = coalesce(library_chemistry, library_chemistry.b),
                 prime = coalesce(prime, prime.b),
                 vdj_kit = coalesce(vdj_kit, vdj_kit.b),
                 prep = coalesce(prep, prep.b)) %>%
          select(-matches("\\.b$"))
        
        # Add cluster column from cell_metadata to temp_meta for matching cell_IDs and sample_IDs
        temp_meta <- temp_meta %>%
          inner_join(cell_metadata %>% 
                       mutate(cell_ID = str_extract(cell_ID, "(?<=_).*")), 
                     by = c("cell_ID" = "cell_ID", "sample_ID"))
        
        #Rename cluster_ID as cluster
        temp_meta <- temp_meta %>%
          rename(cluster = cluster_ID)
        
        # Find cluster IDs from temp_meta
        cluster_IDs <- temp_meta$cluster
        
        # Filter annots_list for cluster IDs matching in cell_metadata
        annots_rows <- annots_list %>% filter(cluster %in% cluster_IDs)
        
        if (nrow(annots_rows) > 0) {
          # Add annotations from annots_list to temp_meta based on matching cluster
          temp_meta <- temp_meta %>%
            left_join(annots_rows, by = "cluster")
        }
      }
    }
    # Add temp_meta to seurat_obj@meta.data for matching cell IDs
    temp_meta <- temp_meta[, !names(temp_meta) %in% c("amp_batch_ID", "sample_ID")]
    seurat_obj@meta.data <- left_join(seurat_obj@meta.data, temp_meta, by = "cell_ID", suffix = c("", ".b")) %>%
      mutate(cluster = coalesce(cluster, cluster.b),
             lineage = coalesce(lineage, lineage.b),
             sub_lineage = coalesce(sub_lineage, sub_lineage.b),
             norm_group = coalesce(norm_group, norm_group.b),
             lig_rec_group = coalesce(lig_rec_group, lig_rec_group.b)) %>%
      select(-cluster.b, -matches("\\.b$"))
  }
  
  # Rename row names, remove cell_ID and update the Seurat object in the list
  rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$cell_ID
  seurat_obj@meta.data <- seurat_obj@meta.data %>% select(-cell_ID)
  seurat_objects[[i]] <- seurat_obj
}

# Remove the ".rds" extension
seu_obj_names <- names(seurat_objects)
seu_obj_names <- gsub("\\.rds$", "", seu_obj_names)
names(seurat_objects) <- seu_obj_names

# Loop through each Seurat object and save it with the modified name
for (i in seq_along(seurat_objects)) {
  # Extract the original name of the Seurat object
  original_name <- names(seurat_objects)[i]
  
  # Create the new file name by appending "_md" before ".rds"
  new_name <- paste0("saved_RDS/",original_name, "_md.rds")
  
  # Save the Seurat object with the new name
  saveRDS(seurat_objects[[i]], file = new_name)
}

#Merging for GSE154826
batch_ids <- sub(".*_(ID_\\d+)_filtered$", "\\1", names(seurat_objects))
GSE154826 <- merge(seurat_objects[[1]], y = seurat_objects[2:77], add.cell.ids = batch_ids)
saveRDS(GSE154826, file = "saved_RDS/GSE154826_filtered_md.rds")

library(Seurat)

#Define the organ/tissue to build the data path
study_subject <- "liver/"
#modify the path according to your data storage location
raw_data_path <- "/parallel_scratch/mp01950/raw_data/"
subject_data_path <- paste(raw_data_path, study_subject, sep = "")
setwd(subject_data_path)

#ONLY RUN THIS PART IF YOU START WITH RAW DATA (e.g. .mtx,.txt)
#GSE98638, SMART seq data, loading the unprocessed matrix
#GSE156625, scRNA seq data, part of GSE156337, loading the unprocessed matrix

# Get the list of all subfolders
folders <- list.dirs(subject_data_path, recursive = FALSE)

# Iterate over each folder
for (folder in folders) {
  # Get the name of the folder
  folder_name <- basename(folder)
  
  # List all gzip files in the folder
  gz_files <- list.files(folder, pattern = "\\.gz$", full.names = TRUE)
  
  # Initialize a flag to check if the folder contains 'mtx' files
  contains_mtx <- any(grepl("mtx\\.gz$", gz_files))
  contains_txt <- any(grepl("txt\\.gz$", gz_files))
  
  if (contains_mtx) {
    # Load the data using Read10X
    dgCMatrix_object <- Read10X(data.dir = folder)
    
    # Create a Seurat object from the dgCMatrix object
    seurat_object <- CreateSeuratObject(counts = dgCMatrix_object)
    
    # Assign the Seurat object to a variable with the folder name
    assign(folder_name, seurat_object)
    
    # Remove the dgCMatrix object from the environment
    rm(dgCMatrix_object)

    # Use the folder_name variable to create the file name & .rds storing path
    file_name <- paste0(folder_name, ".rds")
    file_path <- paste0(subject_data_path,"saved_RDS/", file_name)

    #Save the seurat object as .rds
    saveRDS(seurat_object, file = file_path)
  } else if (contains_txt) {
    # Load the 'txt.gz' file using read.delim
    txt_file <- gz_files[grepl("txt\\.gz$", gz_files)]
    data <- read.delim(gzfile(txt_file), sep = "\t", header = TRUE)
    
    # Find the first column with non-integer data
    non_integer_col <- which(sapply(data, function(x) any(!is.integer(x))))[1]
    if (!is.na(non_integer_col)) {
      #Replace all 'NA' with the values from the -1 column on the same row        
      data[which(is.na(data[[non_integer_col]])),non_integer_col] <- data[which(is.na(data[[non_integer_col]])), non_integer_col-1]
      # Make the first non-integer column as the rownames
      rownames(data) <- data[, non_integer_col]
      
      # Remove the first non-integer column & all columns before the first non-integer column
      data <- data[, (non_integer_col + 1):ncol(data)]
    }
    
    # Convert the dataframe to a matrix
    data_matrix <- as.matrix(data)
    
    # Create a Seurat object from the matrix
    seurat_object <- CreateSeuratObject(counts = data_matrix)
    
    # Assign the Seurat object to a variable with the folder name
    assign(folder_name, seurat_object)

    # Use the folder_name variable to create the file name
    file_name <- paste0(folder_name, ".rds")
    file_path <- paste0(subject_data_path,"saved_RDS/", file_name)

    #Save the seurat object as .rds
    saveRDS(seurat_object, file = file_path)
                                    
    # Remove the dataframe object from the environment
    rm(data)
  }
}
#Remove the seurat object & data_matrix
rm(data_matrix,seurat_object)

# List all objects in the environment
all_objects <- ls()

# Check which of these objects are of class 'Seurat'
seurat_objects <- sapply(all_objects, function(x) {
  obj <- get(x)
  class(obj) == "Seurat"
})

# Extract the names of Seurat objects
seurat_object_names <- names(seurat_objects)[seurat_objects]

# Function to perform quality control on a Seurat object
filter_seurat_object <- function(seurat_obj) {
  # Filter cells with less than 500 genes
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 500)
  
  # Calculate the percentage of mitochondrial reads
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Filter cells with more than 20% mitochondrial reads
  seurat_obj <- subset(seurat_obj, subset = percent.mt <= 20)
  
  # Further filter cells with 0 counts for CD3D, CD3E, or CD3G genes, for Seurat V4
  #seurat_obj <- subset(seurat_obj, subset = CD3D > 0 & CD3E > 0 & CD3G > 0)

  # Further filter cells with 0 counts for CD3D, CD3E, or CD3G genes, for Seurat V5
  # Fetch data for CD3D, CD3E, and CD3G
  gene_data <- FetchData(seurat_obj, vars = c("CD3D", "CD3E", "CD3G"))
  # Identify cells that meet the filtering criteria
  cells_to_keep <- rownames(gene_data)[gene_data$CD3D > 0 & gene_data$CD3E > 0 & gene_data$CD3G > 0]
  # Subset the Seurat object to keep only the selected cells
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    
  return(seurat_obj)
}

# Iterate over each Seurat object and perform quality control
for (obj_name in seurat_object_names) {
  # Get the Seurat object
  seurat_obj <- get(obj_name)
  
  # Perform quality control
  filtered_seurat_obj <- filter_seurat_object(seurat_obj)

  # Use the folder_name variable to create the file name
  file_name <- paste0(obj_name, "_filtered.rds")
  file_path <- paste0(subject_data_path,"saved_RDS/", file_name)
  
  # Save the filtered Seurat object with the specified prefix
  saveRDS(filtered_seurat_obj, file = file_path)
}

# Print completion message
cat("Filtering and saving of Seurat objects completed.\n")













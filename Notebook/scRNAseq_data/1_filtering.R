library(Seurat)
library(Matrix)

#Define the organ/tissue to build the data path
study_subject <- "Lung/"
#modify the path according to your data storage location
raw_data_path <- "/parallel_scratch/mp01950/raw_data/"

subject_data_path <- paste(raw_data_path, study_subject, sep = "")
setwd(subject_data_path)

# Function to read the data files and rename duplicate feature names
read_data_files <- function(barcodes_path, features_path, matrix_path) {
  # Read in the matrix file
  counts <- readMM(matrix_path)
  
  # Read in the barcodes and features files
  barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
  features <- read.table(features_path, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  
  # Check if there are "Custom" strings in the third column
  if (any(features$V3 == "Custom")) {
    
    # Identify duplicated items in the second column
    custom_features <- features$V3 == "Custom"
    
    if (any(custom_features)) {
      message("Modifying feature names for 'Custom' entries in features$V2.")
      
      for (i in which(custom_features)) {
        # Modify the second column to "Custom" if the third column is "Custom"
        features$V2[i] <- "Custom"
        # If the third column is not "Custom", leave the feature name as is
      }
    }
  } else {
    # Check for duplicate feature names in features$V2
    duplicate_features <- duplicated(features$V2) | duplicated(features$V2, fromLast = TRUE)
    
    if (any(duplicate_features)) {
      message("Renaming duplicate feature names in features$V2.")
      
      for (i in which(duplicate_features)) {
        if (startsWith(features$V1[i], "ENSG")) {
          features$V2[i] <- features$V1[i]
        } else {
          features$V2[i] <- paste0(features$V2[i], "_", features$V3[i])
        }
      }
    }
  }
  
  # Add row and column names to the matrix
  colnames(counts) <- barcodes$V1
  rownames(counts) <- features$V2

  # Remove rows where row names are "Custom"
  counts <- counts[rownames(counts) != "Custom", ]
  
  return(counts)
}

# ONLY RUN THIS PART IF YOU START WITH RAW DATA (e.g. .mtx,.txt,.tsv)

# Get the list of all subfolders
folders <- list.dirs(subject_data_path, recursive = FALSE)

# Iterate over each folder
for (folder in folders) {
  # Get the name of the folder
  folder_name <- basename(folder)
  
  # List all gzip files in the folder
  gz_files <- list.files(folder, pattern = "\\.gz$", full.names = TRUE)

  # Initialize a flag to check if the folder contains the corresponding suffix
  contains_mtx <- any(grepl("matrix\\.mtx\\.gz$", gz_files))
  contains_barcodes <- any(grepl("barcodes\\.tsv\\.gz$", gz_files))
  contains_features <- any(grepl("features|genes\\.tsv\\.gz$", gz_files))
  contains_txt_or_tsv <- !contains_mtx && any(grepl("\\.(txt|tsv)\\.gz$", gz_files))

  if (contains_mtx && contains_barcodes && contains_features) {
    # Identify the specific file paths with prefixes
    barcodes_file <- gz_files[grepl("barcodes.tsv.gz$", gz_files)]
    features_file <- gz_files[grepl("(features|genes)\\.tsv\\.gz$", gz_files)]
    matrix_file <- gz_files[grepl("matrix.mtx.gz$", gz_files)]
    
    # Load the data using the custom read_data_files function
    dgCMatrix_object <- read_data_files(barcodes_file, features_file, matrix_file)
    
    # Create a Seurat object from the dgCMatrix object
    seurat_object <- CreateSeuratObject(counts = dgCMatrix_object)
    
    # Assign the Seurat object to a variable with the folder name
    assign(folder_name, seurat_object)
    
    # Remove the dgCMatrix object from the environment
    rm(dgCMatrix_object)

    # Use the folder_name variable to create the file name & .rds storing path
    file_name <- paste0(folder_name, ".rds")
    file_path <- paste0(subject_data_path,"saved_RDS/", file_name)

    # Save the Seurat object as .rds
    saveRDS(seurat_object, file = file_path)
  } else if (contains_txt_or_tsv) {
    # Load the 'txt.gz' or 'tsv.gz' file using read.delim
    txt_file <- gz_files[grepl("\\.(txt|tsv)\\.gz$", gz_files)]
    data <- read.delim(gzfile(txt_file), sep = "\t", header = TRUE)
    
    # Check if the file is a .tsv.gz and transpose the data if necessary
    if (grepl("\\.tsv\\.gz$", txt_file)) {
      data <- t(data)
      colnames(data) <- data[1, ]
      data <- data[-1, ]
    }
    
    # Find the first column with non-integer data
    non_integer_col <- which(sapply(data, function(x) any(!is.integer(x))))[1]
    if (!is.na(non_integer_col) && grepl("\\.txt\\.gz$", txt_file)) {
      # Replace all 'NA' with the values from the -1 column on the same row
      data[which(is.na(data[[non_integer_col]])), non_integer_col] <- data[which(is.na(data[[non_integer_col]])), non_integer_col-1]
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

    # Save the Seurat object as .rds
    saveRDS(seurat_object, file = file_path)
                                    
    # Remove the dataframe object from the environment
    rm(data)
  }
}
# Remove the seurat object & data_matrix
rm(data_matrix, seurat_object)

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

# Iterate over each Seurat object and perform quality control
for (obj_name in seurat_object_names) {
  # Get the Seurat object
  seurat_obj <- get(obj_name)

  # Perform quality control
  filtered_seurat_obj <- filter_seurat_object(seurat_obj, seurat_version)

  # Check if the filtered Seurat object is not NULL
  if (!is.null(filtered_seurat_obj)) {
    # Use the folder_name variable to create the file name
    file_name <- paste0(obj_name, "_filtered.rds")
    file_path <- paste0(subject_data_path,"saved_RDS/", file_name)

    # Save the filtered Seurat object with the specified prefix
    saveRDS(filtered_seurat_obj, file = file_path)
  }
}

# Print completion message
cat("Filtering and saving of Seurat objects completed.\n")



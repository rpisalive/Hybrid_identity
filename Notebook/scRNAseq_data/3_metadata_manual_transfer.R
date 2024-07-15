library(Seurat)

#Set working directory
setwd("/parallel_scratch/mp01950/raw_data/")

# Define the organ directories
organs <- c("Liver", "Colon", "Lung")

# Initialize an empty list to store the loaded data
data_list <- list()

#Loading .rds files into environment
# Loop over each organ directory
for (organ in organs) {
  # Define the path to the saved_RDS directory
  saved_rds_dir <- file.path(organ, "saved_RDS")
  
  # Get the list of all .rds files with the pattern *_md.rds
  rds_files <- list.files(path = saved_rds_dir, pattern = "_md\\.rds$", full.names = TRUE)
  
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

#Rename GSE154826 to avoid updating the tissue column
colnames(data_list[[9]]@meta.data)[colnames(data_list[[9]]@meta.data) == "tissue"] <- "origin"

# Define the tissue mapping
gse_tissue_mapping <- list(
  GSE184362 = "Thyroid",
  GSE191288 = "Thyroid",
  GSE193581 = "Thyroid",
  GSE132465 = "Colon",
  GSE132257 = "Colon",
  GSE144735 = "Colon",
  GSE115469 = "Liver",
  GSE156625 = "Liver",
  GSE98638 = "Liver",
  GSE127465 = "Lung",
  GSE131907 = "Lung",
  GSE154826 = "Lung",
  GSE114727 = "Breast"
)

# Iterate over each Seurat object in data_list to add metadata
for (seurat_name in names(data_list)) {
  # Extract the GSE identifier from the object name
  gse_id <- sub("_filtered_md", "", seurat_name)
  
  # Get the corresponding tissue type from the mapping
  tissue_type <- gse_tissue_mapping[[gse_id]]
  
  # Add 'tissue' and 'GEO.accession' columns to the metadata
  data_list[[seurat_name]]@meta.data$tissue <- tissue_type
  data_list[[seurat_name]]@meta.data$GEO.accession <- gse_id
}

# Define the cancer type mapping
cancer_type_mapping <- list(
  GSE115469 = "Hepatocellular_carcinoma",
  GSE98638 = "Hepatocellular_carcinoma",
  GSE132257 = "Colorectal_cancer",
  GSE132465 = "Colorectal_cancer",
  GSE144735 = "Colorectal_cancer",
  GSE127465 = "NCLC",
  GSE131907 = "Lung_adenocarcinoma"
)

# Iterate over each Seurat object in data_list to add metadata
for (seurat_name in names(data_list)) {
  # Extract the GSE identifier from the object name
  gse_id <- sub("_filtered_md", "", seurat_name)
  
  # Add 'cancer.type' column
  if (gse_id == "GSE154826") {
    data_list[[seurat_name]]@meta.data$cancer.type <- ifelse(
      data_list[[seurat_name]]@meta.data$disease == "LUAD", 
      "Lung_adenocarcinoma",
      ifelse(data_list[[seurat_name]]@meta.data$disease == "LUSC",
             "Lung_squamous_cell_carcinoma",
             NA)
    )
  } else {
    cancer_type <- cancer_type_mapping[[gse_id]]
    data_list[[seurat_name]]@meta.data$cancer.type <- cancer_type
  }
  
  # Add 'sample.origin' column
  if (gse_id == "GSE115469") {
    data_list[[seurat_name]]@meta.data$sample.origin <- "Healthy"
  } else if (gse_id == "GSE156625") {
    data_list[[seurat_name]]@meta.data$sample.origin <- ifelse(
      data_list[[seurat_name]]@meta.data$NormalvsTumor == "N", 
      "Adjacent", 
      ifelse(data_list[[seurat_name]]@meta.data$NormalvsTumor == "T", 
             "Tumor", 
             NA)
    )
  } else if (gse_id == "GSE98638") {
    data_list[[seurat_name]]@meta.data$sample.origin <- ifelse(
      substr(data_list[[seurat_name]]@meta.data$sampleType, 1, 1) == "P", 
      "Blood", 
      ifelse(substr(data_list[[seurat_name]]@meta.data$sampleType, 1, 1) == "T", 
             "Tumor", 
             ifelse(substr(data_list[[seurat_name]]@meta.data$sampleType, 1, 1) == "J", 
                    "Tumor_adj_joint", 
                    ifelse(substr(data_list[[seurat_name]]@meta.data$sampleType, 1, 1) == "N", 
                           "Adjacent",
                           NA)
             )
      )
    )
  } else if (gse_id == "GSE132257" || gse_id == "GSE132465") {
    data_list[[seurat_name]]@meta.data$sample.origin <- ifelse(
      data_list[[seurat_name]]@meta.data$Class == "Tumor", 
      "Tumor", 
      ifelse(data_list[[seurat_name]]@meta.data$Class == "Normal", 
             "Adjacent_mucosa", 
             NA)
    )
  } else if (gse_id == "GSE144735") {
    data_list[[seurat_name]]@meta.data$sample.origin <- ifelse(
      data_list[[seurat_name]]@meta.data$Class == "Tumor", 
      "Tumor_core", 
      ifelse(data_list[[seurat_name]]@meta.data$Class == "Normal", 
             "Adjacent_mucosa", 
             ifelse(data_list[[seurat_name]]@meta.data$Class == "Border", 
                    "Tumor_border", 
                    NA)
      )
    )
  } else if (gse_id == "GSE127465") {
    data_list[[seurat_name]]@meta.data$sample.origin <- ifelse(
      data_list[[seurat_name]]@meta.data$Tissue == "tumor", 
      "Tumor", 
      ifelse(data_list[[seurat_name]]@meta.data$Tissue == "blood", 
             "Blood", 
             NA)
    )
  } else if (gse_id == "GSE131907") {
    data_list[[seurat_name]]@meta.data$sample.origin <- ifelse(
      data_list[[seurat_name]]@meta.data$Sample_Origin == "nLung", 
      "Healthy", 
      ifelse(data_list[[seurat_name]]@meta.data$Sample_Origin == "tLung", 
             "Early_tumor_lung", 
             ifelse(data_list[[seurat_name]]@meta.data$Sample_Origin == "tL/B", 
                    "Advanced_tumor_lung", 
                    ifelse(data_list[[seurat_name]]@meta.data$Sample_Origin == "mLN", 
                           "Lymphnode_metastases", 
                           ifelse(data_list[[seurat_name]]@meta.data$Sample_Origin == "nLN", 
                                  "Healthy", 
                                  ifelse(data_list[[seurat_name]]@meta.data$Sample_Origin == "PE", 
                                         "Pleural_effusion", 
                                         ifelse(data_list[[seurat_name]]@meta.data$Sample_Origin == "mBrain", 
                                                "Brain_metastases", 
                                                NA)
                                  )
                           )
                    )
             )
      )
    )
  } else if (gse_id == "GSE154826") {
    data_list[[seurat_name]]@meta.data$sample.origin <- ifelse(
      data_list[[seurat_name]]@meta.data$origin == "Tumor", 
      "Tumor", 
      ifelse(data_list[[seurat_name]]@meta.data$origin == "PBMC", 
             "Blood", 
             ifelse(data_list[[seurat_name]]@meta.data$origin == "Normal",
                    "Adjacent",
                    NA)
      )
    )
  }
}

# Iterate over each Seurat object in data_list
for (seurat_name in names(data_list)) {
  # Extract the unique tissue information
  tissue <- unique(data_list[[seurat_name]]@meta.data$tissue)
  
  # Make sure there's only one unique value
  if (length(tissue) != 1) {
    warning("More than one unique tissue type found. Using the first one.")
    tissue <- tissue[1]
  }
  
  # Construct the path to the destination directory
  destination_dir <- file.path(tissue, "saved_RDS")
  
  # Create the destination directory if it doesn't exist
  if (!file.exists(destination_dir)) {
    dir.create(destination_dir, recursive = TRUE)
  }
  
  # Construct the file name for saving
  file_name <- paste0(sub("_md$", "_mmd", seurat_name), ".rds")
  
  # Construct the full path to save the .rds file
  save_path <- file.path(destination_dir, file_name)
  
  # Save the Seurat object as .rds
  saveRDS(data_list[[seurat_name]], save_path)
  
  # Print confirmation message
  cat("Saved", seurat_name, "as", file_name, "in", destination_dir, "\n")
}

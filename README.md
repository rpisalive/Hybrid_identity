This is the repository for the hybrid identity project. Scripts of different analysis are stored here.
liver.r is the latest pipeline which you could load .mtx , .txt & .tsv data (must be gziped) and carry out filtering automatically.

Prerequisites:
1. A corresponding file architecture:
   - Each raw dataset must be stored in one directory named by the corresponding accession number. (e.g. sample.mtx.gz in GSE12345 directory) 
   - Directories named by accession codes must be stored in a directory named by its corresponding organ source. (e.g. GSE12345 in liver directory)
   -  Directories named by organ names must be stored in the raw_data_path. 
2. The data loading part of the filtering pipeline does not include .RData loading, please load your .RData object and extract the raw count matrix from the corresponding seurat object in your .RData.
3. In case of 10X format, you do not need to rename the three files as barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz, as long as the name contains "barcodes", "features" and "matrix", the pipeline will handle it automatically. The Read10X is not used in the pipeline.

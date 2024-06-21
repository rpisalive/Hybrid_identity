This is the repository page for the hybrid identity project. Scripts of different analysis are stored here.
liver.r is the latest pipeline which you could load .mtx & .txt data (must be gziped) and carry out filtering automatically.

Prerequisites:
1. A corresponding file architecture:
   - Each raw dataset must be stored in one directory named by the corresponding accession number. (e.g. sample.mtx.gz in GSE12345 directory) 
   - Directories named by accession codes must be stored in a directory named by its corresponding organ source. (e.g. GSE12345 in liver directory)
   -  Directories named by organ names must be stored in the raw_data_path. 
3. For CD3D, CD3E & CD3G filtering, please check your Seurat verion and run the corresponding part.

############### SYNOPSIS ###################
# Download and prepare Mousebrain gene epxression data

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(here)

# ======================================================================= #
# ========================== Download  ================================ #
# ======================================================================= #

file_download <- here("tmp-data/expression/mousebrain-l5_all.loom")

if (!file.exists(file_download)) {
  # Download UMI data  
  download_url <- "https://storage.googleapis.com/linnarsson-lab-loom/l5_all.loom"
  download.file(download_url, destfile=file_download)
}

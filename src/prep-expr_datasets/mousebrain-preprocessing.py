

# Mousebrain

## Description

# Extract Zeisel et al Biorxiv 2018 ('mousebrain') loom file to data and metadat csv files

## Reference 

# http://linnarssonlab.org/loompy/cookbook/index.html

## Setup

### Packages

import loompy
import pandas as pd
import numpy as np

### constants

path_loomfile = "/scratch/data-for_fast_access/pub-others/zeisel-biorxiv-2018/l5_all.loom"

path_fileout_rawData = "/projects/jonatan/pub-perslab/18-mousebrain/RObjects/mousebrain_L5_subset_raw.csv.gz"
path_fileout_metaData = "/projects/jonatan/pub-perslab/18-mousebrain/RObjects/mousebrain_L5_subset_metadata.csv.gz"

cellclasses = ["DEINH3",
"MEGLU11",
"MEINH2",
"MEGLU10",
"DEGLU5",
"MEGLU1",
"MEGLU3",
"TEGLU23",
"DEGLU4",
"TEGLU19",
"HBSER5",
"HBGLU2",
"MEGLU2",
"TEGLU17",
"MEGLU7",
"MEINH3",
"HBGLU5",
"TEINH12",
"TEINH1",
"TEGLU4",
"TEGLU21",
"HBINH8"]

ds = loompy.connect(path_loomfile, mode= "r")

ds.ra.keys()       # Return list of row attribute names
ds.ca.keys()      # Return list of column attribute names

metaData = ds.ca["CellID", "ClusterName"]

metaData = pd.DataFrame(metaData, columns= ["CellID", "ClusterName"])

metaData = metaData[metaData['ClusterName'].isin(cellclasses)]

rawData = ds[:,metaData['ClusterName'].isin(cellclasses)]

rawData = pd.DataFrame(rawData, index =ds.ra["Gene"], columns = metaData["CellID"])

## Write files to disk

metaData.to_csv(path_or_buf=path_fileout_metaData, 
               header=True, 
               index=True, 
               index_label="Gene", 
               compression='infer', 
               quoting=None)

rawData.to_csv(path_or_buf=path_fileout_rawData, 
               header=True, 
               index=True, 
               index_label="Gene", 
               compression='infer', 
               quoting=None)

## wrap i[ 

ds.close()

print("Mousebrain preprocessing done!")

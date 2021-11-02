## #######################################################################################
##
## Get survey geoprecision for Mexico
##
##
## #######################################################################################

library(data.table)
library(dplyr)

working_dir <- '{REDACTED}'
in_fp <- file.path(working_dir, 'mex_bh_surveys.csv')
out_fp <- file.path(working_dir, 'mex_bh_surveys_with_geo.csv')

# Load functions to get spatial codebooks
lbd_core_dir <- '{REDACTED}'
source(file.path(lbd_core_dir, 'data_central/geocodebook_functions.R'))

## Load input NIDs and get all spatial data
bh_mex <- fread(in_fp)
nids <- unique(bh_mex$nid)

geo_cb <- get_geocodebooks(nids = nids, keep_all_columns = TRUE)
geo_cb[!is.na(lat) & !is.na(long), admin_level := 'GPS']
geo_cb <- geo_cb[admin_level == 'GPS' | (!is.na(location_code) & !is.na(shapefile)), ]
geo_cb[, unique_loc := paste(round(lat*1000), round(long*1000), sep='_')]
geo_cb[ admin_level != 'GPS', unique_loc := paste(location_code, shapefile, sep='_')]
geo_agg <- geo_cb[, .(
    cb_admin_level = max(admin_level, na.rm=T),
    cb_n_units = length(unique(unique_loc))
  ), by=nid
]

bh_mex_joined <- merge(bh_mex, geo_agg, by='nid', all.x=TRUE)
fwrite(bh_mex_joined, file=out_fp)

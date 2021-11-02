## #######################################################################################
##
## PREPARE MEXICO NMR INPUT DATA
## AUTHOR: Nat Henry
## CREATED: 18 August 2021
## PURPOSE: Prep script for all (non-simulation) inputs to Mexico NMR model
##
## #######################################################################################

in_libs <- c('haven','data.table','yaml')
invisible(lapply(in_libs, library, character.only=TRUE))

# Set analysis years
a_years <- 2009:2010

# Set input and output paths
data_version <- '20210818'
repo_dir <- '{REDACTED}'
config <- yaml::read_yaml(file.path(repo_dir, 'mexico/config.yaml'))

out_dir <- file.path(config$dirs$prepped_data, data_version)
dir.create(out_dir, showWarnings = FALSE)

vr_dir <- file.path(config$dirs$vr_data, 'mex_adm2')

## LOAD RAW DATA ------------------------------------------------------------------------>

# Load census data (note - this might take a while, uses about 9G of memory)
cen <- haven::read_dta(
  "{REDACTED}/MEX_CENSUS_2010_ALL_VARIABLES_Y2018M08D10.DTA"
)
# NMR: Subset to births in years, excluding the final month of births
cen_nmr <- data.table::as.data.table(
  cen[cen$sex == 2 & cen$age > 12 & (cen$lastbyr %in% a_years), ]
)
# COVARIATES: Subset to women 15-49
cen_cov <- data.table::as.data.table(cen[cen$sex == 2 & cen$age >= 15 & cen$age <= 49, ])
rm(cen)

# Load location metadata
loc_meta <- data.table::fread(file.path(vr_dir, 'location_data/adm_meta_full_current.csv'))

# Load other covariates
covars_dt <- fread(file.path(config$dirs$temp_dir, 'ensanut_covariate_dataset.csv'))
covars_meta <- fread(file.path(config$dirs$temp_dir, 'ensanut_covariate_metadata.csv'))

# Load VR births and deaths for all analysis years
births <- fread(file.path(vr_dir, 'births/mex_adm2_births_current.csv'))
deaths <- fread(file.path(vr_dir, 'u5m/u5m_mex_adm2_current.csv'))


## PREPARE VR BIRTHS AND DEATHS --------------------------------------------------------->

# Aggregate and join into a single dataset
deaths_agg <- deaths[
  (birth_year %in% a_years) & (age_months == 0),
  .(vr_deaths = sum(deaths)),
  by=location_id
]
births_agg <- births[ (birth_year %in% a_years), .(vr_births=sum(births)), by=location_id ]
vr_combined <- merge(x=births_agg, y=deaths_agg, by='location_id', all=TRUE)
vr_combined[is.na(vr_deaths), vr_deaths := 0 ][is.na(vr_births), vr_births := 0]

# Drop the small number of births assigned to the admin0 or admin1 level
vr_combined <- (
  vr_combined[loc_meta, level:=i.level, on='location_id'][level==2,][,level := NULL]
)

## PREPARE CENSUS-BASED NEONATAL MORTALITY ---------------------------------------------->

# Drop observations with missing metadata and in the final month (not yet completed)
cen_nmr <- cen_nmr[ (lastbyr < 2010) | (lastbmo < 5), ]
cen_nmr <- cen_nmr[ (lastbmort %in% 1:2) & (agedeadyr != 998) & !(agedeadmo %in% c(97, 98)), ]

# Define whether a child died before 1 month of age
cen_nmr[, nmr_bool := as.integer((lastbmort==2) & (agedeadyr==997) & (agedeadmo==0))]
cen_nmr[, adm_code := geo2_mx2010 ]
cen_nmr_agg <- cen_nmr[, .(cen_deaths = sum(nmr_bool), cen_births = .N), by=adm_code]


## PREPARE CENSUS-BASED COVARIATES ------------------------------------------------------>

# Define booleans for each relevant indicator
mapfun_1true <- function(x){if(x==1)return(1); if(x==2)return(0); return(as.integer(NA))}
mapfun_2true <- function(x){if(x==2)return(1); if(x==1)return(0); return(as.integer(NA))}
mapfun_water <- function(x){if(x %in% c(99, 0))return(as.integer(NA)); if(x==20)return(0); return(1)}

for(a_cov in c('electric', 'cell', 'indig')){
  cen_cov[[a_cov]] <- sapply(cen_cov[[a_cov]], mapfun_1true)
}
for(a_cov in c('refrig', 'lit')){
  cen_cov[[a_cov]] <- sapply(cen_cov[[a_cov]], mapfun_2true)
}
cen_cov$piped_water <- sapply(cen_cov$watsup, mapfun_water)
cen_cov[yrschool %in% c(98, 99), yrschool := NA ]

collapse_cols <- c('electric','cell','indig','refrig','lit','yrschool','piped_water')
cen_cov[, adm_code := geo2_mx2010 ]
cen_cov_agg <- cen_cov[, lapply(.SD, mean, na.rm=T), .SDcols=collapse_cols, by='adm_code']


## PREPARE ENSANUD COVARIATES ----------------------------------------------------------->

keep_covar_names <- covars_meta[use==1, var_name]
renamed_covars <- covars_meta[use==1, rename_to]
setnames(covars_dt, keep_covar_names, renamed_covars)
covars_sub <- covars_dt[, ..renamed_covars]
covars_sub[, adm_code := ad1_code * 1000 + ad2_code ][, c('ad1_code','ad2_code') := NULL ]


## Merge it all together using the location metadata template --------------------------->

keep_fields <- c('adm_code','location_id','uid','level','parent_code','adm_ascii_name')
loc_parents <- loc_meta[level==1, ..keep_fields]
loc_meta_sub <- loc_meta[ level==2, ..keep_fields]
full_data <- merge(x=loc_meta_sub, y=vr_combined, by='location_id', all=T)
full_data <- merge(x=full_data, y=cen_nmr_agg, by='adm_code', all=T)
full_data <- merge(x=full_data, y=cen_cov_agg, by='adm_code', all=T)
full_data <- merge(x=full_data, y=covars_sub, by='adm_code', all=T)
full_data[is.na(vr_births), vr_births := 0 ][is.na(vr_deaths), vr_deaths := 0 ]

# Drop a few anachronistic admin codes
drop_adm_codes <- c(23010, 23011)
full_data <- full_data[!(adm_code %in% drop_adm_codes), ]
full_data[adm_code == 20047, `:=` (cen_births=vr_births, cen_deaths=vr_deaths)]

# Check for missingness
for(a_var in colnames(full_data)){
  num_missing <- sum(is.na(full_data[[a_var]]))
  if(num_missing > 0) message(a_var, ' missing ', num_missing, ' rows.')
}

# Create a version that is aggregated by UID, for modeling
sum_cols <- c('vr_births','vr_deaths','cen_births','cen_deaths','pop_over20')
wm_cols <- c(
  'electric','cell','indig','refrig','lit','yrschool','lowwage','hcrowd','nohca',
  'propecon','informalwall','dirtfloor','clinicspc','piped_water'
)
grp_cols <- c('uid','parent_code')
full_data_agg <- merge(
  x = full_data[, lapply(.SD, sum), .SDcols=sum_cols, by=grp_cols],
  y = full_data[, lapply(.SD, weighted.mean, w=pop_over20), .SDcols=wm_cols, by=grp_cols],
  by = grp_cols
)

# Create normalized versions of all columns
make_norm <- function(vec, abslim=3){
  normed <- (vec - mean(vec)) / sd(vec)
  normed[normed > abslim] <- abslim
  normed[normed < -abslim] <- -abslim
  return(normed)
}
for(cname in wm_cols){
  full_data_agg[[paste0(cname,'_norm')]] <- make_norm(full_data_agg[[cname]])
}
# Do log for pop density
full_data_agg$popdens_norm <- make_norm(log(full_data_agg$pop_over20))


## Define social exclusion groupings
# Low exclusion (0)
full_data_agg[, excl_group := 0 ]
# Moderate exclusion (1)
full_data_agg[
  (clinicspc==0) & (indig > .5) & ((propecon < .5) | (lit < .9)),
  excl_group := 1
]
# High exclusion (2)
full_data_agg[
  (clinicspc == 0) & (indig > .5) & ((propecon < .25) | (lit < .75)),
  excl_group := 2
]


## Add holdouts
full_data_agg$idx_holdout_bh <- sample(1:5, size = nrow(full_data_agg), replace = TRUE)
full_data_agg$idx_holdout_vr <- sample(1:5, size = nrow(full_data_agg), replace = TRUE)
full_data_agg <- full_data_agg[order(uid)]
full_data_agg[, intercept := 1]
full_data_agg[ vr_deaths > vr_births, `:=` (vr_deaths = 0, vr_births = 0) ]


## Save everything to the output directory ---------------------------------------------->

fwrite(full_data, file = file.path(out_dir, 'prepped_data_full.csv'))
fwrite(full_data_agg, file = file.path(out_dir, 'prepped_data_stable.csv'))
fwrite(loc_meta, file = file.path(out_dir, 'location_metadata.csv'))

# Copy adjacency matrix to data dir
file.copy(
  from = file.path(vr_dir,'shapefile_stable/shapefile_stable_2000_2017_adj_mat_default.RDS'),
  to = file.path(out_dir, 'adjacency_matrix.RDS'),
  overwrite = TRUE
)

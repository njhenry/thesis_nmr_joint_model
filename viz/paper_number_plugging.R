## #######################################################################################
##
## NUMBER PLUG MEXICO CHAPTER
##
## Created: 30 September 2021
## Author: Nat Henry
## Purpose: Automated number-plugging for thesis chapter
##
## #######################################################################################

libs <- c('data.table','matrixStats','glue')
invisible(lapply(libs, library, character.only=T))

# Inputs
model_dir <- '{REDACTED}'
data_version <- '20210818'
runs <- list(
  wide = '20210930_allpois_wide_sds_2',
  narrow = '20210930_allpois_narrow_sds_2',
  sim = '20210930_allpois_sim'
)

admin_all <- fread(file.path(
  '{REDACTED}/mex_nmr/prepped_data/', data_version, 'location_metadata.csv'
))
ad_meta <- admin_all[level==2, .(uid, adm_ascii_name, parent_code)]
ad_meta <- ad_meta[, .(adm_name = paste(adm_ascii_name, collapse='; ')), by=.(uid, parent_code)]
parent_dt <- admin_all[level==1, .(adm_ascii_name, adm_code)]
colnames(parent_dt) <- c('parent_name','parent_code')
ad_meta <- merge(ad_meta, parent_dt, by='parent_code')[, parent_code := NULL ]


## DATA PREP ---------------------------------------------------------------------------->

read_dir <- function(run, fp) fread(glue("{model_dir}/{run}/{fp}.csv"))

preds <- lapply(runs, read_dir, fp='pred_summs')
pred_draws <- lapply(runs, read_dir, fp='pred_draws')
params <- lapply(runs, read_dir, fp='param_summs')
names(preds) <- names(pred_draws) <- names(params) <- names(runs)

for(pred in preds){
  pred[, vr_bias_mean := exp(log_vr_bias_mean)]
  pred[, vr_bias_lower := exp(log_vr_bias_lower)]
  pred[, vr_bias_upper := exp(log_vr_bias_upper)]
}
for(run in names(pred_draws)) pred_draws[[run]][,1] <- NULL

## SIMULATION RESULTS ------------------------------------------------------------------->

sim_data <- copy(preds$sim)
sim_data[, mort_in_ui := as.integer((lower < sim_mort) & (upper > sim_mort))]
sim_data[, bias_in_ui := as.integer((vr_bias_lower < sim_vr_bias) & (vr_bias_upper > sim_vr_bias))]

# Number of mortality results falling within UI bounds
sim_data[, sum(mort_in_ui)]

# Number of completeness results falling within UI bounds (by exclusion group)
sim_data[, .(in_ui = sum(bias_in_ui), total=.N), by=excl_group][order(excl_group)]


## NATIONAL RESULTS --------------------------------------------------------------------->

mort_natl <- copy(pred_draws$wide)
mort_natl$vr_births <- preds$wide$vr_births
draw_cols <- paste0("V",1:1000)
mort_natl <- mort_natl[, lapply(.SD, weighted.mean, w = vr_births), .SDcols=draw_cols]
mort_natl_summ <- data.table(
  mean = rowMeans(mort_natl)*1E3, lower = quantile(unlist(mort_natl), probs=c(0.025))*1E3,
  upper = quantile(unlist(mort_natl), probs=0.975)*1E3
)
mort_natl$mean <- rowMeans(mort_natl)


## MUNICIPALITY RESULTS AND RANKINGS ---------------------------------------------------->

mort_wide <- copy(preds$wide)
mfun <- function(mort) round(mort*1E3, 1)
mort_wide[, death_counts := mean * vr_births]
mort_wide[, `:=` (mean = mfun(mean), lower=mfun(lower), upper=mfun(upper)) ]
mort_wide[, mort_rank := frank(mean) ][, mort_state_rank := frank(mean), by=parent_code]
mort_wide <- merge(x=ad_meta, y=mort_wide, by='uid')[order(mort_rank)]

head(mort_wide)[, .(adm_name, parent_name, mean, lower, upper)]
tail(mort_wide)[, .(adm_name, parent_name, mean, lower, upper)]

mort_wide[, .(prop_low = sum(as.integer(mean < 5))/.N), by=parent_name][order(prop_low)]
mort_wide[, .(prop_high = sum(as.integer(mean > 10))/.N), by=parent_name][order(prop_high)]

# Amount concentrated in top 25?
dc_ordered <- mort_wide[order(-death_counts)][, prop_deaths := cumsum(death_counts)/sum(death_counts)]
num_making_25pct <- dc_ordered[prop_deaths < .25, .N + 1 ]

dc_ordered[1:num_making_25pct, .(adm_name, parent_name, death_counts, prop_deaths)]


## STATE TOTALS - WIDE AND NARROW MODELS ------------------------------------------------>

state_wide <- copy(preds$wide)[
  , lapply(.SD, weighted.mean, w=vr_births), .SDcols=c('mean','lower','upper'), by=parent_code
]
state_narrow <- copy(preds$narrow)[
  , lapply(.SD, weighted.mean, w=vr_births), .SDcols=c('mean','lower','upper'), by=parent_code
]
state_comp <- merge(x=state_wide, y=state_narrow, by='parent_code', suffixes = c('_w','_n'))
state_comp[, diff := mean_w - mean_n ][, diff_formatted := mfun(diff) ]
state_comp[, `:=` (
  mean_w = mfun(mean_w), lower_w = mfun(lower_w), upper_w=mfun(upper_w),
  mean_n=mfun(mean_n), lower_n=mfun(lower_n), upper_n=mfun(upper_n)
)]
state_comp[parent_dt, parent_name := i.parent_name, on='parent_code']
state_comp[order(diff)][c(1,.N),]
state_comp[order(mean_w)][c(1,.N)]


## MUNICIPAL TOTALS - WIDE AND NARROW MODELS -------------------------------------------->

m_wide <- copy(preds$wide)
m_narrow <- copy(preds$narrow)
keep_cols <- c('uid','mean','lower','upper')
m_comp <- merge(
  x=m_wide[, ..keep_cols], y=m_narrow[, ..keep_cols], by='uid', suffixes=c('_w','_n')
)
m_comp <- merge(x=ad_meta, y=m_comp, by='uid')
m_comp[, diff := mean_w - mean_n ][, diff_formatted := mfun(diff) ]
m_comp[, `:=` (
  mean_w = mfun(mean_w), lower_w = mfun(lower_w), upper_w=mfun(upper_w),
  mean_n=mfun(mean_n), lower_n=mfun(lower_n), upper_n=mfun(upper_n)
)]
m_comp[, rank_w := as.integer(frank(-mean_w))]
m_comp[, rank_w_state := as.integer(frank(-mean_w)), by=parent_name]
m_comp[, rank_n := as.integer(frank(-mean_n))]
m_comp[, rank_n_state := as.integer(frank(-mean_n)), by=parent_name]
m_comp[, rank_diff := rank_w - rank_n ]
m_comp[, rank_diff_state := rank_w_state - rank_n_state ]
m_comp[, num_in_state := .N, by=parent_name]

m_comp[abs(diff > 5E-4), .N]
m_comp[abs(diff > 1E-3), .N]
m_comp[order(-abs(rank_diff)), ][1:10, ]
m_comp[order(-abs(rank_diff_state)), ][1:10, ]

cor(m_comp$mean_w, abs(m_comp$diff))

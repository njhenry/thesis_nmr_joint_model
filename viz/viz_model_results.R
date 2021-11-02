## #######################################################################################
##
## VISUALIZE MEXICO MODEL RESULTS
##
## AUTHOR: Nat Henry, github: njhenry
## CREATED: 19 August 2021
## PURPOSE: Vizualize raw data as well as model results
##
## #######################################################################################

required_libs <- c(
  'data.table','dplyr','ggplot2','glue','grid','gridExtra','RColorBrewer','scales','sf'
)
invisible(lapply(required_libs, library, character.only=TRUE))

# Set file paths
repo_dir <- '{REDACTED}'
config <- yaml::read_yaml(file.path(repo_dir, 'mexico/config.yaml'))
data_version <- '20210818'
run_dates <- list(
  narrow = '20210930_allpois_narrow_sds_2',     # Real data, narrow VR bias
  wide = '20210930_allpois_wide_sds_2',         # Real data, wide VR bias
  sim_u = '20210930_allpois_sim' # Simulated data, unbiased BH
)
run_names <- names(run_dates)
# viz_dir <- file.path(config$dirs$viz_dir, gsub('-','',Sys.Date()))
viz_dir <- file.path(config$dirs$viz_dir, '20210930_2')
dir.create(viz_dir, showWarnings=FALSE)

# Load spatial metadata
loc_meta <- fread(file.path(config$dirs$prepped_data, data_version, 'location_metadata.csv'))

# Load shapefiles
ad2_sf <- sf::st_read(file.path(
  config$dirs$vr_data, 'mex_adm2/shapefile_stable/shapefile_stable_2000_2017.shp'
))[, c('uid','geometry')]
ad1_sf <- sf::st_read(file.path(
  config$dirs$vr_data, 'mex_adm2/shapefile_single_year/shapefile_single_year_2017_admin1.shp'
))[, c('GAUL_CODE', 'geometry')]

# Load model summaries
summs <- lapply(run_dates, function(rd){
  mod_dir <- file.path(config$dirs$runs, rd)
  summ_list <- list(
    param = fread(file.path(mod_dir, 'param_summs.csv')),
    pred = fread(file.path(mod_dir, 'pred_summs.csv')),
    fe = fread(file.path(mod_dir, 'fe_summs.csv')),
    sim_args = NULL
  )
  sim_args_fp <- file.path(mod_dir, 'sim_args.RDS')
  if(file.exists(sim_args_fp)) summ_list$sim_args <- readRDS(sim_args_fp)
  return(summ_list)
})

# Minor data prep - create merged admin2 metadata
ad1_meta <- loc_meta[level==1,][, .(location_id, adm_code, adm_ascii_name)]
colnames(ad1_meta) <- c('GAUL_CODE', 'parent_code','parent_name')
ad1_sf <- merge(x=ad1_sf, y=ad1_meta, by='GAUL_CODE', all.x=TRUE)
loc_merge_meta <- (loc_meta
  [level==2, ]
  [, .(parent_code=parent_code[1], adm_name=paste(adm_ascii_name,collapse=', ')), by=uid]
  [ad1_meta, on = 'parent_code']
)
ad2_sf <- merge(x=ad2_sf, y=loc_merge_meta, by='uid')
for(rd in run_names){
  summs[[rd]]$pred <- summs[[rd]]$pred[loc_merge_meta, on='uid'][, i.parent_code := NULL]
}


## FIG 1: DESCRIPTIVE NATIONAL PLOT ----------------------------------------------------->

excl_colors <- c(
  'Less marginalized' = '#b3cde3',
  'Moderately marginalized' = '#8c6bb1',
  'Severely marginalized' = '#6e016b'
)

excl_dt <- copy(summs$narrow$pred)
excl_dt$excl_label <- sapply(excl_dt$excl_group, function(ii) names(excl_colors)[ii+1])
excl_sf <- merge(x=ad2_sf, y=excl_dt[, .(uid, excl_label)])

mex_exclusion_fig <- ggplot() +
  geom_sf(data=excl_sf, aes(fill=excl_label), color='#222222', lwd=0.05) +
  geom_sf(data=ad1_sf, color='#222222', fill=NA, lwd=.25) +
  scale_fill_manual(values = excl_colors) +
  labs(fill = 'Municipality grouping') +
  coord_sf(crs=sf::st_crs(6372)) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
    legend.position = c(0.77, 0.75)
  )
png(file.path(viz_dir, 'excl_groups.png'), height=5.5, width=8, units='in', res=300)
print(mex_exclusion_fig)
dev.off()


## FIG 2: CHARACTERISTICS BY GROUPING --------------------------------------------------->

covar_labs <- data.table(
  cov = c('indig','lit','electric','refrig','lowwage','hcrowd','piped_water'),
  cov_title = c(
    'Identify as indigenous*', 'Literate*','Households electrified*',
    'Own refrigerator*','Low wage workers','Household crowding',
    'Piped water in home*'
  )
)

excl_melted <- melt(
  data = excl_dt,
  id.vars = c('uid', 'excl_label', 'excl_group'),
  measure.vars = covar_labs$cov,
  variable.name = 'cov'
)
excl_agg <- excl_melted[
  , .(val_med=median(value), val_low=quantile(value,0.25), val_high=quantile(value,0.75)),
  by=.(cov, excl_label, excl_group)
][covar_labs, on = 'cov']

exclusion_covars_fig <- ggplot(
    data=excl_agg,
    aes(y=val_med, ymin=val_low, ymax=val_high, x=cov_title, color=excl_label, fill=excl_label)  ) +
  geom_crossbar(position='dodge', color='#222222', width=.3, lwd=.15) +
  scale_fill_manual(values = excl_colors, aesthetics = c('fill')) +
  scale_y_continuous(labels = scales::percent) +
  labs(x='', y='Proportion by municipality', fill='', color='') +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

png(file.path(viz_dir, 'excl_covs.png'), height=5, width=7.5, units='in', res=300)
print(exclusion_covars_fig)
dev.off()


## SIMULATION - MORTALITY WHEN SPECIFIED CORRECTLY -------------------------------------->

sim_u <- copy(summs$sim_u$pred)
sim_u$excl_label <- sapply(sim_u$excl_group, function(ii) names(excl_colors)[ii+1])
sim_u[,sim_mort:=sim_mort*1E3][,mean:=mean*1E3][,lower:=lower*1E3][,upper:=upper*1E3]

plot_max <- max(c(sim_u$upper, sim_u$sim_mort))

sim_fig_a <- ggplot(
    data=sim_u, aes(color=excl_label,x=sim_mort,y=mean,ymin=lower,ymax=upper)
  ) +
  lims(x=c(0, plot_max), y=c(0, plot_max)) +
  labs(
    title="Neonatal mortality per 1,000", x='True (simulated)', y='Model estimate',
    color='Municipality\ngrouping'
  ) +
  geom_point(size=.5) +
  geom_linerange(lwd=.25, alpha=.7) +
  geom_abline(intercept=0, slope=1, linetype=2, color='#888888', alpha=.6) +
  scale_color_manual(values=excl_colors) +
  theme_bw() + theme(legend.position='right')

sim_u$vrb_ratio_mean <- exp(sim_u$log_vr_bias_mean)
sim_u$vrb_ratio_lower <- exp(sim_u$log_vr_bias_lower)
sim_u$vrb_ratio_upper <- exp(sim_u$log_vr_bias_upper)
sim_u$vrb_ratio_true <- copy(sim_u$sim_vr_bias)


sim_u <- sim_u[order(-excl_group), ]

resid_breaks <- c(1/10, 1/5, 1/2, 1, 2, 5, 10)
resid_labels <- c('1:10','1:5','1:2','1','2:1','5:1','10:1')
resid_plot_range <- c(.1, 14)
sim_fig_b <- ggplot(data=sim_u, aes(
    color = excl_label, x = vrb_ratio_true, y = vrb_ratio_mean, ymin = vrb_ratio_lower,
    ymax = vrb_ratio_upper
  )) +
  geom_linerange(lwd=0.25, alpha=.6) +
  geom_point(size=.5) +
  geom_abline(intercept=0, slope=1, linetype=2, color='#888888', alpha=.6) +
  geom_hline(yintercept = 1, linetype = 3, color='#888888', alpha=.6) +
  geom_vline(xintercept = 1, linetype = 3, color='#888888', alpha=.6) +
  labs(title='CRVS bias terms', x='True bias (simulated)', y='Model estimated bias', color='') +
  scale_color_manual(values=excl_colors) +
  scale_y_continuous(
    trans = 'log10', breaks = resid_breaks, limits = resid_plot_range,
    labels = resid_labels, oob=scales::squish
  ) +
  scale_x_continuous(
    trans = 'log10', breaks = resid_breaks, labels = resid_labels, limits = resid_plot_range
  ) +
  theme_bw() + theme(legend.position='none')

png(file.path(viz_dir,'sim_results.png'), height=4.5, width=10, units='in', res=300)
grid.arrange(
  ggplotGrob(sim_fig_a), ggplotGrob(sim_fig_b),
  layout_matrix = matrix(c(1,1,1,1,2,2,2), nrow=1)
)
dev.off()



## PRESENT NMR AND VR BIAS RESULTS FROM FULL MODEL -------------------------------------->

full_est <- copy(summs$wide$pred)

# Set up color scales
nmr_colors <- RColorBrewer::brewer.pal(n=9, name='RdPu')
bias_colors <- RColorBrewer::brewer.pal(n=9, name='BrBG')

# Translate bias into a ratio
full_est[, vrb_mean := exp(log_vr_bias_mean) ]
full_est[, mean_per_1k := mean * 1000 ]
full_sf <- merge(x=ad2_sf, y=full_est[, .(uid, mean_per_1k, vrb_mean)])

## Map NMR
full_nmr_fig <- ggplot() +
  geom_sf(data=full_sf, aes(fill=mean_per_1k), color='#222222', lwd=0.05) +
  geom_sf(data=ad1_sf, color='#222222', fill=NA, lwd=.25) +
  scale_fill_gradientn(
    colors = nmr_colors,
    limits = c(0, 15), breaks=seq(0, 15, by=3),
    oob = scales::squish
  ) +
  labs(fill = 'Neonatal\nMortality Rate') +
  coord_sf(crs=sf::st_crs(6372)) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
    legend.position = c(0.85, 0.75)
  )
png(file.path(viz_dir, 'nmr_wide_model.png'), height=5.5, width=8, units='in', res=300)
print(full_nmr_fig)
dev.off()

## Map VR bias
full_bias_fig <- ggplot() +
  geom_sf(data=full_sf, aes(fill=vrb_mean), color='#222222', lwd=0.05) +
  geom_sf(data=ad1_sf, color='#222222', fill=NA, lwd=.25) +
  scale_fill_gradientn(
    colors = bias_colors,
    limits = c(0.5, 2), breaks=c(0.5, 0.66, 1, 1.5, 2),
    labels=c('1:2','2:3','1:1','3:2','2:1'),
    trans = 'log10',
    oob = scales::squish
  ) +
  labs(fill = 'VR Bias Ratio\n(Mean)') +
  coord_sf(crs=sf::st_crs(6372)) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
    legend.position = c(0.85, 0.75)
  )
bias_fig_sub <- suppressMessages(
  full_bias_fig +
  coord_sf(crs=sf::st_crs(6372), xlim=c(2800000, 3700000), ylim=c(347500, 1078750)) +
  theme(
    legend.position = 'none',
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_rect(colour = NA, fill='white')
  )
)

png(file.path(viz_dir, 'vr_bias_wide_model.png'), height=5.5, width=8, units='in', res=300)
print(full_bias_fig)
# Add inset
npc <- function(x) unit(x, 'npc')
vp <- viewport(x=npc(.19), y=npc(.22), width=npc(.4), height=npc(.48))
grid::pushViewport(vp)
grid.draw(ggplotGrob(bias_fig_sub))
dev.off()


## Map differences between mortality rate at the admin2 and admin1 levels with the narrow
##  model

joined_dt <- merge(
  x = summs$wide$pred,
  y = summs$narrow$pred[, .(uid, mean)],
  by = 'uid',
  suffixes = c('_wide','_narrow')
)
joined_dt[, nmr_diff := (mean_wide - mean_narrow) * 1e3 ]
joined_agg <- joined_dt[, .(
    vr_births = sum(vr_births),
    mean_narrow = weighted.mean(mean_narrow, w=vr_births),
    mean_wide = weighted.mean(mean_wide, w=vr_births)
  ), by=.(parent_code,parent_name)
]
joined_agg[, nmr_diff := (mean_wide - mean_narrow) * 1E3 ]

diff_colors <- rev(RColorBrewer::brewer.pal(name='PiYG',n=9))
diff_breaks <- seq(-2, 2, by=1)
diff_labs <- c('-2', '-1', '0', '+1', '+2')
diff_lims <- range(diff_breaks)

# Top plot: Mortality rate difference at the state level

model_diff_ad1_sf <- merge(x=ad1_sf, y=joined_agg[, .(parent_code, nmr_diff)])

model_diff_ad1_fig <- ggplot() +
  geom_sf(data=model_diff_ad1_sf, aes(fill=nmr_diff), color='#222222', lwd=0.25) +
  scale_fill_gradientn(
    colors = diff_colors,
    limits = diff_lims, breaks=diff_breaks, labels = diff_labs,
    oob = scales::squish
  ) +
  coord_sf(crs=sf::st_crs(6372)) +
  theme_minimal() +
  labs(title = 'A', fill = 'NMR difference\n(per 1,000)') +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
    legend.position = c(0.88, 0.8),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

# Bottom plot: Mortality rate difference at the municipality level

model_diff_ad2_sf <- merge(x=ad2_sf, y=joined_dt[, .(uid, nmr_diff)])

model_diff_ad2_fig <- ggplot() +
  geom_sf(data=model_diff_ad2_sf, aes(fill=nmr_diff), color='#222222', lwd=0.05) +
  geom_sf(data=ad1_sf, color='#222222', fill=NA, lwd=.25) +
  scale_fill_gradientn(
    colors = diff_colors,
    limits = diff_lims, breaks=diff_breaks, labels = diff_labs,
    oob = scales::squish
  ) +
  labs(title = 'B') +
  coord_sf(crs=sf::st_crs(6372), xlim=c(2500000, 4000000), ylim=c(347500, 1322500)) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
    legend.position = 'none',
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

png(file.path(viz_dir, 'crvs_bias_diff.png'), height=11, width=8, units='in', res=300)
grid.arrange(
  ggplotGrob(model_diff_ad1_fig),
  ggplotGrob(model_diff_ad2_fig),
  layout_matrix = matrix(1:2, ncol=1)
)
dev.off()


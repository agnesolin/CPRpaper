# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #
##### CPR paper analysis script #####
# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# Agnes Olin agnes.olin@su.se
# 21 July 2021
# R version 4.1.0

#### initial setup ####


library(boot) 
library(forecast) 
library(ggplot2) 
library(ggpubr) 
library(lubridate) 
library(mgcv) 
library(raster)
library(renv)
library(rworldmap) 
library(rworldxtra) 
library(scales) 
library(signal) 
library(sp) 
library(suncalc) 
library(viridisLite) 

renv::restore() # managing package versions

options(scipen = 999)
par(family = "serif")


#### load and clean up data ####

## load CPR data ##
cpr = read.csv("data/CPR_data_updated.csv", sep = ";")

# sort out date format
cpr$Midpoint_Date_Local = as.Date(cpr$Midpoint_Date_Local)
cpr$Midpoint_Time_Local = as.POSIXlt(cpr$Midpoint_Time_Local, format = "%Y-%m-%d %H:%M", usetz = F)

# adjust sample id (remove numbering)
cpr$Sample_Id = sapply(strsplit(as.character(cpr$Sample_Id), '-'), '[', 1)

# remove Pseudo (included in Pseudo-Para)
cpr$Pseudocalanus.spp..adult.Total = NULL

# remove copepod eggs (not eaten)
cpr$Copepod.eggs = NULL

# save taxa names
taxa = names(cpr)[7:ncol(cpr)]
cpr = cpr[, c(1:6, 6 + order(taxa))]
taxa = sort(taxa)

# replace NAs with 0s (in original NA = absent, 0 = virtually absent but remains present)
cpr[, taxa][is.na(cpr[, taxa])] = 0



## load Stonehaven data ##
stonehaven = read.csv("data/Stonehaven_zooplankton.csv", sep = ";")
stonehaven_nauplii = read.csv("data/Stonehaven_nauplii.csv", sep = ";")

# sort out dates and merge
names(stonehaven)[1] = "date"
stonehaven$date = as.Date(stonehaven$date)
stonehaven_nauplii$date = as.Date(stonehaven_nauplii$date, format = "%Y-%m-%d")
stonehaven = merge(stonehaven,
                   stonehaven_nauplii,
                   by = "date",
                   all.x = T)
rm(stonehaven_nauplii)



## load L4 data ##
L4 = read.csv("data/L4.csv", sep = ";", fileEncoding = "UTF-8-BOM")

# sort out date format
L4$date = as.Date(L4$date)



#### fix filtering volume ####

# extracted from doi.org/10.1093/plankt/24.9.941 Fig 4a (volume filtered as function of greenness index)
c_index = c(0, 1, 2, 7)
volume = c(3.1874999999999996,
           3.119318181818182,
           3.0511363636363638,
           2.7272727272727275)
vol_mod = lm(volume ~ c_index)

par(mfrow = c(1,1))
plot(c_index, volume)
abline(vol_mod)


# adjust so it's per m3 (adjusted for clogging)
real_vol = predict(vol_mod, newdata = data.frame(c_index = cpr$Chlorophyll_Index))
cpr[, taxa] = cpr[, taxa] / real_vol



#### set paras ####

# distance paras
max_dist = 80 # for calculating correction factors
dist_lim = 135 # for aggregating data

# feeding season extent paras
start1 = 80
end1 = 165

start0 = 141
end0 = 212

#### MAP OF SAMPLES (FIG 1) #####

# load locations of sandeel grounds used
loc_sandeels = read.delim("data/locations_sandeel.csv")

# load location of sampling sites
stone_lat = 56.97305556
stone_long = -2.12055556

L4_lat = 50.25
L4_long = -4.21666667


# set up plot
png(
  "figures/sample_map.png",
  width = 26,
  height = 13,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "serif"
)
par(mfrow = c(1, 2), family = "serif", par(mar = c(4.5, 4.5, 1, 1)))

newmap = getMap(resolution = "high")


### map of all samples PLUS CORR FACTOR SITES
plot(
  newmap,
  col = "lightgrey",
  border = "lightgrey",
  bg = NULL,
  xlim = c(-24, 6),
  ylim = c(50, 66),
  xaxt = 'n',
  yaxt = 'n',
  ann = FALSE,
  xlab = "Longitude (?)",
  ylab = "Latitude (?)",
  cex.lab = 1.6,
  cex.axis = 1.3
)


axis(1, cex.axis = 1.3)
axis(2, cex.axis = 1.3)
box(which = "plot", lwd = 0.5)


# calculate distance of CPR samples to Stonehaven
locs = SpatialPoints(
  coords = cbind(cpr$Longitude, cpr$Latitude),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
C = SpatialPoints(
  coords = cbind(stone_long, stone_lat),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
cpr$distance1 = spDistsN1(locs, C, longlat = T)

# calculate distance of CPR samples to L4
C = SpatialPoints(
  coords = cbind(L4_long, L4_lat),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
cpr$distance2 = spDistsN1(locs, C, longlat = T)


# plot out all samples not used to calculate correction factors
points(
  cpr$Longitude[cpr$distance1 > max_dist & cpr$distance2 > max_dist],
  cpr$Latitude[cpr$distance1 > max_dist &
                 cpr$distance2 > max_dist],
  col = alpha("#94B59E", 0.1),
  cex = 0.7,
  pch = 16
)

# plot out samples used to calculate correction factors
points(
  cpr$Longitude[cpr$distance1 <= max_dist |
                  cpr$distance2 <= max_dist],
  cpr$Latitude[cpr$distance1 <= max_dist |
                 cpr$distance2 <= max_dist],
  col = alpha("black", 0.3),
  cex = 0.7,
  pch = 16
)

cpr$distance1 = NULL
cpr$distance2 = NULL

legend("bottomright", "a.", cex = 2.8, bty = "n")


### sandeel sites
plot(
  newmap,
  col = "lightgrey",
  border = "lightgrey",
  bg = NULL,
  xlim = c(-24, 6),
  ylim = c(50, 66),
  xaxt = 'n',
  yaxt = 'n',
  ann = FALSE,
  xlab = "Longitude (?)",
  ylab = "Latitude (?)",
  cex.lab = 1.6,
  cex.axis = 1.3
)

axis(1, cex.axis = 1.3)
axis(2, cex.axis = 1.3)
box(which = "plot", lwd = 0.5)

# loop through and plot out CPR samples for each location + add letters
i = 1
for (l in rev(c(2, 1, 4, 3, 5, 6))) { 
  locs = SpatialPoints(
    coords = cbind(cpr$Longitude, cpr$Latitude),
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  C = SpatialPoints(
    coords = cbind(loc_sandeels$centre_long[l], loc_sandeels$centre_lat[l]),
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  cpr$distance = spDistsN1(locs, C, longlat = T) # calculates distance from colony in km
  
  points(
    cpr$Longitude[cpr$distance <= dist_lim],
    cpr$Latitude[cpr$distance <= dist_lim],
    col = alpha("#94B59E", 0.1),
    cex = 0.7,
    pch = 16
    
  )
  
  text(loc_sandeels$centre_long[l],
       loc_sandeels$centre_lat[l],
       LETTERS[i],
       cex = 1.7)
  i = i + 1
  
}

legend("bottomright", "b.", cex = 2.8, bty = "n")

dev.off()




#### remove taxa that are only present at low levels in study area ####

# create a subset containing all data located in the included sandeel grounds
subset_taxa = data.frame()

for (i in 1:nrow(loc_sandeels)) {
  locs = SpatialPoints(
    coords = cbind(cpr$Longitude, cpr$Latitude),
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  C = SpatialPoints(
    coords = cbind(loc_sandeels$centre_long[i], loc_sandeels$centre_lat[i]),
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  cpr$distance = spDistsN1(locs, C, longlat = T) # calculates distance from colony in km
  sub = cpr[cpr$distance <= dist_lim &
              yday(cpr$Midpoint_Date_Local) >= start1 &
              yday(cpr$Midpoint_Date_Local) <= end0,]
  subset_taxa = rbind(subset_taxa, sub)
}
cpr$distance = NULL

# list of taxa present in less than 5% of samples
rm_taxa = taxa[colSums(subset_taxa[, taxa] != 0) / nrow(subset_taxa[, taxa]) < 0.05]

# adjust so that taxa to be removed do not include taxa which are part of bigger "groups"
rm_taxa = rm_taxa[!rm_taxa %in% c(
  "Calanus.V.VI.unidentified",
  "Centropages.hamatus",
  "Centropages.spp...Unidentified."
)]


cpr = cpr[, names(cpr)[!names(cpr)  %in%  rm_taxa]]
taxa = names(cpr)[7:ncol(cpr)]



#### calculate corr factors ####
source("help_scripts/corr_factors.R")


#### apply to full CPR dataset ####

## sort into day and night samps ##
data = data.frame(
  date = cpr$Midpoint_Date_Local,
  lon = cpr$Longitude,
  lat = cpr$Latitude
)
cpr[, (ncol(cpr) + 1):(ncol(cpr) + 2)] = getSunlightTimes(data = data,
                                                          keep = c("sunrise", "sunset"),
                                                          tz = "UTC")[, 4:5]

cpr$is.day = 0
cpr$is.day[cpr$Midpoint_Time_Local >= cpr$sunrise &
             cpr$Midpoint_Time_Local <= cpr$sunset] = 1



# corrections #

cpr_new = cpr

for (t in taxa) {
  cpr_new[cpr[, t] == 0 & !is.na(cpr[, t]) , t]  = 0
  
  cpr_new[cpr$is.day == 0 & cpr[, t] != 0 & !is.na(cpr[, t]) , t] =
    corr_factors$night_final[corr_factors$taxon == t] *
    cpr[cpr$is.day == 0 & cpr[, t] != 0 & !is.na(cpr[, t]) , t]
  
  cpr_new[cpr$is.day == 1 & cpr[, t] != 0 & !is.na(cpr[, t]) , t] =
    corr_factors$day_final[corr_factors$taxon == t] *
    cpr[cpr$is.day == 1 & cpr[, t] != 0 & !is.na(cpr[, t]) , t]
  
}

cpr_old = cpr
cpr = cpr_new




#### interpolate for sandeel areas ####

# loading in location info for grounds of interest
locations = c("DB", "FoF",  "ECG", "Shetland", "Faroes", "Iceland")
locations_plot_order = locations

# interpolation only during feeding
source("help_scripts/InterpolationFeeding.R")

# interpolation during full year
source("help_scripts/InterpolationFull.R")





#### LOOKING AT TRENDS AND MAKING PLOTS ####

#### PREP ####

# load prey traits
source("help_scripts/prey_info_pub.R")
prey_info = read.csv("data/prey_info.csv")


# sort out format so that taxon names in files match up
prey_info$taxa = gsub(" ", ".", as.character(prey_info$taxa))
prey_info$taxa = gsub("-", ".", as.character(prey_info$taxa))
prey_info$taxa = gsub("[()]", ".", as.character(prey_info$taxa))


# list of most important prey items
key_taxa = c(
  "Acartia.spp...unidentified.",
  "Appendicularia",
  "Calanus.finmarchicus",
  "Calanus.helgolandicus",
  "Calanus.I.IV",
  "Calanus.V.VI.unidentified",
  "Centropages.hamatus",
  "Centropages.spp...Unidentified.",
  "Centropages.typicus",
  "Cirripede.larvae..Total.",
  "Copepod.nauplii",
  "Decapoda.larvae..Total.",
  "Evadne.spp.",
  "Fish.eggs..Total.",
  "Fish.larvae",
  "Metridia.lucens",
  "Oithona.spp.",
  "Para.Pseudocalanus.spp.",
  "Podon.spp.",
  "Temora.longicornis"
)


# calculating total energy per day
df$total_energy = rowSums(t(t(df[, key_taxa]) * prey_info$energy[prey_info$taxa %in% key_taxa]))

# calculating median prey size
df$image_area_median = apply(
  df[, taxa],
  1,
  FUN = function(x)
    median(rep(
      sqrt(prey_info$image_area), as.numeric(round(t(x)))
    ))
)

# adding up abundances of small copepods
df$small_cops = df$Acartia.spp...unidentified. + df$Oithona.spp. + df$Para.Pseudocalanus.spp. + df$Temora.longicornis


# colours for 1+ and 0 group feeding seasons
col1 = viridis(10)[5]
col0 = "goldenrod2"


# set up table for GAM p-values
resp_variables = c("energy_tot", "energyInsideOutside", "size", "cf", "ch", "small_cops")

resp_variables = rep(resp_variables, each = 3)

gamRES = data.frame(
  resp = resp_variables,
  expl = rep(c("age", "smooth0", "smooth1"), length.out = length(resp_variables)),
  DB = rep(NA, length(resp_variables)),
  FoF = rep(NA, length(resp_variables)),
  ECG = rep(NA, length(resp_variables)),
  Shetland = rep(NA, length(resp_variables))
  
  
)


#### FIG 2 ####

source("help_scripts/FIG2.R")


#### FIG 3 ####

source("help_scripts/FIG3.R")

#### FIG 4 ####

source("help_scripts/FIG4.R")


#### FIG 5 ####

source("help_scripts/FIG5.R")


#### REPORTED QUANTITIES ####

# mean and stanrdard error values Iceland & Faroes energy & size
F1 = mean_se(df$total_energy[df$location == "Faroes" &
                               df$doy %in% start1:end1])/1000
c(F1[1], F1[3]-F1[1])

F0 = mean_se(df$total_energy[df$location == "Faroes" &
                               df$doy %in% start0:end0])/1000
c(F0[1], F0[3]-F0[1])


I1 = mean_se(df$total_energy[df$location == "Iceland" &
                               df$doy %in% start1:end1])/1000
c(I1[1], I1[3]-I1[1])

I0 = mean_se(df$total_energy[df$location == "Iceland" &
                               df$doy %in% start0:end0])/1000
c(I0[1], I0[3]-I0[1])



F1 = mean_se(df$image_area_median[df$location == "Faroes" &
                                    df$doy %in% start1:end1])
c(F1[1], F1[3]-F1[1])

F0 = mean_se(df$image_area_median[df$location == "Faroes" &
                                    df$doy %in% start0:end0])
c(F0[1], F0[3]-F0[1])


I1 = mean_se(df$image_area_median[df$location == "Iceland" &
                                    df$doy %in% start1:end1])
c(I1[1], I1[3]-I1[1])

I0 = mean_se(df$image_area_median[df$location == "Iceland" &
                                    df$doy %in% start0:end0])
c(I0[1], I0[3]-I0[1])

# sample sizes
length(unique(df$year[df$location == "Faroes"]))
length(unique(df$year[df$location == "Iceland"]))




### contribution Calanus ###

sub = coll_e[coll_e$group %in% c("Calanus finmarchicus", "Calanus I-IV" ),]

tot = aggregate(sub$energy_prop, list(sub$loc, sub$year), sum)
names(tot) = c("loc", "year", "tot")

aggregate(tot$tot, list(tot$loc), FUN = function(x) mean(x, na.rm = T))



#### SUPPLEMENTARY FIGURES ####
source("help_scripts/SuppPlots.R")



#### VIEW TABLES IN PAPER ####


## corr_factors ##
corr_factors$propSAMPnightX = 
  corr_factors$SAMPnightNOTZEROCPR.x/corr_factors$SAMPnight.x

corr_factors$propSAMPnightY = 
  corr_factors$SAMPnightNOTZEROCPR.y/corr_factors$SAMPnight.y

corr_factors$propSAMPdayX = 
  corr_factors$SAMPdayNOTZEROCPR.x/corr_factors$SAMPday.x

corr_factors$propSAMPdayY = 
  corr_factors$SAMPdayNOTZEROCPR.y/corr_factors$SAMPday.y

corr_factors[, -which(names(corr_factors) =="taxon")] = round(corr_factors[, -which(names(corr_factors) =="taxon")], digits = 2)

corr_factors[, c("night_low.y", "night_high.y",
                 "night_low.x", "night_high.x",
                 "day_low.y", "day_high.y",
                 "day_low.x", "day_high.x")] =
  round(
    corr_factors[, c("night_low.y", "night_high.y",
                     "night_low.x", "night_high.x",
                     "day_low.y", "day_high.y",
                     "day_low.x", "day_high.x")],
    digits = 1
    
  )

View(corr_factors[
  -which(corr_factors$taxon == "Euphausiacea.Total")
  ,
  c("taxon",
    
    "day_final",
    
    # day values Stonehaven
    "day.y",
    "day_low.y",
    "day_high.y",
    "SAMPday.y",
    "propSAMPdayY",
    "day_corr.y",
    
    # day values L4
    "day.x",
    "day_low.x",
    "day_high.x",
    "SAMPday.x",
    "propSAMPdayX",
    "day_corr.x",
    
    "night_final",
    
    # night values Stonehaven
    "night.y",
    "night_low.y",
    "night_high.y",
    "SAMPnight.y",
    "propSAMPnightY",
    "night_corr.y",
    
    # night values L4
    "night.x",
    "night_low.x",
    "night_high.x",
    "SAMPnight.x",
    "propSAMPnightX",
    "night_corr.x"
    
  )
  
])


## prey trait values ##
prey_info$image_area = signif(prey_info$image_area, digits = 2)

View(prey_info[,c("taxa", "size", "weight", "energy_density", "image_area")])



# gam results
gamRES[,c("DB", "FoF", "ECG", "Shetland")] = 
  round(gamRES[,c("DB", "FoF", "ECG", "Shetland")], digits = 2)

View(gamRES)


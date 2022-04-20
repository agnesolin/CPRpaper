#### CORRECTION FACTORS ####

set.seed(1) # for generating boot intervals

#### Stonehaven ####
source("help_scripts/StonehavenComparison.R")

# list of CPR taxa that have a corresponding group in the CPR dataset
poss_taxa = taxa[colSums(!is.na(stonehaven_comparison[, taxa])) != 0]



#### scale of coherence CPR and Stonehaven ####


# set up plot
png(
  "figures/Stonehaven_scale.png",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)

par(
  mfrow = c(5, 4),
  mar = c(3, 3, 2, 0.5) + 0.1,
  family = "sans",
  mgp = c(2, 1, 0)
)

# list of taxon names for plots
taxa_names = c(
  expression(italic("Acartia") ~ "spp."),
  "Appendicularia",
  expression(italic("Calanus finmarchicus")),
  expression(italic("Calanus helgolandicus")),
  expression(italic("Calanus") ~ "I-IV"),
  expression(italic("Centropages hamatus")),
  expression(italic("Centropages typicus")),
  "Cirripede larvae",
  "Copepod nauplii",
  "Decapoda larvae",
  "Fish eggs",
  "Fish larvae",
  "Hyperiidea",
  expression(italic("Metridia lucens")),
  expression(italic("Oithona") ~ "spp."),
  expression(italic("Para-Pseudocalanus") ~ "spp."),
  expression(italic("Temora longicornis"))
)


# calculate distance from all CPR samples to Stonehaven
locs = SpatialPoints(
  coords = cbind(cpr$Longitude, cpr$Latitude),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
C = SpatialPoints(
  coords = cbind(stone_long, stone_lat),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
cpr$distance = spDistsN1(locs, C, longlat = T) # calculates distance from colony in km


# loop through all taxa and make plot looking at corr as function of distance
for (t in poss_taxa) {
  
  # calculate how many non-NA samples exist per year-month combination
  stone_samps = aggregate(stonehaven_comparison[, t],
                          by = (list(
                            year(stonehaven_comparison$date),
                            month(stonehaven_comparison$date)
                          )),
                          function(x)
                            sum(!is.na(x)))
  names(stone_samps) = c("year", "month", "samps")
  
  # calculate mean abundances for each year-month combination
  stone_YM = aggregate(stonehaven_comparison[, t],
                       by = (list(
                         year(stonehaven_comparison$date),
                         month(stonehaven_comparison$date)
                       )),
                       function(x)
                         mean(x, na.rm = T))
  names(stone_YM) = c("year", "month", "stone")
  
  # merge samps and means, and only keep cases where at least 3 samples where available
  # then calculate correlation for year-month mean values calculated from cpr vs stonehaven data 
  stone_YM = merge(stone_samps, stone_YM)
  stone_YM = stone_YM[stone_YM$samps >= 3,]
  
  # vector to store calculated corrs
  ym_res = numeric()
  
  # repeat mean calculations for CPR samples, aggregated over various distances 
  for (d in 50:200) {
    cpr_sub = cpr[cpr$distance <= d,]
    
    cpr_samps = aggregate(cpr_sub[, t], by = (list(
      year(cpr_sub$Midpoint_Date_Local),
      month(cpr_sub$Midpoint_Date_Local)
    )), function(x)
      sum(!is.na(x)))
    names(cpr_samps) = c("year", "month", "samps")
    
    cpr_YM = aggregate(cpr_sub[, t], by = (list(
      year(cpr_sub$Midpoint_Date_Local),
      month(cpr_sub$Midpoint_Date_Local)
    )), function(x)
      mean(x, na.rm = T))
    names(cpr_YM) = c("year", "month", "cpr")
    
    cpr_YM = merge(cpr_samps, cpr_YM)
    cpr_YM = cpr_YM[cpr_YM$samps >= 3,]
    
    YM = merge(stone_YM, cpr_YM, by = c("year", "month"))
    ym_res = c(ym_res,
               cor(YM$stone, YM$cpr, use = "pairwise.complete.obs"))
    
    
  }
  
  # plot correlation as function of distance
  plot(
    50:200,
    ym_res,
    main = taxa_names[which(t == poss_taxa)],
    ylim = c(0, 0.8),
    col = "grey",
    pch = 16,
    cex.axis = 0.8,
    cex.main = 1.2,
    cex.lab = 1,
    xlab = "Radius (km)",
    ylab = expression(italic(r))
  )
  lines(50:200, ma(ym_res, 10), lwd = 2) # moving average, 10 km
  
  # final chosen distances
  abline(v = c(max_dist, dist_lim), lwd = 0.5)
  
  
  
  
}
dev.off()




#### ~~~ calculate correction factors ~~~ ####

# sub-setting
cpr_correction = cpr[cpr$distance <= max_dist &
                       year(cpr$Midpoint_Date_Local) >= min(year(stonehaven$date)) &
                       year(cpr$Midpoint_Date_Local) <= max(year(stonehaven$date)),] # start of time series

# sort into day and night samples
data = data.frame(
  date = cpr_correction$Midpoint_Date_Local,
  lon = cpr_correction$Longitude,
  lat = cpr_correction$Latitude
)
cpr_correction[, (ncol(cpr_correction) + 1):(ncol(cpr_correction) + 2)] = getSunlightTimes(data =
                                                                                             data,
                                                                                           keep = c("sunrise", "sunset"),
                                                                                           tz = "UTC")[, 4:5]

cpr_correction$is.day = 0
cpr_correction$is.day[cpr_correction$Midpoint_Time_Local >= cpr_correction$sunrise &
                        cpr_correction$Midpoint_Time_Local <= cpr_correction$sunset] = 1
cpr$distance = NULL


## considering all Centropages together
poss_taxa = poss_taxa[!poss_taxa %in% c("Centropages.hamatus", "Centropages.typicus")]
poss_taxa = c(poss_taxa, "Centropages")
cpr_correction$Centropages = rowSums(
  cbind(
    cpr_correction$Centropages.typicus,
    cpr_correction$Centropages.hamatus,
    cpr_correction$Centropages.spp...Unidentified.
  )
)
stonehaven_comparison$Centropages = rowSums(
  cbind(
    stonehaven_comparison$Centropages.typicus,
    stonehaven_comparison$Centropages.hamatus,
    stonehaven_comparison$Centropages.spp...Unidentified.
  ),
  na.rm = T
)

poss_taxa = sort(poss_taxa)


## set up empty data frame to store calculated correction factors
corr_factors = data.frame(
  taxon = poss_taxa,
  day = rep(NA, length(poss_taxa)),
  night = rep(NA, length(poss_taxa)),
  day_low = rep(NA, length(poss_taxa)),
  night_low = rep(NA, length(poss_taxa)),
  day_high = rep(NA, length(poss_taxa)),
  night_high = rep(NA, length(poss_taxa)),
  day_corr = rep(NA, length(poss_taxa)),
  night_corr = rep(NA, length(poss_taxa)),
  SAMPday = rep(NA, length(poss_taxa)),
  SAMPnight = rep(NA, length(poss_taxa)),
  SAMPdayNOTZEROCPR = rep(NA, length(poss_taxa)),
  SAMPnightNOTZEROCPR = rep(NA, length(poss_taxa)),
  SAMPdayNOTZEROSTONE = rep(NA, length(poss_taxa)),
  SAMPnightNOTZEROSTONE = rep(NA, length(poss_taxa))
)



for (day in 0:2) {
  # this goes through the cases night, day and day+night (relevant only for fish eggs)
  
  if (day == 2)
    day = c(0, 1) # pooling night and day
  
  # subsetting based on day and assigning samples to months and years
  cpr_sub = cpr_correction[cpr_correction$is.day %in% day, ]
  cpr_sub$month = ceiling(yday(cpr_sub$Midpoint_Date_Local) / ((365 + leap_year(
    year(cpr_sub$Midpoint_Date_Local)
  )) / 12))
  cpr_sub$year = year(cpr_sub$Midpoint_Date_Local)
  
  stonehaven_sub = stonehaven_comparison
  stonehaven_sub$year = year(stonehaven_sub$date)
  stonehaven_sub$month = ceiling(yday(stonehaven_sub$date) / ((365 + leap_year(
    year(stonehaven_sub$date)
  )) / 12))
  
  
  # calculating number of samples per year-month
  samps_stone = table(stonehaven_sub$year, stonehaven_sub$month)
  samps_cpr = table(cpr_sub$year, cpr_sub$month)
  
  
  ### make sure that contribution of each year-month combo is the same
  ### as no of samples are on the same order of magnitude, this works out as the number of samples being identical
  ### in cases where contribution from Stonehaven is greater, keep samples closest in time to CPR samples
  ### in cases where contribution from CPR is greater, keep samples closest in space to Stonehaven samples
  
  # start by removing all month/year combinations where there are more samples at Stonehaven
  more_stone = samps_stone > samps_cpr
  no_rm = samps_stone - samps_cpr
  
  stonehaven_sub$remove = NA
  
  for (i in 1:nrow(stonehaven_sub)) {
    if (more_stone[rownames(more_stone) == stonehaven_sub$year[i], colnames(more_stone) == stonehaven_sub$month[i]]) {
      # if more stonehaven samples
      
      # calculate avg day in cpr
      avg_cpr = mean(yday(cpr_sub$Midpoint_Date_Local[cpr_sub$year == stonehaven_sub$year[i] &
                                                        cpr_sub$month == stonehaven_sub$month[i]]))
      
      # ydays in cpr
      stone_dates = stonehaven_sub$date[stonehaven_sub$year == stonehaven_sub$year[i] &
                                          stonehaven_sub$month == stonehaven_sub$month[i]]
      stone_ydays = yday(stone_dates)
      
      # this is the position of the one(s) to be removed
      day_rm = order(abs(stone_ydays - avg_cpr), decreasing = T)[1:no_rm[rownames(more_stone) == stonehaven_sub$year[i], colnames(more_stone) == stonehaven_sub$month[i]]]
      
      if (stonehaven_sub$date[i] %in% stone_dates[day_rm])
        stonehaven_sub$remove[i] = TRUE
      else
        stonehaven_sub$remove[i] = FALSE
    }
    else
      stonehaven_sub$remove[i] = FALSE
  }
  
  # remove samples
  stonehaven_sub = stonehaven_sub[!stonehaven_sub$remove,]
  
  
  # now remove cpr_samples where there are more cpr samples (based on distance)
  more_cpr = samps_stone < samps_cpr
  no_rm = samps_cpr - samps_stone
  
  cpr_sub$remove = NA
  for (i in 1:nrow(cpr_sub)) {
    if (more_cpr[rownames(more_cpr) == cpr_sub$year[i], colnames(more_cpr) == cpr_sub$month[i]]) {
      # if more stonehaven samples
      
      
      cpr_dists = cpr_sub$distance[cpr_sub$year == cpr_sub$year[i] &
                                     cpr_sub$month == cpr_sub$month[i]]
      
      if (length(unique(cpr_dists)) != length(cpr_dists))
        print("WARNING")
      
      # this is the position of the one(s) to be removed
      day_rm = order(cpr_dists, decreasing = T)[1:no_rm[rownames(more_cpr) == cpr_sub$year[i], colnames(more_cpr) == cpr_sub$month[i]]]
      
      if (cpr_sub$distance[i] %in% cpr_dists[day_rm])
        cpr_sub$remove[i] = TRUE
      else
        cpr_sub$remove[i] = FALSE
    }
    else
      cpr_sub$remove[i] = FALSE
  }
  
  cpr_sub = cpr_sub[!cpr_sub$remove,]
  
  
  
  
  if (length(day) != 2) {
    # if day or night (ie calculating day- and night-specific correction factors - for all except fish eggs)
    
    for (taxon in poss_taxa) {
      # 0 = night, 1 = day
      
      ## mean and CIs
      corr_factors[corr_factors$taxon == taxon , 3 - day] = # mean
        mean(stonehaven_sub[!is.na(cpr_sub[, taxon]) &
                              !is.na(stonehaven_sub[, taxon])  , taxon]) /
        mean(cpr_sub[!is.na(cpr_sub[, taxon]) &
                       !is.na(stonehaven_sub[, taxon])   , taxon])
      
      
      data = data.frame(cpr = cpr_sub[!is.na(cpr_sub[, taxon]) &
                                        !is.na(stonehaven_sub[, taxon]) , taxon], 
                        stonehaven = stonehaven_sub[!is.na(cpr_sub[, taxon]) &
                                                      !is.na(stonehaven_sub[, taxon]) , taxon])
      ratio = function(data, indices) {
        d = data[indices,] # allows boot to select sample
        return(mean(d$stonehaven) / mean(d$cpr))
      }
      
      cis = boot.ci(boot(
        data = data,
        statistic = ratio,
        # bootstrapping with 1000 replications
        R = 1000
      ), type = "bca")$bca[4:5]
      
      corr_factors[corr_factors$taxon == taxon , 5 - day] = cis[1] #lowCI
      corr_factors[corr_factors$taxon == taxon , 7 - day] = cis[2] #highCI
      
      ## monthly corrs
      stone_temp = stonehaven_sub[!is.na(cpr_sub[, taxon]) &
                                    !is.na(stonehaven_sub[, taxon])  ,]
      cpr_temp = cpr_sub[!is.na(cpr_sub[, taxon]) &
                           !is.na(stonehaven_sub[, taxon])   ,]
      
      corr_factors[corr_factors$taxon == taxon , 9 - day] =
        cor(
          aggregate(stone_temp[, taxon], by = list(
            year(stone_temp$date) ,  ceiling(yday(stone_temp$date) / ((
              365 + leap_year(year(stone_temp$date))
            ) / 12))
          ), mean)$x,
          aggregate(cpr_temp[, taxon], by = list(
            year(cpr_temp$Midpoint_Date_Local) , ceiling(yday(cpr_temp$Midpoint_Date_Local) /
                                                           ((
                                                             365 + leap_year(year(cpr_temp$Midpoint_Date_Local))
                                                           ) / 12))
          ), mean)$x
        )
      
      
      ## sample sizes
      corr_factors[corr_factors$taxon == taxon , 11 - day] =
        sum(!is.na(cpr_sub[, taxon]) &
              !is.na(stonehaven_sub[, taxon])) # no of samples
      
      corr_factors[corr_factors$taxon == taxon , 13 - day] =
        sum(!is.na(cpr_sub[, taxon]) &
              !is.na(stonehaven_sub[, taxon])  &
              cpr_sub[, taxon] != 0) # non-zero cpr samples
      
      corr_factors[corr_factors$taxon == taxon , 15 - day] =
        sum(!is.na(cpr_sub[, taxon]) &
              !is.na(stonehaven_sub[, taxon]) &
              stonehaven_sub[, taxon] != 0)  # non-zero stone samples
      
      
    }
  }
  
  
  if (length(day) == 2) {
    # if pooling day and night - only for fish eggs
    taxon = "Fish.eggs..Total."
    
    corr_factors[corr_factors$taxon == taxon , 2:3] =
      mean(stonehaven_sub[!is.na(cpr_sub[, taxon]) &
                            !is.na(stonehaven_sub[, taxon])  , taxon]) /
      mean(cpr_sub[!is.na(cpr_sub[, taxon]) &
                     !is.na(stonehaven_sub[, taxon])   , taxon])
    
    data = data.frame(cpr = cpr_sub[!is.na(cpr_sub[, taxon]) &
                                      !is.na(stonehaven_sub[, taxon]) , taxon], stonehaven = stonehaven_sub[!is.na(cpr_sub[, taxon]) &
                                                                                                              !is.na(stonehaven_sub[, taxon]) , taxon])
    ratio = function(data, indices) {
      d = data[indices,] # allows boot to select sample
      return(mean(d$stone) / mean(d$cpr))
    }
    
    cis = boot.ci(boot(
      data = data,
      statistic = ratio,
      # bootstrapping with 1000 replications
      R = 1000
    ), type = "bca")$bca[4:5]
    
    corr_factors[corr_factors$taxon == taxon , 4:5] = cis[1] #lowCI
    corr_factors[corr_factors$taxon == taxon , 6:7] = cis[2] #highCI
    
    ## monthly corrs
    stone_temp = stonehaven_sub[!is.na(cpr_sub[, taxon]) &
                                  !is.na(stonehaven_sub[, taxon])  ,]
    cpr_temp = cpr_sub[!is.na(cpr_sub[, taxon]) &
                         !is.na(stonehaven_sub[, taxon])   ,]
    
    corr_factors[corr_factors$taxon == taxon , 8:9] =
      cor(
        aggregate(stone_temp[, taxon], by = list(
          year(stone_temp$date) ,  ceiling(yday(stone_temp$date) / ((
            365 + leap_year(year(stone_temp$date))
          ) / 12))
        ), mean)$x,
        aggregate(cpr_temp[, taxon], by = list(
          year(cpr_temp$Midpoint_Date_Local) , ceiling(yday(cpr_temp$Midpoint_Date_Local) /
                                                         ((
                                                           365 + leap_year(year(cpr_temp$Midpoint_Date_Local))
                                                         ) / 12))
        ), mean)$x
      )
    
    
    ## sample sizes
    corr_factors[corr_factors$taxon == taxon , 10:11] =
      sum(!is.na(cpr_sub[, taxon]) &
            !is.na(stonehaven_sub[, taxon])) # no of samples
    
    corr_factors[corr_factors$taxon == taxon , 12:13] =
      sum(!is.na(cpr_sub[, taxon]) &
            !is.na(stonehaven_sub[, taxon])  &
            cpr_sub[, taxon] != 0) # non-zero cpr samples
    
    corr_factors[corr_factors$taxon == taxon , 14:15] =
      sum(!is.na(cpr_sub[, taxon]) &
            !is.na(stonehaven_sub[, taxon]) &
            stonehaven_sub[, taxon] != 0)  # non-zero stone samples
    
    
  }
}







#### fix some formatting ####

corr_factors = merge(corr_factors, data.frame(taxon = taxa), all = T)


# sort format
corr_factors$taxon = as.character(corr_factors$taxon)
corr_factors = corr_factors[order(corr_factors$taxon),]
corr_factors$day = round(corr_factors$day, digits = 1)
corr_factors$night = round(corr_factors$night, digits = 1)

corr_factors_Stonehaven = corr_factors





#### L4 ####

source("help_scripts/L4Comparison.R")

# list of CPR taxa that have a corresponding group in the L4 dataset
poss_taxa = taxa[colSums(!is.na(L4_comparison[, taxa])) != 0]
poss_taxa = poss_taxa[poss_taxa != "Centropages.spp...Unidentified."]


#### scale of coherence CPR and L4 (as above) ####

png(
  "figures/L4_scale.png",
  width = 15,
  height = 15,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)

par(
  mfrow = c(4, 4),
  mar = c(3, 3, 2, 0.5) + 0.1,
  family = "sans",
  mgp = c(2, 1, 0)
)


taxa_names = c(
  expression(italic("Acartia") ~ "spp."),
  "Appendicularia",
  expression(italic("Calanus helgolandicus")),
  expression(italic("Centropages hamatus")),
  expression(italic("Centropages typicus")),
  "Cirripede larvae",
  "Copepod nauplii",
  "Decapoda larvae",
  "Fish eggs",
  "Fish larvae",
  "Hyperiidea",
  expression(italic("Metridia lucens")),
  expression(italic("Oithona") ~ "spp."),
  expression(italic("Para-Pseudocalanus") ~ "spp."),
  expression(italic("Temora longicornis"))
)

locs = SpatialPoints(
  coords = cbind(cpr$Longitude, cpr$Latitude),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
C = SpatialPoints(
  coords = cbind(L4_long, L4_lat),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
cpr$distance = spDistsN1(locs, C, longlat = T) # calculates distance from colony in km



for (t in poss_taxa) {
  L4_samps = aggregate(L4_comparison[, t],
                       by = (list(
                         year(L4_comparison$date), month(L4_comparison$date)
                       )),
                       function(x)
                         sum(!is.na(x)))
  
  L4_YM = aggregate(L4_comparison[, t],
                    by = (list(
                      year(L4_comparison$date), month(L4_comparison$date)
                    )),
                    function(x)
                      mean(x, na.rm = T))
  
  names(L4_samps) = c("year", "month", "samps")
  names(L4_YM) = c("year", "month", "L4")
  L4_YM = merge(L4_samps, L4_YM)
  L4_YM = L4_YM[L4_YM$samps >= 3,]
  
  ym_res = numeric()
  
  for (d in 50:200) {
    cpr_sub = cpr[cpr$distance <= d,]
    cpr_samps = aggregate(cpr_sub[, t], by = (list(
      year(cpr_sub$Midpoint_Date_Local),
      month(cpr_sub$Midpoint_Date_Local)
    )), function(x)
      sum(!is.na(x)))
    cpr_YM = aggregate(cpr_sub[, t], by = (list(
      year(cpr_sub$Midpoint_Date_Local),
      month(cpr_sub$Midpoint_Date_Local)
    )), function(x)
      mean(x, na.rm = T))
    names(cpr_samps) = c("year", "month", "samps")
    names(cpr_YM) = c("year", "month", "cpr")
    cpr_YM = merge(cpr_samps, cpr_YM)
    cpr_YM = cpr_YM[cpr_YM$samps >= 3,]
    
    YM = merge(L4_YM, cpr_YM, by = c("year", "month"))
    ym_res = c(ym_res, cor(YM$L4, YM$cpr, use = "pairwise.complete.obs"))
    
    
  }
  

  
  plot(
    50:200,
    ym_res,
    main = taxa_names[which(t == poss_taxa)],
    ylim = c(0, 0.8),
    col = "grey",
    pch = 16,
    cex.axis = 0.8,
    cex.main = 1.2,
    cex.lab = 1,
    xlab = "Radius (km)",
    ylab = expression(italic(r))
  )
  lines(50:200, ma(ym_res, 10), lwd = 2)

  abline(v = c(max_dist, dist_lim), lwd = 0.5)
  
  
  
  
}
dev.off()


#### ~~~ calculate correction factors (as above) ~~~ ####

# sub-setting
cpr_correction = cpr[cpr$distance <= max_dist &
                       year(cpr$Midpoint_Date_Local) >= min(year(L4$date)) &
                       year(cpr$Midpoint_Date_Local) <= max(year(L4$date)),] # start of time series

# sort into day and night samples
data = data.frame(
  date = cpr_correction$Midpoint_Date_Local,
  lon = cpr_correction$Longitude,
  lat = cpr_correction$Latitude
)
cpr_correction[, (ncol(cpr_correction) + 1):(ncol(cpr_correction) + 2)] = getSunlightTimes(data =
                                                                                             data,
                                                                                           keep = c("sunrise", "sunset"),
                                                                                           tz = "UTC")[, 4:5]

cpr_correction$is.day = 0
cpr_correction$is.day[cpr_correction$Midpoint_Time_Local >= cpr_correction$sunrise &
                        cpr_correction$Midpoint_Time_Local <= cpr_correction$sunset] = 1
cpr$distance = NULL


## developing correction factors ##

## considering all Centropages together
poss_taxa = poss_taxa[!poss_taxa %in% c("Centropages.hamatus", "Centropages.typicus")]
poss_taxa = c(poss_taxa, "Centropages")
cpr_correction$Centropages = rowSums(
  cbind(
    cpr_correction$Centropages.typicus,
    cpr_correction$Centropages.hamatus,
    cpr_correction$Centropages.spp...Unidentified.
  )
)
L4_comparison$Centropages = rowSums(
  cbind(
    L4_comparison$Centropages.typicus,
    L4_comparison$Centropages.hamatus,
    L4_comparison$Centropages.spp...Unidentified.
  ),
  na.rm = T
)
poss_taxa = sort(poss_taxa)


corr_factors = data.frame(
  taxon = poss_taxa,
  day = rep(NA, length(poss_taxa)),
  night = rep(NA, length(poss_taxa)),
  day_low = rep(NA, length(poss_taxa)),
  night_low = rep(NA, length(poss_taxa)),
  day_high = rep(NA, length(poss_taxa)),
  night_high = rep(NA, length(poss_taxa)),
  day_corr = rep(NA, length(poss_taxa)),
  night_corr = rep(NA, length(poss_taxa)),
  SAMPday = rep(NA, length(poss_taxa)),
  SAMPnight = rep(NA, length(poss_taxa)),
  SAMPdayNOTZEROCPR = rep(NA, length(poss_taxa)),
  SAMPnightNOTZEROCPR = rep(NA, length(poss_taxa)),
  SAMPdayNOTZEROSTONE = rep(NA, length(poss_taxa)),
  SAMPnightNOTZEROSTONE = rep(NA, length(poss_taxa))
)




for (day in 0:2) {
  # this goes through the cases night, day and day+night (relevant only for fish eggs)
  
  if (day == 2)
    day = c(0, 1) # pooling night and day
  
  # subsetting based on day and assigning samples to months and years
  cpr_sub = cpr_correction[cpr_correction$is.day %in% day, ]
  cpr_sub$month = ceiling(yday(cpr_sub$Midpoint_Date_Local) / ((365 + leap_year(
    year(cpr_sub$Midpoint_Date_Local)
  )) / 12))
  cpr_sub$year = year(cpr_sub$Midpoint_Date_Local)
  
  L4_sub = L4_comparison
  L4_sub$year = year(L4_sub$date)
  L4_sub$month = ceiling(yday(L4_sub$date) / ((365 + leap_year(
    year(L4_sub$date)
  )) / 12))
  
  # calculating number of samples per year-month
  samps_L4 = table(L4_sub$year, L4_sub$month)
  samps_cpr = table(cpr_sub$year, cpr_sub$month)
  missing_years = as.numeric(row.names(samps_L4))[!as.numeric(row.names(samps_L4)) %in% as.numeric(row.names(samps_cpr))]
  missing_years_repl = matrix(rep(0, 12 * length(missing_years)),
                              ncol = 12,
                              nrow = length(missing_years))
  row.names(missing_years_repl) = missing_years
  samps_cpr = rbind(samps_cpr, missing_years_repl)
  samps_cpr  = samps_cpr[order(row.names(samps_cpr)) , ]
  
  
  # start by removing all month/year combinations where there are more samples at L4
  more_L4 = samps_L4 > samps_cpr
  no_rm = samps_L4 - samps_cpr
  
  L4_sub$remove = NA
  
  for (i in 1:nrow(L4_sub)) {
    if (more_L4[rownames(more_L4) == L4_sub$year[i], colnames(more_L4) == L4_sub$month[i]]) {
      # if more L4 samples
      
      # calculate avg day in cpr
      avg_cpr = mean(yday(cpr_sub$Midpoint_Date_Local[cpr_sub$year == L4_sub$year[i] &
                                                        cpr_sub$month == L4_sub$month[i]]))
      
      # ydays in cpr
      L4_dates = L4_sub$date[L4_sub$year == L4_sub$year[i] &
                               L4_sub$month == L4_sub$month[i]]
      L4_ydays = yday(L4_dates)
      
      # this is the position of the one(s) to be removed
      day_rm = order(abs(L4_ydays - avg_cpr), decreasing = T)[1:no_rm[rownames(more_L4) == L4_sub$year[i], colnames(more_L4) == L4_sub$month[i]]]
      
      if (L4_sub$date[i] %in% L4_dates[day_rm])
        L4_sub$remove[i] = TRUE
      else
        L4_sub$remove[i] = FALSE
    }
    else
      L4_sub$remove[i] = FALSE
  }
  
  # remove samples
  L4_sub = L4_sub[!L4_sub$remove,]
  
  
  # now remove cpr_samples where there are more cpr samples (based on distance)
  more_cpr = samps_L4 < samps_cpr
  no_rm = samps_cpr - samps_L4
  
  cpr_sub$remove = NA
  for (i in 1:nrow(cpr_sub)) {
    if (more_cpr[rownames(more_cpr) == cpr_sub$year[i], colnames(more_cpr) == cpr_sub$month[i]]) {
      # if more L4 samples
      
      
      cpr_dists = cpr_sub$distance[cpr_sub$year == cpr_sub$year[i] &
                                     cpr_sub$month == cpr_sub$month[i]]
      
      if (length(unique(cpr_dists)) != length(cpr_dists))
        print("WARNING")
      
      # this is the position of the one(s) to be removed
      day_rm = order(cpr_dists, decreasing = T)[1:no_rm[rownames(more_cpr) == cpr_sub$year[i], colnames(more_cpr) == cpr_sub$month[i]]]
      
      if (cpr_sub$distance[i] %in% cpr_dists[day_rm])
        cpr_sub$remove[i] = TRUE
      else
        cpr_sub$remove[i] = FALSE
    }
    else
      cpr_sub$remove[i] = FALSE
  }
  
  cpr_sub = cpr_sub[!cpr_sub$remove,]
  
  
  
  
  if (length(day) != 2) {
    # if day or night (ie calculating day- and night-specific correction factors - for all except fish eggs)
    
    for (taxon in poss_taxa) {
      ## mean and CIs
      corr_factors[corr_factors$taxon == taxon , 3 - day] = # mean
        mean(L4_sub[!is.na(cpr_sub[, taxon]) &
                      !is.na(L4_sub[, taxon])  , taxon]) /
        mean(cpr_sub[!is.na(cpr_sub[, taxon]) &
                       !is.na(L4_sub[, taxon])   , taxon])
      
      
      data = data.frame(cpr = cpr_sub[!is.na(cpr_sub[, taxon]) &
                                        !is.na(L4_sub[, taxon]) , taxon], L4 = L4_sub[!is.na(cpr_sub[, taxon]) &
                                                                                        !is.na(L4_sub[, taxon]) , taxon])
      ratio = function(data, indices) {
        d = data[indices,] # allows boot to select sample
        return(mean(d$L4) / mean(d$cpr))
      }
      
      cis = boot.ci(boot(
        data = data,
        statistic = ratio,
        # bootstrapping with 1000 replications
        R = 1000
      ), type = "bca")$bca[4:5]
      
      corr_factors[corr_factors$taxon == taxon , 5 - day] = cis[1] #lowCI
      corr_factors[corr_factors$taxon == taxon , 7 - day] = cis[2] #highCI
      
      ## monthly corrs
      L4_temp = L4_sub[!is.na(cpr_sub[, taxon]) &
                         !is.na(L4_sub[, taxon])  ,]
      cpr_temp = cpr_sub[!is.na(cpr_sub[, taxon]) &
                           !is.na(L4_sub[, taxon])   ,]
      
      corr_factors[corr_factors$taxon == taxon , 9 - day] =
        cor(aggregate(L4_temp[, taxon], by = list(
          year(L4_temp$date) ,  ceiling(yday(L4_temp$date) / ((
            365 + leap_year(year(L4_temp$date))
          ) / 12))
        ), mean)$x,
        aggregate(cpr_temp[, taxon], by = list(
          year(cpr_temp$Midpoint_Date_Local) , ceiling(yday(cpr_temp$Midpoint_Date_Local) /
                                                         ((
                                                           365 + leap_year(year(cpr_temp$Midpoint_Date_Local))
                                                         ) / 12))
        ), mean)$x)
      
      
      
      
      ## sample sizes
      corr_factors[corr_factors$taxon == taxon , 11 - day] =
        sum(!is.na(cpr_sub[, taxon]) &
              !is.na(L4_sub[, taxon])) # no of samples
      
      corr_factors[corr_factors$taxon == taxon , 13 - day] =
        sum(!is.na(cpr_sub[, taxon]) &
              !is.na(L4_sub[, taxon])  &
              cpr_sub[, taxon] != 0) # non-zero cpr samples
      
      corr_factors[corr_factors$taxon == taxon , 15 - day] =
        sum(!is.na(cpr_sub[, taxon]) &
              !is.na(L4_sub[, taxon]) &
              L4_sub[, taxon] != 0)  # non-zero L4 samples
      
      
    }
  }
  
  
  if (length(day) == 2) {
    # if pooling day and night - only for fish eggs
    taxon = "Fish.eggs..Total."
    
    corr_factors[corr_factors$taxon == taxon , 2:3] =
      mean(L4_sub[!is.na(cpr_sub[, taxon]) &
                    !is.na(L4_sub[, taxon])  , taxon]) /
      mean(cpr_sub[!is.na(cpr_sub[, taxon]) &
                     !is.na(L4_sub[, taxon])   , taxon])
    
    data = data.frame(cpr = cpr_sub[!is.na(cpr_sub[, taxon]) &
                                      !is.na(L4_sub[, taxon]) , taxon], L4 = L4_sub[!is.na(cpr_sub[, taxon]) &
                                                                                      !is.na(L4_sub[, taxon]) , taxon])
    ratio = function(data, indices) {
      d = data[indices,] # allows boot to select sample
      return(mean(d$L4) / mean(d$cpr))
    }
    
    cis = boot.ci(boot(
      data = data,
      statistic = ratio,
      # bootstrapping with 1000 replications
      R = 1000
    ), type = "bca")$bca[4:5]
    
    corr_factors[corr_factors$taxon == taxon , 4:5] = cis[1] #lowCI
    corr_factors[corr_factors$taxon == taxon , 6:7] = cis[2] #highCI
    
    ## monthly corrs
    L4_temp = L4_sub[!is.na(cpr_sub[, taxon]) &
                       !is.na(L4_sub[, taxon])  ,]
    cpr_temp = cpr_sub[!is.na(cpr_sub[, taxon]) &
                         !is.na(L4_sub[, taxon])   ,]
    
    corr_factors[corr_factors$taxon == taxon , 8:9] =
      cor(aggregate(L4_temp[, taxon], by = list(
        year(L4_temp$date) ,  ceiling(yday(L4_temp$date) / ((
          365 + leap_year(year(L4_temp$date))
        ) / 12))
      ), mean)$x,
      aggregate(cpr_temp[, taxon], by = list(
        year(cpr_temp$Midpoint_Date_Local) , ceiling(yday(cpr_temp$Midpoint_Date_Local) /
                                                       ((
                                                         365 + leap_year(year(cpr_temp$Midpoint_Date_Local))
                                                       ) / 12))
      ), mean)$x)
    
    ## sample sizes
    corr_factors[corr_factors$taxon == taxon , 10:11] =
      sum(!is.na(cpr_sub[, taxon]) &
            !is.na(L4_sub[, taxon])) # no of samples
    
    corr_factors[corr_factors$taxon == taxon , 12:13] =
      sum(!is.na(cpr_sub[, taxon]) &
            !is.na(L4_sub[, taxon])  &
            cpr_sub[, taxon] != 0) # non-zero cpr samples
    
    corr_factors[corr_factors$taxon == taxon , 14:15] =
      sum(!is.na(cpr_sub[, taxon]) &
            !is.na(L4_sub[, taxon]) &
            L4_sub[, taxon] != 0)  # non-zero L4 samples
    
    
  }
}



# sort format
corr_factors = merge(corr_factors, data.frame(taxon = taxa), all = T)
corr_factors$taxon = as.character(corr_factors$taxon)
corr_factors = corr_factors[order(corr_factors$taxon),]
corr_factors$day = round(corr_factors$day, digits = 1)
corr_factors$night = round(corr_factors$night, digits = 1)

corr_factors_L4 = corr_factors



#### setting up final correction factors ####

corr_factors = merge(corr_factors_L4, corr_factors_Stonehaven, by = "taxon")

corr_factors = corr_factors[, c("taxon", sort(names(corr_factors)[2:ncol(corr_factors)]))]

corr_factors$day_final = NA
corr_factors$night_final = NA


## calculate mean if available for both L4 and Stonehaven ##


# Acartia

corr_factors$day_final[corr_factors$taxon == "Acartia.spp...unidentified."] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Acartia.spp...unidentified."],
         corr_factors$day.y[corr_factors$taxon == "Acartia.spp...unidentified."]))

corr_factors$night_final[corr_factors$taxon == "Acartia.spp...unidentified."] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Acartia.spp...unidentified."],
         corr_factors$night.y[corr_factors$taxon == "Acartia.spp...unidentified."]))



# Appendicularia
corr_factors$day_final[corr_factors$taxon == "Appendicularia"] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Appendicularia"],
         corr_factors$day.y[corr_factors$taxon == "Appendicularia"]))

corr_factors$night_final[corr_factors$taxon == "Appendicularia"] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Appendicularia"],
         corr_factors$night.y[corr_factors$taxon == "Appendicularia"]))


# Calanus finmarchicus
# not available L4 so use Stonehaven
corr_factors$day_final[corr_factors$taxon == "Calanus.finmarchicus"] =
  corr_factors$day.y[corr_factors$taxon == "Calanus.finmarchicus"]


corr_factors$night_final[corr_factors$taxon == "Calanus.finmarchicus"] =
  corr_factors$night.y[corr_factors$taxon == "Calanus.finmarchicus"]


# Calanus helgolandicus
corr_factors$day_final[corr_factors$taxon == "Calanus.helgolandicus"] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Calanus.helgolandicus"],
         corr_factors$day.y[corr_factors$taxon == "Calanus.helgolandicus"]))

corr_factors$night_final[corr_factors$taxon == "Calanus.helgolandicus"] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Calanus.helgolandicus"],
         corr_factors$night.y[corr_factors$taxon == "Calanus.helgolandicus"]))



# Calanus I-IV
# not available L4 so use Stonehaven
corr_factors$day_final[corr_factors$taxon == "Calanus.I.IV"] =
  corr_factors$day.y[corr_factors$taxon == "Calanus.I.IV"]

corr_factors$night_final[corr_factors$taxon == "Calanus.I.IV"] =
  corr_factors$night.y[corr_factors$taxon == "Calanus.I.IV"]


# Calanus V-VI

# "Calanus.V.VI.unidentified"
# use average
corr_factors[corr_factors$taxon == "Calanus.V.VI.unidentified" , c("day_final", "night_final")] =
  c(mean(c(corr_factors[corr_factors$taxon == "Calanus.finmarchicus", "day_final"], corr_factors[corr_factors$taxon == "Calanus.helgolandicus", "day_final"])),
    mean(c(corr_factors[corr_factors$taxon == "Calanus.finmarchicus", "night_final"], corr_factors[corr_factors$taxon == "Calanus.helgolandicus", "night_final"])))




# CENTROPAGES

corr_factors$day_final[corr_factors$taxon == "Centropages"] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Centropages"],
         corr_factors$day.y[corr_factors$taxon == "Centropages"]))

corr_factors$night_final[corr_factors$taxon == "Centropages"] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Centropages"],
         corr_factors$night.y[corr_factors$taxon == "Centropages"]))

corr_factors[corr_factors$taxon == "Centropages.hamatus", c("day_final", "night_final")] =
  corr_factors[corr_factors$taxon == "Centropages", c("day_final", "night_final")]

corr_factors[corr_factors$taxon == "Centropages.typicus", c("day_final", "night_final")] =
  corr_factors[corr_factors$taxon == "Centropages", c("day_final", "night_final")]

corr_factors[corr_factors$taxon == "Centropages.spp...Unidentified.", c("day_final", "night_final")] =
  corr_factors[corr_factors$taxon == "Centropages", c("day_final", "night_final")]

#corr_factors = corr_factors[corr_factors$taxon != "Centropages",]


# cirripede larvae
corr_factors$day_final[corr_factors$taxon == "Cirripede.larvae..Total."] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Cirripede.larvae..Total."],
         corr_factors$day.y[corr_factors$taxon == "Cirripede.larvae..Total."]))

corr_factors$night_final[corr_factors$taxon == "Cirripede.larvae..Total."] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Cirripede.larvae..Total."],
         corr_factors$night.y[corr_factors$taxon == "Cirripede.larvae..Total."]))


# copepod nauplii
corr_factors$day_final[corr_factors$taxon == "Copepod.nauplii"] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Copepod.nauplii"],
         corr_factors$day.y[corr_factors$taxon == "Copepod.nauplii"]))

corr_factors$night_final[corr_factors$taxon == "Copepod.nauplii"] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Copepod.nauplii"],
         corr_factors$night.y[corr_factors$taxon == "Copepod.nauplii"]))


# decapoda larvae
corr_factors$day_final[corr_factors$taxon == "Decapoda.larvae..Total."] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Decapoda.larvae..Total."],
         corr_factors$day.y[corr_factors$taxon == "Decapoda.larvae..Total."]))

corr_factors$night_final[corr_factors$taxon == "Decapoda.larvae..Total."] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Decapoda.larvae..Total."],
         corr_factors$night.y[corr_factors$taxon == "Decapoda.larvae..Total."]))





# fish eggs
corr_factors$day_final[corr_factors$taxon == "Fish.eggs..Total."] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Fish.eggs..Total."],
         corr_factors$day.y[corr_factors$taxon == "Fish.eggs..Total."]))

corr_factors$night_final[corr_factors$taxon == "Fish.eggs..Total."] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Fish.eggs..Total."],
         corr_factors$night.y[corr_factors$taxon == "Fish.eggs..Total."]))

# fish larvae
corr_factors$day_final[corr_factors$taxon == "Fish.larvae"] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Fish.larvae"],
         corr_factors$day.y[corr_factors$taxon == "Fish.larvae"]))

corr_factors$night_final[corr_factors$taxon == "Fish.larvae"] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Fish.larvae"],
         corr_factors$night.y[corr_factors$taxon == "Fish.larvae"]))


# hyperiidea
corr_factors$day_final[corr_factors$taxon == "Hyperiidea..Total."] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Hyperiidea..Total."],
         corr_factors$day.y[corr_factors$taxon == "Hyperiidea..Total."]))

corr_factors$night_final[corr_factors$taxon == "Hyperiidea..Total."] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Hyperiidea..Total."],
         corr_factors$night.y[corr_factors$taxon == "Hyperiidea..Total."]))


# metridia lucens
corr_factors$day_final[corr_factors$taxon == "Metridia.lucens"] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Metridia.lucens"],
         corr_factors$day.y[corr_factors$taxon == "Metridia.lucens"]))

corr_factors$night_final[corr_factors$taxon == "Metridia.lucens"] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Metridia.lucens"],
         corr_factors$night.y[corr_factors$taxon == "Metridia.lucens"]))

# oithona
corr_factors$day_final[corr_factors$taxon == "Oithona.spp."] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Oithona.spp."],
         corr_factors$day.y[corr_factors$taxon == "Oithona.spp."]))

corr_factors$night_final[corr_factors$taxon == "Oithona.spp."] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Oithona.spp."],
         corr_factors$night.y[corr_factors$taxon == "Oithona.spp."]))


# para-pseudo
corr_factors$day_final[corr_factors$taxon == "Para.Pseudocalanus.spp."] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Para.Pseudocalanus.spp."],
         corr_factors$day.y[corr_factors$taxon == "Para.Pseudocalanus.spp."]))

corr_factors$night_final[corr_factors$taxon == "Para.Pseudocalanus.spp."] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Para.Pseudocalanus.spp."],
         corr_factors$night.y[corr_factors$taxon == "Para.Pseudocalanus.spp."]))


# temora
corr_factors$day_final[corr_factors$taxon == "Temora.longicornis"] =
  mean(c(corr_factors$day.x[corr_factors$taxon == "Temora.longicornis"],
         corr_factors$day.y[corr_factors$taxon == "Temora.longicornis"]))

corr_factors$night_final[corr_factors$taxon == "Temora.longicornis"] =
  mean(c(corr_factors$night.x[corr_factors$taxon == "Temora.longicornis"],
         corr_factors$night.y[corr_factors$taxon == "Temora.longicornis"]))

corr_factors$day_final = round(corr_factors$day_final, digits = 1)
corr_factors$night_final = round(corr_factors$night_final, digits = 1)




set.seed(NULL)






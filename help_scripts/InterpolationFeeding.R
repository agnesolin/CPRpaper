#### INTERPOLATION FEEDING SEASON ####


# set up empty data frame to store data
df = data.frame(matrix(nrow = 0, ncol = length(taxa) + 3))
names(df) = c("location", "year", "doy", taxa)

# set up data frame to store some info on variability
var_df = data.frame(taxon = character(), metric = character(), month = numeric(), value = numeric())

# loop through aggregations of sandeel grounds
for (loc in locations) {
  centre_lat = loc_sandeels$centre_lat[loc_sandeels$loc == loc]
  centre_long = loc_sandeels$centre_long[loc_sandeels$loc == loc]
  
  # calculate distance of CPR samples to focal point
  locs = SpatialPoints(
    coords = cbind(cpr$Longitude, cpr$Latitude),
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  C = SpatialPoints(
    coords = cbind(centre_long, centre_lat),
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  cpr$distance = spDistsN1(locs, C, longlat = T) # calculates distance from colony in km
  
  # susbet samples within distance limit
  centre_subset  = cpr[cpr$distance <= dist_lim,]
  
  
  # loop through years in subset
  for (y in min(year(centre_subset$Midpoint_Date_Local)):max(year(centre_subset$Midpoint_Date_Local))) {
    
    
    # subset to year
    cpr_sub = centre_subset[year(centre_subset$Midpoint_Date_Local) == y,]
    
    if (nrow(cpr_sub) != 0) { # check that there are at least some data for this year
      
      # translate dates into monthly sampling periods of standardised length
      cpr_sub$month = ceiling(yday(cpr_sub$Midpoint_Date_Local) / ((365 + leap_year(y)) /
                                                                     12))
      
      # calculate midpoint for each monthly sampling period
      midpoints = data.frame(month = 1:12,
                             mid_point = round((365 + leap_year(y)) / 12 / 2 + (0:11) * (365 + leap_year(y)) /
                                                 12))
      
      # merge midpoints with subset
      cpr_sub = merge(cpr_sub, midpoints, by = "month")
      
      # create empty data frame to store data
      temp = data.frame(matrix(nrow = 365 + leap_year(y), ncol = length(taxa)))
      names(temp) = taxa
      
      
      if (sum(table(factor(cpr_sub$month, levels = 1:12))[3:8] >= 3) == 6) { # check that each standardised monthly sampling period March-August (feeding) contains at least 3 data points
        
        # loop through taxa
        for (t in taxa) {
          
          # for variability estimates
          var_df_temp_1 = data.frame(taxon = t, metric = "SS", 
                                   month = aggregate(cpr_sub[, t], by = list(cpr_sub$month), length)$Group.1,
                                   value = aggregate(cpr_sub[, t], by = list(cpr_sub$month), length)$x)
          var_df_temp_2 = data.frame(taxon = t, metric = "mean", 
                                     month = aggregate(cpr_sub[, t], by = list(cpr_sub$month), mean)$Group.1,
                                     value = aggregate(cpr_sub[, t], by = list(cpr_sub$month), mean)$x)
          var_df_temp_3 = data.frame(taxon = t, metric = "sd", 
                                     month = aggregate(cpr_sub[, t], by = list(cpr_sub$month), sd)$Group.1,
                                     value = aggregate(cpr_sub[, t], by = list(cpr_sub$month), sd)$x)
          
          var_df = rbind(var_df,var_df_temp_1, var_df_temp_2, var_df_temp_3)
          
          # calculate mean abundances for each monthly sampling period
          aggr = aggregate(cpr_sub[, t], by = list(cpr_sub$month), mean)
          names(aggr) = c("month", "abundance")
          aggr = merge(aggr, midpoints)
          
          # interpolate
          xf = 1:(365 + leap_year(y)) # days we want to interpolate for
          xp = c(1, aggr$mid_point, 365 + leap_year(y)) # days we have data for (+ start and end dates of year)
          yp = c(aggr$abundance[1], aggr$abundance, aggr$abundance[nrow(aggr)]) # calculated means (assume constant values first half of Jan, and last half of Dec)
          interp_cpr  = pchip(xp, yp, xf) # hermite interpolation
          
          temp[, t] = interp_cpr
          
          
          
          #### making a plot to illustrate process for a "randomly" chosen taxon, year and location ####
          if (t == "Acartia.spp...unidentified." &
              y == 2014 & loc == "DB") {
            
            # set up plot
            png(
              "figures/interpEX.png",
              width = 15,
              height = 10,
              units = 'cm',
              res = 200,
              pointsize = 9,
              family = "sans"
            )
            par(mar = c(5, 5, 2, 2))
            
            # plot transparent raw (but corrected) cpr data
            plot(
              yday(cpr_sub$Midpoint_Date_Local),
              cpr_sub$Acartia.spp...unidentified.,
              col = alpha("black", 0.1),
              pch = 16,
              cex = 1.3,
              xlim = c(start1, end0),
              ylim = c(0, max(
                cpr_sub$Acartia.spp...unidentified.[yday(cpr_sub$Midpoint_Date_Local) %in% start1:end0]
              )),
              xlab = "Day of year",
              ylab = expression(paste("Abundance m" ^ "-3")),
              cex.lab = 1.6,
              cex.axis = 1.5
            )
            
            # sampling periods
            abline(v =
                     round((365 + leap_year(y)) / 12 + (0:11) * (365 + leap_year(y)) /
                             12),
                   col = "darkgrey",
                   lty = 3)
            
            # raw data but for standardised times
            points(
              cpr_sub$mid_point,
              cpr_sub$Acartia.spp...unidentified.,
              pch = 8,
              cex = 2,
              col = alpha("goldenrod2", 0.3)
              
            )
            
            # calculated means
            points(xp,
                   yp,
                   pch = 1,
                   cex = 3,
                   col = viridis(10)[5])
            
            # linear interpolation
            lines(xp, yp, col = viridis(10)[5])
            
            # hermite interpolation
            lines(xf,
                  interp_cpr,
                  col = viridis(10)[5],
                  lwd = 2.5)
            
            dev.off()
            
          }
          
          
        }
        
        # bind data together
        temp = cbind(loc, y, 1:(365 + leap_year(y)), temp)
        names(temp) = c("location", "year", "doy", taxa)
        df = rbind(df, temp)
      }
    }
  }
  
}


# calculate mean sample size during feeding season
mean(var_df$value[var_df$metric == "SS" & var_df$month %in% 3:8])
sd(var_df$value[var_df$metric == "SS" & var_df$month %in% 3:8])


mean(var_df$value[var_df$metric == "sd" & var_df$month %in% 3:8 & var_df$taxon == "Calanus.finmarchicus"])
mean(var_df$value[var_df$metric == "mean" & var_df$month %in% 3:8 & var_df$taxon == "Calanus.finmarchicus"])





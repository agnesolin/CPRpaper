#### INTERPOLATION FULL YEAR ####

# as interpolation for feeding season, but requires at least 3 samples in all months (not just during the feeding season)

# set up empty data frame to store data
df_full = data.frame(matrix(nrow = 0, ncol = length(taxa) + 3))
names(df_full) = c("location", "year", "doy", taxa)



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
  
  centre_subset  = cpr[cpr$distance <= dist_lim,]
  
  
  for (y in min(year(centre_subset$Midpoint_Date_Local)):max(year(centre_subset$Midpoint_Date_Local))) {
    midpoints = data.frame(month = 1:12,
                           mid_point = round((365 + leap_year(y)) / 12 / 2 + (0:11) * (365 + leap_year(y)) /
                                               12))
    
    cpr_sub = centre_subset[year(centre_subset$Midpoint_Date_Local) == y,]
    
    if (nrow(cpr_sub) != 0) {
      cpr_sub$month = ceiling(yday(cpr_sub$Midpoint_Date_Local) / ((365 + leap_year(y)) /
                                                                     12))
      cpr_sub = merge(cpr_sub, midpoints, by = "month")
      
      temp = data.frame(matrix(nrow = 365 + leap_year(y), ncol = length(taxa)))
      names(temp) = taxa
      
      if (sum(table(factor(cpr_sub$month, levels = 1:12))[1:12] >= 3) == 12) {
        for (t in taxa) {
          aggr = aggregate(cpr_sub[, t], by = list(cpr_sub$month), mean)
          names(aggr) = c("month", "abundance")
          aggr = merge(aggr, midpoints)
          
          xf = 1:(365 + leap_year(y))
          xp = c(1, aggr$mid_point, 365 + leap_year(y))
          yp = c(aggr$abundance[1], aggr$abundance, aggr$abundance[nrow(aggr)])
          interp_cpr  = pchip(xp, yp, xf)
          
          temp[, t] = interp_cpr
          
          
        }
        
        temp = cbind(loc, y, 1:(365 + leap_year(y)), temp)
        names(temp) = c("location", "year", "doy", taxa)
        df_full = rbind(df_full, temp)
      }
    }
  }
  
}

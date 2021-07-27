#### SUPPLEMENTARY PLOTS ####


#### median driver ####


png(
  "~/nonSU/CPR/medianDRIVER.png",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "serif"
)

# all cases where a trend was identified
cases = data.frame(
  loc = c("Shetland", "FoF", "FoF", "DB", "DB"),
  plot_locs = c(
    "Shetland",
    "Firth of Forth",
    "Firth of Forth",
    "Dogger Bank",
    "Dogger Bank"
    
    
  ),
  age = c(1, 1, 0, 1, 0),
  median = c(
    0.3707384, #Shetland 1
    
    0.3648947, # FoF 1
    
    0.3909782,  # FoF 0
    
    0.3299747, # Dogger Bank 1
    
    0.3419371 # Dogger Bank 0
    )
) # predictions from model

# loop through cases
for (i in 1:nrow(cases)) {
  
  # colours
  if(cases$age[i] == 0) col = col0
  if(cases$age[i] == 1) col = col1
  
  
  # abundance of taxa smaller/larger than predicted median at start
  df$abu_small = rowSums(df[,taxa[sqrt(prey_info$image_area) <= cases$median[i]] ])
  df$abu_large = rowSums(df[,taxa[sqrt(prey_info$image_area) >= cases$median[i]] ])
  
  
  # aggregate per year and location and feeding season
  if (cases$age[i] == 1) {
    abuSMALL =
      aggregate(df$abu_small[df$doy >= start1 &
                               df$doy <= end1 & df$location == cases$loc[i]],
                list(df$location[df$doy >= start1 &
                                   df$doy <= end1 &
                                   df$location == cases$loc[i]], df$year[df$doy >= start1 &
                                                                           df$doy <= end1 & df$location == cases$loc[i]]),
                mean)
    names(abuSMALL) = c("loc", "year", "abu_small")
    
    abuLARGE =
      aggregate(df$abu_large[df$doy >= start1 &
                               df$doy <= end1 & df$location == cases$loc[i]],
                list(df$location[df$doy >= start1 &
                                   df$doy <= end1 &
                                   df$location == cases$loc[i]], df$year[df$doy >= start1 &
                                                                           df$doy <= end1 & df$location == cases$loc[i]]),
                mean)
    names(abuLARGE) = c("loc", "year", "abu_large")
    
  }
  
  
  if (cases$age[i] == 0) {
    abuSMALL =
      aggregate(df$abu_small[df$doy >= start0 &
                               df$doy <= end0 & df$location == cases$loc[i]],
                list(df$location[df$doy >= start0 &
                                   df$doy <= end0 &
                                   df$location == cases$loc[i]], df$year[df$doy >= start0 &
                                                                           df$doy <= end0 & df$location == cases$loc[i]]),
                mean)
    names(abuSMALL) = c("loc", "year", "abu_small")
    
    abuLARGE =
      aggregate(df$abu_large[df$doy >= start0 &
                               df$doy <= end0 & df$location == cases$loc[i]],
                list(df$location[df$doy >= start0 &
                                   df$doy <= end0 &
                                   df$location == cases$loc[i]], df$year[df$doy >= start0 &
                                                                           df$doy <= end0 & df$location == cases$loc[i]]),
                mean)
    names(abuLARGE) = c("loc", "year", "abu_large")
    
  }
  
  # calculate max value for setting ylim
  maxVALUE = ceiling(max(c(
    abuLARGE$abu_large, abuSMALL$abu_small
  )) / 1000)
  
  
  # some tricks to get nice plots/legends later
  abuSMALL$age = factor(cases$age[i], levels = c(0,1))
  abuSMALL = rbind(abuSMALL, c(NA, min(abuSMALL$year), NA, abs(cases$age[i]-1)))
  
  abuLARGE$age = factor(cases$age[i], levels = c(0,1))
  abuLARGE = rbind(abuLARGE, c(NA, min(abuLARGE$year), NA, abs(cases$age[i]-1)))
  
  
  ## SMALL ##
  
  # make scatterplot
  fig = ggplot(data = abuSMALL, aes(x=year, y=abu_small/1000, colour = age, shape = age)) +
    
    xlim(1958, 2018) +
    ylim(-(maxVALUE/7), maxVALUE) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_shape_manual(values = c(16,17)) +
    scale_color_manual(values=c(col0, col1))+
    
    geom_point(size = 1) +
    
    labs(x = "Year", y = expression(paste("Abundance (×1000 m" ^ "-3", ")"))) +
    
    theme_bw(base_size = 10) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="serif")) + 
    scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) + 
    scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    theme(legend.position = "bottom")
  
  
  # fit gams
  y = abuSMALL$abu_small / 1000
  x = abuSMALL$year
  gam.mod = gam(y ~ s(x))
  
  
  # check gam
  summary(gam.mod)
  par(mfrow = c(2,2)); gam.check(gam.mod0)
  
  
  # if p < 0.05, add line with prediction intervals
  if (summary(gam.mod)$s.pv < 0.05) {
    p = predict(
      gam.mod,
      data.frame(x = min(x, na.rm = T):max(x, na.rm = T)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(x, na.rm = T):max(x, na.rm = T)
    p = as.data.frame(p)
    p$ymin = p$fit - 1.96 * p$se.fit
    p$ymax = p$fit + 1.96 * p$se.fit
    
    abuSMALL = merge(abuSMALL, p, by = c("year"))
    
    fig = fig +
      geom_line(data = abuSMALL, aes(x = year, y = fit), size = 0.8, col = col) +
      geom_ribbon(data = abuSMALL, aes(x = year, ymin = ymin, ymax = ymax, col = NULL, fill = age),  alpha = .15)   + 
      scale_fill_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+"))
  }
  
  # save plots
  if(i == 1)  fig1small = fig + labs(title = "small")
  if(i == 2)  fig2small = fig
  if(i == 3)  fig3small = fig
  if(i == 4)  fig4small = fig
  if(i == 5)  fig5small = fig
  
  
  
  
  ## LARGE ##
  
  # make scatter plot
  fig =  ggplot(data = abuLARGE, aes(x=year, y=abu_large/1000, colour = age, shape = age)) +
    
    xlim(1958, 2018) +
    ylim(-(maxVALUE/7), maxVALUE) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_shape_manual(values = c(16,17)) +
    scale_color_manual(values=c(col0, col1))+
    
    geom_point(size = 1) +
    
    labs(x = "Year", y = expression(paste("Abundance (×1000 m" ^ "-3", ")"))) +
    
    theme_bw(base_size = 10) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="serif")) + 
    scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) + 
    scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    theme(legend.position = "bottom")
  
  
  # fit gam mod
  y = abuLARGE$abu_large / 1000
  x = abuLARGE$year
  gam.mod = gam(y ~ s(x))
  
  
  # check gam
  summary(gam.mod)
  par(mfrow = c(2,2)); gam.check(gam.mod0)
  
  
  # if p < 0.05, add line with prediction intervals
  if (summary(gam.mod)$s.pv < 0.05) {
    p = predict(
      gam.mod,
      data.frame(x = min(x):max(x)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(x):max(x)
    p = as.data.frame(p)
    p$ymin = p$fit - 1.96 * p$se.fit
    p$ymax = p$fit + 1.96 * p$se.fit
    
    abuLARGE = merge(abuLARGE, p, by = c("year"))
    
    fig = fig +
      geom_line(data = abuLARGE, aes(x = year, y = fit), size = 0.8, col = col) +
      geom_ribbon(data = abuLARGE, aes(x = year, ymin = ymin, ymax = ymax, col = NULL, fill = age),  alpha = .15)  + 
      scale_fill_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+"))
  }
  
  # save plots
  if(i == 1)  fig1large  = fig + labs(title = "large")
  if(i == 2)  fig2large  = fig
  if(i == 3)  fig3large  = fig
  if(i == 4)  fig4large  = fig
  if(i == 5)  fig5large  = fig
  
  
  
  
  
  
}

# arrange plots
print(
  ggarrange(fig1small, fig1large,
            fig2small, fig2large,
            fig3small, fig3large,
            fig4small, fig4large,
            fig5small, fig5large,
            ncol = 2,
            nrow = 5,
            common.legend = TRUE,
            legend = "bottom",
            labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j."),
            font.label = list(size = 15, family = "serif", face = "plain"),
            heights = c(1.1, 1, 1, 1, 1),
            label.x = 0.84,
            label.y = c(0.83, 0.83, rep(0.97,8))
  )
)


dev.off()


#### ICELAND + FAROES ####

# repeat all previous scatterplots for Iceland and Faroes

png(
  "~/nonSU/CPR/IceFar.png",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "serif"
)


for (loc in c("Iceland", "Faroes")) {
  
  #### TOTAL ENERGY ####
  locsub = df[df$loc == loc,]
  age1 = aggregate(locsub$total_energy[locsub$doy >= start1 &
                                         locsub$doy <= end1], list(locsub$year[locsub$doy  >= start1 &
                                                                                 locsub$doy <= end1]), mean)
  names(age1) = c("year", "mean")
  age1$season = as.factor(1)
  
  age0 = aggregate(locsub$total_energy[locsub$doy >= start0 &
                                         locsub$doy <= end0], list(locsub$year[locsub$doy  >= start0 &
                                                                                 locsub$doy <= end0]), mean)
  names(age0) = c("year", "mean")
  age0$season = as.factor(0)
  
  both = rbind(age0, age1)
  both$mean = both$mean/1000
  
  
  fig = ggplot(data = both, aes(x=year, y=mean, shape=season, color=season)) +
    
    xlim(1958, 2018) +
    ylim(0, 4) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_point(size = 1.3) +
    scale_shape_manual(values = c(16,17)) +
    scale_color_manual(values=c(col0, col1))+
    
    labs(x = "Year", y = expression(Total ~ energy ~ (kJ ~ m ^ -3))) +
    
    scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) + 
    scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
    guides(colour = guide_legend(override.aes = list(size = 1.7, fill = c(col0, col1)))) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="serif"))
  
  if(loc == "Iceland") Ice_te = fig 
  if(loc == "Faroes") Far_te = fig
  
  
  #### MEDIAN SIZE ####
  locsub = df[df$loc == loc,]
  age1 = aggregate(locsub$image_area_median[locsub$doy >= start1 &
                                              locsub$doy <= end1], list(locsub$year[locsub$doy  >= start1 &
                                                                                      locsub$doy <= end1]), mean)
  names(age1) = c("year", "mean")
  age1$season = as.factor(1)
  
  age0 = aggregate(locsub$image_area_median[locsub$doy >= start0 &
                                              locsub$doy <= end0], list(locsub$year[locsub$doy  >= start0 &
                                                                                      locsub$doy <= end0]), mean)
  names(age0) = c("year", "mean")
  age0$season = as.factor(0)
  
  both = rbind(age0, age1)
  
  
  fig = ggplot(data = both, aes(x=year, y=mean, shape=season, color=season)) +
    
    xlim(1958, 2018) +
    ylim(0.25, 0.9) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_point(size = 1.3) +
    scale_shape_manual(values = c(16,17)) +
    scale_color_manual(values=c(col0, col1))+
    
    labs(x = "Year", y = "Median prey size (mm)") +
    
    scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) + 
    scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
    guides(colour = guide_legend(override.aes = list(size = 1.7, fill = c(col0, col1)))) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="serif"))
  
  if(loc == "Iceland") Ice_ms = fig 
  if(loc == "Faroes") Far_ms = fig
  
  
  
  
  #### CAL FIN ####
  locsub = df[df$loc == loc,]
  age1 = aggregate(locsub$Calanus.finmarchicus[locsub$doy >= start1 &
                                                 locsub$doy <= end1], list(locsub$year[locsub$doy  >= start1 &
                                                                                         locsub$doy <= end1]), mean)
  names(age1) = c("year", "mean")
  age1$season = as.factor(1)
  
  age0 = aggregate(locsub$Calanus.finmarchicus[locsub$doy >= start0 &
                                                 locsub$doy <= end0], list(locsub$year[locsub$doy  >= start0 &
                                                                                         locsub$doy <= end0]), mean)
  names(age0) = c("year", "mean")
  age0$season = as.factor(0)
  
  both = rbind(age0, age1)
  
  
  fig = ggplot(data = both, aes(x=year, y=mean, shape=season, color=season)) +
    
    xlim(1958, 2018) +
    ylim(0, 600) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_point(size = 1.3) +
    scale_shape_manual(values = c(16,17)) +
    scale_color_manual(values=c(col0, col1))+
    
    labs(x = "Year", y = expression(italic("Calanus finmarchicus") ~ m ^ -3)) +
    
    scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) + 
    scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
    guides(colour = guide_legend(override.aes = list(size = 1.7, fill = c(col0, col1)))) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="serif"))
  
  if(loc == "Iceland") Ice_cf = fig 
  if(loc == "Faroes") Far_cf = fig
  
  
  #### CAL HEL ####
  locsub = df[df$loc == loc,]
  age1 = aggregate(locsub$Calanus.helgolandicus[locsub$doy >= start1 &
                                                  locsub$doy <= end1], list(locsub$year[locsub$doy  >= start1 &
                                                                                          locsub$doy <= end1]), mean)
  names(age1) = c("year", "mean")
  age1$season = as.factor(1)
  
  age0 = aggregate(locsub$Calanus.helgolandicus[locsub$doy >= start0 &
                                                  locsub$doy <= end0], list(locsub$year[locsub$doy  >= start0 &
                                                                                          locsub$doy <= end0]), mean)
  names(age0) = c("year", "mean")
  age0$season = as.factor(0)
  
  both = rbind(age0, age1)
  
  
  fig = ggplot(data = both, aes(x=year, y=mean, shape=season, color=season)) +
    
    xlim(1958, 2018) +
    ylim(0, 9) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_point(size = 1.3) +
    scale_shape_manual(values = c(16,17)) +
    scale_color_manual(values=c(col0, col1))+
    
    labs(x = "Year", y = expression(italic("Calanus helgolandicus") ~ m ^ -3)) +
    
    scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) + 
    scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
    guides(colour = guide_legend(override.aes = list(size = 1.7, fill = c(col0, col1)))) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="serif"))
  
  if(loc == "Iceland") Ice_ch = fig 
  if(loc == "Faroes") Far_ch = fig
  
  
  #### SMALL COPEPODS ####
  locsub = df[df$loc == loc,]
  age1 = aggregate(locsub$small_cops[locsub$doy >= start1 &
                                       locsub$doy <= end1], list(locsub$year[locsub$doy  >= start1 &
                                                                               locsub$doy <= end1]), mean)
  names(age1) = c("year", "mean")
  age1$season = as.factor(1)
  
  age0 = aggregate(locsub$small_cops[locsub$doy >= start0 &
                                       locsub$doy <= end0], list(locsub$year[locsub$doy  >= start0 &
                                                                               locsub$doy <= end0]), mean)
  names(age0) = c("year", "mean")
  age0$season = as.factor(0)
  
  both = rbind(age0, age1)
  both$mean = both$mean/1000
  
  
  fig = ggplot(data = both, aes(x=year, y=mean, shape=season, color=season)) +
    
    xlim(1958, 2018) +
    ylim(0, 9) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_point(size = 1.3) +
    scale_shape_manual(values = c(16,17)) +
    scale_color_manual(values=c(col0, col1))+
    
    labs(x = "Year", y = expression(paste("Small copepods (×1000 m" ^ "-3", ")"))) +
    
    scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) + 
    scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
    guides(colour = guide_legend(override.aes = list(size = 1.7, fill = c(col0, col1)))) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="serif"))
  
  if(loc == "Iceland") Ice_sc = fig 
  if(loc == "Faroes") Far_sc = fig
  
  
  
}

print(
  ggarrange(Ice_te, Far_te,
            Ice_ms, Far_ms,
            Ice_cf, Far_cf,
            Ice_ch, Far_ch,
            Ice_sc, Far_sc,
            ncol = 2,
            nrow = 5,
            common.legend = TRUE,
            legend = "bottom",
            labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j."),
            font.label = list(size = 13, family = "serif", face = "plain"),
            label.x = 0.88,
            label.y = 0.97)
)

dev.off()



#### SMALL COP PHENOLOGY ####



png(
  "~/nonSU/CPR/phenSMALLcop.png",
  width = 16,
  height = 8,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "serif"
)

par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))

## DB ##

# plot small copepods against day in DB for all years
plot(
  df_full$doy[df_full$location == "DB"],
  (
    df_full$Acartia.spp...unidentified.[df_full$location == "DB"] + df_full$Oithona.spp.[df_full$location == "DB"] + df_full$Para.Pseudocalanus.spp.[df_full$location == "DB"] + df_full$Temora.longicornis[df_full$location == "DB"]
  ) / 1000,
  pch = 16,
  col = alpha("grey", 0.1),
  xlab = "Day of year",
  ylab = expression(paste("Small copepods (×1000 m" ^ "-3", ")")),
  cex.lab =  1.2,
  cex.axis = 1.2
)

# plot line for mean values of small copepods per day
lines(aggregate((
  df_full$Acartia.spp...unidentified.[df_full$location == "DB"] + df_full$Oithona.spp.[df_full$location == "DB"] + df_full$Para.Pseudocalanus.spp.[df_full$location == "DB"] + df_full$Temora.longicornis[df_full$location == "DB"]
) / 1000,
list(df_full$doy[df_full$location == "DB"]),
mean
), lwd = 3, col = "grey45")

# vertical lines for extent of feeding season
abline(v = c(start1, end1),
       col = viridis(10)[5],
       lwd = 2)
abline(v = c(start0, end0),
       col = "goldenrod2",
       lwd = 2)

legend("topright", "a.", cex = 2, bty = "n")


## FoF ##

# plot small copepods against day in FoF for all years
plot(
  df_full$doy[df_full$location == "FoF"],
  (
    df_full$Acartia.spp...unidentified.[df_full$location == "FoF"] + df_full$Oithona.spp.[df_full$location == "FoF"] + df_full$Para.Pseudocalanus.spp.[df_full$location == "FoF"] + df_full$Temora.longicornis[df_full$location == "FoF"]
  ) / 1000,
  pch = 16,
  col = alpha("grey", 0.1),
  xlab = "Day of year",
  ylab = expression(paste("Small copepods (×1000 m" ^ "-3", ")")),
  cex.lab =  1.2,
  cex.axis = 1.2
)

# plot line for mean values of small copepods per day
lines(aggregate((
  df_full$Acartia.spp...unidentified.[df_full$location == "FoF"] + df_full$Oithona.spp.[df_full$location == "FoF"] + df_full$Para.Pseudocalanus.spp.[df_full$location == "FoF"] + df_full$Temora.longicornis[df_full$location == "FoF"]
) / 1000,
list(df_full$doy[df_full$location == "FoF"]),
mean
), lwd = 3, col = "grey45")

# vertical lines for extent of feeding season
abline(v = c(start1, end1),
       col = viridis(10)[5],
       lwd = 2)
abline(v = c(start0, end0),
       col = "goldenrod2",
       lwd = 2)

legend("topright", "b.", cex = 2, bty = "n")


dev.off()


#### CALANUS PHENOLOGY ####

png(
  "~/nonSU/CPR/phenCAL.png",
  width = 12,
  height = 6,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "serif"
)

par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))

# plot Cal i-iv in ECG against day
plot(
  df_full$doy[df_full$location == "ECG"],
  df_full$Calanus.I.IV[df_full$location == "ECG"],
  pch = 16,
  col = alpha("grey", 0.1),
  xlab = "Day of year",
  ylab = expression(italic("Calanus") ~ "I-IV" ~ m ^ -3),
  cex.lab =  1.5,
  cex.axis = 1.5
)

# line of means
lines(aggregate(df_full$Calanus.I.IV[df_full$location == "ECG"], list(df_full$doy[df_full$location == "ECG"]), mean),
      lwd = 3,
      col = "grey45")

# indicate extent of feeding season
abline(v = c(start1, end1),
       col = viridis(10)[5],
       lwd = 2)
abline(v = c(start0, end0),
       col = "goldenrod2",
       lwd = 2)

legend("topright", "a.", cex = 2, bty = "n")



# plot Cal fin in ECG against day
plot(
  df_full$doy[df_full$location == "ECG"],
  df_full$Calanus.finmarchicus[df_full$location == "ECG"],
  pch = 16,
  col = alpha("grey", 0.1),
  xlab = "Day of year",
  ylab = expression(italic("Calanus finmarchicus") ~ m ^ -3),
  cex.lab =  1.5,
  cex.axis = 1.5
)

# lines of mean values
lines(
  aggregate(df_full$Calanus.finmarchicus[df_full$location == "ECG"], list(df_full$doy[df_full$location == "ECG"]), mean),
  lwd = 3,
  col = "grey45"
)

# indicate extent of feeding season
abline(v = c(start1, end1),
       col = viridis(10)[5],
       lwd = 2)
abline(v = c(start0, end0),
       col = "goldenrod2",
       lwd = 2)

legend("topright", "b.", cex = 2, bty = "n")



dev.off()







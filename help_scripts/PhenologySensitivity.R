#### EXPLORING IMPACT OF PHENOLOGY #####


#### FIRST: SENSITIVITY CHOICE #####


cols = colorRampPalette(brewer.pal(11, "RdYlBu"))(length(-30:30))

#### energy ####

for (l in c("Shetland", "ECG", "FoF", "DB")) {
  
  
  
  #### age 1 ####
  
  
  preds = data.frame()
  
  for(i in -30:30){
    
    
    
    locsub = df[df$loc == l,]
    
    ## calculate mean values for each location and year ##
    
    age1 = aggregate(locsub$total_energy[locsub$doy >= start1+i &
                                           locsub$doy <= end1+i], 
                     list(locsub$year[locsub$doy  >= start1+i &
                                        locsub$doy <= end1+i]), 
                     mean)
    names(age1) = c("year", "mean")
    age1$mean = age1$mean / 1000
    
    ## fit gams ##
    
    gam.mod = gam(mean ~ s(year), data = age1, method = "REML", family = Gamma(link = "log"))
    
    
    p = predict(
      gam.mod,
      data.frame(year = min(age1$year):max(age1$year)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(age1$year):max(age1$year)
    p = as.data.frame(p)
    p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
    p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
    p$fit = gam.mod$family$linkinv(p$fit)
    p$diff = i
    
    preds = rbind(preds, p)  
    
  }
  
  fig_1 = ggplot(data = preds, aes(col = diff, group = diff)) +
    
    geom_line(aes(x = year, y = fit), size = 0.8) +
    
    xlim(1958, 2018) +
    ylim(0, 5.2) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_colour_gradientn(colours = rev(cols), name = "Shift") +
    
    labs(x = "Year", y = expression(Total ~ energy ~ (kJ ~ m ^ -3))) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_1 = fig_1 + theme(legend.position = "none")
  if(l == "ECG") ECG_1= fig_1 + theme(legend.position = "none")
  if(l == "FoF") FoF_1 = fig_1 + theme(legend.position = "none")
  if(l == "DB") DB_1 = fig_1 + theme(legend.position = "bottom")
  
  
  
  
  
  
  
  
  #### age 0 ####
  
  
  preds = data.frame()
  
  for(i in -30:30){
    
    
    
    locsub = df[df$loc == l,]
    
    ## calculate mean values for each location and year ##
    
    age0 = aggregate(locsub$total_energy[locsub$doy >= start0+i &
                                           locsub$doy <= end0+i], 
                     list(locsub$year[locsub$doy  >= start0+i &
                                        locsub$doy <= end0+i]), 
                     mean)
    names(age0) = c("year", "mean")
    age0$mean = age0$mean / 1000
    
    ## fit gams ##
    
    gam.mod = gam(mean ~ s(year), data = age0, method = "REML", family = Gamma(link = "log"))
    
    
    p = predict(
      gam.mod,
      data.frame(year = min(age0$year):max(age0$year)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(age0$year):max(age0$year)
    p = as.data.frame(p)
    p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
    p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
    p$fit = gam.mod$family$linkinv(p$fit)
    p$diff = i
    
    preds = rbind(preds, p)  
    
  }
  
  fig_0 = ggplot(data = preds, aes(col = diff, group = diff)) +
    
    geom_line(aes(x = year, y = fit), size = 0.8) +
    
    xlim(1958, 2018) +
    ylim(0, 5.2) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_colour_gradientn(colours = rev(cols), name = "Shift") +
    
    labs(x = "Year", y = expression(Total ~ energy ~ (kJ ~ m ^ -3))) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_0 = fig_0 + theme(legend.position = "none")
  if(l == "ECG") ECG_0= fig_0 + theme(legend.position = "none")
  if(l == "FoF") FoF_0 = fig_0 + theme(legend.position = "none")
  if(l == "DB")DB_0 = fig_0 + theme(legend.position = "bottom")
  
  
} 


# save plot #

png(
  "figures/EnergySensitivityChoice.png",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)

ggarrange(Shet_1, Shet_0, ECG_1, ECG_0, FoF_1, FoF_0, DB_1, DB_0,  ncol = 2, nrow = 4,
          font.label = list(size = 15, family = "sans", face = "plain"),
          labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h."), label.x = 0.85, label.y = 0.95,
          heights = c(1,1,1,1.3))

dev.off()







#### size ####


for (l in c("Shetland", "ECG", "FoF", "DB")) {
  
  
  
  #### age 1 ####
  
  
  preds = data.frame()
  
  for(i in -30:30){
    
    
    
    locsub = df[df$loc == l,]
    
    ## calculate mean values for each location and year ##
    
    age1 = aggregate(locsub$image_area_median[locsub$doy >= start1+i &
                                           locsub$doy <= end1+i], 
                     list(locsub$year[locsub$doy  >= start1+i &
                                        locsub$doy <= end1+i]), 
                     mean)
    names(age1) = c("year", "mean")

    
    ## fit gams ##
    
    gam.mod = gam(mean ~ s(year), data = age1, method = "REML", family = Gamma(link = "log"))
    
    
    p = predict(
      gam.mod,
      data.frame(year = min(age1$year):max(age1$year)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(age1$year):max(age1$year)
    p = as.data.frame(p)
    p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
    p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
    p$fit = gam.mod$family$linkinv(p$fit)
    p$diff = i
    
    preds = rbind(preds, p)  
    
  }
  
  fig_1 = ggplot(data = preds, aes(col = diff, group = diff)) +
    
    geom_line(aes(x = year, y = fit), size = 0.8) +
    
    xlim(1958, 2018) +
    ylim(0.3, 0.65) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_colour_gradientn(colours = rev(cols), name = "Shift") +
    
    labs(x = "Year", y = "Median prey size (mm)") +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_1 = fig_1 + theme(legend.position = "none")
  if(l == "ECG") ECG_1= fig_1 + theme(legend.position = "none")
  if(l == "FoF") FoF_1 = fig_1 + theme(legend.position = "none")
  if(l == "DB") DB_1 = fig_1 + theme(legend.position = "bottom")
  
  
  
  
  
  
  
  
  #### age 0 ####
  
  
  preds = data.frame()
  
  for(i in -30:30){
    
    
    
    locsub = df[df$loc == l,]
    
    ## calculate mean values for each location and year ##
    
    age0 = aggregate(locsub$image_area_median[locsub$doy >= start0+i &
                                           locsub$doy <= end0+i], 
                     list(locsub$year[locsub$doy  >= start0+i &
                                        locsub$doy <= end0+i]), 
                     mean)
    names(age0) = c("year", "mean")

    
    ## fit gams ##
    
    gam.mod = gam(mean ~ s(year), data = age0, method = "REML", family = Gamma(link = "log"))
    
    
    p = predict(
      gam.mod,
      data.frame(year = min(age0$year):max(age0$year)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(age0$year):max(age0$year)
    p = as.data.frame(p)
    p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
    p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
    p$fit = gam.mod$family$linkinv(p$fit)
    p$diff = i
    
    preds = rbind(preds, p)  
    
  }
  
  fig_0 = ggplot(data = preds, aes(col = diff, group = diff)) +
    
    geom_line(aes(x = year, y = fit), size = 0.8) +
    
    xlim(1958, 2018) +
    ylim(0.3, 0.65) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_colour_gradientn(colours = rev(cols), name = "Shift") +
    
    labs(x = "Year", y = "Median prey size (mm)") +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_0 = fig_0 + theme(legend.position = "none")
  if(l == "ECG") ECG_0= fig_0 + theme(legend.position = "none")
  if(l == "FoF") FoF_0 = fig_0 + theme(legend.position = "none")
  if(l == "DB")DB_0 = fig_0 + theme(legend.position = "bottom")
  
  
} 


# save plot #

png(
  "figures/SizeSensitivityChoice.png",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)

ggarrange(Shet_1, Shet_0, ECG_1, ECG_0, FoF_1, FoF_0, DB_1, DB_0,  ncol = 2, nrow = 4,
          font.label = list(size = 15, family = "sans", face = "plain"),
          labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h."), label.x = 0.85, label.y = 0.95,
          heights = c(1,1,1,1.3))

dev.off()













#### Calanus finmarchicus ####


for (l in c("Shetland", "ECG", "FoF", "DB")) {
  
  
  
  #### age 1 ####
  
  
  preds = data.frame()
  
  for(i in -30:30){
    
    
    
    locsub = df[df$loc == l,]
    
    ## calculate mean values for each location and year ##
    
    age1 = aggregate(locsub$Calanus.finmarchicus[locsub$doy >= start1+i &
                                                locsub$doy <= end1+i], 
                     list(locsub$year[locsub$doy  >= start1+i &
                                        locsub$doy <= end1+i]), 
                     mean)
    names(age1) = c("year", "mean")
    age1$mean[age1$mean == 0] = 2e-16
    
    ## fit gams ##
    
    gam.mod = gam(mean ~ s(year), data = age1, method = "REML", family = Gamma(link = "log"))
    
    
    p = predict(
      gam.mod,
      data.frame(year = min(age1$year):max(age1$year)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(age1$year):max(age1$year)
    p = as.data.frame(p)
    p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
    p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
    p$fit = gam.mod$family$linkinv(p$fit)
    p$diff = i
    
    preds = rbind(preds, p)  
    
  }
  
  fig_1 = ggplot(data = preds, aes(col = diff, group = diff)) +
    
    geom_line(aes(x = year, y = fit), size = 0.8) +
    
    xlim(1958, 2018) +
    ylim(0, 300) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_colour_gradientn(colours = rev(cols), name = "Shift") +
    
    labs(x = "Year", y = expression(italic("Calanus finmarchicus") ~ m ^ -3)) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_1 = fig_1 + theme(legend.position = "none")
  if(l == "ECG") ECG_1= fig_1 + theme(legend.position = "none")
  if(l == "FoF") FoF_1 = fig_1 + theme(legend.position = "none")
  if(l == "DB") DB_1 = fig_1 + theme(legend.position = "bottom")
  
  
  
  
  
  
  
  
  #### age 0 ####
  
  
  preds = data.frame()
  
  for(i in -30:30){
    
    
    
    locsub = df[df$loc == l,]
    
    ## calculate mean values for each location and year ##
    
    age0 = aggregate(locsub$Calanus.finmarchicus[locsub$doy >= start0+i &
                                                locsub$doy <= end0+i], 
                     list(locsub$year[locsub$doy  >= start0+i &
                                        locsub$doy <= end0+i]), 
                     mean)
    names(age0) = c("year", "mean")
    age0$mean[age0$mean == 0] = 2e-16
    
    
    ## fit gams ##
    
    gam.mod = gam(mean ~ s(year), data = age0, method = "REML", family = Gamma(link = "log"))
    
    
    p = predict(
      gam.mod,
      data.frame(year = min(age0$year):max(age0$year)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(age0$year):max(age0$year)
    p = as.data.frame(p)
    p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
    p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
    p$fit = gam.mod$family$linkinv(p$fit)
    p$diff = i
    
    preds = rbind(preds, p)  
    
  }
  
  fig_0 = ggplot(data = preds, aes(col = diff, group = diff)) +
    
    geom_line(aes(x = year, y = fit), size = 0.8) +
    
    xlim(1958, 2018) +
    ylim(0, 300) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_colour_gradientn(colours = rev(cols), name = "Shift") +
    
    labs(x = "Year", y = expression(italic("Calanus finmarchicus") ~ m ^ -3)) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_0 = fig_0 + theme(legend.position = "none")
  if(l == "ECG") ECG_0= fig_0 + theme(legend.position = "none")
  if(l == "FoF") FoF_0 = fig_0 + theme(legend.position = "none")
  if(l == "DB")DB_0 = fig_0 + theme(legend.position = "bottom")
  
  
} 


# save plot #

png(
  "figures/CfSensitivityChoice.png",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)

ggarrange(Shet_1, Shet_0, ECG_1, ECG_0, FoF_1, FoF_0, DB_1, DB_0,  ncol = 2, nrow = 4,
          font.label = list(size = 15, family = "sans", face = "plain"),
          labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h."), label.x = 0.85, label.y = 0.95,
          heights = c(1,1,1,1.3))

dev.off()














#### Calanus helgolandicus ####


for (l in c("Shetland", "ECG", "FoF", "DB")) {
  
  
  
  #### age 1 ####
  
  
  preds = data.frame()
  
  for(i in -30:30){
    
    
    
    locsub = df[df$loc == l,]
    
    ## calculate mean values for each location and year ##
    
    age1 = aggregate(locsub$Calanus.helgolandicus[locsub$doy >= start1+i &
                                                   locsub$doy <= end1+i], 
                     list(locsub$year[locsub$doy  >= start1+i &
                                        locsub$doy <= end1+i]), 
                     mean)
    names(age1) = c("year", "mean")
    age1$mean[age1$mean == 0] = 2e-16
    
    ## fit gams ##
    
    gam.mod = gam(mean ~ s(year), data = age1, method = "REML", family = Gamma(link = "log"))
    
    
    p = predict(
      gam.mod,
      data.frame(year = min(age1$year):max(age1$year)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(age1$year):max(age1$year)
    p = as.data.frame(p)
    p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
    p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
    p$fit = gam.mod$family$linkinv(p$fit)
    p$diff = i
    
    preds = rbind(preds, p)  
    
  }
  
  fig_1 = ggplot(data = preds, aes(col = diff, group = diff)) +
    
    geom_line(aes(x = year, y = fit), size = 0.8) +
    
    xlim(1958, 2018) +
    ylim(0, 300) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_colour_gradientn(colours = rev(cols), name = "Shift") +
    
    labs(x = "Year", y = expression(italic("Calanus helgolandicus") ~ m ^ -3)) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_1 = fig_1 + theme(legend.position = "none")
  if(l == "ECG") ECG_1= fig_1 + theme(legend.position = "none")
  if(l == "FoF") FoF_1 = fig_1 + theme(legend.position = "none")
  if(l == "DB") DB_1 = fig_1 + theme(legend.position = "bottom")
  
  
  
  
  
  
  
  
  #### age 0 ####
  
  
  preds = data.frame()
  
  for(i in -30:30){
    
    
    
    locsub = df[df$loc == l,]
    
    ## calculate mean values for each location and year ##
    
    age0 = aggregate(locsub$Calanus.helgolandicus[locsub$doy >= start0+i &
                                                   locsub$doy <= end0+i], 
                     list(locsub$year[locsub$doy  >= start0+i &
                                        locsub$doy <= end0+i]), 
                     mean)
    names(age0) = c("year", "mean")
    age0$mean[age0$mean == 0] = 2e-16
    
    
    ## fit gams ##
    
    gam.mod = gam(mean ~ s(year), data = age0, method = "REML", family = Gamma(link = "log"))
    
    
    p = predict(
      gam.mod,
      data.frame(year = min(age0$year):max(age0$year)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(age0$year):max(age0$year)
    p = as.data.frame(p)
    p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
    p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
    p$fit = gam.mod$family$linkinv(p$fit)
    p$diff = i
    
    preds = rbind(preds, p)  
    
  }
  
  fig_0 = ggplot(data = preds, aes(col = diff, group = diff)) +
    
    geom_line(aes(x = year, y = fit), size = 0.8) +
    
    xlim(1958, 2018) +
    ylim(0, 300) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_colour_gradientn(colours = rev(cols), name = "Shift") +
    
    labs(x = "Year", y = expression(italic("Calanus helgolandicus") ~ m ^ -3)) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_0 = fig_0 + theme(legend.position = "none")
  if(l == "ECG") ECG_0= fig_0 + theme(legend.position = "none")
  if(l == "FoF") FoF_0 = fig_0 + theme(legend.position = "none")
  if(l == "DB")DB_0 = fig_0 + theme(legend.position = "bottom")
  
  
} 


# save plot #

png(
  "figures/ChSensitivityChoice.png",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)

ggarrange(Shet_1, Shet_0, ECG_1, ECG_0, FoF_1, FoF_0, DB_1, DB_0,  ncol = 2, nrow = 4,
          font.label = list(size = 15, family = "sans", face = "plain"),
          labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h."), label.x = 0.85, label.y = 0.95,
          heights = c(1,1,1,1.3))

dev.off()














#### Small copepods ####


for (l in c("Shetland", "ECG", "FoF", "DB")) {
  
  
  
  #### age 1 ####
  
  
  preds = data.frame()
  
  for(i in -30:30){
    
    
    
    locsub = df[df$loc == l,]
    
    ## calculate mean values for each location and year ##
    
    age1 = aggregate(locsub$small_cops[locsub$doy >= start1+i &
                                                    locsub$doy <= end1+i], 
                     list(locsub$year[locsub$doy  >= start1+i &
                                        locsub$doy <= end1+i]), 
                     mean)
    names(age1) = c("year", "mean")
    age1$mean[age1$mean == 0] = 2e-16
    age1$mean = age1$mean/1000
    
    ## fit gams ##
    
    gam.mod = gam(mean ~ s(year), data = age1, method = "REML", family = Gamma(link = "log"))
    
    
    p = predict(
      gam.mod,
      data.frame(year = min(age1$year):max(age1$year)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(age1$year):max(age1$year)
    p = as.data.frame(p)
    p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
    p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
    p$fit = gam.mod$family$linkinv(p$fit)
    p$diff = i
    
    preds = rbind(preds, p)  
    
  }
  
  fig_1 = ggplot(data = preds, aes(col = diff, group = diff)) +
    
    geom_line(aes(x = year, y = fit), size = 0.8) +
    
    xlim(1958, 2018) +
    ylim(0, 17) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_colour_gradientn(colours = rev(cols), name = "Shift") +
    
    labs(x = "Year", y = expression(paste("Small copepods ("  %*% "1000 m" ^ "-3", ")"))) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_1 = fig_1 + theme(legend.position = "none")
  if(l == "ECG") ECG_1= fig_1 + theme(legend.position = "none")
  if(l == "FoF") FoF_1 = fig_1 + theme(legend.position = "none")
  if(l == "DB") DB_1 = fig_1 + theme(legend.position = "bottom")
  
  
  
  
  
  
  
  
  #### age 0 ####
  
  
  preds = data.frame()
  
  for(i in -30:30){
    
    
    
    locsub = df[df$loc == l,]
    
    ## calculate mean values for each location and year ##
    
    age0 = aggregate(locsub$small_cops[locsub$doy >= start0+i &
                                                   locsub$doy <= end0+i], 
                     list(locsub$year[locsub$doy  >= start0+i &
                                        locsub$doy <= end0+i]), 
                     mean)
    names(age0) = c("year", "mean")
    age0$mean[age0$mean == 0] = 2e-16
    age0$mean = age0$mean/1000
    
    ## fit gams ##
    
    gam.mod = gam(mean ~ s(year), data = age0, method = "REML", family = Gamma(link = "log"))
    
    
    p = predict(
      gam.mod,
      data.frame(year = min(age0$year):max(age0$year)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(age0$year):max(age0$year)
    p = as.data.frame(p)
    p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
    p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
    p$fit = gam.mod$family$linkinv(p$fit)
    p$diff = i
    
    preds = rbind(preds, p)  
    
  }
  
  fig_0 = ggplot(data = preds, aes(col = diff, group = diff)) +
    
    geom_line(aes(x = year, y = fit), size = 0.8) +
    
    xlim(1958, 2018) +
    ylim(0, 17) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    scale_colour_gradientn(colours = rev(cols), name = "Shift") +
    
    labs(x = "Year", y = expression(paste("Small copepods ("  %*% "1000 m" ^ "-3", ")"))) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_0 = fig_0 + theme(legend.position = "none")
  if(l == "ECG") ECG_0= fig_0 + theme(legend.position = "none")
  if(l == "FoF") FoF_0 = fig_0 + theme(legend.position = "none")
  if(l == "DB")DB_0 = fig_0 + theme(legend.position = "bottom")
  
  
} 


# save plot #

png(
  "figures/scSensitivityChoice.png",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)

ggarrange(Shet_1, Shet_0, ECG_1, ECG_0, FoF_1, FoF_0, DB_1, DB_0,  ncol = 2, nrow = 4,
          font.label = list(size = 15, family = "sans", face = "plain"),
          labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h."), label.x = 0.85, label.y = 0.95,
          heights = c(1,1,1,1.3))

dev.off()













#### optimal shift ####

cols = colorRampPalette(brewer.pal(11, "RdYlBu"))(length(-30:30))


for (l in c("Shetland", "ECG", "FoF", "DB")) {
  
  
  
  #### age 1 ####
  
  
  age1 = data.frame()
  
  locsub = df[df$loc == l,]
  
  for(y in 1958:2018){
    
    ## calculate mean values for each location and year ##
    if(y %in% locsub$year){
      
      # subset whole possible feeding window
      yearsub = locsub[locsub$doy >= start1-30 &
                                      locsub$doy <= end1+30 &
                                      locsub$year == y, ]
      
      starts = (start1-30):(start1+30)
      ends = (end1-30):(end1+30)
      
      res = numeric()
      
      for(i in 1:length(-30:30)){
        res = c(res,
       mean(yearsub$total_energy[yearsub$doy >= starts[i] &
                                        yearsub$doy <= ends[i]]))
         
      }
      
      start_opt = starts[which.max(res)]
      end_opt = ends[which.max(res)]
      
      
      age1_temp = aggregate(locsub$total_energy[locsub$doy >= start_opt &
                                                  locsub$doy <= end_opt  &
                                                  locsub$year == y
      ], 
      list(locsub$year[locsub$doy >= start_opt &
                         locsub$doy <= end_opt  &
                         locsub$year == y
      ]), 
      mean)
      
      age1 = rbind(age1, age1_temp)
    }
    
  }
  
  names(age1) = c("year", "mean")
  age1$mean = age1$mean / 1000
  
  ## fit gams ##
  
  gam.mod = gam(mean ~ s(year), data = age1, method = "REML", family = Gamma(link = "log"))
  
  
  p = predict(
    gam.mod,
    data.frame(year = min(age1$year):max(age1$year)),
    type = "link",
    se.fit = TRUE
    
  )
  p$year = min(age1$year):max(age1$year)
  p = as.data.frame(p)
  p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
  p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
  p$fit = gam.mod$family$linkinv(p$fit)
  
  preds = p
  
  
  
  
  fig_1 = ggplot() +
    
    geom_point(data = age1, aes(x = year, y = mean), col = col1) +
    
    geom_line(data = preds, aes(x = year, y = fit), size = 0.8, col = col1) +
    geom_ribbon(data = preds, aes(x = year, ymin = ymin, ymax = ymax),  alpha = .15, fill = col1) +
    
    xlim(1958, 2018) +
    ylim(0, 6) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    labs(x = "Year", y = expression(Total ~ energy ~ (kJ ~ m ^ -3))) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_1 = fig_1 + theme(legend.position = "none")
  if(l == "ECG") ECG_1= fig_1 + theme(legend.position = "none")
  if(l == "FoF") FoF_1 = fig_1 + theme(legend.position = "none")
  if(l == "DB") DB_1 = fig_1 + theme(legend.position = "none")
  
  
  
  
  
  
  
  

  #### age 0 ####
  
  
  age0 = data.frame()
  
  locsub = df[df$loc == l,]
  
  for(y in 1958:2018){
    
    ## calculate mean values for each location and year ##
    if(y %in% locsub$year){
      
      # subset whole possible feeding window
      yearsub = locsub[locsub$doy >= start0-30 &
                         locsub$doy <= end0+30 &
                         locsub$year == y, ]
      
      starts = (start0-30):(start0+30)
      ends = (end0-30):(end0+30)
      
      res = numeric()
      
      for(i in 1:length(-30:30)){
        res = c(res,
                mean(yearsub$total_energy[yearsub$doy >= starts[i] &
                                            yearsub$doy <= ends[i]]))
        
      }
      
      start_opt = starts[which.max(res)]
      end_opt = ends[which.max(res)]
      
      
      age0_temp = aggregate(locsub$total_energy[locsub$doy >= start_opt &
                                                  locsub$doy <= end_opt  &
                                                  locsub$year == y
      ], 
      list(locsub$year[locsub$doy >= start_opt &
                         locsub$doy <= end_opt  &
                         locsub$year == y
      ]), 
      mean)
      
      age0 = rbind(age0, age0_temp)
    }
    
  }
  
  names(age0) = c("year", "mean")
  age0$mean = age0$mean / 1000
  
  ## fit gams ##
  
  gam.mod = gam(mean ~ s(year), data = age0, method = "REML", family = Gamma(link = "log"))
  
  
  p = predict(
    gam.mod,
    data.frame(year = min(age0$year):max(age0$year)),
    type = "link",
    se.fit = TRUE
    
  )
  p$year = min(age0$year):max(age0$year)
  p = as.data.frame(p)
  p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
  p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
  p$fit = gam.mod$family$linkinv(p$fit)
  
  preds = p
  
  
  
  
  fig_0 = ggplot() +
    
    geom_point(data = age0, aes(x = year, y = mean), col = col0) +
    
    geom_line(data = preds, aes(x = year, y = fit), size = 0.8, col = col0) +
    geom_ribbon(data = preds, aes(x = year, ymin = ymin, ymax = ymax),  alpha = .15,  fill = col0) +
    
    xlim(1958, 2018) +
    ylim(0, 8) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +

    labs(x = "Year", y = expression(Total ~ energy ~ (kJ ~ m ^ -3))) +
    
    theme_bw(base_size = 8) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot
  if(l == "Shetland") Shet_0 = fig_0 + theme(legend.position = "none")
  if(l == "ECG") ECG_0= fig_0 + theme(legend.position = "none")
  if(l == "FoF") FoF_0 = fig_0 + theme(legend.position = "none")
  if(l == "DB") DB_0 = fig_0 + theme(legend.position = "none")
  
  
  
  
  
  
  
  
  
} 


# save plot #

png(
  "figures/EnergySensitivityOptimal.png",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)

ggarrange(Shet_1, Shet_0, ECG_1, ECG_0, FoF_1, FoF_0, DB_1, DB_0,  ncol = 2, nrow = 4,
          font.label = list(size = 15, family = "sans", face = "plain"),
          labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h."), label.x = 0.85, label.y = 0.95)

dev.off()







#### length required ####


for (l in c( "FoF", "DB")) {
  
  #### age 0 ####
  

  
  locsub = df[df$loc == l,]
  
  
  ## calculate mean values for each location and year ##
  
  age0 = aggregate(locsub$total_energy[locsub$doy >= start0 &
                                         locsub$doy <= end0], 
                   list(locsub$year[locsub$doy  >= start0 &
                                      locsub$doy <= end0]), 
                   sum)
  names(age0) = c("year", "sum")
  age0$sum = age0$sum / 1000
  
  ## fit gams ##
  
  gam.mod = gam(sum ~ s(year), data = age0, method = "REML", family = Gamma(link = "log"))
  
  
  goal = as.numeric(gam.mod$family$linkinv(predict(gam.mod, newdata = data.frame(year = 1958))))
  
  min_length = data.frame(year = 1958:2018, length = NA)
  
  for(y in 1958:2018){
    
    ## calculate mean values for each location and year ##
    if(y %in% locsub$year){
      
      # subset whole possible feeding window
      yearsub = locsub[locsub$doy >= start0-30 &
                         locsub$doy <= end0+30 &
                         locsub$year == y, ]
      

      res = data.frame(start = numeric(), end = numeric(), energy = numeric())
      
      for(i in (start0-30):(end0+30-1)){
        for(j in (i+1):end0+30){

          res = rbind(res,
                      data.frame(start = i, end = j, energy = 
                                   sum(yearsub$total_energy[yearsub$doy >= i &
                                                               yearsub$doy <= j])
                                 ))

        }
        
      }
      
      res$energy =res$energy/1000
      res$length = res$end-res$start
      
      min_length$length[min_length$year == y] = min(res$length[res$energy >= goal])

    }
    

    
  }
  min_length$length[is.infinite(min_length$length)] = 0
  min_length$poss = as.factor(min_length$length > 0)
  
  
  
  fig_0 = ggplot() +
    
    geom_point(data = min_length, aes(x = year, y = length, shape = poss)) +
    
    scale_shape_manual(values = c(4,16), labels = c("not possible", "possible"), name = "Window adjustment") +
    
  geom_hline(yintercept = length(start0:end0), linetype = "dashed") +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    labs(x = "Year", y = "Length of window") +
    
    theme_bw(base_size = 10) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="sans"))
  
  
  
  # save each plot

  if(l == "FoF") FoF_0 = fig_0 + theme(legend.position = "none")
  if(l == "DB") DB_0 = fig_0 + theme(legend.position = "none")
  
  
  
  
  
  
  
  
  
} 


# save plot #

png(
  "figures/EnergySensitivityLength.png",
  width = 15,
  height = 7.5,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)

ggarrange(FoF_0, DB_0,  ncol = 2, nrow = 1,
          font.label = list(size = 15, family = "sans", face = "plain"),
          labels = c("a.", "b."), label.x = 0.85, label.y = 0.95,
          common.legend = T,
          legend = "bottom")

dev.off()















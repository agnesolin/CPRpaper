##### FIG 2 ####


#### calculate energy inside/outside season ####

# calculate total energy for full year
df_full$total_energy = rowSums(t(t(df_full[, key_taxa]) * prey_info$energy[prey_info$taxa %in% key_taxa]))

# subset to locations with better data
df_full = df_full[df_full$location %in% c("DB",       "FoF",      "ECG",      "Shetland"), ]

# sum energy within feeding season each year for each location
age1 = aggregate(df_full$total_energy[df_full$doy >= start1 &
                                        df_full$doy <= end1],
                 list(df_full$location[df_full$doy  >= start1 &
                                         df_full$doy <= end1],
                      df_full$year[df_full$doy  >= start1 &
                                     df_full$doy <= end1]),
                 sum)

names(age1) = c("loc", "year", "mean")
age1$season = as.factor(1)


age0 = aggregate(df_full$total_energy[df_full$doy >= start0 &
                                        df_full$doy <= end0],
                 list(df_full$location[df_full$doy  >= start0 &
                                         df_full$doy <= end0],
                      df_full$year[df_full$doy  >= start0 &
                                     df_full$doy <= end0]),
                 sum)
names(age0) = c("loc", "year", "mean")
age0$season = as.factor(0)


# sum energy for full year for each year for each location
full = aggregate(df_full$total_energy,
                 list(df_full$location,
                      df_full$year),
                 sum)
names(full) = c("loc", "year", "meanFULL")

# calculate ratio of energy available within feeding season compared to energy available during full year
age1$ratio = age1$mean/full$meanFULL
age0$ratio = age0$mean/full$meanFULL

all = rbind(age1, age0)



#### FIG 2 ####

png(
  "figures/EnergyPhenology.png",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "serif"
)


for (l in c("Shetland", "ECG", "FoF", "DB")) {
  

  
  #### column 1 ####
  
  locsub = df[df$loc == l,]
  
  ## calculate mean values for each location and year ##
  
  age1 = aggregate(locsub$total_energy[locsub$doy >= start1 &
                                         locsub$doy <= end1], 
                   list(locsub$year[locsub$doy  >= start1 &
                                      locsub$doy <= end1]), 
                   mean)
  names(age1) = c("year", "mean")
  
  age0 = aggregate(locsub$total_energy[locsub$doy >= start0 &
                                         locsub$doy <= end0], 
                   list(locsub$year[locsub$doy  >= start0 &
                                      locsub$doy <= end0]), 
                   mean)
  names(age0) = c("year", "mean")
  
  # kJ rather than J
  age1$mean = age1$mean / 1000
  age0$mean = age0$mean / 1000
  
  both = data.frame(rbind(cbind(age1, data.frame(season = as.numeric(
    1
  ))),  cbind(age0, data.frame(season = (
    0
  )))))
  both$season = as.factor(both$season)
  
  
  ## basic scatter plot ##
  figE = ggplot(data = both, aes(x=year, y=mean, shape=season, color=season)) +
    
    xlim(1958, 2018) +
    ylim(0, 7.5) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_point(size = 1) +
    scale_shape_manual(values = c(16,17)) +
    scale_color_manual(values=c(col0, col1))+
    
    labs(x = "Year", y = expression(Total ~ energy ~ (kJ ~ m ^ -3))) +
    
    theme_bw(base_size = 10) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="serif"))
  
  
  
  ## fit gams ##
  y =  both$mean
  x1 = as.numeric(both$year)
  x2 = as.factor(both$season)
  
  gam.mod = gam(y ~ x2 + s(x1, by = x2))
  
  # save info on p-values for table
  gamsum = summary(gam.mod)
  gamRES[gamRES$resp == "energy_tot" & gamRES$expl == "age", l] = gamsum$p.pv[2]
  gamRES[gamRES$resp == "energy_tot" & gamRES$expl == "smooth0", l] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 0,
      "p-value"]
  gamRES[gamRES$resp == "energy_tot" & gamRES$expl == "smooth1", l] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 1,
      "p-value"]
  
  
  # check gam mods
  gam.mod0 = gam(y[x2 == 0] ~ s(x1[x2 == 0]))
  summary(gam.mod0)
  #par(mfrow = c(2,2)); gam.check(gam.mod0)
  
  gam.mod1 = gam(y[x2 == 1] ~ s(x1[x2 == 1]))
  summary(gam.mod1)
  #par(mfrow = c(2,2)); gam.check(gam.mod1)
  
  
  # if p < 0.05, add line with prediction intervals
  if (summary(gam.mod)$s.pv[1] < 0.05) {
    p = predict(
      gam.mod,
      data.frame(x1 = min(x1):max(x1), x2 = as.factor(0)),
      type = "link",
      se.fit = TRUE
      
    )
    p$year = min(x1):max(x1)
    p = as.data.frame(p)
    p$season = factor(0, levels =c(0,1))
    p$ymin = p$fit - 1.96 * p$se.fit
    p$ymax = p$fit + 1.96 * p$se.fit
    
    both_orig = both
    both = merge(both, p, by = c("year", "season"))
    
    both = rbind(both, # this is done to get a nice legend...
                 data.frame(
                   year = 2000,
                   season = factor(1),
                   mean = NA, fit = NA, se.fit = NA, ymin = NA, ymax = NA
                 ))
    
    figE = figE +
      geom_line(data = both, aes(x = year, y = fit), size = 0.8, col = col0) +
      geom_ribbon(data = both, aes(x = year, ymin = ymin, ymax = ymax, col = NULL, fill = season),  alpha = .15) + 
      scale_fill_manual(values=c(col0, col1))
  }
  
  
  if (summary(gam.mod)$s.pv[2] < 0.05) {
    p = predict(
      gam.mod,
      data.frame(x1 = min(x1):max(x1), x2 = as.factor(1)),
      type = "link",
      se.fit = TRUE
    )
    p$year = min(x1):max(x1)
    p = as.data.frame(p)
    p$season = factor(1, levels = c(0,1))
    p$ymin = p$fit - 1.96 * p$se.fit
    p$ymax = p$fit + 1.96 * p$se.fit
    
    both = merge(both_orig, p, by = c("year", "season"))
   
    both = rbind(both, # this is done to get a nice legend...
          data.frame(
            year = 2000,
            season = factor(0),
            mean = NA, fit = NA, se.fit = NA, ymin = NA, ymax = NA
          ))
    
    figE = figE +
      geom_line(data = both, aes(x = year, y = fit), size = 0.8, col = col1) +
      geom_ribbon(data = both, aes(x = year, ymin = ymin, ymax = ymax, col = NULL, fill = season),  alpha = .15) + 
      scale_fill_manual(values=c(col0, col1))
    
  }
  
  # save each plot
  if(l == "Shetland") ShetA = figE + theme(legend.position = "none")
  if(l == "ECG") ECGA = figE + theme(legend.position = "none")
  if(l == "FoF") FoFA = figE + theme(legend.position = "none")
  if(l == "DB"){ 
    DBA = figE  + # common legend
      scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) + 
      scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
      guides(colour = guide_legend(override.aes = list(size = 1.5))) + 
      scale_fill_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) +
      theme(legend.position = "bottom")}
  
  
  
  
  #### column 2 ####
  
  # subset df with energy ratios
  allSUB = all[all$loc == l,]
  allSUB$season =factor(allSUB$season, levels = c(0,1))
  
  # make basic scatter plot
  figRAT = ggplot(data = allSUB, aes(x=year, y=ratio, shape=season, color=season)) +
    
    xlim(1958, 2018) +
    ylim(min(all$ratio), max(all$ratio)) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_point(size = 1) +
    scale_shape_manual(values = c(16,17)) +
    scale_color_manual(values=c(col0, col1))+
    
    labs(x = "Year", y = "Inside/outside season") +
    
    theme_bw(base_size = 10) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="serif"))
  
  
  
  # fit gams
  gam.mod = gam(ratio ~ season + s(year, by = season), data = allSUB)
  summary(gam.mod)
  
  # save p-values
  gamsum = summary(gam.mod)
  gamRES[gamRES$resp == "energyInsideOutside" & gamRES$expl == "age", l] = gamsum$p.pv[2]
  gamRES[gamRES$resp == "energyInsideOutside" & gamRES$expl == "smooth0", l] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 0,
      "p-value"]
  gamRES[gamRES$resp == "energyInsideOutside" & gamRES$expl == "smooth1", l] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 1,
      "p-value"]
  
  
  # check gams
  gam.mod0 = gam(ratio ~ s(year), data = allSUB)
  summary(gam.mod0)
  #par(mfrow = c(2,2)); gam.check(gam.mod0)
  
  gam.mod1 = gam(ratio ~ s(year), data = allSUB)
  summary(gam.mod1)
  #par(mfrow = c(2,2)); gam.check(gam.mod1)
  
  
  # if p < 0.05, add line with prediction intervals
  if (summary(gam.mod)$s.pv[1] < 0.05) {
    p = predict(gam.mod,
                data.frame(year = min(allSUB$year):max(allSUB$year), season = 0),
                type = "link",
                se.fit = TRUE)
    
    p$year = min(allSUB$year):max(allSUB$year)
    p = as.data.frame(p)
    p$season = factor(0, levels = c(0,1))
    p$ymin = p$fit - 1.96 * p$se.fit
    p$ymax = p$fit + 1.96 * p$se.fit
    
    allSUB_orig = allSUB
    allSUB = merge(allSUB, p, by = c("year", "season"))
    
    allSUB = rbind(allSUB, # this is done to get a nice legend...
                 data.frame(
                   year = 2000,
                   season = factor(1),
                   loc = l,
                   mean = NA, ratio = NA, fit = NA, se.fit = NA, ymin = NA, ymax = NA
                 ))
    
    figRAT = figRAT +
      geom_line(data = allSUB, aes(x = year, y = fit), size = 0.8, col = col0) +
      geom_ribbon(data = allSUB, aes(x = year, ymin = ymin, ymax = ymax, col = NULL, fill = season),  alpha = .15) + 
      scale_fill_manual(values=c(col0, col1))
  }
  
  
  
  if (summary(gam.mod)$s.pv[2] < 0.05) {
    p = predict(gam.mod,
                data.frame(year = min(allSUB$year):max(allSUB$year), season = 1),
                type = "link",
                se.fit = TRUE)
    
    
    p$year = min(allSUB$year):max(allSUB$year)
    p = as.data.frame(p)
    p$season = factor(1, levels = c(0,1))
    p$ymin = p$fit - 1.96 * p$se.fit
    p$ymax = p$fit + 1.96 * p$se.fit
    
    allSUB = merge(allSUB_orig, p, by = c("year", "season"))
    
    allSUB = rbind(allSUB, # this is done to get a nice legend...
                   data.frame(
                     year = 2000,
                     season = factor(0),
                     loc = l,
                     mean = NA, ratio = NA, fit = NA, se.fit = NA, ymin = NA, ymax = NA
                   ))
    
    figRAT = figRAT +
      geom_line(data = allSUB, aes(x = year, y = fit), size = 0.8, col = col1) +
      geom_ribbon(data = allSUB, aes(x = year, ymin = ymin, ymax = ymax, col = NULL, fill = season), alpha = .15) + 
      scale_fill_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+"))
  }
  
  # save plots
  if(l == "Shetland") ShetB = figRAT + theme(legend.position = "none")
  if(l == "ECG") ECGB = figRAT + theme(legend.position = "none")
  if(l == "FoF") FoFB = figRAT + theme(legend.position = "none")
  if(l == "DB"){ 
    DBB = figRAT  + 
      scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) + 
      scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) + 
      scale_fill_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) +
      guides(colour = guide_legend(override.aes = list(size = 1.5))) +
      theme(legend.position = "bottom")}
  
  
  
  
  
  #### column 3 ####
  
  # subset to location
  sub = df_full[df_full$location == l, ]
  sub = sub[order(sub$doy),]
  
  
  # standardise within years
  max_values = aggregate(sub$total_energy, list(sub$year), max)
  names(max_values) = c("year", "max")
  sub = merge(sub, max_values, by = c("year"))
  sub$rel_energy = sub$total_energy/sub$max
  
  
  # make heat map of peak energy availability within each year
  figHM = ggplot(sub, aes(x = year, y = doy, fill = rel_energy)) + 
    
    xlim(1958, 2018) +
    ylim(0, 366) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_tile() + 
    scale_fill_viridis_c(name = "Available\nenergy", limits = c(0,1), breaks = c(0, 0.5, 1)) +
    
    geom_hline(yintercept = c(start1, end1, start0, end0), linetype = 1, color = c(col1, col1, col0,col0), size = 1) +
    
    labs(x = "Year", y = "Day of year") +
    
    theme_bw(base_size = 10) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          text=element_text(family="serif"))
  
  
  # save plots
  if(l == "Shetland") ShetC = figHM + theme(legend.position = "none")
  if(l == "ECG") ECGC = figHM + theme(legend.position = "none")
  if(l == "FoF") FoFC = figHM + theme(legend.position = "none")
  if(l == "DB"){ 
    DBC = figHM  +  
      theme(legend.position = "bottom") +
      guides(fill = guide_colourbar(barwidth = 4, barheight = 0.35))}
  
  
  
}

# arrange plots
print(
  ggarrange(ShetA, ShetB, ShetC, 
            ECGA, ECGB, ECGC,
            FoFA, FoFB, FoFC,
            DBA, DBB,DBC,
            ncol = 3,
            nrow = 4,
            labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l."),
            font.label = list(size = 15, family = "serif", face = "plain"),
            label.x = c(
              0.085,
              0.095,
              0.11,
              0.085,
              0.095,
              0.13,
              0.085,
              0.093,
              0.13,
              0.12,
              0.093,
              0.13
              
              
            ),
            common.legend = FALSE,
            heights = c(1,1,1,1.3))
)

dev.off()

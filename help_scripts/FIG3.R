#### FIG 3 ####


png(
  "figures/IA.png",
  width = 8.5,
  height = 24,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "serif"
)


for (loc in c("Shetland", "ECG", "FoF", "DB")) {
  
  # calculate mean median image area per year and location
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
  
  # make scatter plot
  fig = ggplot(data = both, aes(x=year, y=mean, shape=season, color=season)) +
    
    xlim(1958, 2018) +
    ylim(0.27, 0.82) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_point(size = 1.3) +
    scale_shape_manual(values = c(16,17)) +
    scale_color_manual(values=c(col0, col1))+
    
    labs(x = "Year", y = "Median prey size (mm)") +
    
    scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+")) + 
    scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
    guides(colour = guide_legend(override.aes = list(size = 1.7, fill = c(col0, col1)))) +
    
    theme_bw(base_size = 12) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(family="serif"))
  
  
  # fit gams
  gam.mod = gam(mean ~ season + s(year, by = season), data = both)
  summary(gam.mod)
  
  
  # save gam p-values
  gamsum = summary(gam.mod)
  gamRES[gamRES$resp == "size" & gamRES$expl == "age", loc] = gamsum$p.pv[2]
  gamRES[gamRES$resp == "size" & gamRES$expl == "smooth0", loc] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 0,
      "p-value"]
  gamRES[gamRES$resp == "size" & gamRES$expl == "smooth1", loc] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 1,
      "p-value"]
  
  # check gams
  gam.mod0 = gam(y[x2 == 0] ~ s(x1[x2 == 0]), method = "REML")
  summary(gam.mod0)
  #par(mfrow = c(2,2)); gam.check(gam.mod0)
  
  gam.mod1 = gam(y[x2 == 1] ~ s(x1[x2 == 1]), method = "REML")
  summary(gam.mod1)
  #par(mfrow = c(2,2)); gam.check(gam.mod1)
  
  
  # if p < 0.05, add line with prediction intervals
  if (summary(gam.mod)$s.pv[1] < 0.05) {
    p = predict(
      gam.mod,
      data.frame(year = min(both$year):max(both$year), season = as.factor(0)),
      type = "link",
      se.fit = TRUE
    )
    
    p$year = min(both$year):max(both$year)
    p = as.data.frame(p)
    p$season = as.factor(0)
    p$ymin = p$fit - 1.96 * p$se.fit
    p$ymax = p$fit + 1.96 * p$se.fit
    
    both_orig = both
    both = merge(both, p, by = c("year", "season"))
    
    fig = fig +
      geom_line(data = both, aes(x = year, y = fit), size = 1, col = col0) +
      geom_ribbon(data = both, aes(x = year, ymin = ymin, ymax = ymax, col = NULL), fill = col0, alpha = .15)  +
      scale_fill_manual(values = c(col0, col1), name = "Age group", labels = c("0", "1+")) 
    
    
    
  }
  
  if (summary(gam.mod)$s.pv[2] < 0.05) {
    p = predict(
      gam.mod,
      data.frame(year = min(both$year):max(both$year), season = as.factor(1)),
      type = "link",
      se.fit = TRUE
    )
    
    p$year = min(both$year):max(both$year)
    p = as.data.frame(p)
    p$season = as.factor(1)
    p$ymin = p$fit - 1.96 * p$se.fit
    p$ymax = p$fit + 1.96 * p$se.fit
    
    both = merge(both_orig, p, by = c("year", "season"))
    
    fig = fig +
      geom_line(data = both, aes(x = year, y = fit), size = 1, col = col1) +
      geom_ribbon(data = both, aes(x = year, ymin = ymin, ymax = ymax, col = NULL), fill = col1, alpha = .15) +
      scale_fill_manual(values = c(col0, col1), name = "Age group", labels = c("0", "1+")) 
    
    
    
  }
  
  # save each plot
  if(loc == "Shetland") Shet_f = fig + geom_text(x = 1975, y = 0.72, label = "Increasing abundances\nof large prey", col = "black")
  if(loc == "ECG") ECG_f = fig
  if(loc == "FoF") FoF_f = fig + geom_text(x = 1975, y = 0.72, label = "Decreasing abundances\nof small prey", col = "black")
  if(loc == "DB") DB_f = fig + geom_text(x = 1975, y = 0.72, label = "Decreasing abundances\nof small prey", col = "black")
  
}

# arrange plots
print(
  ggarrange(Shet_f, ECG_f, FoF_f, DB_f,
            ncol = 1,
            nrow = 4,
            labels = c("a.", "b.", "c.", "d."),
            font.label = list(size = 18, family = "serif", face = "plain"),
            label.x = 0.88,
            label.y = 0.98,
            legend = "bottom",
            common.legend = TRUE)
)


dev.off()
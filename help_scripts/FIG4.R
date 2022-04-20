#### FIG 4 ####




for (loc in c("Shetland", "ECG", "FoF", "DB")) {
  
  #### CAL FIN ####
  
  # calculate means for each location
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
  both$mean[both$mean == 0] = 2e-16
  
  # create scatter plot
  fig = ggplot(data = both, aes(x=year, y=mean, shape=season, color=season)) +
    
    ylim (0, 1000) +
    xlim (1958, 2018) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_point(size = 1.1) +
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
          text=element_text(family="sans"))
  
  
  # fit gam mod
  gam.mod = gam(mean ~ season + s(year, by = season), data = both, method = "REML",  family = Gamma(link = "log"))
  summary(gam.mod)
  
  
  # save p-values
  gamsum = summary(gam.mod)
  gamRES[gamRES$resp == "cf" & gamRES$expl == "age", loc] = gamsum$p.pv[2]
  gamRES[gamRES$resp == "cf" & gamRES$expl == "smooth0", loc] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 0,
      "p-value"]
  gamRES[gamRES$resp == "cf" & gamRES$expl == "smooth1", loc] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 1,
      "p-value"]
  
  
  # qqplot 
  qq_data = data.frame(
    theo_quant = qq.gam(gam.mod),
    resid =  sort(residuals.gam(gam.mod, type = "deviance"))
  )
  
  qq = ggplot(data = qq_data, aes(x = theo_quant, y = resid)) +
    geom_segment(aes(x = min(qq_data), y = min(qq_data), xend = max(qq_data), yend = max(qq_data)), colour = "red") +
    geom_point() +
    labs(x = "Theoretical quantiles", y = "Deviance residuals", title = "QQ plot of residuals", subtitle = "Method: uniform") +
    theme_bw(base_size = 7) + theme(text=element_text(family="sans"))
  
  
  # linear predictor vs residuals
  res_linpred = residuals_linpred_plot(gam.mod) + theme_bw(base_size = 7) + theme(text=element_text(family="sans"))
  
  # year vs residuals
  res_predictor = ggplot(data = both) +
    geom_point(aes(x = year, y = residuals(gam.mod), group = season, colour = season, shape = season), size = 0.8) +
    geom_smooth(aes(x = year, y = residuals(gam.mod), group = season, colour = season, fill = season), method = "loess") +
    
    scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
    scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+"))+
    scale_fill_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+"))+
    
    labs(x = "Year", y = "Deviance residuals", title = "Residuals vs year", subtitle = "") +
    
    theme_bw(base_size = 7) + theme(text=element_text(family="sans"))
  
  
  
  
  p = predict(
    gam.mod,
    data.frame(year = min(both$year):max(both$year), season = as.factor(0)),
    type = "link",
    se.fit = TRUE
  )
  
  p$year = min(both$year):max(both$year)
  p = as.data.frame(p)
  p$season = as.factor(0)
  p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
  p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
  p$fit = gam.mod$family$linkinv(p$fit)
  
  
  both_orig = both
  both = merge(both, p, by = c("year", "season"), all.y = T)
  
  fig = fig +
    geom_line(data = both, aes(x = year, y = fit), size = 0.8, col = col0) +
    geom_ribbon(data = both, aes(x = year, ymin = ymin, ymax = ymax, col = NULL), fill = col0, alpha = .15)  +
    scale_fill_manual(values = c(col0, col1), name = "Age group", labels = c("0", "1+")) 
  
  
  
  
  
  
  p = predict(
    gam.mod,
    data.frame(year = min(both$year):max(both$year), season = as.factor(1)),
    type = "link",
    se.fit = TRUE
  )
  
  p$year = min(both$year):max(both$year)
  p = as.data.frame(p)
  p$season = as.factor(1)
  p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
  p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
  p$fit = gam.mod$family$linkinv(p$fit)
  
  
  both = merge(both_orig, p, by = c("year", "season"), all.y = TRUE)
  
  fig = fig +
    geom_line(data = both, aes(x = year, y = fit), size = 0.8, col = col1) +
    geom_ribbon(data = both, aes(x = year, ymin = ymin, ymax = ymax, col = NULL), fill = col1, alpha = .15) +
    scale_fill_manual(values = c(col0, col1), name = "Age group", labels = c("0", "1+")) 
  
  
  
  # save plots
  if(loc == "Shetland"){ Shet_cf = fig 
  Shet_diag_cf = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("a.", "b.", "c."), font.label = list(size = 15, family = "sans", face = "plain"))}
  if(loc == "ECG"){ ECG_cf = fig
  ECG_diag_cf = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("d.", "e.", "f."), font.label = list(size = 15, family = "sans", face = "plain"))}
  if(loc == "FoF"){ FoF_cf = fig 
  FoF_diag_cf = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("g.", "h.", "i."), font.label = list(size = 15, family = "sans", face = "plain"))}
  if(loc == "DB"){ DB_cf = fig 
  DB_diag_cf = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("j.", "k.", "l."), font.label = list(size = 15, family = "sans", face = "plain"))}
  
  
  
  
  
  
  
  #### CAL HEL ####
  
  
  # calculate means for each year and location
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
  both$mean[both$mean == 0] = 2e-16
  
  
  # make scatter plot
  fig = ggplot(data = both, aes(x=year, y=mean, shape=season, color=season)) +
    
    ylim(0, 670) + 
    xlim(1958, 2018) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_point(size = 1.1) +
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
          text=element_text(family="sans"))
  
  
  # fit gams
  gam.mod = gam(mean ~ season + s(year, by = season), data = both, method = "REML", family = Gamma(link = "log"))
  summary(gam.mod)
  
  
  # save p-values
  gamsum = summary(gam.mod)
  gamRES[gamRES$resp == "ch" & gamRES$expl == "age", loc] = gamsum$p.pv[2]
  gamRES[gamRES$resp == "ch" & gamRES$expl == "smooth0", loc] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 0,
      "p-value"]
  gamRES[gamRES$resp == "ch" & gamRES$expl == "smooth1", loc] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 1,
      "p-value"]
  
  
  # qqplot 
  qq_data = data.frame(
    theo_quant = qq.gam(gam.mod),
    resid =  sort(residuals.gam(gam.mod, type = "deviance"))
  )
  
  qq = ggplot(data = qq_data, aes(x = theo_quant, y = resid)) +
    geom_segment(aes(x = min(qq_data), y = min(qq_data), xend = max(qq_data), yend = max(qq_data)), colour = "red") +
    geom_point() +
    labs(x = "Theoretical quantiles", y = "Deviance residuals", title = "QQ plot of residuals", subtitle = "Method: uniform") +
    theme_bw(base_size = 7) + theme(text=element_text(family="sans"))
  
  # linear predictor vs residuals
  res_linpred = residuals_linpred_plot(gam.mod) + theme_bw(base_size = 7) + theme(text=element_text(family="sans"))
  
  # year vs residuals
  res_predictor = ggplot(data = both) +
    geom_point(aes(x = year, y = residuals(gam.mod), group = season, colour = season, shape = season), size = 0.8) +
    geom_smooth(aes(x = year, y = residuals(gam.mod), group = season, colour = season, fill = season), method = "loess") +
    
    scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
    scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+"))+
    scale_fill_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+"))+
    
    labs(x = "Year", y = "Deviance residuals", title = "Residuals vs year", subtitle = "") +
    
    theme_bw(base_size = 7) + theme(text=element_text(family="sans"))
  
  
  
  
  p = predict(
    gam.mod,
    data.frame(year = min(both$year):max(both$year), season = as.factor(0)),
    type = "link",
    se.fit = TRUE
  )
  
  p$year = min(both$year):max(both$year)
  p = as.data.frame(p)
  p$season = as.factor(0)
  p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
  p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
  p$fit = gam.mod$family$linkinv(p$fit)
  
  
  both_orig = both
  both = merge(both, p, by = c("year", "season"), all.y = T)
  
  fig = fig +
    geom_line(data = both, aes(x = year, y = fit), size = 0.8, col = col0) +
    geom_ribbon(data = both, aes(x = year, ymin = ymin, ymax = ymax, col = NULL), fill = col0, alpha = .15)  +
    scale_fill_manual(values = c(col0, col1), name = "Age group", labels = c("0", "1+")) 
  
  
  
  
  p = predict(
    gam.mod,
    data.frame(year = min(both$year):max(both$year), season = as.factor(1)),
    type = "link",
    se.fit = TRUE
  )
  
  p$year = min(both$year):max(both$year)
  p = as.data.frame(p)
  p$season = as.factor(1)
  p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
  p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
  p$fit = gam.mod$family$linkinv(p$fit)
  
  
  both = merge(both_orig, p, by = c("year", "season"), all.y = TRUE)
  
  fig = fig +
    geom_line(data = both, aes(x = year, y = fit), size = 0.8, col = col1) +
    geom_ribbon(data = both, aes(x = year, ymin = ymin, ymax = ymax, col = NULL), fill = col1, alpha = .15) +
    scale_fill_manual(values = c(col0, col1), name = "Age group", labels = c("0", "1+")) 
  
  
  
  
  
  
  # save plots
  if(loc == "Shetland"){ Shet_ch = fig 
  Shet_diag_ch = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("a.", "b.", "c."), font.label = list(size = 15, family = "sans", face = "plain"))}
  if(loc == "ECG"){ ECG_ch = fig
  ECG_diag_ch = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("d.", "e.", "f."), font.label = list(size = 15, family = "sans", face = "plain"))}
  if(loc == "FoF"){ FoF_ch = fig 
  FoF_diag_ch = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("g.", "h.", "i."), font.label = list(size = 15, family = "sans", face = "plain"))}
  if(loc == "DB"){ DB_ch = fig 
  DB_diag_ch = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("j.", "k.", "l."), font.label = list(size = 15, family = "sans", face = "plain"))}
  
  
  
  
  #### SMALL COPS ####
  
  
  # calculate means for each year and location
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
  
  both$mean[both$mean == 0] = 2e-16
  
  
  # make scatter plot
  fig = ggplot(data = both, aes(x=year, y=mean, shape=season, color=season)) +
    
    ylim (0, 30) + 
    xlim(1958, 2018) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    
    geom_point(size = 1.1) +
    scale_shape_manual(values = c(16,17)) +
    scale_color_manual(values=c(col0, col1))+
    
    labs(x = "Year", y = expression(paste("Small copepods ("  %*% "1000 m" ^ "-3", ")"))) +
    
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
          text=element_text(family="sans"))
  
  
  # fit gam mod
  gam.mod = gam(mean ~ season + s(year, by = season), data = both, method = "REML",  family = Gamma(link = "log"))
  summary(gam.mod)
  
  
  # save p-values
  gamsum = summary(gam.mod)
  gamRES[gamRES$resp == "small_cops" & gamRES$expl == "age", loc] = gamsum$p.pv[2]
  gamRES[gamRES$resp == "small_cops" & gamRES$expl == "smooth0", loc] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 0,
      "p-value"]
  gamRES[gamRES$resp == "small_cops" & gamRES$expl == "smooth1", loc] = 
    gamsum$s.table[
      substr(rownames(gamsum$s.table), nchar(rownames(gamsum$s.table)),  nchar(rownames(gamsum$s.table))) == 1,
      "p-value"]
  
  
  # qqplot 
  qq_data = data.frame(
    theo_quant = qq.gam(gam.mod),
    resid =  sort(residuals.gam(gam.mod, type = "deviance"))
  )
  
  qq = ggplot(data = qq_data, aes(x = theo_quant, y = resid)) +
    geom_segment(aes(x = min(qq_data), y = min(qq_data), xend = max(qq_data), yend = max(qq_data)), colour = "red") +
    geom_point() +
    labs(x = "Theoretical quantiles", y = "Deviance residuals", title = "QQ plot of residuals", subtitle = "Method: uniform") +
    theme_bw(base_size = 7) + theme(text=element_text(family="sans"))
  
  # linear predictor vs residuals
  res_linpred = residuals_linpred_plot(gam.mod) + theme_bw(base_size = 7) + theme(text=element_text(family="sans"))
  
  # year vs residuals
  res_predictor = ggplot(data = both) +
    geom_point(aes(x = year, y = residuals(gam.mod), group = season, colour = season, shape = season), size = 0.8) +
    geom_smooth(aes(x = year, y = residuals(gam.mod), group = season, colour = season, fill = season), method = "loess") +
    
    scale_shape_manual(values = c(16,17), name = "Age group", labels = c("0", "1+")) +
    scale_color_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+"))+
    scale_fill_manual(values=c(col0, col1), name = "Age group", labels = c("0", "1+"))+
    
    labs(x = "Year", y = "Deviance residuals", title = "Residuals vs year", subtitle = "") +
    
    theme_bw(base_size = 7) + theme(text=element_text(family="sans"))
  
  
  
  
  
  p = predict(
    gam.mod,
    data.frame(year = min(both$year):max(both$year), season = as.factor(0)),
    type = "link",
    se.fit = TRUE
  )
  
  p$year = min(both$year):max(both$year)
  p = as.data.frame(p)
  p$season = as.factor(0)
  p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
  p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
  p$fit = gam.mod$family$linkinv(p$fit)
  
  
  both_orig = both
  both = merge(both, p, by = c("year", "season"), all.y = T)
  
  fig = fig +
    geom_line(data = both, aes(x = year, y = fit), size = 0.8, col = col0) +
    geom_ribbon(data = both, aes(x = year, ymin = ymin, ymax = ymax, col = NULL), fill = col0, alpha = .15)  +
    scale_fill_manual(values = c(col0, col1), name = "Age group", labels = c("0", "1+")) 
  
  
  
  p = predict(
    gam.mod,
    data.frame(year = min(both$year):max(both$year), season = as.factor(1)),
    type = "link",
    se.fit = TRUE
  )
  
  p$year = min(both$year):max(both$year)
  p = as.data.frame(p)
  p$season = as.factor(1)
  p$ymin = gam.mod$family$linkinv(p$fit - 1.96 * p$se.fit)
  p$ymax = gam.mod$family$linkinv(p$fit + 1.96 * p$se.fit)
  p$fit = gam.mod$family$linkinv(p$fit)
  
  
  both = merge(both_orig, p, by = c("year", "season"), all.y = TRUE)
  
  fig = fig +
    geom_line(data = both, aes(x = year, y = fit), size = 0.8, col = col1) +
    geom_ribbon(data = both, aes(x = year, ymin = ymin, ymax = ymax, col = NULL), fill = col1, alpha = .15) +
    scale_fill_manual(values = c(col0, col1), name = "Age group", labels = c("0", "1+")) 
  
  
  
  # save plots
  if(loc == "Shetland"){ Shet_sc = fig 
  Shet_diag_sc = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("a.", "b.", "c."), font.label = list(size = 15, family = "sans", face = "plain"))}
  if(loc == "ECG"){ ECG_sc = fig
  ECG_diag_sc = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("d.", "e.", "f."), font.label = list(size = 15, family = "sans", face = "plain"))}
  if(loc == "FoF"){ FoF_sc = fig 
  FoF_diag_sc = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("g.", "h.", "i."), font.label = list(size = 15, family = "sans", face = "plain"))}
  if(loc == "DB"){ DB_sc = fig 
  DB_diag_sc = ggarrange(qq,  res_linpred, res_predictor, ncol = 3, nrow = 1, widths = c(1,1,1.26), labels = c("j.", "k.", "l."), font.label = list(size = 15, family = "sans", face = "plain"))}
  
  
  
}



jpeg(
  "figures/taxTRENDS.jpeg",
  width = 17,
  height = (17/15)*20,
  units = 'cm',
  res = 600,
  pointsize = 9,
  family = "sans"
)



# arrange plots
print(
  ggarrange(Shet_cf, Shet_ch, Shet_sc, 
            ECG_cf, ECG_ch, ECG_sc,
            FoF_cf, FoF_ch, FoF_sc,
            DB_cf, DB_ch, DB_sc,
            ncol = 3,
            nrow = 4,
            common.legend = TRUE,
            legend = "bottom",
            labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l."),
            font.label = list(size = 12, family = "sans", face = "plain"),
            label.x = 0.84,
            label.y = 0.97)
)



dev.off()




# open plot
jpeg(
  "figures/DiagPlots_cf.jpeg",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)


ggarrange(Shet_diag_cf, ECG_diag_cf, FoF_diag_cf, DB_diag_cf, ncol = 1, nrow = 4)

dev.off()


# open plot
jpeg(
  "figures/DiagPlots_ch.jpeg",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)


ggarrange(Shet_diag_ch, ECG_diag_ch, FoF_diag_ch, DB_diag_ch, ncol = 1, nrow = 4)

dev.off()


# open plot
jpeg(
  "figures/DiagPlots_sc.jpeg",
  width = 15,
  height = 20,
  units = 'cm',
  res = 200,
  pointsize = 9,
  family = "sans"
)


ggarrange(Shet_diag_sc, ECG_diag_sc, FoF_diag_sc, DB_diag_sc, ncol = 1, nrow = 4)

dev.off()

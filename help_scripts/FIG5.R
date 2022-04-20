#### FIG 5 ####

plot_locs = c("DB", "FoF", "ECG", "Shetland")

# group taxa together 
groups = list(
  c(
    "Copepod.nauplii",
    "Cirripede.larvae..Total.",
    "Decapoda.larvae..Total."
  ),
  c("Appendicularia"),
  c("Fish.eggs..Total."),
  c(
    "Acartia.spp...unidentified.",
    "Para.Pseudocalanus.spp.",
    "Temora.longicornis",
    "Oithona.spp."
  ),
  c(
    "Centropages.hamatus",
    "Centropages.typicus",
    "Centropages.spp...Unidentified."
  ),
  c("Calanus.I.IV"),
  c("Metridia.lucens"),
  c("Calanus.helgolandicus"),
  c("Calanus.finmarchicus"),
  c("Fish.larvae"),
  c("Hyperiidea..Total.")
)

# names for groups 
names(groups) = c(
  "Crustacean larvae",
  "Appendicularia",
  "Fish eggs",
  "Small copepods",
  "Centropages spp.",
  "Calanus I-IV",
  "Metridia lucens",
  "Calanus helgolandicus",
  "Calanus finmarchicus",
  "Fish larvae",
  "Hyperiids"
)



#### standardise values ####


## abundance (left-hand column) ##

# create empty dataframe for standardised abundances
coll = data.frame(
  loc = character(),
  year = numeric(),
  taxon = character(),
  abundance = numeric()
)


for (l in plot_locs) {
  for (t in 1:length(groups)) {
    x = aggregate(
      apply(as.data.frame(df[df$doy >= start1 & # sum values for each group for each day within location and over total feeding season
                               df$doy <= end0 &
                               df$location == l, groups[[t]]]), 1, sum),
      
      list(df$year[df$doy >= start1 & # aggregate per year
                     df$doy <= end0 &
                     df$location == l]),
      mean) # calculate mean
    names(x) = c("year", "abundance")
    
    
    x = merge(x, data.frame(year = 1958:2018), all = T) # to make sure NAs are included for years with no data
    x = data.frame(
      loc = l,
      year = x$year,
      group = names(groups)[t],
      abundance = x$abundance
    )
    
    coll = rbind(coll, x)
    
  }
  
}


# recalculate so each is divided by maximum abundance for that taxon
coll$abundance = log10(as.numeric(coll$abundance) + 1)
coll$abu_prop = NA

for (t in names(groups)) {
  coll[coll$group == t, "abu_prop"] =
    coll[coll$group == t, "abundance"] /
    max(coll[coll$group == t, "abundance"], na.rm = T)
  
}





## energy ##

# create empty data frame 
coll_e = data.frame(
  loc = character(),
  year = numeric(),
  taxon = character(),
  energy = numeric()
)

# calculate mean energy per year for each group (see above)
for (l in c(plot_locs, "Faroes", "Iceland")) {
  for (t in 1:length(groups)) {
    x = aggregate(rowSums(t(
      t(as.data.frame(df[df$doy >= start1 &
                           df$doy <= end0 &
                           df$location == l, groups[[t]]])) *
        prey_info$energy[prey_info$taxa %in% groups[[t]]]
    ))
    
    ,
    list(df$year[df$doy >= start1 &
                   df$doy <= end0 &
                   df$location == l]),
    mean)
    names(x) = c("year", "energy")
    x = merge(x, data.frame(year = 1958:2018), all = T)
    x = data.frame(
      loc = l,
      year = x$year,
      group = names(groups)[t],
      energy = x$energy
    )
    
    
    coll_e = rbind(coll_e, x)
    
  }
  
}


# recalculate so each is divided by total energy for year and location

coll_e$energy_prop = NA
for (l in c(plot_locs, "Faroes", "Iceland")) {
  for (t in 1958:2018) {
    coll_e[coll_e$year == t & coll_e$loc == l, "energy_prop"] =
      coll_e[coll_e$year == t & coll_e$loc == l, "energy"] /
      sum(coll_e[coll_e$year == t &
                   coll_e$loc == l, "energy"], na.rm = T)
    
  }
}


# calculate averages over all years for each location for barplot
res = expand.grid(group = unique(coll_e$group),
                  loc = plot_locs)
res$prop = NA

for (l in plot_locs) {
  for (g in unique(coll_e$group)) {
    res$prop[res$group == g & res$loc == l] =
      sum(coll_e$energy[coll_e$group == g &
                          coll_e$loc == l], na.rm = T) / sum(coll_e$energy[coll_e$loc == l], na.rm = T)
  }
}


# labels for groups
plot_names = c(
  "Crustacean larvae",
  "Appendicularia",
  "Fish eggs",
  "Small copepods",
  expression(italic("Centropages") ~ "spp."),
  expression(italic("Calanus") ~ "I-IV"),
  expression(italic("Metridia lucens")),
  expression(italic("C. helgolandicus")),
  expression(italic("C. finmarchicus")),
  "Fish larvae",
  "Hyperiids"
)


#### FIG 5 ####

jpeg(
  "figures/indTAXA.jpeg",
  width = 17,
  height = (17/15)*20,
  units = 'cm',
  res = 600,
  pointsize = 9,
  family = "sans"
)


for (l in plot_locs) {
  
  #### heatmap 1 ####
  
  # subset to location
  sub = coll[coll$loc == l , c("group", "year", "abu_prop")]
  
  
  # create heatmap
  fig = ggplot(sub, aes(x = year, y = group, fill = abu_prop)) + 
    
    xlim(1958, 2018) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    scale_y_discrete(limits = rev(names(groups)), labels = rev(plot_names)) +
    
    geom_tile() + 
    scale_fill_viridis_c(name = "Relative abundance/Proportional energy contribution", limits = c(0,1), breaks = c(0, 0.5, 1)) +
    
    geom_hline(yintercept = c(start1, end1, start0, end0), linetype = 1, color = c(col1, col1, col0,col0), size = 1) +
    
    labs(x = "Year", y = "") +
    
    theme_bw(base_size = 7) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          plot.margin=unit(c(0.5,0.1,0.1,0.1),"cm"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          text=element_text(family="sans")) +
    guides(fill = guide_colourbar(barwidth = 3.5, barheight = 0.35, title.position = "left", title.hjust = 0, title.vjust = 1, title.theme = element_text(size = 7)))
  
  # save plots
  if(l == "Shetland") Shet_1 = fig 
  if(l == "ECG") ECG_1 = fig
  if(l == "FoF") FoF_1 = fig 
  if(l == "DB") DB_1 = fig 
  
  
  
  
  #### heatmap 2 ####
  
  # subset to location
  sub = coll_e[coll_e$loc == l , c("group", "year", "energy_prop")]
  
  # create heatmap
  fig = ggplot(sub, aes(x = year, y = group, fill = energy_prop)) + 
    
    xlim(1958, 2018) +
    
    scale_x_continuous(name = "Year", breaks = c(1970,1990,2010), labels = c("1970","1990","2010")) +
    scale_y_discrete(limits = rev(names(groups))) +
    
    geom_tile() + 
    scale_fill_viridis_c(name = "Relative abundances/Proportional energy contribution", limits = c(0,1), breaks = c(0, 0.5, 1)) +
    
    geom_hline(yintercept = c(start1, end1, start0, end0), linetype = 1, color = c(col1, col1, col0,col0), size = 1) +
    
    labs(x = "Year", y = "") +
    
    theme_bw(base_size = 7) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          plot.margin=unit(c(0.5,0.1,0.1,0.1),"cm"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          text=element_text(family="sans")) +
    guides(fill = guide_colourbar(barwidth = 3.5, barheight = 0.35, title.position = "left", title.hjust = 0,  title.vjust = 1, title.theme = element_text(size = 7)))
  
  # save plots
  if(l == "Shetland") Shet_2 = fig 
  if(l == "ECG") ECG_2 = fig
  if(l == "FoF") FoF_2 = fig 
  if(l == "DB") DB_2 = fig 
  
  #### barplot ####
  
  # subset to location
  sub = res[res$loc == l, ]
  
  
  # create barplot
  fig = ggplot(data=sub, aes(x = group, y = prop)) +
    geom_bar(stat="identity") + 
    coord_flip() +
    scale_x_discrete(limits = rev(levels(sub$group))) +
    ylim(c(0,0.7)) +
    labs(y = "Proportional contribution") +
    
    theme_bw(base_size = 7) +
    theme(panel.border = 
            element_rect(
              fill = NA,
              size = 1), 
          plot.margin=unit(c(0.5,0.1,0.1,0.1),"cm"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          text=element_text(family="sans"))
  
  # save plots
  if(l == "Shetland") Shet_3 = fig 
  if(l == "ECG") ECG_3 = fig
  if(l == "FoF") FoF_3 = fig 
  if(l == "DB") DB_3 = fig 
  
  
  
}

# arrange pltos
print(
  ggarrange(Shet_1, Shet_2, Shet_3, 
            ECG_1, ECG_2, ECG_3,
            FoF_1, FoF_2, FoF_3,
            DB_1, DB_2, DB_3,
            ncol = 3,
            nrow = 4,
            common.legend = TRUE,
            legend = "bottom",
            widths = c(1.55, 1,1),
            labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l."),
            font.label = list(size = 11, family = "sans", face = "plain"),
            label.x = rep(c(0.65,0.45,0.45),4),
            label.y = 1.02)
)

dev.off()













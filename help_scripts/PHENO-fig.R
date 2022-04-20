

#### PHENOLOGY PLOTS ####

# plots of average values per day across all years

## energy ##
energy_pheno = aggregate(total_energy ~ location + doy, data = df_full, mean)
energy_pheno$location = factor(energy_pheno$location, levels = c("DB", "FoF", "ECG", "Shetland"))

p1 = ggplot() +
  
  geom_rect(aes(xmin=start1, xmax=end1, ymin=-5, ymax=10000), fill = viridis(10)[5], alpha = 0.3, lwd = 0) +
  geom_rect(aes(xmin=start0, xmax=end0, ymin=-5, ymax=10000), fill = "goldenrod", alpha = 0.3, lwd = 0) +
  
  geom_line(data = energy_pheno, aes(x = doy, y = total_energy, group = location, colour = location), lwd =  0.65) +
  
  scale_colour_manual(values = brewer.pal(5, "Oranges")[2:5], name = "Location", labels = c("Dogger Bank", "Firth of Forth", "ECG", "Shetland"))  + 
  
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  
  coord_cartesian(xlim = c(50, 320), ylim = c(0, 2700)) +
  
  xlab("Day of year") +
  ylab(expression(Total ~ energy ~ (kJ ~ m ^ -3))) +
  
  theme_bw(base_size = 8) +
  theme(panel.border = 
          element_rect(
            fill = NA,
            size = 1), 
        panel.grid = element_blank(),
        legend.position = "none",
        text=element_text(family="sans"))



df_full$image_area_median = apply(
  df_full[, taxa],
  1,
  FUN = function(x)
    median(rep(
      sqrt(prey_info$image_area), as.numeric(round(t(x)))
    ))
)

## size ##
size_pheno = aggregate(image_area_median ~ location + doy, data = df_full, mean)
size_pheno$location = factor(size_pheno$location, levels = c("DB", "FoF", "ECG", "Shetland"))
length(size_pheno$image_area_median[size_pheno$location == "DB"])

size_pheno$image_area_median[size_pheno$location == "Shetland"] = 
  c(rep(NA, 3), rollmean(size_pheno$image_area_median[size_pheno$location == "Shetland"], k = 7), rep(NA, 3))


size_pheno$image_area_median[size_pheno$location == "ECG"] = 
  c(rep(NA, 3), rollmean(size_pheno$image_area_median[size_pheno$location == "ECG"], k = 7), rep(NA, 3))

size_pheno$image_area_median[size_pheno$location == "FoF"] = 
  c(rep(NA, 3), rollmean(size_pheno$image_area_median[size_pheno$location == "FoF"], k = 7), rep(NA, 3))


size_pheno$image_area_median[size_pheno$location == "DB"] = 
c(rep(NA, 3), rollmean(size_pheno$image_area_median[size_pheno$location == "DB"], k = 7), rep(NA, 3))



p2 = ggplot() +
  
  geom_rect(aes(xmin=start1, xmax=end1, ymin=-5, ymax=10000), fill = viridis(10)[5], alpha = 0.3, lwd = 0) +
  geom_rect(aes(xmin=start0, xmax=end0, ymin=-5, ymax=10000), fill = "goldenrod", alpha = 0.3, lwd = 0) +
  
  geom_line(data = size_pheno, aes(x = doy, y = image_area_median, group = location, colour = location), lwd =  0.65) +
  
  scale_colour_manual(values = brewer.pal(5, "Oranges")[2:5], name = "Location", labels = c("Dogger Bank", "Firth of Forth", "ECG", "Shetland"))  + 
  
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  
  coord_cartesian(xlim = c(50, 320), ylim = c(0.2, 0.8)) +
  
  xlab("Day of year") +
  ylab("Median prey size (mm)") +
  
  theme_bw(base_size = 8) +
  theme(panel.border = 
          element_rect(
            fill = NA,
            size = 1), 
        legend.position = "none",
        panel.grid = element_blank(),
        text=element_text(family="sans"))







## C fin ##
cf_pheno = aggregate(Calanus.finmarchicus ~ location + doy, data = df_full, mean)
cf_pheno$location = factor(cf_pheno$location, levels = c("DB", "FoF", "ECG", "Shetland"))


p3 = ggplot() +
  
  geom_rect(aes(xmin=start1, xmax=end1, ymin=-5, ymax=1000), fill = viridis(10)[5], alpha = 0.3, lwd = 0) +
  geom_rect(aes(xmin=start0, xmax=end0, ymin=-5, ymax=1000), fill = "goldenrod", alpha = 0.3, lwd = 0) +
  
  geom_line(data = cf_pheno, aes(x = doy, y = Calanus.finmarchicus, group = location, colour = location), lwd =  0.65) +
  
  scale_colour_manual(values = brewer.pal(5, "Oranges")[2:5], name = "Location", labels = c("Dogger Bank", "Firth of Forth", "ECG", "Shetland"))  + 
  
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  
  coord_cartesian(xlim = c(50, 320), ylim = c(0, 350)) +
  
  xlab("Day of year") +
  ylab( expression(italic("Calanus finmarchicus") ~ m ^ -3)) +
  
  theme_bw(base_size = 8) +
  theme(panel.border = 
          element_rect(
            fill = NA,
            size = 1), 
        panel.grid = element_blank(),
        legend.text=element_text(size=8),
        text=element_text(family="sans"))


## C hel ##
ch_pheno = aggregate(Calanus.helgolandicus ~ location + doy, data = df_full, mean)
ch_pheno$location = factor(ch_pheno$location, levels = c("DB", "FoF", "ECG", "Shetland"))


p4 = ggplot() +
  
  geom_rect(aes(xmin=start1, xmax=end1, ymin=-5, ymax=100), fill = viridis(10)[5], alpha = 0.3, lwd = 0) +
  geom_rect(aes(xmin=start0, xmax=end0, ymin=-5, ymax=100), fill = "goldenrod", alpha = 0.3, lwd = 0) +
  
  geom_line(data = ch_pheno, aes(x = doy, y = Calanus.helgolandicus, group = location, colour = location), lwd =  0.65) +
  
  scale_colour_manual(values = brewer.pal(5, "Oranges")[2:5], name = "Location", labels = c("Dogger Bank", "Firth of Forth", "ECG", "Shetland"))  + 
  
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  
  coord_cartesian(xlim = c(50, 320), ylim = c(0, 80)) +
  
  xlab("Day of year") +
  ylab( expression(italic("Calanus helgolandicus") ~ m ^ -3)) +
  
  theme_bw(base_size = 8) +
  theme(panel.border = 
          element_rect(
            fill = NA,
            size = 1), 
        panel.grid = element_blank(),
        legend.text=element_text(size=8),
        text=element_text(family="sans"))


## small cops ##
df_full$small_copepods = (df_full$Acartia.spp...unidentified. +df_full$Oithona.spp. + df_full$Para.Pseudocalanus.spp. + df_full$Temora.longicornis)/1000
small_cop_pheno = aggregate(small_copepods ~ location + doy, data = df_full, mean)
small_cop_pheno$location = factor(small_cop_pheno$location, levels = c("DB", "FoF", "ECG", "Shetland"))


p5 = ggplot() +
  
  geom_rect(aes(xmin=start1, xmax=end1, ymin=-5, ymax=11), fill = viridis(10)[5], alpha = 0.3, lwd = 0) +
  geom_rect(aes(xmin=start0, xmax=end0, ymin=-5, ymax=11), fill = "goldenrod", alpha = 0.3, lwd = 0) +
  
  geom_line(data = small_cop_pheno, aes(x = doy, y = small_copepods, group = location, colour = location), lwd =  0.65) +
  
  scale_colour_manual(values = brewer.pal(5, "Oranges")[2:5], name = "Location", labels = c("Dogger Bank", "Firth of Forth", "ECG", "Shetland"))  + 
  
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  
  coord_cartesian(xlim = c(50, 320), ylim = c(0, 8.7)) +
  
  xlab("Day of year") +
  ylab(expression(paste("Small copepods ("  %*% "1000 m" ^ "-3", ")"))) +
  
  theme_bw(base_size = 8) +
  theme(panel.border = 
          element_rect(
            fill = NA,
            size = 1), 
        panel.grid = element_blank(),
        legend.text=element_text(size=8),
        text=element_text(family="sans"))


## save fig ##

jpeg(
  "figures/phenoFIG.jpeg",
  width = 17,
  height = (12/15)*17,
  units = 'cm',
  res = 600,
  pointsize = 9,
  family = "sans"
)


p1 = ggarrange(p1, p2,  ncol = 2, nrow = 1, 
          labels = c("a.", "b."),
          font.label = list(size = 12, family = "sans", face = "plain"),
          label.x = 0.90,
          label.y = 0.97
          )

p2 = ggarrange(p3, p4, p5, ncol = 3, nrow = 1, common.legend = T, legend = "bottom",
          labels = c("c.", "d.", "e."),
          font.label = list(size = 12, family = "sans", face = "plain"),
          label.x = 0.85,
          label.y = 0.97
          )


ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1,1.3))

dev.off()


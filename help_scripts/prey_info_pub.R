# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~  prey info for CPR DATA ~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# see supplementary table for sources

#### load empty data frame ####

setwd("~/nonSU/CPR")
prey_info = read.csv("prey_info_empty.csv", sep = ";")

#### add sizes ####
prey_info$size = c(1.15,
                   1,
                   2.7,
                   2.68,
                   1.65,
                   2.48,
                   1.3,
                   1.63,
                   1.55,
                   0.68, 
                   0.19,
                   0.9, # assume shorter values more representative
                   17,
                   0.5,
                   1,
                   12,
                   16,
                   2.27,
                   0.68,
                   0.7,
                   1,
                   1)
prey_info$size = signif(prey_info$size, digits = 2)

 
#### add weights ####
prey_info$weight =c(signif(0.08*prey_info$size[1]^2.1, digits = 2),
                    signif(((10^(2.6270*log10(prey_info$size[2]*1000)-7.1348))/ ((1-0.054)/100) )/1000, digits = 2), # King + tunicate energy density
                    signif(0.08*prey_info$size[3]^2.1, digits = 2),
                    signif(0.08*prey_info$size[4]^2.1, digits = 2),
                    signif(0.08*prey_info$size[5]^2.1, digits = 2),
                    signif(0.08*prey_info$size[6]^2.1, digits = 2),
                    signif(0.08*prey_info$size[7]^2.1, digits = 2),
                    signif(0.08*prey_info$size[8]^2.1, digits = 2),
                    signif(0.08*prey_info$size[9]^2.1, digits = 2),
                    signif((0.00201)/0.0961, digits = 2), # Muxagata
                    signif((0.096*prey_info$size[11]^3.21)/0.162, digits = 2),
                    signif(((10^(2.58*log10(prey_info$size[12])+2.04))/0.183)/1000, digits = 2),
                    signif(0.012* prey_info$size[13]^2.98 , digits = 2),
                    signif((2/1000)/0.0961, digits = 2),
                    signif(0.11/0.07, digits = 2) , ## riis
                    2,
                    signif((0.0064*prey_info$size[17]^2.4614)/0.183, digits = 2)  ,
                    signif(0.08*prey_info$size[18]^2.1, digits = 2),
                    signif(0.08*prey_info$size[19]^2.1, digits = 2),
                    signif(0.08*prey_info$size[20]^2.1, digits = 2),
                    signif(((10^( 3.9*log10(1000*prey_info$size[21])-10.12))/0.183)/1000, digits = 2),
                    signif(0.08*prey_info$size[22]^2.1, digits = 2))


#### add image area ####


# assume copepods are ellipse with length l and width l/2 (van Deurs)
# assume crustacean larvae are circular simplification
# copepods ~ same shape all direction
# for adult crustaceans we assume that appendages (legs/tail) are not visible
# for krill and hyperiidea - is length equivalent
# assume 0.75 for all calanoid copepods based on Conway

prey_info$image_area = c(pi * (prey_info$size[1]*(0.75))/2 * (prey_info$size[1]*(0.75)/2)/2, 
                         pi *prey_info$size[2]/2 * (prey_info$size[2]/2)/2, # trunk length ~ 1mm (https://www.researchgate.net/figure/Growth-trunk-length-of-Oikopleura-dioica-expressed-as-a-function-of-a-age-h-post_fig3_5764742) - assumme cannot see tail
                         pi * (prey_info$size[3]*(0.75))/2 * (prey_info$size[3]*(0.75)/2)/2,  
                         pi * (prey_info$size[4]*(0.75))/2 * (prey_info$size[4]*(0.75)/2)/2, 
                         pi * (prey_info$size[5]*(0.75))/2 * (prey_info$size[5]*(0.75)/2)/2,
                         pi * (prey_info$size[6]*(0.75))/2 * (prey_info$size[6]*(0.75)/2)/2,
                         pi * (prey_info$size[7]*(0.75))/2 * (prey_info$size[7]*(0.75)/2)/2,
                         pi * (prey_info$size[8]*(0.75))/2 * (prey_info$size[8]*(0.75)/2)/2,
                         pi * (prey_info$size[9]*(0.75))/2 * (prey_info$size[9]*(0.75)/2)/2,
                         pi * (prey_info$size[10]/2)^2, # assume crustacean larvae are circular
                         pi * (prey_info$size[11]/2)^2, # assume crustacean larvae are circular
                         pi * (prey_info$size[12]/2)^2, # assume crustacean larvae are circular
                         pi * (prey_info$size[13])/2 * (prey_info$size[13]/8)/2, # assume krill width is ~1/8 of length, assume ellipse shape (from side) # https://swfsc.noaa.gov/uploadedFiles/Operating_units/FRD/Survey_Technology/Broad-bandwidth%20TTS%20and%20absorption%20from%20m%20norvegica,%20mysids,%20and%20crangon.pdf
                         pi * (prey_info$size[14]/2)^2, # assume evadne is circular
                         pi * (prey_info$size[15]/2)^2, # egg is circular
                         10^(2.62 *log10( prey_info$size[16] ) -2.01), # fish larvae from https://www.cambridge.org/core/services/aop-cambridge-core/content/view/00790979F633729954C59C8B938BFB53/S0025315400032756a.pdf/div-class-title-developmental-changes-in-the-opacity-of-larval-herring-span-class-italic-clupea-harengus-span-and-their-implications-for-vulnerability-to-predation-div.pdf
                         pi * (prey_info$size[17])/2 * (prey_info$size[17]*(3/11))/2, # measuremet from imageJ to get ratio https://www.researchgate.net/figure/Themisto-compressa-Ovigerous-individuals-with-fertilized-eggs-and-recently-hatched_fig2_273269311
                         pi * (prey_info$size[18]*(0.75))/2 * (prey_info$size[18]*(0.75)/2)/2,
                         pi * (prey_info$size[19]*0.75)/2* (prey_info$size[19]*(0.75)/2)/2, # varies a ot but using same os calanoid probs reasonable https://copepodes.obs-banyuls.fr/en/fichesp.php?sp=1803
                         pi * (prey_info$size[20]*(0.75))/2* (prey_info$size[20]*(0.75)/2)/2,
                         pi * (prey_info$size[21]/2)^2, # assume podon is circular
                         pi * (prey_info$size[22]*(0.75))/2 * (prey_info$size[22]*(0.75)/2)/2)


#### add energy densities ####


cal2J = 4.184 # conversion factor calories/J
copd2w = 0.162 # percentage water cops


prey_info$energy_density = c(signif(5160.0 * cal2J * copd2w, digits = 2), # Laurence
                             signif(3176, digits = 2) , # Davis et al. 1998
                             signif(6425.1 * cal2J * copd2w, digits = 2), # Laurence
                             signif(6425.1 * cal2J * copd2w, digits = 2), # Laurence
                             signif((((-0.276*4+9.39)/1000000)*46000)/(0.08* ((-11.8*4+1157)/1000) ^2.1 /1000), digits = 2), #Campbell
                             signif(6425.1 * cal2J * copd2w, digits = 2), # Laurence
                             signif(4998.6 * cal2J * copd2w, digits = 2), # Laurence
                             signif(mean(c( 4998.6 * cal2J * copd2w , 5244.7 * cal2J * copd2w )), digits = 2),
                             signif(5244.7 * cal2J * copd2w, digits = 2), # Laurence
                             signif(((0.00201/1000)*46000)/ (prey_info$weight[10]/1000), digits = 2), # Muxagata
                             signif(((4.36*0.19^2.3)/1000000)*46000 / (prey_info$weight[11]/1000), digits = 2), # Tanskanen
                             signif((((10^(0.983 *log10(  10^(2.58*log10(prey_info$size[12])+2.04) ) -0.38))/1000000)*46000)  / (prey_info$weight[12]/1000), digits = 2), # Lindley
                             signif( 4700*cal2J*0.228614, digits = 2), # Kulka & Corey & Kiorboe
                             signif(((2/1000000) * 46000 ) / (prey_info$weight[14]/1000), digits = 2), # Rodhouse
                             signif(24000/(1/0.07), digits = 2), # Riis for water content and Paul for energy density
                             signif(2000, digits = 2), # Arrhenius
                             signif((-2.1303 + 4.0982*(0.0064*prey_info$size[17]^2.46))*cal2J  / (prey_info$weight[17]/1000), digits = 2),
                             signif((((prey_info$weight[18]*0.162*0.421)/1000)*46000)/(prey_info$weight[18]/1000), digits = 2), # Lnidley 97 & Kiorboe
                             signif((((10^(1.45*log10(1000*prey_info$size[19]*0.75)-4.25))/1000000) *46000)/(  10^(1.84*log10(1000*prey_info$size[19]*0.75)-4.84)/1000000 /0.162), digits = 2), # Uye & Kiorboe
                             signif(5070.9 * cal2J * copd2w, digits = 2), # Laurence
                             signif(((10^(4.15*log10(1000*prey_info$size[21])-11.15))/1000000) *46000 /  (prey_info$weight[21]/1000), digits = 2), # Uye
                             signif(4466.3 * cal2J * copd2w, digits = 2)) # Laurence


# calculate total energy
prey_info$energy = prey_info$energy_density*(prey_info$weight/1000)



#### save data ####
write.csv(prey_info, "prey_info.csv")


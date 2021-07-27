# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
#### HELP SCRIPT THAT MATCHES UP TAXA ACROSS L4/STONEHAVEN ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# set up empty data frame that matches L4 data in length, and no of taxa in cpr


L4_comparison = data.frame(
  matrix(NA, 
         nrow = nrow(L4),
         ncol = (1+length(taxa)))
)
names(L4_comparison) = c("date", taxa)
L4_comparison$date = L4$date


#### goes through each taxon in turn and matches up across datasets ####


L4_comparison$Acartia.spp...unidentified. = L4$Total.Acartia.clausi

L4_comparison$Appendicularia = L4$Appendicularia

L4_comparison$Calanus.helgolandicus = L4$Total.Calanus.helgolandicus # not sure if age groups align now

L4_comparison$Centropages.hamatus = L4$Centropages.hamatus

L4_comparison$Centropages.typicus = L4$Total.Centropages.typicus

L4_comparison$Centropages.spp...Unidentified. = L4$Centropages.chierchiae

L4_comparison$Cirripede.larvae..Total. = rowSums(cbind(L4$Cirripede.nauplii,
                                                       L4$Cirripede.cyprid))

L4_comparison$Copepod.nauplii = L4$Copepod.nauplii


L4_comparison$Decapoda.larvae..Total. = L4$Total.Decapoda


L4_comparison$Euphausiacea.Total = L4$Euphausiid.adult

L4_comparison$Evadne.spp. = L4$Evadne.spp.

L4_comparison$Fish.eggs..Total. = L4$Total.Fish.Eggs

L4_comparison$Fish.larvae = L4$Fish.larvae


L4_comparison$Hyperiidea..Total. = L4$Hyperiida 



L4_comparison$Metridia.lucens = L4$Metridia.lucens


L4_comparison$Oithona.spp. = L4$Oithona.spp.

L4_comparison$Para.Pseudocalanus.spp. =   rowSums(cbind(L4$Total.Paracalanus.parvus..Calculated.,
                                                        L4$Total.Pseudocalanus.elongatus..Calculated.,
                                                        L4$Total.Ctenocalanus.vanus..Calculated.,
                                                        L4$Total.Clausocalanus.spp...Calculated.,
                                                        L4$Microcalanus.spp.
))




L4_comparison$Podon.spp. = L4$Podon.spp.

L4_comparison$Temora.longicornis = L4$Total.Temora.longicornis



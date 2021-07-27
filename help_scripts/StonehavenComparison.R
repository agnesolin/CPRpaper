# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
#### HELP SCRIPT THAT MATCHES UP TAXA ACROSS CPR/STONEHAVEN ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# set up empty data frame that matches Stonehaven data in length, and no of taxa in cpr
stonehaven_comparison = data.frame(
matrix(NA, 
       nrow = nrow(stonehaven),
       ncol = (1+length(taxa)))
)

names(stonehaven_comparison) = c("date", taxa)
stonehaven_comparison$date = stonehaven$date



#### goes through each taxon in turn and matches up across datasets ####

stonehaven_comparison$Acartia.spp...unidentified. = rowSums(cbind(stonehaven$Acartia..Acartiura..clausi.C5,
                                                                  stonehaven$Acartia..Acartiura..clausi.C5.F, 
                                                                  stonehaven$Acartia..Acartiura..clausi.C6.F, 
                                                                  stonehaven$Acartia..Acartiura..clausi.C5.M, 
                                                                  stonehaven$Acartia..Acartiura..clausi.C6.M,
                                                                  stonehaven$Acartia..Acartiura..longiremis.C6.F,
                                                                  stonehaven$Acartia..Acartiura..longiremis.C6.M))
  


stonehaven_comparison$Appendicularia = stonehaven$Appendicularia

stonehaven_comparison$Calanus.finmarchicus = rowSums(cbind(stonehaven$Calanus.finmarchicus.C5, 
                                                           stonehaven$Calanus.finmarchicus.C6.F, 
                                                           stonehaven$Calanus.finmarchicus.C6.M))

stonehaven_comparison$Calanus.helgolandicus = rowSums(cbind(stonehaven$Calanus.helgolandicus.C5,
                                                            stonehaven$Calanus.helgolandicus.C6.F,
                                                            stonehaven$Calanus.helgolandicus.C6.M))

stonehaven_comparison$Calanus.I.IV = rowSums(cbind(stonehaven$Calanus.spp..C1, 
                                                   stonehaven$Calanus.spp..C2, 
                                                   stonehaven$Calanus.spp..C3,
                                                   stonehaven$Calanus.spp..C4))



stonehaven_comparison$Centropages.hamatus = rowSums(cbind(stonehaven$Centropages.hamatus.C5,
                                       stonehaven$Centropages.hamatus.C5.F,
                                       stonehaven$Centropages.hamatus.C6.F,
                                       stonehaven$Centropages.hamatus.C5.M,
                                       stonehaven$Centropages.hamatus.C6.M))

stonehaven_comparison$Centropages.typicus = rowSums(cbind(stonehaven$Centropages.typicus.C5,
                                                          stonehaven$Centropages.typicus.C6.F,
                                                          stonehaven$Centropages.typicus.C6.M))

stonehaven_comparison$Cirripede.larvae..Total. = rowSums(cbind(stonehaven$Cirripedia.nauplius,
                                                               stonehaven$Cirripedia.cypris))

stonehaven_comparison$Copepod.nauplii = stonehaven$copepod_nauplii


stonehaven_comparison$Decapoda.larvae..Total. = rowSums(cbind(stonehaven$Decapoda.larva ,
                                                              stonehaven$Nephrops.norvegicus.larva), na.rm = T)


stonehaven_comparison$Euphausiacea.Total = rowSums(cbind(stonehaven$Nyctiphanes.couchii.adult,
                                       stonehaven$Nyctiphanes.couchii.juvenile,
                                       stonehaven$Thysanoessa.inermis.adult,
                                       stonehaven$Thysanoessa.longicaudata.adult))

stonehaven_comparison$Evadne.spp. = stonehaven$Evadne.nordmanni

stonehaven_comparison$Fish.eggs..Total. = stonehaven$Pisces.egg

stonehaven_comparison$Fish.larvae = rowSums(cbind(stonehaven$Ammodytidae.larva, 
                                       stonehaven$Clupeidae.larva,
                                       stonehaven$Gadiformes.larva,
                                       stonehaven$Pisces.larva), na.rm = T)


stonehaven_comparison$Hyperiidea..Total. = rowSums(cbind(stonehaven$Hyperia.spp.,
                                      stonehaven$Themisto.spp.))



stonehaven_comparison$Metridia.lucens = rowSums(cbind(stonehaven$Metridia.lucens.C5,
                                           stonehaven$Metridia.lucens.C6.F,
                                           stonehaven$Metridia.lucens.C6.M))


stonehaven_comparison$Oithona.spp. = rowSums(cbind(stonehaven$Oithona.spp..C4.5,
                                   stonehaven$Oithona.spp..C6.M,
                                   stonehaven$Oithona.spp..C6.F))

stonehaven_comparison$Para.Pseudocalanus.spp.= rowSums(cbind(stonehaven$Paracalanus.parvus.C5.F,
                                                             stonehaven$Paracalanus.parvus.C5.M,
                                                             stonehaven$Paracalanus.parvus.C6.F,
                                                             stonehaven$Paracalanus.parvus.C6.M,
                                                             stonehaven$Pseudocalanus.minutus.elongatus.C5.F,
                                                             stonehaven$Pseudocalanus.minutus.elongatus.C5.M,
                                                             stonehaven$Pseudocalanus.minutus.elongatus.C6.F,
                                                             stonehaven$Pseudocalanus.minutus.elongatus.C6.M,
                                                             stonehaven$Ctenocalanus.vanus.C5,
                                                             stonehaven$Ctenocalanus.vanus.C6.F,
                                                             stonehaven$Ctenocalanus.vanus.C6.M,
                                                             stonehaven$Microcalanus.pusillus.C1.5,
                                                             stonehaven$Microcalanus.pusillus.C3,
                                                             stonehaven$Microcalanus.pusillus.C4,
                                                             stonehaven$Microcalanus.pusillus.C5,
                                                             stonehaven$Microcalanus.pusillus.C6.F,
                                                             stonehaven$Microcalanus.pusillus.C6.M,
                                                             stonehaven$Para.pseudo.cteno.clausocalanus.C1.5), na.rm = T)



stonehaven_comparison$Podon.spp. = rowSums(cbind(stonehaven$Podon.spp,
                                 stonehaven$Podon.leuckartii,
                                 stonehaven$Podon.intermedius), na.rm = T)

stonehaven_comparison$Temora.longicornis = rowSums(cbind(stonehaven$Temora.longicornis.C5,
                                              stonehaven$Temora.longicornis.C5.F,
                                              stonehaven$Temora.longicornis.C5.M,
                                              stonehaven$Temora.longicornis.C6.F,
                                              stonehaven$Temora.longicornis.C6.M))





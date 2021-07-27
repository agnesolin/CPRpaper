
# at end of InterpolationFeeding

## save as input for model ##
all_food = df
names(all_food)[1:3] = c("loc", "year", "jd")
write.csv(all_food,
          "~/nonSU/sandeel_model/EnvironmentalDrivers/Food/all_food.csv")


# at end of prey info




prey_info$mode = c("A",
                   "A",
                   "B",
                   "B",
                   "B",
                   "B",
                   "B",
                   "B",
                   "B",
                   "A",
                   "A",
                   "A",
                   "C",
                   "A",
                   "A",
                   "C",
                   "C",
                   "B",
                   "A",
                   "A",
                   "A",
                   "A")




write.csv(prey_info, "~/nonSU/sandeel_model/EnvironmentalDrivers/Food/prey_info.csv")



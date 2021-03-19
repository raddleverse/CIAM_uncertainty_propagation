### Make RefA

refA <- read_csv("/Users/catherineledna/.julia/dev/MimiCIAM/data/input/refA.csv")

refA_H <- refA %>% mutate(value=ifelse(variable=="H",value,0),variable="H")
refA_R <- refA %>% mutate(value=ifelse(variable=="R",value,0),variable="R")

write_csv(refA_H,"/Users/catherineledna/.julia/dev/MimiCIAM/data/input/refA_H.csv")
write_csv(refA_R,"/Users/catherineledna/.julia/dev/MimiCIAM/data/input/refA_R.csv")

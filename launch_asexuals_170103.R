

library(compiler)



source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/WITH ASEXUALS/master_asexuals_ded_170103.R")

master <- cmpfun(master)

source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/WITH ASEXUALS/initialization_asexuals_170103.R")

genesis <- cmpfun(genesis)

source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/WITH ASEXUALS/competition_170103.R")

dogeatdog <- cmpfun(dogeatdog)

source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/WITH ASEXUALS/faircompet_170103.R")

faircompetition <- cmpfun(faircompetition)

source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/WITH ASEXUALS/diaspora_customkernel_asexuals_170103.R")

diaspora <- cmpfun(diaspora)

source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/WITH ASEXUALS/meetic_asexuals_170103.R")

meetic <- cmpfun(meetic)

source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/WITH ASEXUALS/stillborns_asexuals_170103.R")

stillborn <- cmpfun(stillborn)

source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/WITH ASEXUALS/plotdensity_asexuals.R")

plotdensity <- cmpfun(plotdensity)

source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/plottruc.R")

plottruc <- cmpfun(plottruc)

source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/plotGRbabies.R")

plotGR <- cmpfun(plotGR)


 competition <- dogeatdog
 # competition <- faircompetition
 
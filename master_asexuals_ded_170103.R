
master <- function( s, 
                    K, 
                    fec,
                    fecasex,
                    G,
                    probamating,
                    dispkernel,
                    bsline,
                    Ka,
                    B,
                    M,
                    pmut,
                    tps, 
                    run
                    )
  {

  
  Ninit <- K * fec
  
  # Initialization
  
  # input: K, Ninit, G
  # output: popgenome, sex, newbabyR, newbabyB
   
  TEMP <- genesis( K = K, fec = fec, Ninit = Ninit, G = G) 
  
  popgenome <- TEMP$popgenome
  sex <- TEMP$sex
  newbabyR <- TEMP$newbabyR
  newbabyB <- TEMP$newbabyB
  repro  <- TEMP$repro
  
  
  # 0. Initial dispersal
  # 
  # input: newbabyR, newbabyB, 
  # output: newpopR, newpopB, newpopRB, coordalive, haspartner
  
  TEMP <- diaspora( newbabyR = newbabyR, newbabyB = newbabyB, dispkernel = dispkernel, sex = sex, repro = repro )
  
  popRB <- TEMP$newpopRB
  popB <- TEMP$newpopB
  popR <- TEMP$newpopR
  coordalive <- TEMP$coordalive
  haspartner <- TEMP$haspartner
  
#   # saves the reproductive population size per patch, after the first ever round of dispersion
#   popsize_past <- tapply( !is.na(coordalive), coordalive, sum )
#   popRB_past <- popRB
#   popR_past <- popR
#   popB_past <- popB
  
  # preparing data gathering
  
  INVASION_tot <- numeric( length = tps ) # pop size gathered once adults have settled and survived competition, before reproduction
  INVASION_asex <- numeric( length = tps )
  
  ######### --- enter loop --- ########
  
  for ( t in 1: tps ){
    
    
    # 1. Competition
    # 
    # input: haspartner, coordalive, popgenome, K, G
    # output: popsurvival
    
    ## to switch on genotype competition
    popsurvival <- competition( haspartner = haspartner, coordalive = coordalive, popgenome = popgenome, K = K, G = G)
    

    # numer of patches with surviving individuals; and with more than half of their population asexual
    INVASION_tot[ t ] <- length( unique( coordalive[ popsurvival == 1 ] ) )
    temptable <- table( coordalive[ popsurvival == 1 ], factor( repro[ popsurvival == 1 ], levels = c('a', 's') ) )
    INVASION_asex[ t ] <- sum( temptable[,1]/temptable[,2] > .5 )
    
    
    # 2. Reproduction: Allee effects and offspring production - meetic
    # 
    # input: sex, coordalive, parameters Allee effect, popgenome, popsurvival
    # output: babyR, babyB, babysex, babygenome
    
    TEMP <- meetic( sex = sex, coordalive = coordalive, K = K, probamating = probamating, fec = fec, fecasex = fecasex, popgenome = popgenome, popsurvival = popsurvival, repro = repro, popRB = popRB, popR = popR, popB = popB, G = G )
    
    if( !is.list( TEMP )) break
    
    babysex <- TEMP$babysex
    babyR <- TEMP$babyR
    babyB <- TEMP$babyB
    babygenome <- TEMP$babygenome
    babyrepro <- TEMP$babyrepro
    
    
    
    # 4. Survival of inbreeding to juvenile stage - stillbornn
    # 
    # input: babygenome, babysex, babyR, babyB
    # output: newbabygenome, newbabysex, newbabyR, newbabyB
    # 
    
    TEMP <- stillborn( babygenome = babygenome, babysex = babysex, babyR = babyR, babyB = babyB, B = B, M = M, Ka = Ka, bsline = bsline, pmut = pmut, babyrepro = babyrepro, G = G )
    
    sex <- TEMP$newbabysex
    popgenome <- as.matrix( TEMP$newbabygenome, drop = FALSE )
    newbabyR <- TEMP$newbabyR
    newbabyB <- TEMP$newbabyB
    
    repro <- TEMP$newbabyrepro
    
    

    
    
    # Plot density of total, and asexual populations.
    
#     if ( round( t/1 ) == t/1){
#     popsize <- tapply( popsurvival, coordalive, sum )
#     asexpopsize <- tapply( popsurvival[ repro == "a"], coordalive[ repro == "a"], sum ) 
#     plotdensity( popsize = popsize, asexpopsize = asexpopsize, plotname = 'Density per patch' )
#     }
    
    
    
    
    
    # 5. Dispersal - voyage
    # 
    # input: newbabyR, newbabyB, dispersal kernel, sex, repro
    # output: newpopR, newpopB, newpopRB, coordalive, haspartner
    
    # TEMP <- diaspora( newbabyR = newbabyR, newbabyB = newbabyB, short_dist_sd = short_dist_sd, plongdist = plongdist, long_dist_sd = long_dist_sd, long_dist_mean = long_dist_mean)
    
    TEMP <- diaspora( newbabyR = newbabyR, newbabyB = newbabyB, dispkernel = dispkernel, sex = sex, repro = repro )
    
    popRB <- TEMP$newpopRB
    popB <- TEMP$newpopB
    popR <- TEMP$newpopR
    coordalive <- TEMP$coordalive
    haspartner <- TEMP$haspartner
    
    # Limit cases
    if( length( popRB ) == 0 | sum ( is.na(coordalive) ) == length(coordalive) ){
      popsurvival <- NA
      break
    } 
    if( sum( repro == "s" ) / sum( repro == "a") < 0.1 ){
      popsurvival <- competition( haspartner = haspartner, coordalive = coordalive, popgenome = popgenome, K = K, G = G)
      break
    } 
    
    ######### --- loop --- ########
    
    print( paste( run, ':', t ) )
    
  }
  
  if( t == tps ) popsurvival <- competition( haspartner = haspartner, coordalive = coordalive, popgenome = popgenome, K = K, G = G)
  
  
  
 namerun <- paste( .vcompet[run], "K", .vK[run], '_fs', .vfec[run], '_fa', .vfecasex[run], '_', .vprobamating[run],
                   '_G', .vG[run], '_b', .vbsline[run], '_Ka', .vKa[run], '_B', .vB[run],'M_', .vM[run], 
                   '_pm', .vpmut[run], '_c', .vc[run], '_a', .va[run], '.RData', sep = "" )
 # backup <- paste( 'C:/Users/Anais Tilquin/Desktop/SimulLaCiotat/', namerun, sep = "" ) 
  
 save( list = ls( all.names = TRUE ), file = namerun, envir = environment() )
 # save( list = ls( all.names = TRUE ), file = backup, envir = environment() )
  
  rm( list = setdiff( ls(), lsf.str() ) ) # removes all variables but functions
}


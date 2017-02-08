rm( list = ls( all = T ))

setwd("/home/ubuntu/simulations")

#############################################################

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
                    run)
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
                    '_G', .vG[run], '_b', .vbsline[run], '_Ka', .vKa[run], '_B', .vB[run],'_M', .vM[run], 
                    '_pm', .vpmut[run], '_c', .vc[run], '_md', .vmean_distance[run], '.RData', sep = "" )
  # backup <- paste( 'C:/Users/Anais Tilquin/Desktop/SimulLaCiotat/', namerun, sep = "" ) 
  
  save( list = ls( all.names = TRUE ), file = namerun, envir = environment() )
  # save( list = ls( all.names = TRUE ), file = backup, envir = environment() )
  
  rm( list = ls() ) # removes all variables but functions
}


#############################################################





##### creating the initial population kernel

# INPUT: K, Ninit, fec, G
# 
# OUTPUT: sex, popgenome, newbabyR, newbabyB, repro

genesis <- function(K, 
                    fec, 
                    Ninit, G){
  
  # pop starts half female, half male, all sexuals, all on patch (0;0)
  sex <- c( rep( "fem", each = ceiling( Ninit/2 )), rep( "mal", each= floor( Ninit/2 )))       
  repro <- rep("s", Ninit)
  newbabyR <- rep( 0, Ninit )
  newbabyB <- rep( 0, Ninit )
  
  # creating a collection of Ninit genotypes for one locus at equilibrium proportions (1/4 00, 1/4 11, 1/2 01)
  col1 <- head( sample( c( rep(1, ceiling( Ninit/2 )), rep(0, ceiling( Ninit/4 )), rep(2, ceiling( Ninit/4 ))) ) , Ninit )
  
  # matrix with population genotype for 1 gene
  popgenome <- as.matrix( col1, drop = FALSE )
  
  # extended to G genes if G > 1
  if( G != 1 ){
    for ( g in 1:(G-1) ){
      popgenome <- cbind( popgenome, sample( col1, Ninit ) )
    }}
  
  
  return( list( sex = sex, popgenome = popgenome, newbabyR = newbabyR, newbabyB = newbabyB, repro = repro ))
}


#############################################################


### where different genotypes can exploid different (limited) resources

# input: haspartner, coordalive, popgenome, G, K
# output : popsurvival

# each patch can hold at max K individuals, as it has K units of each resource type (ie for each gene)

dogeatdog <- function( coordalive, 
                       haspartner, 
                       popgenome, 
                       G, 
                       K){
  
  occupiedpatches <- sort( unique( coordalive ) )
  
  # prepare for computing the number of copies of each allele present in each patch
  patch_compet <- matrix( nrow = length(occupiedpatches) , ncol= 2*G )
  rownames( patch_compet ) = occupiedpatches # /!/ names are sorted
  
  for (g in 1:G){
    Tgenepatch <- table( factor( popgenome[,g], levels = 0:2), coordalive ) # table with allelic sums per patch for gene g
    patch_compet[,(g*2-1)] <- Tgenepatch[1,]*2 + Tgenepatch[2,] # vector with number of copies of the 0 allele of gene g on the patch
    patch_compet[,(g*2)] <- Tgenepatch[3,]*2 + Tgenepatch[2,] # vector with number of copies of the 1 allele of gene g on the patch
  }
  
  # table with reward of possessing each allele (col: al0Gen1 al1Gen1 al0Gen2 al1Gen2...) on each patch (row)
  patch_compet <- K / patch_compet # total resource of each type available in the patch, divided by number of corresponding alleles in the patch
  patch_compet[ patch_compet > 1 ] <- 1 # on patches with excess food, individual share is bounded at 1; Inf values also become 1 but that doesn't matter as they stand for alleles no-one has
  
  popgenome_haspartner <- popgenome[haspartner==T,]
  reward_haspartner <-  patch_compet[ coordalive[ haspartner == T ], ] # matrix of the reward got by each individual (row) from each of his alleles (g1a0, g1a1, g2a0 etc)
  
  popfitness_haspartner <- vector( length = sum(haspartner) ) 
  for( g in 1:G ){
    fitness_g <- popgenome_haspartner[,g] * reward_haspartner[, 2*g] + ( 2- popgenome_haspartner[,g] ) * reward_haspartner[, 2*g-1]
    popfitness_haspartner <- popfitness_haspartner + fitness_g
  }
  popfitness_haspartner <- popfitness_haspartner / (2*G) # the survival probability resulting from resource acquisition, of each individual with reproductive prospects
  
  # each individual of the population now stochastically survives (1) or not (0). Without reproductive prospects, the individual dies.
  popsurvival <- numeric( length(haspartner) )
  popsurvival[ haspartner==T ] <- mapply( FUN = rbinom, prob = popfitness_haspartner, size = 1, n = 1)
  
  return(popsurvival)
}

#############################################################

### where everybody shares equally the resource of the patch

# input: coordalive, K ... (should accept popgenome, G, haspartner, from the other competition function)
# output : popsurvival

faircompetition <- function( coordalive, 
                             K, 
                             ... ){
  
  pop_patch <- table( coordalive )
  popfitness <- K / pop_patch[ coordalive ]
  popfitness[ popfitness > 1 ] <- 1
  condition <- !is.na( popfitness )
  popsurvival <- rep( 0, length(coordalive) )
  popsurvival[ condition ] <- mapply( FUN = rbinom, prob = popfitness[ condition ], size = 1, n = 1)
  
  return( popsurvival )
}





#############################################################


### when the juveniles go see further if grass is greener there

# input: newbabyR, newbabyB, dispersal kernel, sex, repro
# output: newpopR, newpopB, newRB, coordalive

diaspora <- function( newbabyR, 
                      newbabyB, 
                      dispkernel, 
                      sex, 
                      repro ) {
  
  
  distance <- sample( x = dispkernel, size = length(newbabyR), replace = T )
  
  angle <- runif( n = length(newbabyR), min=0, max = 2*pi )
  
  # translate movement on Red and Blue hexagonal grid axes
  movR <- round( distance * cos( angle ) )
  movB <- round( distance * cos( angle + (4/3)*pi) )
  
  newpopB <- newbabyB + movB
  newpopR <- newbabyR + movR
  
  ## do they end up in a patch without potential mating partners?
  newpopRB <- paste( newpopR,newpopB,sep="," )
  
  # create vector for females / males saying if there are mating partners in their patch 
  # (aim: save time as no mating partner = as good as dead, so no need to calculate further survival)
  
  haspartner <- vector( length = length( newpopRB ))
  haspartner[ sex == "fem" & repro == "s" ] <- newpopRB[ sex == "fem" & repro == "s" ] %in% newpopRB[ sex == "mal" ] 
  haspartner[ sex == "fem" & repro == "a" ] <- TRUE
  haspartner[ sex == "mal" ] <- newpopRB[ sex =="mal" ] %in% newpopRB[ sex =="fem" & repro == "s" ] 
  
  coordalive <- newpopRB
  coordalive[!haspartner] <- NA # coordinates of individuals, lonely and doomed ones are marked NA
  
  return( list( newpopR = newpopR, newpopB = newpopB, newpopRB = newpopRB, coordalive = coordalive, haspartner = haspartner ))
  
}



#############################################################


### where they find a soul mate

# INPUT: sex, coordalive, K, parameters Allee effect: probamating, theta, fecundity f, popRB, popR, popB, G
# OUTPUT: babysex, babyR, babyB, babygenome, babyrepro


meetic <- function( sex, 
                    coordalive, 
                    K, 
                    probamating, 
                    theta, 
                    fec, 
                    fecasex, 
                    popgenome, 
                    popsurvival, 
                    repro, 
                    popRB, 
                    popR, 
                    popB, 
                    G ){
  
  ## mate-finding Allee-effect
  
  if( sum( sex[popsurvival==1]=='mal' ) == 0 ) return( 'nomales' )
  
  # table of number of males per patch
  m_per_patch <- table( sex[popsurvival==1], coordalive[popsurvival==1])['mal', , drop=F] 
  if ( length( m_per_patch ) == 1) names(m_per_patch) <- unique( coordalive[popsurvival==1 & sex=='mal'] )
  
  if ( is.numeric(probamating) ) {
    pmating <- probamating * popsurvival
    pmating[ coordalive %in% colnames(m_per_patch)[m_per_patch==0] ] <- 0 # only females with males on their patch have pmating, otherwise 0
  }
  if ( probamating == "allee") {
    if ( missing('theta') ) {
      theta <- 3 / ( 0.5 * K ) # this theta gives such a curve for mating probability:
      # curve(expr= 1 - exp ( - 6 * x ), from=0, to=1, xlab='percent of K that are males', ylab='female mating proba')
    }
    pmating <-  (1 - exp ( - theta * m_per_patch[popRB] )) 
    pmating <- pmating * popsurvival # removing the females that didn't survive the competition earlier
    pmating[is.na(pmating)] <- 0 # otherwise rbinom errors over the NA, see if can't code them as 0 from the start
  }
  
  # Only females get a 1 or 0 (so the vector is the size of the female pop), and males are picked by females randomly; saves some random draw for the males
  realized_mating <- mapply( FUN = rbinom, prob = pmating[ sex =="fem" & repro == "s" ], size = 1, n = 1 ) 
  if( sum(realized_mating) == 0 ) return( 'nomating')
  
  # coordinates of the sexual mothers
  mum_patch <- coordalive[ sex == "fem" & repro == "s" ][ realized_mating == 1 ]
  
  # total number of sexually + asexually produced offspring
  newpopsize <- fec * sum( realized_mating ) + fecasex * sum( repro == 'a' & popsurvival == 1 )
  
  # each sexual female chooses one mate on her patch
  potential_dads <-lapply( X = unique(popRB[which( popsurvival == 1 & sex == "mal" )]), 
                           FUN = function(x) which( popRB[ which( popsurvival == 1 & sex == "mal" )] == x ) )
  names(potential_dads) <-  unique(popRB[which( popsurvival == 1 & sex == "mal" )])
  
  dad_id <- unlist(lapply( X = mum_patch,
                           FUN = function(x) sample(potential_dads[[which(names(potential_dads)== x)]], 1) ))
  
  # genotypes of babies produced sexually
  mumgenes <- popgenome[ sex=="fem" & repro == "s" , , drop = F][ realized_mating == 1 , , drop = F ]
  dadgenes <- popgenome[ dad_id, ]
  mendel_parents <- mumgenes + dadgenes 
  mendel <- mendel_parents[rep(1:nrow(mendel_parents), times = fec), ]
  
  mendel0 <- mendel == 0
  mendel4 <- mendel == 4
  mendel1 <- mendel == 1
  mendel3 <- mendel == 3
  mendel2 <- mendel == 2
  mendelmom <- mumgenes[rep(1:nrow(mumgenes), times = fec), ] == 1
  
  # REMARK: ( 1 baby per female ), fec times ; NOT ( fec baby ), nb of female times
  sex_babygenome <- matrix( nrow = nrow(mendel) , ncol = G )
  sex_babygenome[ mendel0 ] <- 0
  sex_babygenome[ mendel4 ] <- 2
  sex_babygenome[ mendel1 ] <- sample(0:1, size = sum( mendel1 ), replace = TRUE)
  sex_babygenome[ mendel3 ] <- sample(1:2, size = sum( mendel3 ), replace = TRUE)
  sex_babygenome[ mendel2 & mendelmom ] <- sample( c(0,1,1,2) , size = sum( mendel2 & mendelmom ), replace = TRUE )
  sex_babygenome[ mendel2 & !mendelmom ] <- 1
  
  # genotypes of babies produced asexually
  asex_babygenome <- matrix( nrow = sum( repro == "a" & popsurvival == 1 ) * fecasex , ncol = G )
  if( sum( repro == "a" & popsurvival == 1 ) > 0 ){
    asex_babygenome_unique <- popgenome[ repro == "a" & popsurvival == 1, , drop = F ]  
    asex_babygenome <- asex_babygenome_unique[rep(1:nrow(asex_babygenome_unique), times = fecasex), ] # this line repeats the ENTIRE matrix fecasex fois
  }
  
  babygenome <- rbind(sex_babygenome, asex_babygenome)
  
  # sex of babies, sexuals then asexuals, in one vector
  babysex <- sample( c("fem", "mal"), replace = T, size = fec * sum( realized_mating ) )
  babysex <- c( babysex, rep("fem", fecasex * sum( repro == 'a' & popsurvival == 1 )) ) 
  
  # reproductive mode of babies, sexuals then asexuals, in one vector
  babyrepro <- c( rep("s", fec*sum( realized_mating ) ), rep( "a", fecasex * sum( repro =="a" & popsurvival == 1 )) )
  
  # babies birthplace, sexuals then asexuals, in one vector
  babyR <- c( rep( popR[ sex=="fem" & repro=="s"][ realized_mating ==1], fec ), rep( popR[ repro == "a" & popsurvival == 1 ], fecasex ) )
  babyB <- c( rep( popB[ sex=="fem" & repro=="s"][ realized_mating ==1], fec ), rep( popB[ repro == "a" & popsurvival == 1 ], fecasex ) )
  
  return( list( babysex = babysex, babyR = babyR, babyB = babyB, babygenome = babygenome, babyrepro = babyrepro ) )
  
}


#############################################################


### where viability selection occurs depending on the levels of homozygosity

# INPUT: babygenome, babysex, babyR, babyB, G
# OUTPUT: newbabygenome, newbabysex, newbabyR, newbabyB 

# the relationship between heterozygosity and survival is:
# curve(expr = Ka /( 1 + 1*exp( -B*(x-M) ) ), from=0, to=1)

stillborn <- function( babygenome, 
                       babysex, 
                       babyR, 
                       babyB, 
                       B, 
                       M, 
                       Ka, 
                       bsline, 
                       pmut, 
                       babyrepro, 
                       G ){
  
  # each juvenile's proportion of heterozygous loci
  heteroz <- apply( babygenome == 1, 1 , sum ) / G
  
  # each juvenile's survival probability to inbreeding depression
  inbreeding <- Ka /( 1 + 1*exp( -B*(heteroz-M) )) + bsline
  
  # each juvenile's actual survival
  babysurvival <- as.logical( mapply( FUN = rbinom, prob = inbreeding, size = 1, n = 1 ) )
  
  newbabygenome <- babygenome[ babysurvival, ]
  newbabysex <- babysex[ babysurvival ]
  newbabyR <- babyR[ babysurvival ]
  newbabyB <- babyB[ babysurvival ]
  newbabyrepro <- babyrepro[ babysurvival ]
  
  nmutants <- rbinom(n = 1, size = sum( newbabysex == "fem" & newbabyrepro == "s" ), prob = pmut)
  idmutants <- sample( 1: sum( newbabysex == "fem" & newbabyrepro == "s" ), nmutants)
  
  newbabyrepro[ newbabysex == 'fem' & newbabyrepro == "s"][ idmutants ] <- "a"
  
  return( list( newbabygenome = newbabygenome, newbabysex = newbabysex, newbabyR = newbabyR, newbabyB = newbabyB, newbabyrepro = newbabyrepro ))
  
}



#############################################################

library(compiler)


.genesis <- cmpfun(genesis)

.dogeatdog <- cmpfun(dogeatdog)

.faircompetition <- cmpfun(faircompetition)

.diaspora <- cmpfun(diaspora)

.meetic <- cmpfun(meetic)

.stillborn <- cmpfun(stillborn)

.master <- cmpfun(master)


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


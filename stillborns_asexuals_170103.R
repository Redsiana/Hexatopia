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


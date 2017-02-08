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
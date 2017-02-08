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





# to homogeneize initial conditions between trials, pop always start with exactly half males and half females, all heterozygous



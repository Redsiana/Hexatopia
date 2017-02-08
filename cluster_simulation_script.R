.nsimul <- 100

.vK <- sample( x = 10:30, size = .nsimul, replace = T )
.vfec <- sample( x = 3:10, size = .nsimul, replace = T )
.vfecasex <- numeric( .nsimul )
for( i in 1:.nsimul){
  .vfecasex[i] <- sample( x = 2:.vfec[i], size = 1, replace = T )
} 
.vG <- sample( x = 2:6, size = .nsimul, replace = T )
# .vprobamating <- sample( x = c(1, 'allee'), size = .nsimul, replace = T )
.vprobamating <- round( runif( n = .nsimul, min = 0.8, max = 1 ), digits = 2 )

.vcompet <- sample( x = c('ded_', 'fc_'), size = .nsimul, replace = T )
.vc <- round( runif( n = .nsimul, min = 0.5, max = 3 ), digits = 2 )
.vmean_distance <- round( runif( n = .nsimul, min = 0.5, max = 4 ), digits = 1 )

.vbsline <- round( runif( n = .nsimul, min = 0, max = 1 ), digits = 2 )
.vKa <- 1 - .vbsline
.vB <- round( runif( n = .nsimul, min = 0, max = 50 ), digits = 2 )
.vM <- round( runif( n = .nsimul, min = 0, max = 1 ), digits = 2 )
.vpmut <- round( runif( n = .nsimul, min = 0.0005, max = 0.005), digits = 5 )

.tps <- 30

.s <- 2/3




for ( .run in 1: .nsimul )
{
  
  run <- .run
  
  master <- .master
  genesis <- .genesis
  dogeatdog <- .dogeatdog
  faircompetition <- .faircompetition
  diaspora <- .diaspora
  meetic <- .meetic
  stillborn <- .stillborn
  
  if( .vcompet[run] == 'ded_') {
    competition <- .dogeatdog
  } else {
    competition <- .faircompetition
  }
  
  a <- .vmean_distance[run] / ( gamma( 2/.vc[run] ) / gamma( 1/.vc[run] ) )
  nsamples <- 1E5
  x <- runif(nsamples, min = 0.5001, max = 1)
  Gamma <- function(a, z) pgamma(z, a, lower = FALSE) * gamma(a) # define the same incomplete gamma function as Mathematica, which is gamma(a) from z to infinity
  f2 <- function(x, u) 1 - ( x* (abs(x/a)^.vc[run]) ^(-1/.vc[run]) * Gamma( 1/.vc[run] , abs(x/a)^.vc[run] ) ) / ( 2*a*gamma(1/.vc[run]) ) - u
  my.uniroot2 <- function(x) uniroot(f2, interval = c(0.000001, 100), u = x, tol = 0.0001)$root
  my.uniroot2 <- cmpfun( my.uniroot2 )
  dispkernel <- vapply(x, my.uniroot2, numeric(1))
  
  
  
  if ( .vprobamating[run] != 'allee' ) probamating <- as.numeric( .vprobamating[run] ) else probamating <- .vprobamating[run] 
  
  master( s = .s,
          K = .vK[run],
          fec = .vfec[run],
          fecasex = .vfecasex[run],
          G = .vG[run],
          probamating = probamating,
          dispkernel = dispkernel,
          bsline = .vbsline[run],
          Ka = .vKa[run],
          B = .vB[run],
          M = .vM[run],
          pmut = .vpmut[run],
          tps = .tps,
          run = run ) 
  print(run)
}


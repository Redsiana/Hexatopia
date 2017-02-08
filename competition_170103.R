### where different genotypes can exploid different (limited) resources

# input: haspartner, coordalive, popgenome, G, K
# output : popsurvival

# each patch can hold at max K individuals, as it has K units of each resource type (ie for each gene)

dogeatdog <- function( coordalive, haspartner, popgenome, G , K){
  
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





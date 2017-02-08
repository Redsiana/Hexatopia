### when the juveniles go see further if grass is greener there
## et donc, on test git-hub
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
  
  ################### TO TEST - PLOT IT ALL
  # 
  # ############# Plot Hexagonia
  # y = 3/2 * s * Blueworld # get cartesian coordinates of centers on axis y
  # x = sqrt(3) * s * ( Blueworld/2 + Redworld) # get cartesian coordinates of centers on axis x
  # S <- SpatialPoints(cbind(x,y))
  # hex_grid <- HexPoints2SpatialPolygons(S)
  # plot(hex_grid, main="Hexatopia", border="white")
  # abline(0,-0.5773503, col="red")
  # abline(v=0, col="blue")
  
  # ############ Plot where the babies start
  # y = 3/2 * s * newbabyB # get cartesian coordinates of centers on axis y
  # x = sqrt(3) * s * ( newbabyB/2 + newbabyR) # get cartesian coordinates of centers on axis x
  # S <- SpatialPoints(cbind(x,y))
  # hex_grid <- HexPoints2SpatialPolygons(S, dx=1)
  # plot(hex_grid, main="Hexatopia", col="cadetblue3", add=T)
  
  
  
  ############# Plot trajectories on Hexagonia
  # 
  # initX <- sqrt(3) * s * ( newbabyB/2 + newbabyR)
  # initY <- 3/2 * s * newbabyB
  # finX <- sqrt(3) * s * ( newpopB/2 + newpopR)
  # finY <- 3/2 * s * newpopB
  # arrows(x0=initX, y0=initY, x1=finX, y1=finY, length=0.1, col="red")
  # 
  # nameplot <- paste('hexatopia', 2*t-1,'.png', sep="")
  # dev.copy(png, nameplot)
  # dev.off()
  
  # ############# Highlight newly colonized patches
  # y = 3/2 * s * newpopB # get cartesian coordinates of centers on axis y
  # x = sqrt(3) * s * ( newpopB/2 + newpopR) # get cartesian coordinates of centers on axis x
  # S <- SpatialPoints(cbind(x,y))
  # hex_grid <- HexPoints2SpatialPolygons(S, dx=sqrt(3)*2/3) # if don't specify a dx, plot will screw up when only one patch is plotted
  # plot(hex_grid, main= paste("Hexatopia", t), col="red", add=T)
  # # arrows(x0=initX, y0=initY, x1=finX, y1=finY, length=0.1, col="red")
  # 
  # nameplot <- paste('hexatopia', 2*t,'.png', sep="")
  # dev.copy(png, nameplot)
  # dev.off()
  
}



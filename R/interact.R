interact <-
function(x, y, z = NULL, numPerm = 100, numFDR = 1000, method = "Pearson", verbose = TRUE){
  
  this.call <- match.call()

  if(method == "Pearson"){
    method = "pearson"
  }
  if(method == "Spearman"){
    method = "spearman"
  }

  if(is.null(colnames(x))){
    colnames(x) <- 1:ncol(x)
  }
  if(!is.null(z)){
    if(is.null(colnames(z))){
      colnames(z) <- 1:ncol(z)
    }
  }
  
  if(is.null(z)){
    numFDR <- min(numFDR, ncol(x) * (ncol(x) - 1)/2)
  }

  if(!is.null(z)){
    numFDR <- min(numFDR, ncol(z) * ncol(x))
  }

  
  if(length(colnames(x)) == 0){
    colnames(x) <- 1:ncol(x)
  }
  
  groups <- unique(y)

  ind1 <- which(y == groups[1])
  ind2 <- which(y == groups[2])

  x1 <- x[ind1,,drop=F]
  x2 <- x[ind2,,drop=F]

  x1 <- t((t(x1) - apply(x1, 2, mean))/apply(x1,2,sd))
  x2 <- t((t(x2) - apply(x2, 2, mean))/apply(x2,2,sd))

  if(!is.null(z)){
    z1 <- z[ind1,,drop=F]
    z2 <- z[ind2,,drop=F]
    
    z1 <- t((t(z1) - apply(z1, 2, mean))/apply(z1,2,sd))
    z2 <- t((t(z2) - apply(z2, 2, mean))/apply(z2,2,sd))
    }

  if(is.null(z)){
    cors1 <- cor(x1, x1, method = method)
    cors2 <- cor(x2, x2, method = method)
  }
  if(!is.null(z)){
    cors1 <- cor(x1, z1, method = method)
    cors2 <- cor(x2, z2, method = method)
  }
    
  stats <- atanh(cors1) - atanh(cors2)
  rm(cors1)
  rm(cors2)

  ## Making a matrix to figure out which interactions are important ##

  indOrder <- matrix(1:length(stats), ncol = ncol(stats), nrow = nrow(stats))

  if(is.null(z)){
    indOfOrder <- indOrder[upper.tri(indOrder)]
    rm(indOrder)
    stats <- stats[upper.tri(stats)]
    num.greater <- rep(0, numFDR)
    ordInd <- order(abs(stats), decreasing = TRUE)
    sortIndOfOrder <- indOfOrder[ordInd]
    rm(indOfOrder)
    G1 <- sortIndOfOrder %% ncol(x)
    G2 <- floor(sortIndOfOrder / ncol(x)) + 1
    interactionOrder <- cbind(colnames(x)[G1],colnames(x)[G2])
    rm(sortIndOfOrder)
  }

    if(!is.null(z)){
    indOfOrder <- indOrder
    rm(indOrder)
    num.greater <- rep(0, numFDR)
    ordInd <- order(abs(stats), decreasing = TRUE)
    sortIndOfOrder <- indOfOrder[ordInd]
    rm(indOfOrder)
    G1 <- sortIndOfOrder %% ncol(x)
    G2 <- floor(sortIndOfOrder / ncol(x)) + 1
    ind0 <- which(G1 == 0)
    G1[ind0] <- ncol(x)
    G2[ind0] <- G2[ind0] - 1
    interactionOrder <- cbind(colnames(x)[G1],colnames(z)[G2])
    rm(sortIndOfOrder)
  }

  ## Done ##

  sortStats <- abs(stats)[ordInd]  

  x.new <- rbind(x1,x2)
  
  if(!is.null(z)){
    z.new <- rbind(z1,z2)   
  }
##################

  for(i in 1:numPerm){    if(verbose == TRUE){
      write(paste("Permutation #",i, sep = ""),"")
    }
    
    permInd <- sample(1:nrow(x.new), replace = FALSE)

    ind1.perm <- which(y[permInd] == groups[1])
    ind2.perm <- which(y[permInd] == groups[2])

    x1.perm <- x.new[ind1.perm,,drop=F]
    x2.perm <- x.new[ind2.perm,,drop=F]

    if(!is.null(z)){
      z1.perm <- z.new[ind1.perm,,drop = F]
      z2.perm <- z.new[ind2.perm,,drop = F]
    }
    
    if(is.null(z)){
      perm.cors1 <- cor(x1.perm, x1.perm, method = method)
      perm.cors2 <- cor(x2.perm, x2.perm, method = method)
    }

    if(!is.null(z)){
      perm.cors1 <- cor(x1.perm, z1.perm, method = method)
      perm.cors2 <- cor(x2.perm, z2.perm, method = method)
    }
    
    perm.stats <- atanh(perm.cors1) - atanh(perm.cors2)

    if(is.null(z)){
      perm.stats <- perm.stats[upper.tri(perm.stats)]
    }
    
    rm(perm.cors1)
    rm(perm.cors2)

    perm.sortStats <- sort(abs(perm.stats), decreasing = TRUE)

    rm(perm.stats)


    count <- 1

    junk <- .C("calcFDR", numFDR = as.integer(numFDR), stats = as.double(sortStats[1:numFDR]), permStats = as.double(perm.sortStats), numGreater = as.integer(rep(0,numFDR)), numPermStats = as.integer(length(perm.sortStats)))

    num.greater <- junk$numGreater + num.greater
  }
  FDR <- pmin(num.greater/ (numPerm * (1:numFDR)),1)

  ord.FDR <- FDR
  junk <- .C("cumMax", stats = as.single(FDR), monotonicStats = as.single(ord.FDR), length = as.integer(length(FDR))) ## This makes FDR nondecreasing

  FDR <- junk$monotonicStats

  interactionOrder <- data.frame(interactionOrder[1:numFDR,], FDR)
  colnames(interactionOrder) <- c("feature1", "feature2", "qval")

  stuff <- list(interaction.ordered = interactionOrder, internals = list(FDR = FDR, stats = stats[ordInd], call = this.call))
  class(stuff) = "interact"
  
  return(stuff)
}

print.interact <- function(x, ...){
  cat("Call:\n")
  dput((x$internals)$call)
  cat("\n")


  mat = x$interaction.ordered[1:(min(10,nrow(x$interaction.ordered))),]
  cat("Most significant marginal interactions:\n")
  print(mat, quote = FALSE)
  invisible()
}

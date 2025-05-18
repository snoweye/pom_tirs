library(MASS)

### Lunceford and Davidian (2004): Z.given.X = TRUE, fixed.n = FALSE
### Modified version: Z.given.X = FALSE, fixed.n = TRUE (default)
generate.continuous.LD04m <- function(n.ctr = 2000, n.trt = 100, param,
    Z.given.X = FALSE, fixed.n = TRUE){
  nu <- param$nu
  a <- param$a
  xi <- param$xi
  p.X.3 <- param$p.X.3
  w.V.3 <- param$w.V.3
  tau.1 <- param$tau.1
  tau.0 <- param$tau.0
  Sigma.1 <- param$Sigma.1
  Sigma.0 <- param$Sigma.0

  err.sd <- 1
  if(exists("err.sd", param)){
    err.sd <- param$err.sd
  }

  ### Total n
  n <- n.ctr + n.trt

  n.ctr.tmp <- n.ctr
  n.trt.tmp <- n.trt
  X.tmp <- NULL
  V.tmp <- NULL
  Z.tmp <- NULL

  repeat{
    ### X.3
    X.3 <- rbinom(n, 1, p.X.3)

    ### V.3
    p.V.3 <- w.V.3 * X.3 + (1 - w.V.3) * (1 - X.3)
    V.3 <- rbinom(n, 1, p.V.3)

    ### (X.1, V.1, X.2, V.2)
    n.X.3 <- sum(X.3)
    X.V.1.2 <- matrix(0.0, nrow = n, ncol = 4)
    if(n.X.3 > 0){
      X.V.1.2[X.3 == 1,] <- mvrnorm(n.X.3, tau.1, Sigma.1)
    }
    if(n.X.3 < n){
      X.V.1.2[X.3 == 0,] <- mvrnorm(n - n.X.3, tau.0, Sigma.0)
    }

    ### X and V
    X <- cbind(X.V.1.2[, c(1, 3)], X.3)
    V <- cbind(X.V.1.2[, c(2, 4)], V.3)
    attr(X, "dimnames") <- attr(V, "dimnames") <- NULL

    ### Z
    if(!Z.given.X){
      ### Z is not dependent on X
      Z <- rep(c(1, 0), c(n.trt, n.ctr))
      break
    } else{
      ### Z is dependent on X
      e.X.a <- 1.0 / (1.0 + exp(-(cbind(1, X) %*% a)))
      Z <- rbinom(n, 1, e.X.a)

      if(!fixed.n){
        ### n.trt and n.ctr need not to be as requested
        break
      } else{
        ### n.trt and n.ctr need to be as requested
        id.tmp <- Z == 0
        n.ctr.Z <- sum(id.tmp)
        n.trt.Z <- n - n.ctr.Z

        ### For control arm
        if(n.ctr.tmp > 0){
          if(n.ctr.Z > n.ctr.tmp){
            X.tmp <- rbind(X.tmp, X[which(id.tmp)[1:n.ctr.tmp],])
            V.tmp <- rbind(V.tmp, V[which(id.tmp)[1:n.ctr.tmp],])
            Z.tmp <- c(Z.tmp, Z[which(id.tmp)[1:n.ctr.tmp]])
            n.ctr.tmp <- 0
          } else{
            X.tmp <- rbind(X.tmp, X[id.tmp,])
            V.tmp <- rbind(V.tmp, V[id.tmp,])
            Z.tmp <- c(Z.tmp, Z[id.tmp])
            n.ctr.tmp <- n.ctr.tmp - n.ctr.Z
          }
        }

        ### For treatment arm
        id.tmp <- !id.tmp
        if(n.trt.tmp > 0){
          if(n.trt.Z > n.trt.tmp){
            X.tmp <- rbind(X.tmp, X[which(id.tmp)[1:n.trt.tmp],])
            V.tmp <- rbind(V.tmp, V[which(id.tmp)[1:n.trt.tmp],])
            Z.tmp <- c(Z.tmp, Z[which(id.tmp)[1:n.trt.tmp]])
            n.trt.tmp <- 0
          } else{
            X.tmp <- rbind(X.tmp, X[id.tmp,])
            V.tmp <- rbind(V.tmp, V[id.tmp,])
            Z.tmp <- c(Z.tmp, Z[id.tmp])
            n.trt.tmp <- n.trt.tmp - n.trt.Z
          }
        }

        if(n.ctr.tmp == 0 & n.trt.tmp == 0){
          ### n.trt and n.ctr are reached to the numbers as requested
          X <- X.tmp
          V <- V.tmp
          Z <- Z.tmp
          break
        }
      } # End of fixed.n
    } # End of Z.given.X
  } # End of repeat
  
  ### Y
  err <- rnorm(n, sd = err.sd)
  Y <- cbind(1, X, Z) %*% nu + V %*% xi + err 

  ### Return ctr and trt
  X.combine <- cbind(X[, 3], V[, 3], X[, 1:2], V[, 1:2])
  id <- Z == 0
  da.ctr <- list(y = Y[id], X = X.combine[id,], err = err[id], p.bin = 2)
  da.trt <- list(y = Y[!id], X = X.combine[!id,], err = err[!id], p.bin = 2)
  ret <- list(da.ctr = da.ctr, da.trt = da.trt)

  ret 
} # End of generate.continuous.LD04m().


drop.continuous <- function(da, n.drop.bin = 1, n.drop.cont = 1,
    id.drop.bin = NULL, id.drop.cont = NULL){
  id.drop <- NULL

  ### Check binary covariates
  if(da$p.bin > 0){
    if(is.null(id.drop.bin)){
      if(da$p.bin < n.drop.bin){
        stop("p.bin is too small")
      } else{
        if(n.drop.bin > 0){
          id.drop <- 1:n.drop.bin
          p.bin <- da$p.bin - n.drop.bin
        }
      }
    } else{
      n.drop.bin <- length(id.drop.bin)
      if(da$p.bin < n.drop.bin){
        stop("p.bin is too small")
      } else{
        ### TODO: need to check the range of id.drop.bin
        id.drop <- id.drop.bin
        p.bin <- da$p.bin - n.drop.bin
      }
    }
  }

  ### Check continuous covariates
  p.cont <- ncol(da$X) - da$p.bin
  if(is.null(id.drop.cont)){
    if(p.cont < n.drop.cont){
      stop("p.cont is too small")
    } else{
      if(n.drop.cont > 0){
        id.drop <- c(id.drop, (1:n.drop.cont) + da$p.bin)
      }
    }
  } else{
    ### TODO: need to check the range of id.drop.cont
    id.drop <- c(id.drop, id.drop.cont + da$p.bin)
  }

  ### Drop covariates and return
  X <- as.matrix(da$X[, -id.drop])
  ret <- list(y = da$y, X = X, err = da$err, p.bin = p.bin)
  ret
} # End of drop.continuous.LD04m().


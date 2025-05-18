library(MASS)

fit.continuous.control.mixed <- function(da.ctr.mix){
  X <- data.frame(da.ctr.mix$X)
  if(da.ctr.mix$p.bin > 0){
    for(i.col in 1:da.ctr.mix$p.bin){
      X[, i.col] <- factor(X[, i.col], c(0, 1), c(0, 1))
    }
  }
  da.combine <- cbind(y = da.ctr.mix$y, X)
  # fm <- paste("y~", paste("X", 1:ncol(X), sep = "", collapse="+"))
  # fm <- as.formula(fm)
  # ret <- lm(fm, data = da.combine)
  ret <- lm(y ~ ., data = da.combine)
  ret
} # End of fit.continuous.control.mixed().

generate.continuous.po.mixed <- function(n.simu.po,
    sm.ctr.mix, da.trt.mix,
    method.beta = c("common", "individual", "fixed",
                    "common.noerr", "individual.noerr", "fixed.noerr"),
    method.Sigma = c("vcov", "independent")){
  sm <- sm.ctr.mix

  po.X <- da.trt.mix$X
  n <- nrow(po.X)
  p <- ncol(po.X) + 1
  po.X.intercept <- cbind(rep(1, n), po.X)

  mu.beta <- sm$mu.beta
  RSS <- sm$RSS
  df.sd <- sm$df.sd
  Sigma <- sm$Sigma
  sd <- sm$sd

  if(length(mu.beta) != p){
    stop("length(mu.beta) != p")
  }

  if(method.Sigma[1] == "vcov"){
    Sigma.beta <- Sigma
  } else if(method.Sigma[1] == "independent"){
    Sigma.beta <- diag(p)
    diag(Sigma.beta) <- diag(Sigma)
  } else{
    stop("method.Sigma is not implemented.")
  }

  ret <- list()
  for(i in 1:n.simu.po){
    if(method.beta[1] == "common"){
      po.beta <- mvrnorm(1, mu = mu.beta, Sigma = Sigma.beta)
      po.main <- as.vector(po.X.intercept %*% po.beta)

      ### 1/sigma^2 * sum_{i=1}^n r_i^2 ~ chi^2_{n - p - 1}
      po.sd <- sqrt(RSS / rchisq(1, df = df.sd))
      po.err <- rnorm(n, mean = 0, sd = po.sd)
    } else if(method.beta[1] == "individual"){
      po.beta <- mvrnorm(n, mu = mu.beta, Sigma = Sigma.beta)
      po.main <- as.vector(rowSums(po.X.intercept * po.beta))

      po.sd <- sqrt(RSS / rchisq(n, df = df.sd))
      po.err <- rnorm(n, mean = 0, sd = po.sd)
    } else if(method.beta[1] == "fixed"){
      po.beta <- mu.beta
      po.main <- as.vector(po.X.intercept %*% po.beta)

      po.err <- rnorm(n, mean = 0, sd = sd)
    } else if(method.beta[1] == "common.noerr"){
      po.beta <- mvrnorm(1, mu = mu.beta, Sigma = Sigma.beta)
      po.main <- as.vector(po.X.intercept %*% po.beta)

      po.err <- rep(0, n)
    } else if(method.beta[1] == "individual.noerr"){
      po.beta <- mvrnorm(n, mu = mu.beta, Sigma = Sigma.beta)
      po.main <- as.vector(rowSums(po.X.intercept * po.beta))

      po.err <- rep(0, n)
    } else if(method.beta[1] == "fixed.noerr"){
      po.beta <- mu.beta
      po.main <- as.vector(po.X.intercept %*% po.beta)

      po.err <- rep(0, n)
    } else{
      stop("method.beta is not implemented.")
    }

    po.y <- as.vector(po.main + po.err)

    ret[[i]] <- list(y = po.y, beta = po.beta, err = po.err)

    if(method.beta[1] == "fixed.noerr"){
      break
    }
  }

  ret
} # End of generate.continuous.po.mixed().

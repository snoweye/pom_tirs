# setwd("C:/Users/.../")
rm(list = ls())
i.gen.da <- 1    # 1 for same covariate, 2 for different
i.xi <- 1        # 1 for strong effect to outcome, 2 for none
i.a <- 1         # 1 for strong association with Z assignment
i.npo <- 3       # 1,2,3 for # of potential outcome
set.seed(5234567 + i.gen.da * 1000 + i.xi * 100 + i.a * 10 + i.npo)

# .Random.seed.new <- .Random.seed
source("./code/continuous_outcome_ld04m.r")
source("./code/get_po_est_mcvar.r")

### Parameter
param <- list()
param$nu <- c(0.0, -1.0, 1.0, -1.0, NA)
param$a <- rep(0.0, 4)
param$xi <- rep(0.0, 3)
param$p.X.3 <- 0.2
param$w.V.3 <- 0.75
param$tau.1 <- c(1.0, 1.0, -1.0, -1.0)
param$tau.0 <- c(-1.0, -1.0, 1.0, 1.0)
param$Sigma.1 <- matrix(-0.5, nrow = 4, ncol = 4)
param$Sigma.1[2, 1] <- param$Sigma.1[1, 2] <-
  param$Sigma.1[3, 4] <- param$Sigma.1[4, 3] <- 0.5
diag(param$Sigma.1) <- 1.0
param$sd <- 0.06
param$Sigma.1 <- param$Sigma.1 * param$sd^2
param$Sigma.0 <- param$Sigma.1

### Parameter
xi <- list()
xi[["str"]] <- c(-1.0, 1.0, 1.0)
xi[["no"]] <- c(0.0, 0.0, 0.0)
a <- list()
a[["str"]] <- c(0.0, 0.6, -0.6, 0.6)

### OC parameter
oc.trt.eff <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
alpha.level <- 0.05
simu.npo <- c(5, 10, 30)
n.simu.po <- simu.npo[i.npo]
n.mc.stderr <- 100

### Generate data
n.sim <- 3000
n.ctr <- 2000
n.trt <- 200

### Main loop
# for(i.gen.da in 1:2){
  ret.oc <- NULL
  ret.bias <- NULL
  ret.mse <- NULL
  ret.pval <- NULL
  ret.trt.eff <- NULL
  ret.cicr <- NULL

  # for(i.xi in 1:length(xi)){
    param$xi <- xi[[i.xi]]

    # for(i.a in 1:length(a)){
      param$a <- a[[i.a]]

      ### Skip the RCT data generation cases because a does not have
      ### any effect on treatment assignments.
      if(i.gen.da == 1 && i.a == 2){
        next
      }

      for(i.trt.eff in 1:length(oc.trt.eff)){
        trt.eff <- oc.trt.eff[i.trt.eff]
        param$nu[5] <- trt.eff

        cat(paste(i.xi, i.a, i.trt.eff, date(), "\n"))
        ret.n.var <- 8
        ret.sim <- matrix(0, nrow = n.sim, ncol = ret.n.var)
        ret.delta <- matrix(0, nrow = n.sim, ncol = ret.n.var)
        ret.ci.lb <- matrix(0, nrow = n.sim, ncol = ret.n.var)
        ret.ci.ub <- matrix(0, nrow = n.sim, ncol = ret.n.var)

        for(i.sim in 1:n.sim){
          # .Random.seed <- .Random.seed.new

          if(i.gen.da == 1){
            ### Modified version: Z.given.X = FALSE, fixed.n = TRUE (default)
            ### Treatment assignment has nothing to do with X.
            tmp.da <- generate.continuous.LD04m(n.ctr, n.trt, param)
          } else if(i.gen.da == 2){
            ### Modified version: Z.given.X = TRUE, fixed.n = TRUE
            ### The sample size is pre-specified but the baseline distribution
            ### would not be what it should be.
            tmp.da <- generate.continuous.LD04m(n.ctr, n.trt, param,
                                                Z.given.X = TRUE)
          } else if(i.gen.da == 3){
            ### Lunceford and Davidian (2004): Z.given.X = TRUE, fixed.n = FALSE
            ### The sample size can not be pre-specified in this case.
            tmp.da <- generate.continuous.LD04m(n.ctr, n.trt, param,
                                                Z.given.X = TRUE,
                                                fixed.n = FALSE)
          } else{
            stop("method.gen.da is not found.")
          }

          # .Random.seed.new <- .Random.seed

          ### Get ctr and trt datasets
          da.ctr <- tmp.da$da.ctr
          da.trt <- tmp.da$da.trt

          ### Get po est
          sm.ctr <- get.po.sm.ctr(da.ctr)
          ret.po.vcov <- get.po.pv(da.trt, sm.ctr, n.simu.po = n.simu.po,
                                   method.beta = "common",
                                   method.Sigma = "vcov",
                                   alternative = "two.sided",
                                   n.mc.stderr = n.mc.stderr)
          ret.po.indep <- get.po.pv(da.trt, sm.ctr, n.simu.po = n.simu.po,
                                    method.beta = "common",
                                    method.Sigma = "independent",
                                    alternative = "two.sided",
                                    n.mc.stderr = n.mc.stderr)

          ### Dump out returns
          tmp.sim <- c(ret.po.vcov$trt.pv,
                       ret.po.vcov$trt.pv.1st,
                       ret.po.vcov$trt.pv.mean,
                       ret.po.vcov$trt.pv.mcvar,
                       ret.po.indep$trt.pv,
                       ret.po.indep$trt.pv.1st,
                       ret.po.indep$trt.pv.mean,
                       ret.po.indep$trt.pv.mcvar)
          ret.sim[i.sim,] <- tmp.sim

          tmp.delta <- c(ret.po.vcov$trt.eff,
                         ret.po.vcov$trt.eff.1st,
                         ret.po.vcov$trt.eff.mean,
                         ret.po.vcov$trt.eff.mcvar,
                         ret.po.indep$trt.eff,
                         ret.po.indep$trt.eff.1st,
                         ret.po.indep$trt.eff.mean,
                         ret.po.indep$trt.eff.mcvar)
          ret.delta[i.sim,] <- tmp.delta

          tmp.ci <- rbind(ret.po.vcov$trt.ci,
                          ret.po.vcov$trt.ci.1st,
                          ret.po.vcov$trt.ci.mean,
                          ret.po.vcov$trt.ci.mcvar,
                          ret.po.indep$trt.ci,
                          ret.po.indep$trt.ci.1st,
                          ret.po.indep$trt.ci.mean,
                          ret.po.indep$trt.ci.mcvar)
          ret.ci.lb[i.sim,] <- tmp.ci[, 1]
          ret.ci.ub[i.sim,] <- tmp.ci[, 2]
        } # End of i.sim

        alpha.level.adj <- rep(alpha.level, 8)
        tmp.oc <- colMeans(ret.sim < alpha.level.adj)
        ret.oc <- rbind(ret.oc, c(i.xi, i.a, trt.eff, tmp.oc))

        tmp.bias <- colMeans(ret.delta - trt.eff)
        ret.bias <- rbind(ret.bias, c(i.xi, i.a, trt.eff, tmp.bias))

        tmp.mse <- colMeans((ret.delta - trt.eff)^2)
        ret.mse <- rbind(ret.mse, c(i.xi, i.a, trt.eff, tmp.mse))

        ret.pval <- rbind(ret.pval, cbind(i.xi, i.a, trt.eff, ret.sim))
        ret.trt.eff <- rbind(ret.trt.eff, cbind(i.xi, i.a, trt.eff, ret.delta))

        tmp.cicr <- (ret.ci.lb < trt.eff) & (ret.ci.ub > trt.eff)
        tmp.cicr <- colMeans(tmp.cicr)
        ret.cicr <- rbind(ret.cicr, c(i.xi, i.a, trt.eff, tmp.cicr))

      } # End of i.trt.eff
    # } # End of i.a
  # } # End of i.xi

  tmp.name <- c("xi", "a", "trt.eff",
                "vcov.med", "vcov.1st", "vcov.mean", "vcov.mcvar",
                "indep.med", "indep.1st", "indep.mean","indep.mcvar")
  tmp.name.d <- c("xi", "a", "trt.eff",
                  "vcov.med", "vcov.1st", "vcov.mean", "vcov.mcvar",
                  "indep.med", "indep.1st", "indep.mean", "indep.mcvar")
  colnames(ret.oc) <- tmp.name
  colnames(ret.bias) <- tmp.name.d
  colnames(ret.mse) <- tmp.name.d
  colnames(ret.cicr) <- tmp.name.d
  cat("OC:\n")
  print(ret.oc)
  cat("Bias:\n")
  print(ret.bias)
  cat("MSE:\n")
  print(ret.mse)
  cat("CI coverage rate:\n")
  print(ret.cicr)

  colnames(ret.pval) <- tmp.name
  colnames(ret.trt.eff) <- tmp.name.d

  fn <- paste("./results_oc_pval/08s0-oc_",
              "genda_", i.gen.da, "_xi_", i.xi, "_a_", i.a,
              "_po_", i.npo, ".rda",
              sep = "")
  save(list = c("ret.oc", "ret.bias", "ret.mse", "ret.cicr",
                "ret.pval", "ret.trt.eff",
                "ret.sim", "ret.delta", "ret.ci.lb", "ret.ci.ub",
                "i.gen.da", "i.xi", "i.a", "i.npo",
                "param", "xi", "a", "oc.trt.eff", "simu.npo",
                "n.sim", "n.ctr", "n.trt", "i.gen.da", "n.simu.po",
                "n.mc.stderr"),
       file = fn)

# } # End i.gen.da


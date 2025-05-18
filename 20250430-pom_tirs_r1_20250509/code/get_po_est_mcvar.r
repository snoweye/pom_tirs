### https://www.stata-journal.com/sjpdf.html?articlenum=st0235

source("./code/continuous_outcome_mixed.r")

get.po.sm.ctr <- function(da.ctr){
  ### Build model from ctr dataset
  m.ctr <- fit.continuous.control.mixed(da.ctr)
  sm <- summary(m.ctr)
  sm.ctr <- list(mu.beta = as.vector(sm$coefficients[, 1]),
                 RSS = sum(sm$residuals^2),
                 df.sd = m.ctr$df.residual,
                 Sigma = vcov(sm),
                 sd = sm$sigma)
  class(sm.ctr) <- "sm.ctr"
  sm.ctr
} # End of get.po.sm.ctr().


get.po.pv <- function(da.trt, sm.ctr,
    n.simu.po = 10, method.beta = "common", method.Sigma = "vcov",
    alternative = c("greater", "two.sided", "less"),
    mu.null = 0, var.equal = FALSE,
    conf.level = 0.95, n.mc.stderr = 100, ...){
  ### Get potential outcome from the model
  da.po <- generate.continuous.po.mixed(n.simu.po, sm.ctr,
                                        da.trt,
                                        method.beta = method.beta[1],
                                        method.Sigma = method.Sigma[1])

  ### Get po est
  tmp <- lapply(1:n.simu.po,
                function(i.po){
                  mean(da.trt$y - da.po[[i.po]]$y, na.rm = TRUE)
                })
  tmp <- do.call("c", tmp)
  trt.eff <- median(tmp)
  trt.eff.1st <- tmp[1]
  trt.eff.mean <- mean(tmp)

  ### Get po pv via paired t-test
  tmp.pt <- lapply(1:n.simu.po,
                   function(i.po){
                     tmp <- t.test(da.trt$y, da.po[[i.po]]$y,
                                   alternative = alternative[1],
                                   mu = mu.null, paired = TRUE, var.equal = var.equal,
                                   conf.level = conf.level, ...)
                     c(tmp$p.value, tmp$conf.int)      
                   })
  tmp.pt <- do.call("rbind", tmp.pt)
  tmp.conf.int <- tmp.pt[, 2:3]

  tmp <- tmp.pt[, 1]
  trt.pv <- median(tmp)
  trt.pv.1st <- tmp[1]
  trt.pv.mean <- mean(tmp)

  trt.ci <- c(median(tmp.conf.int[, 1]), median(tmp.conf.int[, 2]))
  trt.ci.1st <- tmp.conf.int[1, ]
  trt.ci.mean <- c(mean(tmp.conf.int[, 1]), mean(tmp.conf.int[, 2]))

  ### Get po pv via mcvar
  tmp <- lapply(1:n.simu.po,
                function(i.po){
                  da.trt$y - da.po[[i.po]]$y
                })
  tmp.trt.eff <- do.call("c", tmp)

  ### Get estimated stderr
  tmp.y.delta <- NULL
  for(i in 1:n.mc.stderr){
    da.null <- generate.continuous.po.mixed(1, sm.ctr,
                                            da.trt,
                                            method.beta = method.beta[1],
                                            method.Sigma = method.Sigma[1])
    da.po.null <- generate.continuous.po.mixed(n.simu.po, sm.ctr,
                                               da.trt,
                                               method.beta = method.beta[1],
                                               method.Sigma = method.Sigma[1])
    tmp <- lapply(1:n.simu.po,
                  function(i.po){
                    mean(da.null[[1]]$y - da.po.null[[i.po]]$y, na.rm = TRUE)
                  })
    y.delta <- do.call("c", tmp)

    tmp.y.delta <- c(tmp.y.delta, mean(y.delta))
  }
  stderr.mcvar <- sqrt(sum((tmp.y.delta - mu.null)^2) / (n.mc.stderr - 1))

  ### Get pval
  trt.eff.mcvar <- mean(tmp.trt.eff)
  stat <- (trt.eff.mcvar - mu.null) / stderr.mcvar
  df <- length(da.trt$y) - 1
  if(alternative[1] == "greater"){
    pv <- pt(stat, df, lower.tail = FALSE)
  } else if(alternative[1] == "two.sided"){
    pv <- pt(abs(stat), df, lower.tail = FALSE)
    pv <- 2 * pv
  } else if(alternative[1] == "less"){
    pv <- pt(stat, df)
  } else{
    stop("Alternative is not found.")
  }
  trt.pv.mcvar <- pv

  ### Get CI
  qt.df.alpha.2 <- qt((1 - conf.level) / 2, df, lower.tail = FALSE) * c(-1, 1)
  trt.ci.mcvar <- trt.eff.mcvar + qt.df.alpha.2 * stderr.mcvar

  ### For return
  ret <- list(trt.pv = trt.pv,                # median of K paired t-tests
              trt.pv.1st = trt.pv.1st,        # single paired t-test
              trt.pv.mean = trt.pv.mean,      # mean of K paired t-tests
              trt.pv.mcvar = trt.pv.mcvar,    # test with bootstrap se.
              trt.eff = trt.eff,
              trt.eff.1st = trt.eff.1st,
              trt.eff.mean = trt.eff.mean,
              trt.eff.mcvar = trt.eff.mcvar,
              trt.ci = trt.ci,
              trt.ci.1st = trt.ci.1st,
              trt.ci.mean = trt.ci.mean,
              trt.ci.mcvar = trt.ci.mcvar)
  ret
} # End of get.po.pv().



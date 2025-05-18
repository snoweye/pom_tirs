# setwd("C:/Users/.../")
rm(list = ls())
library(ggplot2)

pn.o <- "./"
pn.d <- "./"

npo.all <- 3
npo.all.list <- c(1, 5, 10, 30)

n.trt.all <- 3
n.trt.all.list <- c(50, 100, 200)

da.org <- NULL
for(i.genda in 1:2){
  for(i.xi in 1:2){
    for(i.a in 1){
      for(i.npo in 1:npo.all){
        for(i.n.trt in 1:n.trt.all){
          if(i.genda == 1 && i.a == 2){
            next
          }

          n.trt <- n.trt.all.list[i.n.trt]

          if(n.trt == 200){
            fn <- paste(pn.o, "./results_oc_pval/08s0-oc_genda_", i.genda,
                        "_xi_", i.xi, "_a_", i.a, "_po_", i.npo,
                        ".rda", sep = "")
          } else{
            fn <- paste(pn.o, "./results_oc_pval_ntrt/08s0-oc_genda_", i.genda,
                        "_xi_", i.xi, "_a_", i.a, "_po_", i.npo,
                        "_ntrt_", n.trt, ".rda", sep = "")
          }
          if(!file.exists(fn)){
            next
          }
          load(fn)

          if(i.genda == 1 && i.a == 1){
            ret.bias[, 2] <- 0  ### a
          }

          tmp.bias <- cbind(rep(i.genda, nrow(ret.bias)),
                            ret.bias,
                            n.simu.po,
                            n.trt)
          da.org <- rbind(da.org, tmp.bias)
        }
      }
    }
  }
}
colnames(da.org) <- c("genda", colnames(ret.bias), "npo", "n.trt")

da <- as.data.frame(da.org)
da$genda <- factor(da$genda, labels = c("Z random", "Z|U,a"))
da$xi <- factor(da$xi, labels = c("xi str"))
da <- da[da$xi == "xi str",]


da.g.vcov <- NULL
for(i.genda in levels(da$genda)){
  for(i.xi in levels(da$xi)){
    for(i.npo in npo.all.list){
      for(i.n.trt in 1:n.trt.all){

        n.trt <- n.trt.all.list[i.n.trt]

        if(i.npo == 1){
          id <- da$genda == i.genda &
                da$xi == i.xi &
                da$npo == 5 &
                da$n.trt == n.trt
          da.g.vcov <- rbind(da.g.vcov,
                             cbind(da[id, 1:4],
                                   val = da$vcov.1st[id],
                                   npo = rep(1, sum(id)),
                                   n.trt = rep(n.trt, sum(id)),
                                   method = rep("1st", sum(id))))
        } else{
          id <- da$genda == i.genda &
                da$xi == i.xi &
                da$npo == i.npo &
                da$n.trt == n.trt
          da.g.vcov <- rbind(da.g.vcov,
                             cbind(da[id, 1:4],
                                   val = da$vcov.med[id],
                                   npo = rep(i.npo, sum(id)),
                                   n.trt = rep(n.trt, sum(id)),
                                   method = rep("median", sum(id))))
          da.g.vcov <- rbind(da.g.vcov,
                             cbind(da[id, 1:4],
                                   val = da$vcov.mean[id],
                                   npo = rep(i.npo, sum(id)),
                                   n.trt = rep(n.trt, sum(id)),
                                   method = rep("mean", sum(id))))
          da.g.vcov <- rbind(da.g.vcov,
                             cbind(da[id, 1:4],
                                   val = da$vcov.mcvar[id],
                                   npo = rep(i.npo, sum(id)),
                                   n.trt = rep(n.trt, sum(id)),
                                   method = rep("mcvar", sum(id))))
        }
      }
    }
  }
}
da.g.vcov$val <- da.g.vcov$val * 100


width <- 2400 * 1.5
height <- 1200 * 1.0 * 3
res <- 400 * 1.5

### Plot
da.g1.vcov <- da.g.vcov[da.g.vcov$method %in% c("1st", "mcvar"),]
da.g1.vcov$K <- factor(da.g1.vcov$npo)
da.g1.vcov$method.var <- factor("vcov", levels = c("vcov"))
da.g1.vcov$n.trt <- paste("n.trt = ", da.g1.vcov$n.trt, sep = "")
da.g1.vcov$n.trt <- factor(da.g1.vcov$n.trt,
                           levels = paste("n.trt = ", n.trt.all.list,
                                          sep = ""))
da.g1 <- rbind(da.g1.vcov)

xlim <- range(da.g1$trt.eff)
ylim <- range(da.g1$val)

g1 <- ggplot(data = da.g1,
             aes(x = trt.eff, y = val, group = K)) +
        geom_line(aes(linetype = K, color = K)) +
        geom_point(aes(shape = K, color = K)) +
        geom_hline(yintercept = 0.0, linetype = "dashed") +
        scale_shape_manual(values = c(0, 4, 2, 5, 3, 6, 8, 9)) +
        ylim(ylim[1], ylim[2]) +
        xlim(xlim[1], xlim[2]) +
        xlab("Treatment Effect") +
        ylab("Bias * 100") +
        theme_bw() +
        # theme(legend.position = "bottom") +
        # facet_grid(method.var ~ genda + xi)
        facet_grid(n.trt ~ genda, scales = "free_y")
 png(paste(pn.d, "./paper_supp_plots_v4/Rplot_ntrt_cont_s0_comb_bias_mcvar.png", sep = ""),
     width = width, height = height, res = res)
print(g1)
 dev.off()


# setwd("C:/Users/.../")
rm(list = ls())
library(ggplot2)

pn.o <- "./"
pn.d <- "./"

npo.all <- 3
npo.all.list <- c(1, 5, 10, 30)

da.org <- NULL
for(i.genda in 1:2){
  for(i.xi in 1:2){
    for(i.a in 1){
      for(i.npo in 1:npo.all){
        if(i.genda == 1 && i.a == 2){
          next
        }

        fn <- paste(pn.o, "./results_oc_pval_noise/08s0-oc_genda_", i.genda,
                    "_xi_", i.xi, "_a_", i.a, "_po_", i.npo, ".rda", sep = "")
        if(!file.exists(fn)){
          next
        }
        load(fn)

        if(i.genda == 1 && i.a == 1){
          ret.cicr[, 2] <- 0  ### a
        }

        tmp.cicr <- cbind(rep(i.genda, nrow(ret.cicr)),
                          ret.cicr,
                          n.simu.po)
        da.org <- rbind(da.org, tmp.cicr)
      }
    }
  }
}
colnames(da.org) <- c("genda", colnames(ret.cicr), "npo")

da <- as.data.frame(da.org)
da$genda <- factor(da$genda, labels = c("Z random", "Z|U,a"))
da$xi <- factor(da$xi, labels = c("xi str"))
da <- da[da$xi == "xi str",]


da.g.vcov <- NULL
for(i.genda in levels(da$genda)){
  for(i.xi in levels(da$xi)){
    for(i.npo in npo.all.list){
      if(i.npo == 1){
        id <- da$genda == i.genda &
              da$xi == i.xi &
              da$npo == 5
        da.g.vcov <- rbind(da.g.vcov,
                           cbind(da[id, 1:4],
                                 val = da$vcov.1st[id],
                                 npo = rep(1, sum(id)),
                                 method = rep("1st", sum(id))))
      } else{
        id <- da$genda == i.genda &
              da$xi == i.xi &
              da$npo == i.npo 
        da.g.vcov <- rbind(da.g.vcov,
                           cbind(da[id, 1:4],
                                 val = da$vcov.med[id],
                                 npo = rep(i.npo, sum(id)),
                                 method = rep("median", sum(id))))
        da.g.vcov <- rbind(da.g.vcov,
                           cbind(da[id, 1:4],
                                 val = da$vcov.mean[id],
                                 npo = rep(i.npo, sum(id)),
                                 method = rep("mean", sum(id))))
        da.g.vcov <- rbind(da.g.vcov,
                           cbind(da[id, 1:4],
                                 val = da$vcov.mcvar[id],
                                 npo = rep(i.npo, sum(id)),
                                 method = rep("mcvar", sum(id))))
      }
    }
  }
}
da.g.vcov$val <- da.g.vcov$val


width <- 2400 * 1.5
height <- 1200 * 1.5
res <- 400 * 1.5

### Plot
da.g1.vcov <- da.g.vcov[da.g.vcov$method %in% c("1st", "mcvar"),]
da.g1.vcov$K <- factor(da.g1.vcov$npo)
da.g1.vcov$method.var <- factor("vcov", levels = c("vcov"))
da.g1 <- rbind(da.g1.vcov)

xlim <- range(da.g1$trt.eff)
ylim <- range(da.g1$val)

g1 <- ggplot(data = da.g1,
             aes(x = trt.eff, y = val, group = K)) +
        geom_line(aes(linetype = K, color = K)) +
        geom_point(aes(shape = K, color = K)) +
        geom_hline(yintercept = 0.0, linetype = "dashed") +
        scale_shape_manual(values = c(0, 4, 2, 5, 3, 6, 8, 9)) +
        ylim(0.85, 1.0) +
        xlim(xlim[1], xlim[2]) +
        xlab("Treatment Effect") +
        ylab("95% CI Coverage Rate") +
        theme_bw() +
        # theme(legend.position = "bottom") +
        # facet_grid(method.var ~ genda + xi)
        facet_grid( ~ genda, scales = "free_y")
 png(paste(pn.d, "./paper_supp_plots_v4/Rplot_noise_cont_s0_comb_cicr_mcvar.png", sep = ""),
     width = width, height = height, res = res)
print(g1)
 dev.off()


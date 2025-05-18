# setwd("C:/Users/.../")
rm(list = ls())
library(ggplot2)
library(qqplotr)

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

        fn <- paste(pn.o, "./results_oc_pval/08s0-oc_genda_", i.genda,
                    "_xi_", i.xi, "_a_", i.a, "_po_", i.npo, ".rda", sep = "")
        if(!file.exists(fn)){
          next
        }
        load(fn)

        if(i.genda == 1 && i.a == 1){
          ret.pval[, 2] <- 0  ### a
        }

        ret.pval <- ret.pval[ret.pval[, 3] == 0,]   ### only the null

        tmp.pval <- cbind(rep(i.genda, nrow(ret.pval)),
                          ret.pval,
                          n.simu.po)
        da.org <- rbind(da.org, tmp.pval)
      }
    }
  }
}
colnames(da.org) <- c("genda", colnames(ret.pval), "npo")

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


width <- 2400 * 1.5
height <- 1200 * 1.0 * 4
res <- 400 * 1.5

### Plot
da.g1.vcov <- da.g.vcov[da.g.vcov$method %in% c("1st", "mcvar"),]
da.g1.vcov$K <- factor(da.g1.vcov$npo)
da.g1.vcov$method.var <- factor("vcov", levels = c("vcov"))
da.g1 <- rbind(da.g1.vcov)

g1 <- ggplot(data = da.g1,
             mapping = aes(sample = val,
                           group = K,
                           # shape = K,
                           color = K)) +
        stat_qq_point(distribution = "unif") +
        geom_abline(slope = 1, linetype = 3) +
        scale_shape_manual(values = c(0, 4, 2, 5, 3, 6, 8, 9)) +
        ylim(0, 1) +
        xlim(0, 1) +
        xlab("Theoretical Quantiles (Uinform(0, 1))") +
        ylab("Sample Quantiles (p-values)") +
        theme_bw() +
        # theme(legend.position = "bottom") +
        # facet_grid(method.var ~ genda + xi)
        facet_grid(K ~ genda, scales = "free_y")
 png(paste(pn.d, "./paper_plots_v4/Rplot_cont_s0_comb_qqplot_mcvar.png", sep = ""),
     width = width, height = height, res = res)
print(g1)
 dev.off()


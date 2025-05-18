.pn <- "./paper_supp_scripts_v4/"
# .pn <- "./"


.fn.all <- list.files(.pn, "08s0-ntrt-plot")

for(i.fn in 1:length(.fn.all)){
  source(paste(.pn, .fn.all[i.fn], sep = ""))
}


.fn.all <- list.files(.pn, "08s0-noise-plot")

for(i.fn in 1:length(.fn.all)){
  source(paste(.pn, .fn.all[i.fn], sep = ""))
}


.fn.all <- list.files(.pn, "08s0-ps-plot")

for(i.fn in 1:length(.fn.all)){
  source(paste(.pn, .fn.all[i.fn], sep = ""))
}


.pn <- "./paper_scripts_v4/"
# .pn <- "./"


.fn.all <- list.files(.pn, "08s0-plot")

for(i.fn in 1:length(.fn.all)){
  source(paste(.pn, .fn.all[i.fn], sep = ""))
}


.fn.all <- list.files(.pn, "08s0-supp-plot")

for(i.fn in 1:length(.fn.all)){
  source(paste(.pn, .fn.all[i.fn], sep = ""))
}

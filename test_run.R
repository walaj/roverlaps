library(roverlaps)
library(gUtils)

k=1e6
verbose=FALSE
o1 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)
k=1e8
o2 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)

for (i in 1:5) {
  print(i)
  system.time(o <- roverlaps(o1, o2, robust=robust, cores=cores, verbose=verbose, index_only=index_only))
  gc()
  Sys.sleep(15)
  system.time(o <- gr.findoverlaps(o1, o2, verbose=verbose, max.chunk=1e16))
  gc()
  Sys.sleep(15)
}



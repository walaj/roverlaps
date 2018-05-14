library(roverlaps)
library(gUtils)

cores=1
verbose=FALSE
index_only=FALSE

k=1e6
o1 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)

k=2e7
o2 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)
gc()

for (i in 1:5) {
  #print(paste("...running roverlaps #", i))
  write(system.time(o <- roverlaps(o1, o2, cores=cores, verbose=verbose, index_only=index_only)), stdout())
  gc()
  #print("...")
  #Sys.sleep(15)
  #print(paste("...running gr.findoverlaps #",i))
  #write(system.time(o <- gr.findoverlaps(o1, o2, verbose=verbose, max.chunk=1e16)), stdout())
  #gc()
  #print("...")
  #Sys.sleep(15)
}

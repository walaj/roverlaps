context("test-roverlaps.R")
library(data.table)
test_that("unit test 1", {
  k=1e6
  verbose=TRUE
  cores=1
  index_only=FALSE
  o1 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)
  k=1e6
  o2 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)
  system.time(o <- roverlaps(o1, o2, robust=robust, cores=cores, verbose=verbose, index_only=index_only))
  system.time(o <- gr.findoverlaps(o1, o2, verbose=TRUE, max.chunk=1e16))
  rm(d2)
  expect_equal(as.character(o$seqnames), as.character("1"))
  expect_equal(nrow(o), 1)
  expect_equal(o$start, 5)
  expect_equal(o$end, 9)
  expect_equal(o$query.id, 1)
  expect_equal(o$subject.id, 2)
})

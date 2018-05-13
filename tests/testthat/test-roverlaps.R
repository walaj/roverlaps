context("test-roverlaps.R")

test_that("unit test 1", {
  d1 <- data.table(seqnames=c(1, "X"), start=c(1,10), end=c(10,20))
  d2 <- data.table(seqnames=c("Y",1), start=c(10,5), end=c(20,9))
  o <- roverlaps(d1, d2)
  expect_equal(as.character(o$seqnames), as.character("1"))
  expect_equal(nrow(o), 1)
  expect_equal(o$start, 5)
  expect_equal(o$end, 9)
  expect_equal(o$query.id, 1)
  expect_equal(o$subject.id, 2)
})

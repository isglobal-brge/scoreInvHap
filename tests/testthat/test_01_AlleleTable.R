context("Create Allele Table")

map <- data.frame(position = 1:3, chromosome = rep("chr7", 3), allele.1 = c("A", "C", "A"),
                  allele.2 = c("T", "A", "G"), names = c("SNP1", "SNP2", "SNP3"))
rownames(map) <- map$names

test_that("check function", {
  alleleTable <- getAlleleTable(map)
  expect_equal(ncol(alleleTable), 4)
  expect_equal(colnames(alleleTable), c("AA", "AB", "BB", "Fail"))
  expect_equal(alleleTable[, 1], c("AA", "CC", "AA"))
  expect_equal(alleleTable[, 2], c("AT", "AC", "AG"))
  expect_equal(alleleTable[, 3], c("TT", "AA", "GG"))
  expect_equal(alleleTable[, 4], c("NN", "NN", "NN"))
})

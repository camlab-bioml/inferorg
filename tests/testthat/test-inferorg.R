test_that("inferorg handles input correctly", {

  expect_error(inferorg(NA))

})

test_that("organism guessed correctly", {
  genes_human_symbol <- c("HLA-A", "HLA-B", "HLA-C")
  genes_mouse_ensembl <- c("ENSMUSG00000026171", "ENSMUSG00000033276", "ENSMUSG00000033257", "ENSMUSG00000026170", "ENSMUSG00000100827", "ENSMUSG00000100553", "ENSMUSG00000098611", "ENSMUSG00000100534")

  io1 <- inferorg(genes_human_symbol)
  io2 <- inferorg(genes_mouse_ensembl)

  expect_equal(io1$organism, 'human')
  expect_equal(io1$format, 'symbol')

  expect_equal(io2$organism, 'mouse')
  expect_equal(io2$format, 'ensgene')
})

test_that("autoconvert works correctly", {
  expect_equal(autoconvert(c("madeup", "gene", "name")), c(NA, NA, NA))
})

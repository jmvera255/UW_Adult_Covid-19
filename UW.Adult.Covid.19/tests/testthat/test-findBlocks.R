test_that("findBlocks works", {
    protein.df = data.frame(Protein = c("A", "A"), Pos = c(1,3), Pvalue = c(1,1), logFC = c(0,0), Mean.Signal = c(1,1))
    blocks = findBlocks(protein.df)
    eids = getEpitopeID(blocks$Protein, blocks$Start, blocks$Stop)
  expect_equal(eids, c("A_1_1", "A_3_3"))
})

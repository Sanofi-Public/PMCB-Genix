test_that("test genix", {
  # get data
  data("genixDemoData")
  set.seed(123)
  # test constructNets
  rslts <- constructNets(gex=genixDemoData$gex, 
                         glasso.rho=1.0, 
                         metadata=genixDemoData$metadata, 
                         field="group",
                         scale.gex=TRUE, qc.gex=c(0.01,0.5), maxit=1e3)
  from_name_expected <- c("RPL22", "ENO1", "EFHD2")
  to_name_expected <- c("RPS8", "EEF1A1", "RPL5")
  genix_observed <- igraph::as_long_data_frame(rslts$grp1@icm_graph)[c(2,28,32), c("from_name", "to_name")]
  expect_equal(from_name_expected, genix_observed$from_name)
  expect_equal(to_name_expected, genix_observed$to_name)
  # test compileNets
  cmpl_rslts <- compileNets(grph=rslts, degree.th=3, hubs.th=0.01, 
                            minModuleSize=7, deepSplit=TRUE)
  hubs_exp <- c("FTH1", "MT-ND5", "RPLP0")
  hubs_obs <- names(cmpl_rslts$grp1@hubs)[c(2,28,32)]
  expect_equal(hubs_exp, hubs_obs)
  # test compareNets
  cmpr_rslts <- compareNets(grph_1=cmpl_rslts$grp1, 
                            grph_2=cmpl_rslts$grp2, 
                            n.perm=25)
  expect_equal(0.2179, round(cmpr_rslts@ji, 4))
  expect_equal(91, cmpr_rslts@inAB_mtx[1,2])
  expect_equal(47, cmpr_rslts@inAB_mtx[2,1])
  expect_equal(583, cmpr_rslts@inAB_mtx[2,2])
})
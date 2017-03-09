library(SNSequate)
context("Equated values for different models")

tol <- 0.7

test_that("EG design", {
  data(Math20EG)
  mod <- ker.eq(scores=Math20EG,kert="gauss",degree=c(2,3),design="EG")
  
  expect_equal(mod$nx, 1453)
  expect_equal(mod$ny, 1455)
  
  # 6 first and 3 last
  rj <- c(3.30, 6.44, 11.77, 20.17, 32.43, 48.89, 28.93, 17.71, 10.16)
  expect_equal(mod$rj[c(1:6, 19:21)]*mod$nx, rj, tolerance=tol, check.attributes=F)
  
  sk <- c(1.71, 3.77, 7.65, 14.24, 24.44, 38.75, 47.73, 34.54, 24.18)
  expect_equal(mod$sk[c(1:6, 19:21)]*mod$ny, sk, tolerance=tol, check.attributes=F)
  
  pre10_yx <- c(0.01,0.01,0.02,0.04,0.07,0.13,0.22,0.33,0.48,0.67)
  pre10_xy <- c(-0.01,-0.02,-0.05,-0.11,-0.21,-0.37,-0.58,-0.88,-1.25,-1.71)
  prep10 <- SNSequate::PREp(mod, 10)
  expect_equal(prep10$preXy, pre10_xy, tolerance=tol, check.attributes=F, scale=pre10_xy)
  expect_equal(prep10$preYx, pre10_yx, tolerance=tol, check.attributes=F, scale=pre10_yx)
  
  see_xy <- c(0.145,0.225,0.275,0.279,0.199,0.210,0.205,0.144)
  expect_equal(mod$SEEXy[c(1:4, 18:21)], see_xy, tolerance=tol, check.attributes=F, scale=see_xy)
  
  see_yx <- c(0.220,0.289,0.287,0.266,0.199,0.170,0.119,0.070)
  expect_equal(mod$SEEYx[c(1:4, 18:21)], see_yx, tolerance=tol, check.attributes=F, scale=see_yx)
})


test_that("SG design", {
  data(Math20SG)
  mod <- ker.eq(scores=Math20SG,kert="gauss",degree=c(3,3,1,1),design="SG")
  
  expect_equal(mod$nx, 1453)
  expect_equal(mod$ny, 1453)
  
  rj <- c(2.3, 5.17, 10.47, 19.22,  44.93, 30.39, 19.42, 11.67)
  expect_equal(mod$rj[c(1:4,18:21)]*mod$nx, rj, tolerance=tol, scale=rj)
  
  sk <- c(2.29, 5.27, 10.86, 20.31,  32.25, 18.83, 9.93, 4.67)
  expect_equal(mod$sk[c(1:4,18:21)]*mod$ny, sk, tolerance=tol, scale=sk)
  
  pre10_yx <- c(-0.0031, -0.0133, -0.0332, -0.0701, -0.1333, -0.2330, -0.3793, -0.5817, -0.8481, -1.1851)
  pre10_xy <- c(0.0007, 0.0059, 0.0148, 0.0309, 0.0590, 0.1042, 0.1714, 0.2654, 0.3900, 0.5485)
  prep10 <- SNSequate::PREp(mod, 10)
  expect_equal(prep10$preXy, pre10_xy, tolerance=tol, check.attributes=F)
  expect_equal(prep10$preYx, pre10_yx, tolerance=tol, check.attributes=F, scale=pre10_yx)
  
  see_xy <- c(0.1617, 0.2208, 0.2208, 0.1931,  0.1404, 0.1674, 0.1854, 0.1581)
  expect_equal(mod$SEEXy[c(1:4, 18:21)], see_xy, tolerance=tol, check.attributes=F, scale=see_xy)
  
  see_yx <- c(0.1579, 0.2236, 0.2254, 0.1970,  0.1670, 0.1675, 0.1338, 0.0885)
  expect_equal(mod$SEEYx[c(1:4, 18:21)], see_yx, tolerance=tol, check.attributes=F, scale=see_yx)
})


test_that("CB design", {
#   data(CBdata)
#   mod <- ker.eq(scores=CBdata$datx1y2,scores2=CBdata$datx2y1,kert="gauss",degree=c(3,3,1,1,3,3,1,1),design="CB",
#                 J=75,K=76,wx=0.5,wy=0.5)
#   
#   expect_equal(mod$N12, 143)
#   expect_equal(mod$N21, 140)
  
  ## ... no more data in the book to compare 
})


test_that("NEAT/CE design", {  
  
})


test_that("NEAT/PSE design", {  
  load("NEAT_DATA.RData")
  mod = ker.eq(scores=XA, design="NEAT_PSE", kert="gauss", 
                scores2=YA, degreeXA=c(5,4,2,2), degreeYA = c(2,4,2,2), 
                J=21, K=21, L=11, w=0.5)
  
  #   Hardcoded directly from the output of ETS software
  eqETS <- c(0.092015296,1.005279303,1.710671782,2.403847694,3.168253899,
             4.02464056,4.956838131,5.933265686,6.923579693,7.905854702,
             8.866832733,9.800272942,10.70591831,11.58936405,12.4622221,
             13.34186268,14.25089836,15.21767712,16.28080368,17.50959015,
             19.06346893)
  
  seeETS <- c(0.194235444,0.29360047,0.300062269,0.278455943,0.248480558,
              0.217080116,0.187951937,0.164718613,0.149325773,0.141054183,
              0.137982562,0.139190942,0.145542458,0.158268914,0.176969305,
              0.199965313,0.227453128,0.264277279,0.315503955,0.365412503,
              0.311620802)
  
  expect_equal(mod$eqYx, eqETS, tolerance=tol, check.attributes=F, scale=eqETS)
  expect_equal(mod$SEEYx, seeETS, tolerance=tol, check.attributes=F, scale=seeETS)
  
})
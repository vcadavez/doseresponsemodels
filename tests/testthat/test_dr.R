test_that("Risk computation by DR function", {
  risk <-  DR(5:10, "Pouillot", 1:11)
  expect_true(is.vector(risk)==FALSE)
})

# test_that("DR for all subpopulations", {
#   # Very long test. Should not be put in the package
#   risk <-NULL
#   
#   for(Poisson in c(TRUE, FALSE)){
#     
#     if(Poisson) dose <- c(0,10^seq(-5,12,length=10))
#     else dose <- c(0,10^seq(0,12,length=10))
#     
#     risk <-  cbind(risk, DR(dose, "JEMRA", 0, Poisson = Poisson, method="cubature"))
#     risk <- cbind(risk, DR(dose, "Pouillot", 0, Poisson = Poisson, method="cubature"))
#     risk <- cbind(risk, DR(dose, "Fritsch", 1:3, Poisson = Poisson, method="cubature"))
#     risk <- cbind(risk, DR(dose, "EFSA", 0, Poisson = Poisson, method="cubature"))
#     risk <- cbind(risk, DR(dose, "EFSALV", 0, Poisson = Poisson, method="cubature"))
#     risk <- cbind(risk, DR(dose, "EFSAV", 0, Poisson = Poisson, method="cubature"))
#     risk <- cbind(risk, DR(dose, "EFSAMV", 0, Poisson = Poisson, method="cubature"))
#   }
# 
#   expect_true(!any(is.na(risk)))
# })

test_that("DRQuick for all subpopulations", {
  risk <-NULL

  for(Poisson in c(TRUE, FALSE)){
    
    if(Poisson) dose <- c(0,10^seq(-5,12,length=100)) else dose <- 10^seq(0,12,length=101)
                          
    risk <-  cbind(risk, DRQuick(dose, "JEMRA", 0:2, Poisson = Poisson))
    risk <- cbind(risk, DRQuick(dose, "Pouillot", 0:11, Poisson = Poisson))
    risk <- cbind(risk, DRQuick(dose, "Fritsch", 0:3, Poisson = Poisson))
    risk <- cbind(risk, DRQuick(dose, "EFSA", 0:14, Poisson = Poisson))
    risk <- cbind(risk, DRQuick(dose, "EFSALV", 0:14, Poisson = Poisson))
    risk <- cbind(risk, DRQuick(dose, "EFSAV", 0:14, Poisson = Poisson))
    risk <- cbind(risk, DRQuick(dose, "EFSAMV", 0:14, Poisson = Poisson))
  }
  
  expect_true(!any(is.na(risk)) & all(dim(risk) == c(101,158)))
})

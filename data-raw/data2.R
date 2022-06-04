## code to prepare `DATASET` dataset goes here

usethis::use_data(data2, overwrite = TRUE)

library(sampling)
library(dplyr)

set.seed(369)

simbaru <- function(out=0, m=30, s=100, c1=1, c2=2){

  #set.seed(369)

  Ni = rep(s, m)
  N = sum(Ni)

  # variabel penyerta
  z <- rbinom(N,1,out)
  aux.var <- (1-z)*rnorm(N,5,sqrt(1)) + z*(rnorm(N,5,sqrt(1))+10)

  n <- rbinom(N,1,0.1)

  # membangkitkan data efek acak area
  v <- (1-n)*rnorm(N,0,sqrt(3)) + n*rnorm(N,0,sqrt(30))

  # membangkitkan data efek acak individu
  e <- (1-n)*rnorm(N,0,sqrt(6)) + n*rnorm(N,0,sqrt(150))

  # membentuk code area
  code <- rep(1:m, each = s)

  # membentuk variable of interest
  y = c1 + c2*aux.var + v + e

  # menggabungkan data menjadi dataframe populasi
  pop.matrix <- cbind(y,aux.var,v,e,code)
  pop <- as.data.frame(pop.matrix)
  names(pop) <- c("y", "x","v","e", "area")

  return(pop)
}

sim1 <- simbaru(out=0.01)
sim1$area <- as.factor(sim1$area)

# skenario 1
strata1 = strata(sim1,c("area"),size=c(rep(5,15),rep(10,15)),method = "srswor")
stratified_data1 = getdata(sim1,strata1)

sim1_id <- sim1
sim1_id$ID_unit <- c(1:3000)
outs1 <- anti_join(sim1_id,stratified_data1, by="ID_unit")


data2 <- as.data.frame(cbind(x.out = outs1$x, reg.r = outs1$area))
data2[,2] <- as.factor(data2[,2])

use_data(data2, overwrite = TRUE)

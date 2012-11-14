library(MASS)
library(MuMIn)
library(nnet)

Jdata<- read.delim("http://r.789695.n4.nabble.com/file/n4636105/20120709_JLittle_data_file.txt")

model1 <- multinom(JVeg5 ~ Elevation + Lat_Y_pos + Subregion +  Wind_310 + TPI + Landform, data=Jdata, maxit=600)

dredge1 <- dredge(model1) 


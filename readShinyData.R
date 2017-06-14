#defines sliders based on most recent 18-lm configuration of 4-cell embryos
#

curvepts.2<-matrix(c(
  2,5,6,
  5,6,1,
  1,7,8,
  7,8,3,
  3,9,10,
  9,10,4,
  4,11,12,
  11,12,2
), ncol=3, byrow=T)
colnames(curvepts.2) <- c("before", "slide", "after")

curvepts.4<-matrix(c(
  4,11,1,
  1,12,2,
  2,13,3,
  3,14,7,
  7,15,10,
  10,16,9,
  9,17,8,
  8,18,4,
  12,2,13,
  14,7,15,
  16,9,17,
  18,4,11
), ncol=3, byrow=T)
colnames(curvepts.4) <- c("before", "slide", "after")


cell.link.2 <- matrix(c(2,5,
                        5,6,
                        6,1,
                        1,7,
                        7,8,
                        8,3,
                        3,9,
                        9,10,
                        10,4,
                        4,11,
                        11,12,
                        12,2,
                        1,4                        
),ncol=2,byrow=T)

cell.link.4 <- matrix(c(4,11,
                        11,1,
                        1,12,
                        12,2,
                        2,13,
                        13,3,
                        3,14,
                        14,7,
                        7,15,
                        15,10,
                        10,16,
                        16,9,
                        9,17,
                        17,8,
                        8,18,
                        18,4,
                        1,5,
                        5,8,
                        5,6,
                        6,3,
                        6,10
),ncol=2,byrow=T)

embryo2 <- read.csv("data/embryo2.csv")
embryo4 <- read.csv("data/embryo4.csv")

embryo2.3D <- arrayspecs(embryo2[,c(1:24)], 12, 2)
embryo4.3D <- arrayspecs(embryo4[,c(1:36)], 18, 2)
ID <- embryo2[,25]
devTemp <- embryo2[,26]

TA2 <- trajectory.analysis(embryo2.3D ~ ID * devTemp)
TA4 <- trajectory.analysis(embryo4.3D ~ ID * devTemp)

cells2 <- read.csv("data/cells2.csv")
cells4 <- read.csv("data/cells4.csv")

ID2 <- as.factor(paste(ID, devTemp, sep="-"))
col.gp <- rainbow(length(levels(ID2)))
names(col.gp) <- levels(ID2)
col.gp <- col.gp[match(ID2, names(col.gp))] # col.gp must NOT be a factor

GP3 <- gridPar(pt.bg = "blue", link.col = "blue", pt.size = 1, tar.pt.bg = "limegreen", tar.link.col = "limegreen") 

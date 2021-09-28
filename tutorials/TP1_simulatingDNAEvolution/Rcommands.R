## Exponential distribution for waiting times
simulatedValues=read.table("simulatedExpForR.txt", h=F)$V1

hist(simulatedValues, nc=50, freq=F)
x<-(1:600)/10
y<-dexp(x, rate=0.1)
lines(x, y, t="l", col="red", lwd=2)

## Probabilities of drawing states
simulatedValues=read.table("simulatedUnifForR.txt", h=F)$V1
table(simulatedValues) / length(simulatedValues)

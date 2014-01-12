
f <- read.table("f2.txt")
ff = (f$V1)^2 + (f$V2)^2
plot(ff[1:40], type='b')

r = (1:40)*0.2
th = 4*r^2*exp(-2*r)
lines(th,col="green")
# plot(ff[1:40], col="green", add=TRUE)

f_ini = read.table("f1.txt")
ff_ini = (f_ini$V1)^2 + (f_ini$V2)^2
# plot(ff_ini[1:40], col="red", add=TRUE)
lines(ff_ini[1:40], col="red")
points(ff_ini[1:40], col="red")

abline(v=1:40, col="gray60")

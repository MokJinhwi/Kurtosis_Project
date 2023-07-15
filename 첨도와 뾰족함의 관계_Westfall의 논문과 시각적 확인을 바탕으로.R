###### Illustrating example - cauchy ######

set.seed(12344) 
rc <- rcauchy(1e3) 
mean(rc); sd(rc) # 표본평균과 표본표준편차 구하기 
sort(rc)[1:7] # Outlier 확인
hist(rc, breaks = 1e2, xlab="x", main="", ylim=c(0,120)) # Histogram
abline(v=-sd(rc), col="blue", lty=3 ); abline(v=sd(rc), col="blue", lty=3)

Z.rc <- (rc-mean(rc))/sd(rc)
mean(Z.rc^4) # 첨도
Z.rc.outlier <- Z.rc[abs(Z.rc)>1] # Without 1 sd
Z.rc.within <- Z.rc[abs(Z.rc)<=1] # Within 1 sd
sum(Z.rc.within^4)/(1e3)
(sum(Z.rc.within^4)/(sum(Z.rc.outlier^4)+sum(Z.rc.within^4)))*100 #비율


-----------------------------------------------------------------------------------


###### Table of the central part to the kurtosis within 1 sd ######

# unif
set.seed(412)
x.1 <- runif(1e6,0,1)
Z.1 <- (x.1-mean(x.1))/sd(x.1)
mean(Z.1^4)
Z.1.outlier <- Z.1[abs(Z.1)>1]
Z.1.within <- Z.1[abs(Z.1)<=1]
sum(Z.1.within^4)/(1e6)
sum(Z.1.within^4)/(sum(Z.1.outlier^4)+sum(Z.1.within^4))

# normal
set.seed(412)
x.2 <- rnorm(1e6,0,1)
Z.2 <- (x.2-mean(x.2))/sd(x.2)
mean(Z.2^4)
Z.2.outlier <- Z.2[abs(Z.2)>1]
Z.2.within <- Z.2[abs(Z.2)<=1]
sum(Z.2.within^4)/(1e6)
sum(Z.2.within^4)/(sum(Z.2.outlier^4)+sum(Z.2.within^4))

# exp
set.seed(412)
x.3 <- rexp(1e6)
Z.3 <- (x.3-mean(x.3))/sd(x.3)
mean(Z.3^4)
Z.3.outlier <- Z.3[abs(Z.3)>1]
Z.3.within <- Z.3[abs(Z.3)<=1]
sum(Z.3.within^4)/(1e6)
sum(Z.3.within^4)/(sum(Z.3.outlier^4)+sum(Z.3.within^4))

# cauchy
set.seed(412)
x.4 <- rcauchy(1e6)
Z.4 <- (x.4-mean(x.4))/sd(x.4)
mean(Z.4^4)
Z.4.outlier <- Z.4[abs(Z.4)>1]
Z.4.within <- Z.4[abs(Z.4)<=1]
sum(Z.4.within^4)/(1e6)
sum(Z.4.within^4)/(sum(Z.4.outlier^4)+sum(Z.4.within^4))

# poisson
set.seed(412)
x.5 <- rpois(1e6,5)
Z.5 <- (x.5-mean(x.5))/sd(x.5)
mean(Z.5^4)
Z.5.outlier <- Z.5[abs(Z.5)>1]
Z.5.within <- Z.5[abs(Z.5)<=1]
sum(Z.5.within^4)/(1e6)
sum(Z.5.within^4)/(sum(Z.5.outlier^4)+sum(Z.5.within^4))

# Bimodal
set.seed(412)
n <- 10^6
rand.beta<-rbeta(n,2,2)
x.6 <- c(rand.beta[1:(n/2)]*2+1, rand.beta[(n/2+1):n]*2+3)
Z.6 <- (x.6-mean(x.6))/sd(x.6)
mean(Z.6^4)
Z.6.outlier <- Z.6[abs(Z.6)>1]
Z.6.within <- Z.6[abs(Z.6)<=1]
sum(Z.6.within^4)/(1e6)
sum(Z.6.within^4)/(sum(Z.6.outlier^4)+sum(Z.6.within^4))



-----------------------------------------------------------------------------------



####### Kaplansky Example #######

#### P(x) ####
support_px <- seq(-5, 5, 0.01)
Px <- function(x){
  1/(3*sqrt(pi))*(9/4+x^4)*(exp(-x^2))
}
probs_px <- Px(support_px)
plot(support_px, probs_px, type = "l", xlab = "x", ylab = "P(x)")  # Plot pdf

# Sampling
samp_fun_px <- function(n) sample(
  support_px, 
  n, # n draws
  TRUE,  # with replacement
  probs_px # using these probabilities
)
simulated_px <- samp_fun_px(1e6)

# Plot simulated values
hist(
  simulated_px, 
  breaks = seq(-5, 5, 0.25),
  freq = FALSE,
  xlab = "Support"
)

p_mean <- mean(simulated_px); p_sd <- sd(simulated_px)
mean(((simulated_px-p_mean)/p_sd)^4) # Sample Kurtosis of P(x)



#### S(x) ####
support_sx <- seq(-5, 5, 0.01)
Sx <- function(x){
  (3*sqrt(3))/(16*sqrt(pi))*(2+x^2)*(exp(-3*x^2/4))
}
probs_sx <- Sx(support_sx)
plot(support_sx, probs_sx, type = "l")

# Sampling
samp_fun_sx <- function(n) sample(
  support_sx, 
  n, # n draws
  TRUE,  # with replacement
  probs_sx # using these probabilities
) 
simulated_sx <- samp_fun_sx(1e6)


# Plot simulated values
hist(
  simulated_sx,
  breaks = seq(-5, 5, 0.25),
  freq = FALSE,
  xlab = "Support"
)

s_mean <- mean(simulated_sx); s_sd <- sd(simulated_sx)
mean(((simulated_sx-s_mean)/s_sd)^4) # Sample Kurtosis of S(x)


#### Q(x) ####
support_qx <- seq(-5, 5, 0.01)
Qx <- function(x){
  3/(2*sqrt(2*pi))*exp(-x^2/2)-1/(6*sqrt(pi))*(9/4+x^4)*exp(-x^2)
}
probs_qx <- Qx(support_qx)
plot(support_qx, probs_qx, type = "l")

# Sampling
samp_fun_qx <- function(n) sample(
  support_qx, 
  n, # n draws
  TRUE,  # with replacement
  probs_qx # using these probabilities
) 
simulated_qx <- samp_fun_qx(1e6)

# Plot simulated values
hist(
  simulated_qx,
  breaks = seq(-5, 5, 0.25),
  freq = FALSE,
  xlab = "Support"
)
q_mean <- mean(simulated_qx); q_sd <- sd(simulated_qx)
mean(((simulated_qx-q_mean)/q_sd)^4) # Sample Kurtosis of Q(x)



#### R(x) ####
support_rx <- seq(-5, 5, 0.01)
Rx <- function(x){
  1/(6*sqrt(pi))*(exp(-x^2/4)+4*exp(-x^2))
}
probs_rx <- Rx(support_rx)
plot(support_rx, probs_rx, type = "l")

# Sampling
samp_fun_rx <- function(n) sample(
  support_rx, 
  n, # n draws
  TRUE,  # with replacement
  probs_rx # using these probabilities
) 
simulated_rx <- samp_fun_rx(1e6)

# Plot simulated values
hist(
  simulated_rx,
  breaks = seq(-5, 5, 0.25),
  freq = FALSE,
  xlab = "Support"
)
r_mean <- mean(simulated_rx); r_sd <- sd(simulated_rx)
mean(((simulated_rx-r_mean)/r_sd)^4) # Sample Kurtosis of R(x)

#### Plot with Normal Distribution ####
par(mfrow=c(1,2))
x <- seq(-5, 5, 0.01)
y <- dnorm(x, 0, 1)
plot(x, y, type="l", ylim = c(0, 0.6), col="black", ylab = "Prob", main = "P(x), S(x), Normal")
lines(support_px, probs_px, type = "l", col="blue")
lines(support_sx, probs_sx, type = "l", col="red")
legend("topleft",
       c("N(0,1)","P(x)", "S(x)"),
       col=c("black","blue","red"),lwd=2)

plot(x, y, type="l", ylim = c(0, 0.6), col="black", ylab = "Prob", main = "Q(x), R(x), Normal")
lines(support_qx, probs_qx, type = "l", col="blue")
lines(support_rx, probs_rx, type = "l", col="red")
legend("topleft",
       c("N(0,1)", "Q(x)", "R(x)"),
       col=c("black","blue","red"),lwd=2)




-----------------------------------------------------------------------------------



####### Westfall Example #######

library(EnvStats); par(mfrow=c(3,2))
kt <- function(x) sum((x-mean(x))^4/(sd(x)^4))/length(x) #kurtosis

### Tringular dist / mixed with T(df=1) (prob=0.0001)
set.seed(4)
x <- rtri(10^6,-2.4495,2.4495,0)
kt(x)
x1 <- rtri(10^6*0.9999,-2.4495,2.4495,0) 
x2 <- rt(10^6*0.0001,1)
y <- c(x1,x2)
kt(y)
# different kurtosis with same shape
hist(x, nclass = 100, xlim = c(-5,5), freq = FALSE, main = "Histogram of Triangular (k=2.4)")
hist(y, nclass = 10000, xlim = c(-5,5), freq = FALSE, main = "Histogram of Mixed Traingular (k=578808)")



### Slip-dress dist / mixed with T(df=1) (prob=0.0001)
n <- 10^6; set.seed(123); rand.b<-rbeta(n,0.5,1)
r.1 <- sample(rand.b, n/4, replace = F ); r.2 <- sample(rand.b, n/4, replace = F )
r.3 <- sample(rand.b, n/4, replace = F ); r.4 <- sample(rand.b, n/4, replace = F )
x <- c(0.7241+1.5423*r.1, 0.7241-1.5423*r.2, 
       -0.7241+1.5423*r.3, -0.7241-1.5423*r.4)
x1 <- sample(x, 10^6*0.9999, replace = F)
x2 <- rt(10^6*0.0001, 1)
y <- c(x1,x2)
kt(x);kt(y)
# different kurtosis with same shape
hist(x, nclass = 100, xlim = c(-2,2), ylim = c(0,1), freq = FALSE, main = "Histogram of Slip-dress (k=2.4)")
hist(y, nclass = 10000, xlim = c(-2,2), ylim = c(0,1),freq = FALSE, main = "Histogram of Mixed Slip-dress (k=183320)")



### Bimodal dist / mixed with T(df=1) (prob=0.0001)
n <- 10^6
set.seed(32)
x <- rbeta(n,0.5,0.5)
x1 <- sample(x, 10^6*0.9999, replace = F)
x2 <- rt(10^6*0.0001, 1)
y <- c(x1,x2)
kt(x);kt(y)
# different kurtosis with same shape
hist(x, nclass = 100, xlim = c(0,1), freq = FALSE, main = "Histogram of Bimodal (k= 1.58)")
hist(y, nclass = 100000, xlim = c(0,1), freq = FALSE, main = "Histogram of Mixed Bimodal (k= 753658)")




### Tringular dist / mixed with T(df=1) (prob=0.00005)
set.seed(95)
x <- rtri(10^6,-2.4495,2.4495,0)
kt(x)
x1 <- rtri(10^6*0.99995,-2.4495,2.4495,0) 
x2 <- rt(10^6*0.00005, 1)
y <- c(x1,x2)
kt(y)
# different kurtosis with same shape
hist(x, nclass = 100, xlim = c(-5,5), freq = FALSE, main = "Histogram of Triangular (k=2.4)")
hist(y, nclass = 100000, xlim = c(-5,5), freq = FALSE, main = "Histogram of Mixed Traingular (k=1684)")



### Slip-dress dist / mixed with cauchy(df=1) (prob=0.00005)
n <- 10^6; set.seed(123); rand.b<-rbeta(n,0.5,1)
r.1 <- sample(rand.b, n/4, replace = F ); r.2 <- sample(rand.b, n/4, replace = F )
r.3 <- sample(rand.b, n/4, replace = F ); r.4 <- sample(rand.b, n/4, replace = F )
x <- c(0.7241+1.5423*r.1, 0.7241-1.5423*r.2, 
       -0.7241+1.5423*r.3, -0.7241-1.5423*r.4)
set.seed(9)
x1 <- sample(x, 10^6*0.99995, replace = F)
x2 <- rt(10^6*0.00005, 1)
y <- c(x1,x2)
kt(x);kt(y)

# different kurtosis with same shape
hist(x, nclass = 100, xlim = c(-2,2), freq = FALSE, main = "Histogram of Slip-dress (k=2.4)")
hist(y, nclass = 10000, xlim = c(-2,2), freq = FALSE, main = "Histogram of Mixed Slip-dress (k=604905)")




# Bimodal dist / mixed with cauchy(df=1) (prob=0.00005)
n <- 10^6
set.seed(30)
x <- rbeta(n,0.5,0.5)
x1 <- sample(x, 10^6*0.99995, replace = F)
x2 <- rt(10^6*0.00005, 1)
y <- c(x1,x2)
kt(x);kt(y)
# different kurtosis with same shape
hist(x, nclass = 100, xlim = c(0,1), freq = FALSE, main = "Histogram of Bimodal (k= 1.49)")
hist(y, nclass = 10000, xlim = c(0,1), freq = FALSE, main = "Histogram of Mixed Bimodal (k= 22451)")







### Tringular dist / mixed with T(df=1) (prob=0.001)
set.seed(14)
x <- rtri(10^6,-2.4495,2.4495,0)
kt(x)
x1 <- rtri(10^6*0.999,-2.4495,2.4495,0) 
x2 <- rt(10^6*0.001, 1)
y <- c(x1,x2)
kt(y)
# different kurtosis with same shape
hist(x, nclass = 100, xlim = c(-5,5), freq = FALSE, main = "Histogram of Triangular (k=2.4)")
hist(y, nclass = 10000, xlim = c(-5,5), freq = FALSE, main = "Histogram of Mixed Traingular (k=27261)")



### Slip-dress dist / mixed with T(df=1) (prob=0.001)
n <- 10^6; set.seed(3); rand.b<-rbeta(n,0.5,1)
r.1 <- sample(rand.b, n/4, replace = F ); r.2 <- sample(rand.b, n/4, replace = F )

r.3 <- sample(rand.b, n/4, replace = F ); r.4 <- sample(rand.b, n/4, replace = F )
x <- c(0.7241+1.5423*r.1, 0.7241-1.5423*r.2, 
       -0.7241+1.5423*r.3, -0.7241-1.5423*r.4)
x1 <- sample(x, 10^6*0.999, replace = F)
x2 <- rt(10^6*0.001, 1)
y <- c(x1,x2)
kt(x);kt(y)
# different kurtosis with same shape
hist(x, nclass = 100, xlim = c(-2,2), ylim=c(0,1), freq = FALSE, main = "Histogram of Slip-dress (k=2.4)")
hist(y, nclass = 10000, xlim = c(-2,2), ylim=c(0,1), freq = FALSE, main = "Histogram of Mixed Slip-dress (k=73377)")



### Bimodal dist / mixed with T(df=1) (prob=0.001)
n <- 10^6
set.seed(4)
x <- rbeta(n,0.5,0.5)
x1 <- sample(x, 10^6*0.999, replace = F)
x2 <- rt(10^6*0.001, 1)
y <- c(x1,x2)
kt(x);kt(y)
# different kurtosis with same shape
hist(x, nclass = 100, xlim = c(0,1), freq = FALSE, main = "Histogram of Bimodal (k= 1.5)")
hist(y, nclass = 100000, xlim = c(0,1), freq = FALSE, main = "Histogram of Mixed Bimodal (k= 742185)")






######## Asymmetric case ########



x <- seq(0, 10, length.out = 101)
y8_0.5 <- dgamma(x, 8, scale = 0.5)
y10_0.5 <- dgamma(x, 10, scale = scale)
plot(x,y8_0.5, ylim=c(0,0.4),type="l", col = "blue", ylab = "density", main = "Plot of two Gamma Distributions")
lines(x,y10_0.5, type="l", col = "red")
legend("topleft", 
       legend = c("Gamma(8, 0.5)","Gamma(10, 0.5)"),
       col =c("blue","red"), lwd = 2)


x <- seq(-3, 20, length.out = 101); y <- dnorm(x, 0, 1)
y8_0.5 <- dgamma(x, 8, scale = 0.5)
y10_0.5 <- dgamma(x, 10, scale = scale)
plot(x, y, type="l", ylim = c(0, 0.55), xlim = c(-3,6), ylab = "density", col="black", main = "Compare with Normal")
lines(x-3.5, y8_0.5, col = "blue")
lines(x-4.5, y10_0.5, col = "red")
legend("topleft", 
       legend = c("Standard Normal", "Gamma(8, 0.5) with x-axis translation","Gamma(10, 0.5) with x-axis translation"),
       col =c("black", "blue","red"), lwd = 2)


kt(rgamma(1e6, 8, scale=0.5)) #Kurtosis
kt(rgamma(1e6, 10, scale=0.5)) # Kurtosis



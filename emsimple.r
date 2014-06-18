#### test the optim function
## testopti = function(theta) {
##     (theta[1]^2 + theta[2]^2)
## }
## optim(c(.5, .2), method="L-BFGS-B", fn=testopti, lower=c(-1000, -1000), upper = c(1000, 1000))

nToss = 50
thetaTrue1 = .2
thetaTrue2 = .9
thetaTrue = c(thetaTrue1, thetaTrue2)
x1 = rbinom(1, nToss, thetaTrue1)
x2 = rbinom(1, nToss, thetaTrue2)
print(x1)
print(x2)
x = c(x1, x2)

# four possible configurations,
# coin1 coin1; coin1, coin2; ...
z1 = c(1, 1)
z2 = c(1, 2)
z3 = c(2, 1)
z4 = c(2, 2)
z = list(z1, z2, z3, z4)
print(z)


# likelihood function
lhfun = function(x, zi, theta) {
    x1 = x[1]
    x2 = x[2]
    theta1 = theta[zi[1]]
    theta2 = theta[zi[2]]
    theta1^x1 * (1-theta1)^(nToss-x1) * theta2^x2 * (1-theta2)^(nToss-x2)
}

lhsum = function(x, z, theta) {
    lhlist = lapply(z, function(zi) {
        lhfun(x, zi, theta)
    })
    do.call(sum, lhlist)
}


pzi = function(x, z, zi, thetaHat) {
    likelihoodSum = lhsum(x, z, thetaHat)
    likelihoodZi = lhfun(x, zi, thetaHat)
    likelihoodZi / likelihoodSum
}
sapply(1:length(z), function(i) {
    pzi(x, z, z[[i]], c(.2, .8))
})

pzi(x, z, z[[2]], c(.3, .8))

## estimation function
eFunc = function(x, z, thetaHat) {
    toMaximize = function(theta) {
        glist = lapply(1:length(z), function(i) {
            pzires = pzi(x, z, z[[i]], thetaHat)
            pzires * log(lhfun(x, z[[i]], theta) / pzires)
        })
        -do.call(sum, glist)
    }
}
eFuncToMax = eFunc(x, z, c(.3, .7))
eFuncToMax(c(.3, .8))

optim(c(.2, .3), method="L-BFGS-B", fn=eFuncToMax, lower=c(0.001, 0.001), upper = c(.999, .999))

emSimple = function(x, z, thetaHatInit) {
    thetaHat0 = optim(
          # initial guess for the optimization procedure, not so important.
          c(.2, .3),
          method="L-BFGS-B",
          fn=eFunc(x, z, thetaHatInit),
          lower=c(0.001, 0.001),
          upper = c(.999, .999)
          )$par
    thetaHat1 = optim(
          c(.2, .3),
          method="L-BFGS-B",
          fn=eFunc(x, z, thetaHat0),
          lower=c(0.001, 0.001),
          upper = c(.999, .999)
          )$par
    nIters = 1
    while(max(abs(thetaHat1 - thetaHat0)) > 0.00001) {
        nIters = nIters + 1
        thetaHat0 = thetaHat1
        thetaHat1 = optim(
              c(.2, .3),
              method="L-BFGS-B",
              fn=eFunc(x, z, thetaHat0),
              lower=c(0.001, 0.001),
              upper = c(.999, .999)
              )$par
    }
    list(EMestimate=thetaHat1, nIters)
}
emSimple(x, z, c(.5, .6))

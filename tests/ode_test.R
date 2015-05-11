#!/usr/bin/env Rscript

library(deSolve)

y0 <- c(I=1, S=4999)
times <- seq(0, 100000, 100000/100)
beta <- 0.011
gamma <- 0.005
mu <- 0.005

dx <- function (t, y, parms, ...) {
    I <- y["I"]
    S <- y["S"]
    N <- I + S

    list(c(I = parms$beta*I*S/N - (parms$mu+parms$gamma)*I,
           S = -parms$beta*I*S/N + (parms$mu+parms$gamma)*I))
}

parms <- list(beta=beta, gamma=gamma, mu=mu)
system.time(res <- ode(y0, times, dx, parms, method="rk4"))
res

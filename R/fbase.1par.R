fbase.bern.ilogit <- function(u,y) {
  eu <- exp(u)
  f <- -log(1+1/eu)-(1-y)*u
  g <- 1/(1+eu)-(1-y)
  h <- -eu/(1+eu)^2
  return (list(f=f,g=g,h=h))
}

fbase.pois.exp <- function(u,y) {
  eu <- exp(u)
  f <- y*u-eu - lfactorial(y)
  g <- y-eu
  h <- -eu
  return (list(f=f,g=g,h=h))
}

fbase.exp.exp <- function(u,y) {
  eu <- exp(u)
  f <- u-y*eu
  g <- 1-y*eu
  h <- -y*eu
  return (list(f=f,g=g,h=h))
}

fbase.geo.ilogit <- function(u,y) {
  eu <- exp(u)
  f <- -(y*u+(1+y)*log(1+1/eu))
  g <- -y+(1+y)/(1+eu)
  h <- -(1+y)*eu/(1+eu)^2
  return (list(f=f,g=g,h=h))
}

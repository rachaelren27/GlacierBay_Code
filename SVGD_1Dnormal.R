log_grad <- function(x, mu, sigma){
  (mu - x)/(sigma^2)
}

SVGD_update <- function(x.particles, epsilon, mu, sigma){
  n.particles <- length(x.particles)
  
  pairwise.dist <- as.matrix(dist(x.particles))
  h <- (median(pairwise.dist)^2) / log(n.particles)
  kernel.matrix <- exp(-pairwise.dist^2 / h)
  
  log.grads <- log_grad(x.particles, mu, sigma)
  
  x.updates <- rep(0, n.particles)
  
  for (i in 1:n.particles) {
    x.i <- x.particles[i]
    
    kernel.grads <- (x.particles - x.i) * (-2 / h) * kernel.matrix[i, ]
    
    x.updates[ i] <- sum(
      kernel.matrix[i,]*log.grads + kernel.grads
    )
  }
  
  x.particles <- x.particles + (epsilon / n.particles) * x.updates
  
  return(x.particles)
}

SVGD <- function(x.init, epsilon, n.iter, mu, sigma){
  SVGD.out <- list()
  SVGD.out[[1]] <- x.init
  for(k in 2:n.iter){
    x.particles <- SVGD.out[[k-1]]
    SVGD.out[[k]] <- SVGD_update(x.particles, epsilon, mu, beta)
    if(k %% 100 == 0){
      print(k)
    }
  }
  
  return(SVGD.out)
}

n.particles <- 100
x.init <- rnorm(n.particles)
epsilon <- 0.01
n.iter <- 10000
mu <- 5
sigma <- 1.5

SVGD.out <- SVGD(x.init, epsilon, n.iter, mu, sigma)

SVGD.save <- c()
for(k in 1:n.iter){
  SVGD.save[k] <- sd(SVGD.out[[k]])
}
plot(SVGD.save)

mean(SVGD.out[[n.iter]])
sd(SVGD.out[[n.iter]])

x <- seq(-5, 10, length.out = 500)
y <- dnorm(x, mean = 5, sd = 1.5)

plot(x, y, type = "l", lwd = 2, 
     xlab = "X", ylab = "Density")
for(i in 1:n.particles){
  points(x = SVGD.out[[n.iter]][i], y = dnorm(SVGD.out[[n.iter]][i], mu, beta),
         pch = 16, cex = 1,col = "red")
}
lines(density(SVGD.out[[n.iter]]), col = "green")



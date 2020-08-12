dhalf = function(x, sigma){
  return(exp(-x^2 / (2 * sigma^2)))
}

plt = sapply(1:2000, dhalf, sigma = 600)
plot(plt, col = "GREEN", main = "Half-normal detection function for varying sigma values (g0 = 1)")
points(sapply(1:2000, dhalf, sigma = 700), col = "BLUE")
points(sapply(1:2000, dhalf, sigma = 800), col = "RED")

legend(x = 1700, y = 0.95, c("600", "700", "800"), lty=c(1,1), lwd=c(2.5,2.5), col=c("GREEN", "BLUE", "RED"))
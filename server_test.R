x <- 5
y <- 3

z <- x + y

pdf("server_test.pdf")
plot(x = x, y = y)
dev.off()

save(x,y,z, file = "test_output.RData")


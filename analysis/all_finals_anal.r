sims <- getfile()
sim <- read.table(sims$files, header = TRUE, sep = " ")
sim$Generation <- as.numeric(vapply(strsplit(as.character(sim$File), "_"), function(x) strsplit(x[10], "\\.")[[1]][1], "butt"))
sim$Samp.Size <- factor(sim$Samp.Size, levels = c(10,25,50,100))

simnarm <- sim[-which(is.na(sim$rbarD)), ]

Ia.bar <- ggplot(simnarm) + geom_boxplot(aes(x = factor(Generation*1000), y = Ia, fill = Samp.Size), outlier.shape = 20, outlier.size = 1, notch = TRUE) + facet_wrap(~Sex.Rate, scales = "free_y", nrow = 2) + theme_bw() + scale_fill_grey(start = 0.4, end = 1)
rbarD.bar <- ggplot(simnarm) + geom_boxplot(aes(x = factor(Generation*1000), y = rbarD, fill = Samp.Size), outlier.shape = 20, outlier.size = 1, notch = TRUE) + facet_wrap(~Sex.Rate, scales = "free_y", nrow = 2) + theme_bw() + scale_fill_grey(start = 0.4, end = 1)


Ia.point <- ggplot(simnarm, aes(x = factor(Generation*1000), y = Ia, color = Samp.Size, group = Samp.Size, alpha = 0.5)) + geom_point(position = "jitter") + stat_smooth() + facet_grid(Sex.Rate ~ Samp.Size) + theme_bw()
rbarD.point <- ggplot(simnarm, aes(x = factor(Generation*1000), y = rbarD, color = Samp.Size, group = Samp.Size, alpha = 0.5)) + geom_point(position = "jitter") + stat_smooth() + facet_grid(Sex.Rate ~ Samp.Size) + theme_bw()


# Idea for generating alternative hypothesis
data(nancycats)
nan9 <- popsub(nancycats, 9)
nan9null <- poppr:::.sampling(pop=seploc(nan9), iterations=999, quiet=TRUE, type="codom")
nan9alt <- sapply(1:1000, function(x) ia(nan9[sample(1:9, 5), ]))
library(reshape)
nan9null <- melt(nan9null)
nan9null$dist <- "null"
nan9alt <- melt(nan9alt)
nan9alt$dist <- "alt"
nan9test <- rbind(nan9null, nan9alt)
ggplot(nan9test, aes(x = value, fill = dist)) + geom_histogram(alpha = 0.5, position = "identity") + geom_rug(alpha = 0.5, aes(color = dist)) + facet_wrap(~variable, scales = "free") + theme_classic()


data(Aeut)
Athena <- popsub(Aeut, "Athena")
Anull <- .sampling(Athena, iterations= 999, quiet=TRUE, type="PA")
Aalt <- sapply(1:1000, function(x) ia(Athena[sample(1:nInd(Athena), nInd(Athena)/2), ]))
Anull$dist <- "null"
Anull <- melt(Anull)
Aalt <- melt(data.frame(t(Aalt)))
Aalt$dist <- "alt"
Atest <- rbind(Anull, Aalt)
ggplot(Atest, aes(x = value, fill = dist)) + geom_histogram(alpha = 0.5, position = "identity") + geom_rug(alpha = 0.5, aes(color = dist)) + facet_wrap(~variable, scales = "free") + theme_classic()




nannull <- .sampling(seploc(nancycats), iterations= 999, quiet=TRUE, type="codom")
nanalt <- sapply(1:1000, function(x) ia(nancycats[sample(1:nInd(nancycats), nInd(nancycats)/2), ]))
nannull$dist <- "null"
nannull <- melt(nannull)
nanalt <- melt(data.frame(t(nanalt)))
nanalt$dist <- "alt"
nantest <- rbind(nannull, nanalt)
ggplot(nantest, aes(x = value, fill = dist)) + geom_histogram(alpha = 0.5, position = "identity") + geom_rug(alpha = 0.5, aes(color = dist)) + facet_wrap(~variable, scales = "free") + theme_classic()




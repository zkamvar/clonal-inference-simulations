library(poppr)
library(reshape)
sims <- getfile()
sim <- read.table(sims$files, header = TRUE, sep = " ")
sim$Generation <- as.numeric(vapply(strsplit(as.character(sim$File), "_"), function(x) strsplit(x[10], "\\.")[[1]][1], "butt"))
sim$Sex.Rate <- round(sim$Sex.Rate, 5)
sim$Sex.Rate <- factor(sim$Sex.Rate)
sim$Samp.Size <- factor(paste("n =",sim$Samp.Size))
sim$Samp.Size <- factor(sim$Samp.Size, levels(sim$Samp.Size)[c(1,3:4,2)])
simnarm <- sim[-which(is.na(sim$rbarD)), ]
themes <-  theme_bw() + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1, family = "Helvetica"), 
        axis.text.y = element_text(family = "Helvetica"),
        axis.title = element_text(size = rel(2)),
        panel.grid.major.x = element_blank(), 
        strip.background = element_rect(color = "black", fill = "white"), 
        strip.text = element_text(face = "bold", family = "Helvetica"), 
        plot.title = element_text(face = "bold", size = rel(2), vjust = 1, 
                                  family = "Helvetica"))

theme2 <- theme_bw() + 
  theme(legend.position = c(1, 0), 
        legend.justification = c(1, 0), 
        legend.direction = "horizontal",
        legend.title = element_text(size = rel(1), family = "Helvetica"),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1, family = "Helvetica"), 
        axis.text.y = element_text(family = "Helvetica"), 
        axis.title.y = element_text(angle = 0, face = "bold", size = rel(2)), 
        axis.title.x = element_text(angle = 0, face = "bold", size = rel(2)),
        panel.grid.major.x = element_blank(), 
        plot.title = element_text(face = "bold", size = rel(2), vjust = 1, 
                                  family = "Helvetica")
  )

theme3 <- theme_bw() + 
  theme(#legend.position = c(0, 1), 
    #legend.justification = c(0, 1), 
    legend.position = "top",
    legend.title = element_text(size = rel(1), family = "Helvetica"),
    axis.text.x = element_text(angle = 45, hjust=1, vjust=1, family = "Helvetica"), 
    axis.text.y = element_text(family = "Helvetica"), 
    axis.title.y = element_text(angle = 90, face = "bold", size = rel(2)), 
    axis.title.x = element_text(angle = 0, face = "bold", size = rel(2)),
    panel.grid.major.x = element_blank(), 
    plot.title = element_text(face = "bold", size = rel(2), vjust = 1, 
                              family = "Helvetica")
  )
Ialab <- labs(list(y = expression(I[A]), x = "Generation (x1000)", fill = "Sample Size"))
rbarDlab <- labs(list(y = expression(bar(r)[d]), x = "Generation (x1000)", fill = "Sample Size"))
Ia.bar <- ggplot(simnarm) + geom_boxplot(aes(x = factor(Generation), y = Ia, fill = Samp.Size), outlier.shape = 20, outlier.size = 1, notch = TRUE) + facet_wrap(~Sex.Rate, scales = "free_y") + theme_bw() + scale_fill_grey(start = 0.4, end = 1)
rbarD.bar <- ggplot(simnarm) + geom_boxplot(aes(x = factor(Generation), y = rbarD, fill = Samp.Size), outlier.shape = 20, outlier.size = 1, notch = TRUE) + facet_wrap(~Sex.Rate, scales = "free_y") + theme_bw() + scale_fill_grey(start = 0.4, end = 1)
Ia.bar + Ialab + theme3
rbarD.bar + rbarDlab + theme3

#Ia.point <- ggplot(simnarm, aes(x = factor(Generation), y = Ia, color = Samp.Size, group = Samp.Size, alpha = 0.0001)) + geom_point(position = "jitter") + stat_smooth(alpha = 1) + facet_wrap(~Sex.Rate, scales = "free_y") + theme_bw()
#rbarD.point <- ggplot(simnarm, aes(x = factor(Generation), y = rbarD, color = Samp.Size, group = Samp.Size, alpha = 0.5)) + geom_point(position = "jitter") + stat_smooth() + facet_grid(Sex.Rate ~ Samp.Size) + theme_bw()

# Missing stats:

Samp.Sex <- table(sim[is.na(sim$rbarD) & sim$Generation == 10, c("Sex.Rate", "Samp.Size")])
MLG.Sex <- table(sim[is.na(sim$rbarD) & sim$Generation == 10, c("Sex.Rate", "MLG")])
Samp.Sex.df <- data.frame(t(Samp.Sex[1:2, ]))
names(Samp.Sex.df) <- c(0, 1e-04)
Samp.Sex.df$Samp.Size <- factor(rownames(Samp.Sex.df), levels = rownames(Samp.Sex.df))
ggplot(melt(Samp.Sex.df)) + 
  geom_bar(aes(x = Samp.Size, y = value, group = variable, fill = variable), stat = "identity", position = "dodge") +
  labs(list(x = "Sample Size", y = "Number of Missing Values", fill = "Sex Rate")) +
  theme3 + scale_fill_grey(start = 0, end = 0.5)

MLG.Sex.df <- data.frame(t(MLG.Sex[1:2, ]))
names(MLG.Sex.df) <- c(0, 1e-04)
MLG.Sex.df$MLG <- factor(rownames(MLG.Sex.df), levels = rownames(MLG.Sex.df))
ggplot(melt(MLG.Sex.df)) + 
  geom_bar(aes(x = MLG, y = value, group = variable, fill = variable), stat = "identity", position = "dodge") +
  labs(list(x = "Number of MLGs", y = "Number of Missing Values", fill = "Sex Rate")) +
  theme3 + scale_fill_grey(start = 0, end = 0.5)

MLG.Sex.Samp <- table(sim[is.na(sim$rbarD) & sim$Generation == 10, c("MLG", "Samp.Size", "Sex.Rate")])[, , 1:2]
MLG.Sex.Samp.df <- melt(MLG.Sex.Samp)

ggplot(MLG.Sex.Samp.df) + geom_bar(aes(x = factor(MLG), y = value, color = Samp.Size, fill = Samp.Size, alpha = rev(factor(Sex.Rate))), stat = "identity", position = "stack") + scale_alpha_manual(values=c(0.3, 0.9)) + theme_bw()

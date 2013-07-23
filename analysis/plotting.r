# Creating the functions
showdensity2 <- function(pop, sample=999){
  if(pop@type=="PA"){
    pop_meth <- lapply(1:4, function(x) poppr:::.sampling(pop, ifelse(x !=4, sample, sample*2), type=pop@type, quiet=F, method=x))
  }
  else{
    pop_meth <- lapply(1:4, function(x) poppr:::.sampling(seploc(pop), ifelse(x !=4, sample, sample*2), type=pop@type, quiet=F, method=x))
  }
  df <- NULL
  methlev <- c("Original (multilocus)", "Permutataion", "Parametric Bootstrap", "Non-Parametric Bootstrap")
  invisible(lapply(1:4, function(x) df <<- rbind(df, data.frame(list(Value=unlist(pop_meth[[x]]), Index=rep(c("Ia","rbarD"), each=ifelse(x!=4,sample,sample*2)), Method=rep(methlev[x], sample*ifelse(x!=4,2,4)))))))
  method_plot <- ggplot(data=df, aes(Value, color=Method)) + 
    geom_density(aes(fill=Method), alpha=0.1) + 
    facet_wrap(. ~ Index, scales="free") + theme_classic()
  print(method_plot)
  return(list(DataFrame=df, Plot=method_plot))
}

showdensity <- function(pop, sample=999){
  Indexfac <- factor(1:2, levels=1:2, labels=c("I[A]","bar(r)[D]"))
  if(pop@type=="PA"){
    pop_meth <- lapply(1:4, function(x) poppr:::.sampling(pop, sample, type=pop@type, quiet=F, method=x))
  }
  else{
    pop_meth <- lapply(1:4, function(x) poppr:::.sampling(seploc(pop), sample, type=pop@type, quiet=F, method=x))
  }
  df <- NULL
  methlev <- c("Original (multilocus)", "Permutataion", "Parametric Bootstrap", "Non-Parametric Bootstrap")
  invisible(lapply(1:4, function(x) df <<- rbind(df, data.frame(list(Value=unlist(pop_meth[[x]]), Index=rep(Indexfac, each=sample), Method=rep(methlev[x], sample*2))))))
  method_plot <- ggplot(data=df, aes(Value, color=Method)) + 
    geom_density(aes(fill=Method), alpha=0.1) + 
    geom_rug() +
    facet_wrap(~ Index, nrow=1, scales="free") + theme_classic()
  print(method_plot)
  return(list(DataFrame=df, Plot=method_plot))
}



showplots <- function(x, name){
  cat("|", name,"\n")
  themed <- theme(axis.text.x=element_text(size = 10, angle=45, hjust=1, vjust = 1), legend.position = "bottom")
  thickens <- x + labs(title=paste(name, "over 1000 data sets")) + xlab("Sex Rate") + ylab(name) + themed
#   if(name %in% IndexNames[c(2,10)]){
#   	thickens <- thickens + geom_hline(aes(yintercept = 0.05, color = "red")) + annotate("text", x = 0, y = 0.08, label = "p = 0.05", color = "red") + coord_trans(y = "log2")
#   }
  ggsave(filename=paste(name,"pdf",sep="."), plot=thickens, width = 16, height = 8, units="in", dpi=300)
  return
}
pvalprop.rd <- function(df, alpha=0.05){
  sam <- levels(df$Samp.Size)
  sex <- levels(df$Sex.Rate)
  #cat(sam, "\n")
  #cat(sex, "\n")
  derp <- sapply(sex, function(x, alpha){
    #cat("x:",x,"\n")
    temp <- df[df$Sex.Rate == x, ]
    head(temp)
    vapply(sam, function(y, alpha){
      #cat("y:",y,"\n")
      #res <- sum(temp[temp$Samp.Size == y, ]$p.Ia <= 0.05, na.rm=TRUE)
      res <- sum(temp[temp$Samp.Size == y, ]$p.rD <= alpha, na.rm=TRUE)
      return(res)
    }, 1, alpha)
  }, alpha)
  rownames(derp) <- sam
  return(derp)
}
pvalprop.Ia <- function(df, alpha=0.05){
  sam <- levels(df$Samp.Size)
  sex <- levels(df$Sex.Rate)
  #cat(sam, "\n")
  #cat(sex, "\n")
  derp <- sapply(sex, function(x, alpha){
    #cat("x:",x,"\n")
    temp <- df[df$Sex.Rate == x, ]
    head(temp)
    vapply(sam, function(y, alpha){
      #cat("y:",y,"\n")
      res <- sum(temp[temp$Samp.Size == y, ]$p.Ia <= alpha, na.rm=TRUE)
      #res[2] <- sum(temp[temp$Samp.Size == y, ]$p.rD <= 0.05, na.rm=TRUE)
      return(res)
    }, 1, alpha)
  }, alpha)
  rownames(derp) <- sam
  return(derp)
}

splitbymethod <- function(meth, index, alpha=0.05){
  mlist <- NULL
  mlist$ml <- meth[meth$Method == levels(meth$Method)[1], ]
  mlist$perm <- meth[meth$Method == levels(meth$Method)[2], ]
  mlist$pb <- meth[meth$Method == levels(meth$Method)[3], ]
  mlist$npb <- meth[meth$Method == levels(meth$Method)[4], ]
  if(index=="Ia"){
    lapply(mlist, pvalprop.Ia, alpha)
  }
  else{
    lapply(mlist, pvalprop.rd, alpha) 
  }
}

# Starting the script

library(poppr)
x <- getfile(mult=TRUE, pattern="^final.+?csv$")
setwd(x$path)
an <- lapply(x$files, read.table, header=TRUE)

#==============================================================================#
# Edit this for the number of methods you use.
#
methan <- lapply(1:4, function(y) cbind(an[[y]], list(Method=rep(y, 40000))))
meth <- rbind(methan[[1]], methan[[2]], methan[[3]], methan[[4]])
#
# Got it?
#==============================================================================#




meth$Sex.Rate <- round(meth$Sex.Rate, 5)
meth$Sex.Rate <- factor(meth$Sex.Rate)
meth$Method <- factor(meth$Method)
meth$Samp.Size <- factor(paste("n =",meth$Samp.Size))
meth$Samp.Size <- factor(meth$Samp.Size, levels(meth$Samp.Size)[c(1,3:4,2)])




#==============================================================================#
# Hey, this is where you should make sure that you have the correct methods.
#
methlev <- c("Original (multilocus)", "Permutataion", "Parametric Bootstrap", "Non-Parametric Bootstrap")
samlev <- levels(meth$Samp.Size)
sexlev <- levels(meth$Sex.Rate)
#
#
#==============================================================================#

levels(meth$Method) <- methlev
IndexNames <- c("I[A]", "I[A] p-value", "Resampled I[A] min", "Resampled Ia max", "Resampled Ia median", "Resampled Ia mean", "Resampled Ia variance", 
                "Resampled Ia standard deviation", "rbarD", "rbarD p-value", "Resampled rbarD min", "Resampled rbarD max", "Resampled rbarD median", "Resampled rbarD mean", 
                "Resampled rbarD variance", "Resampled rbarD standard deviation")
names(IndexNames) <- names(meth)[10:25]
Plots <- NULL
cat("\n\nCreating Plots\n\n")

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
  theme(legend.position = c(1, 1), 
        legend.justification = c(1, 1), 
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

meth2 <- meth[-which(is.na(meth$rbarD)), ]

Iagrey <- ggplot(meth2, aes_string(x = "Sex.Rate", y = "Ia", fill = "Samp.Size")) + geom_boxplot(outlier.shape = 20, outlier.size = 1, notch = TRUE) + theme_bw() + labs(y = expression(I[A])) + xlab("Sex Rate") + theme2 + labs(fill = "Sample Size") + scale_fill_grey(start = 0.4, end = 1) #+ labs(title = expression(paste(I[A], " over 1000 data sets")), size = rel(2), face = "bold") 

rbarDgrey <- ggplot(meth2, aes_string(x = "Sex.Rate", y = "rbarD", fill = "Samp.Size")) + geom_boxplot(outlier.shape = 20, outlier.size = 1, notch = TRUE) + labs(y = expression(bar(r)[d])) + xlab("Sex Rate") + theme2 + labs(fill = "Sample Size") + scale_fill_grey(start = 0.4, end = 1) #+ labs(title = expression(paste(bar(r)[d], " over 1000 data sets")), size = rel(2), face = "bold") 

Iapvalgrey <- ggplot(meth2[meth2$Method == levels(meth2$Method)[2], ], aes_string(x = "Sex.Rate", y = "p.Ia", fill = "Samp.Size")) + geom_boxplot(outlier.shape = 20, outlier.size = 1) + labs(y = expression(paste(I[A], " p-value (log scale)"))) + xlab("Sex Rate") + theme3 + labs(fill = "Sample Size") + scale_fill_grey(start = 0.4, end = 1) + coord_trans(y = "log2") + geom_hline(aes(yintercept = 0.05), linetype = 2) + annotate("text", x = 0, y = 0.04, label = "p = 0.05")

rbarDpvalgrey <- ggplot(meth2[meth2$Method == levels(meth2$Method)[2], ], aes_string(x = "Sex.Rate", y = "p.rD", fill = "Samp.Size")) + geom_boxplot(outlier.shape = 20, outlier.size = 1) + labs(y = expression(paste(bar(r)[d], " p-value (log scale)"))) + xlab("Sex Rate") + theme3 + labs(fill = "Sample Size") + scale_fill_grey(start = 0.4, end = 1) + coord_trans(y = "log2") + geom_hline(aes(yintercept = 0.05), linetype = 2) + annotate("text", x = 0, y = 0.04, label = "p = 0.05")

E5grey <- ggplot(meth2, aes_string(x = "Sex.Rate", y = "E.5", fill = "Samp.Size")) + geom_boxplot(outlier.shape = 20, outlier.size = 1, notch = TRUE) + theme_bw() + labs(y = expression(E[5])) + xlab("Sex Rate") + theme2 + labs(fill = "Sample Size") + scale_fill_grey(start = 0.4, end = 1) #+ labs(title = expression(paste(E[5], " over 1000 data sets")), size = rel(2), face = "bold") 

Hexpgrey <- ggplot(meth2, aes_string(x = "Sex.Rate", y = "Hexp", fill = "Samp.Size")) + geom_boxplot(outlier.shape = 20, outlier.size = 1, notch = TRUE) + theme_bw() + labs(y = "Expected Heterozygosity") + xlab("Sex Rate") + theme3 + labs(fill = "Sample Size") + scale_fill_grey(start = 0.4, end = 1) #+ labs(title = expression(paste(E[5], " over 1000 data sets")), size = rel(2), face = "bold") 

Hgrey <- ggplot(meth2, aes_string(x = "Sex.Rate", y = "H", fill = "Samp.Size")) + geom_boxplot(outlier.shape = 20, outlier.size = 1, notch = TRUE) + theme_bw() + labs(y = "Shannon-Wiener Index") + xlab("Sex Rate") + theme3 + labs(fill = "Sample Size") + scale_fill_grey(start = 0.4, end = 1) #+ labs(title = expression(paste(E[5], " over 1000 data sets")), size = rel(2), face = "bold") 

Ggrey <- ggplot(meth2, aes_string(x = "Sex.Rate", y = "G", fill = "Samp.Size")) + geom_boxplot(outlier.shape = 20, outlier.size = 1, notch = TRUE) + theme_bw() + labs(y = "Stoddart and Taylor's Index") + xlab("Sex Rate") + theme3 + labs(fill = "Sample Size") + scale_fill_grey(start = 0.4, end = 1) #+ labs(title = expression(paste(E[5], " over 1000 data sets")), size = rel(2), face = "bold") 

Mgrey <- ggplot(meth2, aes_string(x = "Sex.Rate", y = "MLG", fill = "Samp.Size")) + geom_boxplot(outlier.shape = 20, outlier.size = 1, notch = TRUE) + theme_bw() + labs(y = "Multilocus genotypes") + xlab("Sex Rate") + theme3 + labs(fill = "Sample Size") + scale_fill_grey(start = 0.4, end = 1) #+ labs(title = expression(paste(E[5], " over 1000 data sets")), size = rel(2), face = "bold") 


quantsrd <- lapply(levels(meth2$Sex.Rate), function(x) lapply(levels(meth2$Samp.Size), function(y){
  p <- shapiro.test(meth2[meth2$Sex.Rate == x & meth2$Samp.Size == y & meth2$Method == methlev[1], ]$rbarD)$p.value
  derp <- qqnorm(meth2[meth2$Sex.Rate == x & meth2$Samp.Size == y & meth2$Method == methlev[1], ]$rbarD, main = paste("rbarD, p-value:",p,"\nSample Size:",y,"Sex Rate:",x))
  qqline(derp$y)
}))


kruskalrd <- sapply(levels(meth2$Sex.Rate), function(x){
  kruskal.test(meth2[meth2$Sex.Rate == x & 
                       meth2$Method == methlev[1], ]$p.rD, 
               meth2[meth2$Sex.Rate == x & 
                       meth2$Method == methlev[1], ]$Samp.Size)$p.value
})

kruskalia <- sapply(levels(meth2$Sex.Rate), function(x){
  kruskal.test(meth2[meth2$Sex.Rate == x & 
                       meth2$Method == methlev[1], ]$p.Ia, 
               meth2[meth2$Sex.Rate == x & 
                       meth2$Method == methlev[1], ]$Samp.Size)$p.value
})
kruskal <- data.frame(list(Ia = kruskalia,rbarD = kruskalrd))


lapply(names(meth)[10:25], function(derp){
  cat("|", derp,"\n")
  Plots[[derp]] <<- ggplot(meth, aes_string(x = "Sex.Rate", y = derp, fill = "Method")) + 
    #geom_violin(trim=FALSE) + 
    geom_boxplot(outlier.shape = 20, outlier.size = 1, notch = TRUE) +
    facet_wrap(~ Samp.Size, nrow = 2) + 
  	themes + 
  	ylab(IndexNames[derp]) +
  	xlab("Sex Rate")
  return(0)
  })
names(Plots) <- IndexNames


cat("\n\nPrinting Plots (May take a little while)\n\n")
invisible(lapply(IndexNames, function(x) showplots(Plots[[x]], x)))


cat("\n\nCreating the Probability of Rejection plots...rbarD\n\n")
rd.tot <- splitbymethod(meth, "rbarD", alpha=0.05)
rd.comb <- NULL
Ia.comb <- NULL
lapply(rd.tot, function(x) rd.comb <<- rbind(rd.comb, x))
rd.comb <- as.vector(rd.comb/1000)

rd.df <- data.frame(list(Percent_Reject = rd.comb, 
                         Samp.Size = rep(levels(meth$Samp.Size), 10*length(methlev)), 
                         Sex.Rate = rep(levels(meth$Sex.Rate), each=4*length(methlev)), 
                         Method = rep(rep(levels(meth$Method), each=4), 10))
                    )
rd.df$Sex.Rate <- factor(rd.df$Sex.Rate, levels(rd.df$Sex.Rate)[c(1,9:10,2:8)])
rd.df$Samp.Size <- factor(rd.df$Samp.Size, levels(rd.df$Samp.Size)[c(1,3:4,2)])
rd.df$Method <- factor(rd.df$Method, methlev)

cat("\n\nCreating the Probability of Rejection plots...Ia\n\n")

Ia.tot <- splitbymethod(meth, "Ia", alpha=0.05)
lapply(Ia.tot, function(x) Ia.comb <<- rbind(Ia.comb, x))
Ia.comb <- as.vector(Ia.comb/1000)
Ia.df <- data.frame(list(Percent_Reject = Ia.comb, 
                         Samp.Size = rep(levels(meth$Samp.Size), 10*length(methlev)), 
                         Sex.Rate = rep(levels(meth$Sex.Rate), each=4*length(methlev)), 
                         Method = rep(rep(levels(meth$Method), each=4), 10))
                    )
Ia.df$Sex.Rate <- factor(Ia.df$Sex.Rate, levels(Ia.df$Sex.Rate)[c(1,9:10,2:8)])
Ia.df$Samp.Size <- factor(Ia.df$Samp.Size, levels(Ia.df$Samp.Size)[c(1,3:4,2)])
Ia.df$Method <- factor(Ia.df$Method, methlev)

png("Probability_of_Rejection%03d.png", width=1024, height=1225)
rd.plot <- ggplot(rd.df, aes(x = Sex.Rate, y = Percent_Reject, color=Method))
rd.plot <- rd.plot + geom_point() + geom_hline(aes(yintercept=0.050), color="red") + 
  geom_hline(aes(yintercept=0.950), color="blue") + facet_grid(Samp.Size ~ .) +
  labs(title=expression(paste("Probability of Rejecting the Null Hypothesis with ", bar(r)[D]))) + 
  xlab("Sex Rate") + ylab("Probability of Null Hypothesis Rejection") + 
  annotate("text", x = 1, y = 0.10, label="5% rejection", color="red") +
  annotate("text", x = 10, y = 0.90, label="95% rejection", color="blue")
print(rd.plot + geom_line(aes(group = Method)) + theme_bw())

Ia.plot <- ggplot(Ia.df, aes(x = Sex.Rate, y = Percent_Reject, color=Method))
Ia.plot <- Ia.plot + geom_point() + geom_hline(aes(yintercept=0.050), color="red") + 
  geom_hline(aes(yintercept=0.950), color="blue") + facet_grid(Samp.Size ~ .) +
  labs(title=expression(paste("Probability of Rejecting the Null Hypothesis with ", I[A]))) + 
  xlab("Sex Rate") + ylab("Probability of Null Hypothesis Rejection") + 
	annotate("text", x = 1, y = 0.10, label="5% rejection", color="red") +
  annotate("text", x = 10, y = 0.90, label="95% rejection", color="blue")
print(Ia.plot + geom_line(aes(group = Method)) + theme_bw()) 


rd.plot.method <- ggplot(rd.df, aes(x = Sex.Rate, y = Percent_Reject, color=Samp.Size))
rd.plot.method <- rd.plot.method + geom_point() + geom_hline(aes(yintercept=0.050), color="red") + 
  geom_hline(aes(yintercept=0.950), color="blue") + facet_grid(Method ~ .) +
  labs(title=expression(paste("Probability of Rejecting the Null Hypothesis with ", bar(r)[D]))) + 
  xlab("Sex Rate") + ylab("Probability of Null Hypothesis Rejection") + 
  annotate("text", x = 1, y = 0.10, label="5% rejection", color="red") +
  annotate("text", x = 10, y = 0.90, label="95% rejection", color="blue")
print(rd.plot.method + geom_line(aes(group = Samp.Size)) + theme_bw())

Ia.plot.method <- ggplot(Ia.df, aes(x = Sex.Rate, y = Percent_Reject, color=Samp.Size))
Ia.plot.method <- Ia.plot.method + geom_point() + geom_hline(aes(yintercept=0.050), color="red") + 
  geom_hline(aes(yintercept=0.950), color="blue") + facet_grid(Method ~ .) +
  labs(title=expression(paste("Probability of Rejecting the Null Hypothesis with ", I[A]))) + 
  xlab("Sex Rate") + ylab("Probability of Null Hypothesis Rejection") + 
  annotate("text", x = 1, y = 0.10, label="5% rejection", color="red") +
  annotate("text", x = 10, y = 0.90, label="95% rejection", color="blue")
print(Ia.plot.method + geom_line(aes(group = Samp.Size)) + theme_bw())
dev.off()

#==============================================================================#
#=====================================ROC PLOTS================================#
#alphacut <- c(seq(0,0.2,0.001), seq(0.3,1,0.1))
# alphacut <- seq(0,1,0.001)
# nalpha <- length(alphacut)
# 
# 
# cat("Constructing Ia ROC data. This will take a while....\n")
# print(system.time(Ia.tot2 <- lapply(alphacut, function(x) splitbymethod(meth, "Ia", alpha=x))))
# cat("Constructing rbarD ROC data. Again, this will take a while...\n")
# print(system.time(rd.tot2 <- lapply(alphacut, function(x) splitbymethod(meth, "rbarD", alpha=x))))
# Ia.comb2 <- NULL
# rd.comb2 <- NULL
# cat("Combining data sets...\n")
# print(system.time(lapply(Ia.tot2, function(y) lapply(y, function(x) Ia.comb2 <<- rbind(Ia.comb2, x)))))
# print(system.time(lapply(rd.tot2, function(y) lapply(y, function(x) rd.comb2 <<- rbind(rd.comb2, x)))))
# 
# 
# 
# Ia.comb2 <- as.vector(Ia.comb2/1000)
# Ia.df2 <- data.frame(list(Percent_Reject = Ia.comb2, 
#                          Samp.Size = rep(rep(levels(meth$Samp.Size), 10*length(methlev)), nalpha), 
#                          Sex.Rate = rep(rep(levels(meth$Sex.Rate), each=4*length(methlev)), each=nalpha), 
#                          Method = rep(rep(rep(levels(meth$Method), each=4), 10), nalpha))
# )
# Ia.df2$Sex.Rate <- factor(Ia.df2$Sex.Rate, levels(Ia.df2$Sex.Rate)[c(1,9:10,2:8)])
# Ia.df2$Samp.Size <- factor(Ia.df2$Samp.Size, levels(Ia.df2$Samp.Size)[c(1,3:4,2)])
# Ia.df2$Method <- factor(Ia.df2$Method, methlev)
# Ia.df2$alpha <- rep(rep(alphacut, each = 16), 10)
# Ia.df2n <- Ia.df2[Ia.df2$Sex.Rate == 1, ]
# Ia.df2 <- Ia.df2[Ia.df2$Sex.Rate != 1, ]
# #nullidx <- unlist(lapply(unique(Ia.df2n$alpha), function(x) rep(which(Ia.df2n$alpha == x), 9)))
# #Ia.df2$Null_Percent <- Ia.df2n$Percent_Reject[nullidx]
# Ia.df2$Null_Percent <- rep(Ia.df2n$Percent_Reject, 9)
# write.csv(Ia.df2, file="Ia_ROC.csv", row.names=FALSE)
# 
# rd.comb2 <- as.vector(rd.comb2/1000)
# rd.df2 <- data.frame(list(Percent_Reject = rd.comb2, 
#                           Samp.Size = rep(rep(levels(meth$Samp.Size), 10*length(methlev)), nalpha), 
#                           Sex.Rate = rep(rep(levels(meth$Sex.Rate), each=4*length(methlev)), each=nalpha), 
#                           Method = rep(rep(rep(levels(meth$Method), each=4), 10), nalpha))
# )
# rd.df2$Sex.Rate <- factor(rd.df2$Sex.Rate, levels(rd.df2$Sex.Rate)[c(1,9:10,2:8)])
# rd.df2$Samp.Size <- factor(rd.df2$Samp.Size, levels(rd.df2$Samp.Size)[c(1,3:4,2)])
# rd.df2$Method <- factor(rd.df2$Method, methlev)
# rd.df2$alpha <- rep(rep(alphacut, each = 16), 10)
# rd.df2n <- rd.df2[rd.df2$Sex.Rate == 1, ]
# rd.df2 <- rd.df2[rd.df2$Sex.Rate != 1, ]
# #nullidx <- unlist(lapply(unique(rd.df2n$alpha), function(x) rep(which(rd.df2n$alpha == x), 9)))
# #rd.df2$Null_Percent <- rd.df2n$Percent_Reject[nullidx]
# rd.df2$Null_Percent <- rep(rd.df2n$Percent_Reject, 9)
# write.csv(rd.df2, file="rd_ROC.csv", row.names=FALSE)


#==============================================================================#
# Reading in and printing ROC curves
#==============================================================================#

Ia.df2 <- read.csv("Ia_ROC.csv", header=TRUE)
Ia.df2$Sex.Rate <- factor(Ia.df2$Sex.Rate, levels(meth$Sex.Rate))
Ia.df2$Samp.Size <- factor(Ia.df2$Samp.Size, levels(meth$Samp.Size))
Ia.df2$Method <- factor(Ia.df2$Method, methlev)

library(pracma)
samlev <- levels(Ia.df2$Samp.Size)

Ia.plot <- ggplot(Ia.df2, aes(x = Null_Percent, y = Percent_Reject, xlim=1))
Ia.plot <- Ia.plot + 
  # geom_point(aes(color=Method, alpha=alpha)) + 
  geom_line(aes(color=Method, linetype=Samp.Size)) + 
  facet_wrap(~ Sex.Rate, nrow=3) +
  # facet_grid(Sex.Rate ~ Samp.Size) +
  labs(title=expression(paste("ROC analysis of ", I[A]))) + 
  xlab("False Positive Fraction") + ylab("True Positive Fraction")
print(Ia.plot + theme_classic())
ggsave("ROCIa.pn", width=12, height=6, units="in")
rd.df2 <- read.csv("rd_ROC.csv", header=TRUE)
rd.df2$Sex.Rate <- factor(rd.df2$Sex.Rate, levels(meth$Sex.Rate))
rd.df2$Samp.Size <- factor(rd.df2$Samp.Size, levels(meth$Samp.Size))
rd.df2$Method <- factor(rd.df2$Method, methlev)

rd.plot <- ggplot(rd.df2, aes(x = Null_Percent, y = Percent_Reject, xlim=1))
rd.plot <- rd.plot + 
#  geom_point(aes(color=Method, alpha=alpha)) + 
  geom_line(aes(color=Method, linetype=Samp.Size)) + 
  facet_wrap(~ Sex.Rate, nrow=3) +
#  facet_grid(Sex.Rate ~ Samp.Size) +
  labs(title=expression(paste("ROC analysis of ", bar(r)[D]))) + 
  xlab("False Positive Fraction") + ylab("True Positive Fraction")
print(rd.plot + theme_classic())
ggsave("ROCrd.png", width=12, height=6, units="in")

#==============================================================================#
# Constructing data for ANOVA of the area under the ROC curve
#==============================================================================#

AUC.rd <- lapply(levels(rd.df2$Sex.Rate), function(w) sapply(methlev, function(x) sapply(samlev, function(y) trapz(rd.df2[rd.df2$Samp.Size == y & rd.df2$Method == x & rd.df2$Sex.Rate == w, ]$Null_Percent, rd.df2[rd.df2$Samp.Size == y & rd.df2$Method == x & rd.df2$Sex.Rate == w, ]$Percent_Reject))))

names(AUC.rd) <- levels(rd.df2$Sex.Rate)
AUC.rd.df <- data.frame(list(AUC = unlist(AUC.rd), Method = rep(rep(methlev, each=4), 10), Sex.Rate=rep(names(AUC.rd), each=16), Samp.Size = rep(rep(samlev, 4), 10)))

AUC.Ia <- lapply(levels(Ia.df2$Sex.Rate), function(w) sapply(methlev, function(x) sapply(samlev, function(y) trapz(Ia.df2[Ia.df2$Samp.Size == y & Ia.df2$Method == x & Ia.df2$Sex.Rate == w, ]$Null_Percent, Ia.df2[Ia.df2$Samp.Size == y & Ia.df2$Method == x & Ia.df2$Sex.Rate == w, ]$Percent_Reject))))

names(AUC.Ia) <- levels(Ia.df2$Sex.Rate)
AUC.Ia.df <- data.frame(list(AUC = unlist(AUC.Ia), Method = rep(rep(methlev, each=4), 10), Sex.Rate=rep(names(AUC.Ia), each=16), Samp.Size = rep(rep(samlev, 4), 10)))

AUC.rd.df <- read.table("rd_AUROCC.csv", sep = ",", head = TRUE)
AUC.Ia.df <- read.table("Ia_AUROCC.csv", sep = ",", head = TRUE)

Ia.TP <- aov(Percent_Reject ~ Sex.Rate + Samp.Size * Method, data=Ia.df2)
Ia.FP <- aov(Null_Percent ~ Sex.Rate + Samp.Size*Method, data=Ia.df2)
Ia.ROC.aov <- aov(AUC ~ Sex.Rate + Samp.Size + Method, data=AUC.Ia.df)
rd.TP <- aov(Percent_Reject ~ Sex.Rate + Samp.Size + Method, data=rd.df2)
rd.FP <- aov(Null_Percent ~ Sex.Rate + Samp.Size + Method, data=rd.df2)
rd.ROC.aov <- aov(AUC ~ Sex.Rate + Samp.Size + Method, data=AUC.rd.df)


#==============================================================================#
# Comparing the different methods using the original data set with the kruskal
# wallace test.
#

kruskalia.samp <- sapply(levels(meth2$Sex.Rate), function(x) sapply(samlev, function(y){
  kruskal.test(meth2[meth2$Sex.Rate == x & meth2$Samp.Size == y, ]$p.Ia, 
               meth2[meth2$Sex.Rate == x & meth2$Samp.Size == y, ]$Method)$p.value
}))

kruskalrd.samp <- sapply(levels(meth2$Sex.Rate), function(x) sapply(samlev, function(y){
  kruskal.test(meth2[meth2$Sex.Rate == x & meth2$Samp.Size == y, ]$p.rD, 
               meth2[meth2$Sex.Rate == x & meth2$Samp.Size == y, ]$Method)$p.value
}))

library(lawstat)
levene.boot <- c(
  `AUC for Ia` = levene.test(AUC.Ia.df$AUC, AUC.Ia.df$Method)$p.value,
  `AUC for rbarD` = levene.test(AUC.rd.df$AUC, AUC.rd.df$Method)$p.value,
  `False Positive rbarD` = levene.test(rd.df2$Null_Percent, rd.df2$Method)$p.value,
  `True Positive rbarD` = levene.test(rd.df2$Percent_Reject, rd.df2$Method)$p.value,
  `False Positive Ia` = levene.test(Ia.df2$Null_Percent, Ia.df2$Method)$p.value,
  `True Positive Ia` = levene.test(Ia.df2$Percent_Reject, Ia.df2$Method)$p.value,
  `p-value Ia` = levene.test(meth2$p.Ia, meth2$Method)$p.value,
  `p-value rbarD` = levene.test(meth2$p.rD, meth2$Method)$p.value
)

kruskal.boot <- c(
  `AUC for Ia` = kruskal.test(AUC.Ia.df$AUC, AUC.Ia.df$Method)$p.value,
  `AUC for rbarD` = kruskal.test(AUC.rd.df$AUC, AUC.rd.df$Method)$p.value,
  `False Positive rbarD` = kruskal.test(rd.df2$Null_Percent, rd.df2$Method)$p.value,
  `True Positive rbarD` = kruskal.test(rd.df2$Percent_Reject, rd.df2$Method)$p.value,
  `False Positive Ia` = kruskal.test(Ia.df2$Null_Percent, Ia.df2$Method)$p.value,
  `True Positive Ia` = kruskal.test(Ia.df2$Percent_Reject, Ia.df2$Method)$p.value,
  `p-value Ia` = kruskal.test(meth2$p.Ia, meth2$Method)$p.value,
  `p-value rbarD` = kruskal.test(meth2$p.rD, meth2$Method)$p.value
)


levene.no.boot <- c(
  `AUC for Ia` = levene.test(AUC.Ia.df[AUC.Ia.df$Method != methlev[4], ]$AUC, AUC.Ia.df[AUC.Ia.df$Method != methlev[4], ]$Method)$p.value,
  `AUC for rbarD` = levene.test(AUC.rd.df[AUC.rd.df$Method != methlev[4], ]$AUC, AUC.rd.df[AUC.rd.df$Method != methlev[4], ]$Method)$p.value,
  `False Positive rbarD` = levene.test(rd.df2[rd.df2$Method != methlev[4], ]$Null_Percent, rd.df2[rd.df2$Method != methlev[4], ]$Method)$p.value,
  `True Positive rbarD` = levene.test(rd.df2[rd.df2$Method != methlev[4], ]$Percent_Reject, rd.df2[rd.df2$Method != methlev[4], ]$Method)$p.value,
  `False Positive Ia` = levene.test(Ia.df2[Ia.df2$Method != methlev[4], ]$Null_Percent, Ia.df2[Ia.df2$Method != methlev[4], ]$Method)$p.value,
  `True Positive Ia` = levene.test(Ia.df2[Ia.df2$Method != methlev[4], ]$Percent_Reject, Ia.df2[Ia.df2$Method != methlev[4], ]$Method)$p.value,
  `p-value Ia` = levene.test(meth2[meth2$Method != methlev[4], ]$p.Ia, meth2[meth2$Method != methlev[4], ]$Method)$p.value,
  `p-value rbarD` = levene.test(meth2[meth2$Method != methlev[4], ]$p.rD, meth2[meth2$Method != methlev[4], ]$Method)$p.value
)

kruskal.no.boot <- c(
`AUC for Ia` = kruskal.test(AUC.Ia.df[AUC.Ia.df$Method != methlev[4], ]$AUC, AUC.Ia.df[AUC.Ia.df$Method != methlev[4], ]$Method)$p.value,
`AUC for rbarD` = kruskal.test(AUC.rd.df[AUC.rd.df$Method != methlev[4], ]$AUC, AUC.rd.df[AUC.rd.df$Method != methlev[4], ]$Method)$p.value,
`False Positive rbarD` = kruskal.test(rd.df2[rd.df2$Method != methlev[4], ]$Null_Percent, rd.df2[rd.df2$Method != methlev[4], ]$Method)$p.value,
`True Positive rbarD` = kruskal.test(rd.df2[rd.df2$Method != methlev[4], ]$Percent_Reject, rd.df2[rd.df2$Method != methlev[4], ]$Method)$p.value,
`False Positive Ia` = kruskal.test(Ia.df2[Ia.df2$Method != methlev[4], ]$Null_Percent, Ia.df2[Ia.df2$Method != methlev[4], ]$Method)$p.value,
`True Positive Ia` = kruskal.test(Ia.df2[Ia.df2$Method != methlev[4], ]$Percent_Reject, Ia.df2[Ia.df2$Method != methlev[4], ]$Method)$p.value,
`p-value Ia` = kruskal.test(meth2[meth2$Method != methlev[4], ]$p.Ia, meth2[meth2$Method != methlev[4], ]$Method)$p.value,
`p-value rbarD` = kruskal.test(meth2[meth2$Method != methlev[4], ]$p.rD, meth2[meth2$Method != methlev[4], ]$Method)$p.value
)

levene.df <- data.frame(list(`With Non-Parametric Bootstrap` = levene.boot, `Without Non-Parametric Bootstrap` = levene.no.boot))

kruskal.df <- data.frame(list(`With Non-Parametric Bootstrap` = kruskal.boot, `Without Non-Parametric Bootstrap` = kruskal.no.boot))






summary(Ia.TP)
summary(Ia.FP)
summary(Ia.ROC.aov)
summary(rd.TP)
summary(rd.FP)
summary(rd.ROC.aov)
#==============================================================================#

#==============================================================================#





#==============================================================================#
#==============================================================================#
#==============================================================================#



# Creating plots looking at the differences between max values of distribution
# and observed values.

rddiff <- ggplot(meth, aes(x = (rbarD - Rd.max), y = rbarD, color = Sex.Rate))

Iadiff <- ggplot(meth, aes(x = (Ia - Ia.max), y = Ia, color = Sex.Rate))

cat("Printing Observed vs. Distance from maximum value...\n")

png("Distances%d.png", width=1024, height=1225)
print(rddiff + geom_point() + facet_grid(Method ~ Samp.Size) + labs(title="Distance of observed value from theoretical distribution maximum") + ylab(expression(bar(r)[D])) + xlab(expression(bar(r)[D] - paste(scriptstyle(max), bar(r)[D]))) + theme_classic())

print(rddiff + geom_point() + facet_grid(Sex.Rate ~ Method + Samp.Size) + labs(title="Distance of observed value from theoretical distribution maximum") + ylab(expression(bar(r)[D])) + xlab(expression(bar(r)[D] - paste(scriptstyle(max), ,bar(r)[D]))) + theme_bw())


print(Iadiff + geom_point() + facet_grid(Method ~ Samp.Size) + labs(title="Distance of observed value from theoretical distribution maximum") + ylab(expression(I[A])) + xlab(expression(I[A] - paste(scriptstyle(max), I[A]))) + theme_bw())

print(Iadiff + geom_point() + facet_grid(Sex.Rate ~ Method + Samp.Size) + labs(title="Distance of observed value from theoretical distribution maximum") + ylab(expression(I[A])) + xlab(expression(I[A] - paste(scriptstyle(max), I[A]))) + theme_bw())
dev.off()

# Standard Deviations from the mean
rddiff_sd <- ggplot(meth, aes(x = (rbarD / Rd.sd), y = (rbarD), color=Sex.Rate))
Iadiff_sd <- ggplot(meth, aes(x = (Ia / Ia.sd), y = (Ia), color=Sex.Rate))


cat("Printing plots of Observed vs. SD ratio...\n")

png("SDplots%d.png", width=1024, height=1225)
print(
rddiff_sd + geom_point() + facet_grid(Method ~ Samp.Size) + labs(title="Observed/SD(Observed) 
 vs. 
 Observed") + ylab(expression(bar(r)[D])) + xlab(expression(frac(bar(r)[D],SD(bar(r)[D])))) + theme_bw()
)

print(
rddiff_sd + geom_point() + facet_grid(Sex.Rate ~ Method + Samp.Size) + labs(title="Observed/SD(Observed) 
 vs. 
 Observed") + ylab(expression(bar(r)[D])) + xlab(expression(frac(bar(r)[D],SD(bar(r)[D])))) + theme_bw()
)

print(
Iadiff_sd + geom_point() + facet_grid(Method ~ Samp.Size) + labs(title="Observed/SD(Observed) 
 vs. 
 Observed") + ylab(expression(I[A])) + xlab(expression(frac(I[A],SD(I[A])))) + theme_bw()
)

print(
Iadiff_sd + geom_point() + facet_grid(Sex.Rate ~ Method + Samp.Size) + labs(title="Observed/SD(Observed) 
 vs. 
 Observed") + ylab(expression(I[A])) + xlab(expression(frac(I[A],SD(I[A])))) + theme_bw()
)
dev.off()


# Standard Deviations from the SD
rddiff_mn <- ggplot(meth, aes(x = (rbarD / Rd.sd), y = (Rd.sd), color=Sex.Rate))
Iadiff_mn <- ggplot(meth, aes(x = (Ia / Ia.sd), y = (Ia.sd), color=Sex.Rate))


cat("Printing SD vs SD Ratio...\n")

png("SDplots_2_%d.png", width=1024, height=1225)

print(rddiff_mn + geom_point() + facet_grid(Method ~ Samp.Size) + labs(title="Observed (units = standard deviations) \n vs. \n Standard Deviation") + ylab(expression(SD(bar(r)[D]))) + xlab(expression(frac(bar(r)[D],SD(bar(r)[D])))) + theme_bw())

print(rddiff_mn + geom_point() + facet_grid(Sex.Rate ~ Method + Samp.Size) + labs(title="Observed (units = standard deviations) \n vs. \n Standard Deviation") + ylab(expression(SD(bar(r)[D]))) + xlab(expression(frac(bar(r)[D],SD(bar(r)[D])))) + theme_bw())

print(Iadiff_mn + geom_point() + facet_grid(Method ~ Samp.Size) + labs(title="Observed (units = standard deviations) \n vs. \n Standard Deviation") + ylab(expression(SD(I[A]))) + xlab(expression(frac(I[A],SD(I[A])))) + theme_bw())

print(Iadiff_mn + geom_point() + facet_grid(Sex.Rate ~ Method + Samp.Size) + labs(title="Observed (units = standard deviations) \n vs. \n Standard Deviation") + ylab(expression(SD(I[A]))) + xlab(expression(frac(I[A],SD(I[A])))) + theme_bw())
dev.off()


# Analyzing at alpha == 0.05

meth95 <- meth[which(meth$p.rD == 0.05 | meth$p.Ia == 0.05), ]

rddiff95 <- ggplot(meth95, aes(x = (rbarD - Rd.max), y = rbarD, color = Sex.Rate))
Iadiff95 <- ggplot(meth95, aes(x = (Ia - Ia.max), y = Ia, color = Sex.Rate))
rddiff95_sd <- ggplot(meth95, aes(x = (rbarD / Rd.sd), y = (Rd.sd), color=Sex.Rate))
Iadiff95_sd <- ggplot(meth95, aes(x = (Ia / Ia.sd), y = (Ia.sd), color=Sex.Rate))


ggsave(rddiff95 + geom_point() + facet_grid(Method ~ Samp.Size) + labs(title="Distance of observed value from theoretical distribution maximum \n at p = 0.05") + ylab(expression(bar(r)[D])) + xlab(expression(bar(r)[D] - paste(scriptstyle(max), ,bar(r)[D]))) + theme_bw(), file="rddiff95.png", width=8, height=9.57)
ggsave(
Iadiff95 + geom_point() + facet_grid(Method ~ Samp.Size) + labs(title="Distance of observed value from theoretical distribution maximum \n at p = 0.05") + ylab(expression(I[A])) + xlab(expression(I[A] - paste(scriptstyle(max), I[A]))) + theme_bw(), file="Iadiff95.png", width=8, height=9.57)

print(
rddiff95_sd + geom_point() + facet_grid(Method ~ Samp.Size) + labs(title="Observed/SD(Observed) 
 vs. 
 Observed") + ylab(expression(bar(r)[D])) + xlab(expression(frac(bar(r)[D],SD(bar(r)[D])))) + theme_bw()
)
print(
rddiff95_sd + geom_point() + facet_grid(Sex.Rate ~ Method + Samp.Size) + labs(title="Observed/SD(Observed) 
 vs. 
 Observed") + ylab(expression(bar(r)[D])) + xlab(expression(frac(bar(r)[D],SD(bar(r)[D])))) + theme_bw()
)
print(
Iadiff95_sd + geom_point() + facet_grid(Method ~ Samp.Size) + labs(title="Observed/SD(Observed) 
 vs. 
 Observed") + ylab(expression(I[A])) + xlab(expression(frac(I[A],SD(I[A])))) + theme_bw()
)
print(
Iadiff95_sd + geom_point() + facet_grid(Sex.Rate ~ Method + Samp.Size) + labs(title="Observed/SD(Observed) 
 vs. 
 Observed") + ylab(expression(I[A])) + xlab(expression(frac(I[A],SD(I[A])))) + theme_bw()
)


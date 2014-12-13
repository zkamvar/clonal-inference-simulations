library(poppr)

x <- read.table("deleteme.csv", sep = ",", head = TRUE)

xgc <- poppr:::read_allele_columns(x[-1], ploidy = 2)
xgc <- as.genclone(df2genind(xgc, sep = "/"))
sethierarchy(xgc) <- x[1]
addhierarchy(xgc) <- data.frame(logx = floor(log(x[[1]])))
setpop(xgc) <- ~logx
print(poppr(xgc))
poppr:::jackbootcomp(xgc)
print(locus_table(xgc))
replen <- as.numeric(sub("^Locus.{3,4}([[:digit:]]).+?$", "\\1", xgc@loc.names))
bruvo.boot(xgc, replen, tree = "nj", cutoff = 75)
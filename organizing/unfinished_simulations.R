x <- readLines("analysis/MISSING_DATASETS.txt")
out <- matrix(unlist(strsplit(x, "10000")), ncol = 2, byrow = TRUE)
expected <- unlist(lapply(c(paste0("0", 0:9), 10), function(i) paste0("../zhian_ssr_data/", out[, 1], i, "000", out[, 2])))
dirs <- sapply(paste0("../zhian_ssr_data/twenty_loci", 1:16), dir, full.names = TRUE)
dirs <- grep("pop", unlist(dirs, use.names = FALSE), value = TRUE)

have <- expected[expected %in% dirs]
gens <- as.integer(gsub("^.+?gen_([0-9]{5})_rep.+?$", "\\1", have))
nogen <- factor(gsub("gen_[0-9]{5}", "", have))


finals <- vapply(levels(nogen), function(x) max(which(nogen == x)), integer(1), USE.NAMES = FALSE)
have[finals]

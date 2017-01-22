library("magrittr")
runs    <- as.character(1:8)

twelve <- 5:8
thirteen <- 1:4
sexrate <- c("0.0000",
             "0.0001",
             "0.0005",
             "0.0010",
             "0.0050",
             "0.0100",
             "0.0500"
             )
twelve_seed <- expand.grid("more_alleles_high_mutation", twelve,
                           "/seed_", as.character(0:11), # twelve seeds
                           "_sex_", sexrate,
                           "_gen_10000",
                           "_rep_", paste0("0", as.character(0:9)), ".pop") %>%
  apply(1, paste, collapse = "")
thirteen_seed <- expand.grid("more_alleles_high_mutation", thirteen,
                             "/seed_", as.character(0:12), # thirteen seeds
                             "_sex_", sexrate,
                             "_gen_10000",
                             "_rep_", paste0("0", as.character(0:9)), ".pop") %>%
  apply(1, paste, collapse = "")
expected <- sort(c(twelve_seed, thirteen_seed)) # Expect 40,000


dirs <- sapply(paste0("../zhian_ssr_data/more_alleles_high_mutation", 1:8), dir, full.names = TRUE)
dirs <- grep("pop", unlist(dirs, use.names = FALSE), value = TRUE)
dirs <- grep(paste0("(", sexrate[1], "|", sexrate[2], ")"), dirs, value = TRUE)

expected <- paste0("../zhian_ssr_data/", expected)
missing <- expected[!expected %in% dirs]

out <- matrix(unlist(strsplit(missing, "10000")), ncol = 2, byrow = TRUE)
expected <- unlist(lapply(c(paste0("0", 0:9), 10), function(i) paste0(out[, 1], i, "000", out[, 2])))

have <- expected[expected %in% dirs]
gens <- as.integer(gsub("^.+?gen_([0-9]{5})_rep.+?$", "\\1", have))
nogen <- factor(gsub("gen_[0-9]{5}", "", have))


finals <- vapply(levels(nogen), function(x) max(which(nogen == x)), integer(1), USE.NAMES = FALSE)
mutations <- c("1e-3", rep("1e-5", 20))
finalruns <- normalizePath(have[finals])
run <- paste(normalizePath("simulations/recover_ssr.py"),
             "--murate", paste(mutations, collapse = " "),
             "--popfile", finalruns)

cat(run, sep = "\n", file = "analysis/more_alleles_high_mutation_rerun.txt")

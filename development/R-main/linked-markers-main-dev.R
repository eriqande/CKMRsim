library(dplyr)
library(CKMRsim)



# make some data with 20 chromosomes, each of 100 Mb and positions
# on each random
set.seed(1)
linked_mhaps <- mhaps %>%
  group_by(Chrom, Locus) %>%
  mutate(NewChrom = sample(1:20, 1),
         NewPos = floor(runif(1, 1, 10^8))
  ) %>%
  ungroup() %>%
  mutate(Chrom = NewChrom,
         Pos = NewPos) %>%
  select(-starts_with("New")) %>%
  CKMRsim::reindex_markers()



# this is what I am using while working on the linked markers
mhlist <- long_markers_to_X_l_list(D = linked_mhaps,
                                   kappa_matrix = kappas[c("PO", "FS", "U"), ])

# add the matrices that account for genotyping error
mhlist2 <- insert_C_l_matrices(mhlist,
                               snp_err_rates = 0.005,
                               scale_by_num_snps = TRUE)

# do the matrix multiplication that gives the Y_l matrices
mhlist3 <- insert_Y_l_matrices(mhlist2)


# now we do the linked simulation stuff.  For that we are going to need the number
# of alleles in the proper order.
alle_nums <- linked_mhaps %>%
  mutate(list_name = paste(Chrom, Locus, Pos, sep = ".")) %>%
  mutate(list_name = factor(list_name, levels = unique(list_name))) %>%
  group_by(list_name) %>%
  tally()


# we are doing do to linked markers for PO and FS
relats <- c("HS", "FS")
names(relats) <- relats

true_genos <- lapply(relats, function(x) {
  tmpDir = tempdir()
  write_all_mendel_files("mendel-example", 1000, floor(runif(1, min = 100, max = 100000)), linked_mhaps, pedigrees[[x]], Dir = tmpDir)
  run_mendel(tmpDir, "mendel-example-Control.in")
  outgenos <- read_mendel_outped(file.path(tmpDir, "mendel-example-Ped.out"), alle_nums$n) %>%
    lapply(function(x) {
      dimnames(x) = list(NULL, locus = as.character(alle_nums$list_name))
      x})
})

# at this juncture, we have the genotypes of each individual at all the loci, but we
# need to apply the true genotyping errors to them.  This we do by, for either HS or FS,
# lapplying over the loci. (as the names of mhlist2).
G1 <- true_genos$FS$indiv1
G2 <- true_genos$FS$indiv2
locs <- names(mhlist3)
names(locs) <- locs
results <- lapply(locs, function(n) {
  g1 <- G1[, n]
  g2 <- G2[, n]
  Clt <- mhlist2[[n]]$C_l_true
  g1e <- samp_from_mat(Clt[g1,])
  g2e <- samp_from_mat(Clt[g2,])

  # here is the number of genos
  nG <- nrow(Clt)

  # so, the position of the genotype of the two individuals in a matrix of
  # joint probabilities will be  nG * (g2 - 1) + g1
  jg_idx <- nG * (g2e - 1) + g1e

  spoogie <- list(g1=g1, g2=g2, Clt=Clt, g1e = g1e, g2e = g2e, jg_idx = jg_idx)



  # check the observed geno freqs:
  obsy <- as.data.frame(table(spoogie$jg_idx)/length(spoogie$jg_idx)) %>%
    tbl_df() %>%
    mutate(alle = as.character(Var1)) %>%
    select(alle, Freq)

  # and compare against the joint probs
  tmp <- as.numeric(mhlist3[[n]]$Y_l_true$FS)
  names(tmp) <- 1:length(tmp)
  predi <- data_frame(alle = names(tmp), Freq_pred = tmp)

  full_join(predi, obsy)
}) %>%
  bind_rows(.id = "Locus")


# then plot that to make sure our simulated values are coming out like
# we have predicted.
ggplot(results, aes(x = Freq_pred, y = Freq)) +
  geom_point(colour = "blue") +
  geom_abline(intercept = 0, slope = 1)

library(dplyr)
library(CKMRsim)
library(ggplot2)





# this is what I am using while working on the linked markers
mhlist <- long_markers_to_X_l_list(D = linked_mhaps,
                                   kappa_matrix = kappas[c("PO", "FS", "U"), ])

# add the matrices that account for genotyping error
mhlist2 <- insert_C_l_matrices(mhlist,
                               snp_err_rates = 0.005,
                               scale_by_num_snps = TRUE)

# do the matrix multiplication that gives the Y_l matrices
mhlist3 <- insert_Y_l_matrices(mhlist2)



# now we want to simulate the linked genos from Mendel, let's say for Full siblings
FSs <- sample_linked_genotype_pairs(df = linked_mhaps, ped = pedigrees$FS, C = mhlist3, num = 1000)


# that is really all there is to it.  Now, it would be cool to confirm that
# simulated genotype frequencies are marginally what we expect.
fs_probs <- lapply(mhlist3, function(x) {y <- as.numeric(x$Y_l_true$FS); data_frame(geno_idx = 1:length(y), prob = y)}) %>%
  bind_rows(.id = "locus")

fs_freqs <- lapply(FSs, function(x) data_frame(geno_idx = x)) %>% bind_rows(.id = "locus") %>%
  group_by(locus, geno_idx) %>%
  tally() %>%
  mutate(freq = n / sum(n))

# now stick those together and see what we get:
tmp <- full_join(fs_probs, fs_freqs) %>%
  mutate(freq = ifelse(is.na(freq), 0, freq))

ggplot(tmp, aes(x = prob, y = freq)) +
  geom_point(colour = "blue") +
  geom_abline(intercept = 0, slope = 1)



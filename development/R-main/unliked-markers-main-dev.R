library(dplyr)
library(CKMRsim)
library(ggplot2)





# read in the "linked mhaps" but treat them as unlinked
mhlist <- long_markers_to_X_l_list(D = linked_mhaps,
                                   kappa_matrix = kappas[c("MZ", "PO", "FS", "HS", "U"), ])

# add the matrices that account for genotyping error
mhlist2 <- insert_C_l_matrices(mhlist,
                               snp_err_rates = 0.005,
                               scale_by_num_snps = TRUE)

# do the matrix multiplication that gives the Y_l matrices
mhlist3 <- insert_Y_l_matrices(mhlist2)


# simulate values from each relationship, assuming unlinked. This takes
# about 12 seconds
Qvals <- simulate_and_calc_Q(mhlist3, reps = 10^4)


# collect the MC averages
mc_sample_simple(Qvals, nu = c("PO", "FS", "HS"), method = "both")

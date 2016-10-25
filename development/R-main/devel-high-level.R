
# workflows while developing high-level S3-oriented functions

library(CKMRsim)

CK <- create_ckmr(linked_mhaps)



# could do it like this for simple unlinked
QU <- simulate_Qij(CK, froms = c("PO", "FS", "HS", "U"), tos = c("PO", "FS", "U"), reps = 500)

# or, to save typing and do all of the relationships:
QU1 <- simulate_Qij(CK, froms = ckmr_relats(CK), tos = ckmr_relats(CK), reps = 1000)

# or do it linked
QU2 <- simulate_Qij(CK, froms = "FS", tos = ckmr_relats(CK), reps = 100, unlinked = FALSE, pedigree_list = pedigrees)


# or linked with both HS and FS (and U and PO)
QU3 <- simulate_Qij(CK, froms = c("HS", "FS", "U", "PO"), tos = ckmr_relats(CK), reps = 750, unlinked = FALSE, pedigree_list = pedigrees)

# Same as above, but force PO to be simulated from MENDEL
QU4 <- simulate_Qij(CK, froms = c("HS", "FS", "U", "PO"), tos = ckmr_relats(CK), reps = 750, unlinked = FALSE,
                    pedigree_list = pedigrees, forceLinkagePO = TRUE)


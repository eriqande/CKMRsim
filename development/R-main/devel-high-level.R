
# workflows while developing high-level S3-oriented functions

library(CKMRsim)
library(ggplot2)

CK <- create_ckmr(linked_mhaps)



# could do it like this for simple unlinked
QU <- simulate_Qij(CK, froms = c("PO", "FS", "HS", "U"), tos = c("PO", "FS", "U"), reps = 500)

# or, to save typing and do all of the relationships:
QU1 <- simulate_Qij(CK, froms = ckmr_relats(CK), tos = ckmr_relats(CK), reps = 1000)

# or do it linked
QU2 <- simulate_Qij(CK, froms = "FS", tos = ckmr_relats(CK), reps = 1250, unlinked = FALSE, pedigree_list = pedigrees)


# or linked with both HS and FS (and U and PO)
QU3 <- simulate_Qij(CK, froms = c("HS", "FS", "U", "PO"), tos = ckmr_relats(CK), reps = 750, unlinked = FALSE, pedigree_list = pedigrees)

# Same as above, but force PO to be simulated from MENDEL
QU4 <- simulate_Qij(CK, froms = c("HS", "FS", "U", "PO"), tos = ckmr_relats(CK), reps = 750, unlinked = FALSE,
                    pedigree_list = pedigrees, forceLinkagePO = TRUE)



## Or, if we wanted to go after the this for parentage and full sibs only for rockfish:
RF <- simulate_Qij(CK, froms = c("PO", "U", "FS"), tos = c("PO", "U", "FS"), reps = 10000)

RF_linked <- simulate_Qij(CK, froms = c("PO", "U", "FS"), tos = c("PO", "U", "FS"), reps = 10000, unlinked = FALSE, forceLinkagePO = TRUE, pedigree_list = pedigrees)

# then extract and plot those
df <- bind_rows(
  extract_logls(RF, numer = c(PO=1), denom = c(U=1)),
  extract_logls(RF_linked, numer = c(PO=1), denom = c(U=1)),
  extract_logls(RF, numer = c(FS=1), denom = c(PO=1)),
  extract_logls(RF_linked, numer = c(FS=1), denom = c(PO=1)),
  extract_logls(RF, numer = c(FS=1), denom = c(U=1)),
  extract_logls(RF_linked, numer = c(FS=1), denom = c(U=1))
)


ggplot(df, aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7) +
  facet_grid(numer_wts + simtype ~ denom_wts)


# that is pretty cool.  It shows that full sibs can be hard to distinguish from PO, surprisingly...

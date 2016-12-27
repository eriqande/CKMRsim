
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



## If we wanted to simulate a certain amount of missing data we could do:
QU_miss1 <- simulate_Qij(CK, froms = c("PO", "FS", "HS", "U"), tos = c("PO", "FS", "U"), reps = 500,  rando_miss_wts = 1/(1:4), rando_miss_n = 10)




## Or, if we wanted to go after the this for parentage and full sibs only for rockfish:
RF <- simulate_Qij(CK, froms = c("PO", "U", "FS"), tos = c("PO", "U", "FS"), reps = 10000)

RF_linked <- simulate_Qij(CK, froms = c("PO", "U", "FS"), tos = c("PO", "U", "FS"), reps = 10000, unlinked = FALSE, forceLinkagePO = TRUE, pedigree_list = pedigrees)

# then extract and plot those
df <- dplyr::bind_rows(
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

# now, let's just check that we are getting reasonable results for the missing data part
RFmiss10 <- simulate_Qij(CK, froms = c("PO", "U", "FS"), tos = c("PO", "U", "FS"), reps = 10000, rando_miss_wts = 1/(1:4), rando_miss_n = 10)
RFmiss20 <- simulate_Qij(CK, froms = c("PO", "U", "FS"), tos = c("PO", "U", "FS"), reps = 10000, rando_miss_wts = 1/(1:4), rando_miss_n = 20)
RFmiss100 <- simulate_Qij(CK, froms = c("PO", "U", "FS"), tos = c("PO", "U", "FS"), reps = 10000, rando_miss_wts = 1/(1:4), rando_miss_n = 100)

df2 <- dplyr::bind_rows(
  extract_logls(RF, numer = c(PO=1), denom = c(U=1)),
  extract_logls(RFmiss10, numer = c(PO=1), denom = c(U=1)),
  extract_logls(RFmiss20, numer = c(PO=1), denom = c(U=1)),
  extract_logls(RFmiss100, numer = c(PO=1), denom = c(U=1))
)

ggplot(df2, aes(x = logl_ratio, fill = as.factor(rando_miss_n))) +
  geom_density(alpha = 0.5)

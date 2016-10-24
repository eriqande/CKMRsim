
# workflows while developing high-level S3-oriented functions

library(CKMRsim)

CK <- create_ckmr(microhaps)

# could do it like this for simple unlinked
QU <- simulate_Qij(CK, froms = c("PO", "FS", "HS", "U"), tos = c("PO", "FS", "U"))


# or, to save typing and do all of the relationships:
QU <- simulate_Qij(CK, froms = ckmr_relats(CK), tos = ckmr_relats(CK), reps = 100)

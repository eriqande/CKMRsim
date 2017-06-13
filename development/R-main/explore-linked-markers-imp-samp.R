
library(tidyverse)
library(CKMRsim)
library(ggplot2)


# this is how we go about doing importance sampling with linked markers:
# we simulate unlinked, and then compute the FNRs with the linked markers.
CK <- create_ckmr(linked_mhaps)


QU_unlinked <- simulate_Qij(CK, froms = c("U", "FS"), tos = ckmr_relats(CK), reps = 5e4)
QU_linked <- simulate_Qij(CK, froms = c("FS"), tos = ckmr_relats(CK), reps = 5e04, unlinked = FALSE, pedigree_list = pedigrees)


result <- mc_sample_simple(Q = QU_unlinked, nu = "FS", de = "U", method = "both", Q_for_fnrs = QU_linked)



# let's just plot the logls for fun, too
logls <- list(
unlinked = extract_logls(QU_unlinked, numer = c(FS = 1), denom = c(U = 1)),
linked = extract_logls(QU_linked, numer = c(FS = 1), denom = c(U = 1))
) %>%
  bind_rows(.id = "linkage") %>%
  filter(true_relat == "FS")


ggplot(logls, aes(x = logl_ratio, colour = linkage)) +
  geom_density()



# This is some code for testing out the indiv_comparison funcs.
library(readr)



# get the allele freqs and stuff
CK <- create_ckmr(linked_mhaps)

# flatten those out
po_flat <- flatten_ckmr(CK, "PO")
unrel_flat <- flatten_ckmr(CK, "U")

# then compute the log-likelihoods for the parent offspring relationship
po_logl_flat <- po_flat
po_logl_flat$probs <- log(po_flat$probs / unrel_flat$probs)


## choose which spip output to use:
#GENO <- "development/data/spip400_geno.txt"
GENO <- "development/data/spip4K_geno.txt.gz"
PED <- "development/data/spip4K_ped.txt.gz"

#### now prepare the genotypes in the right order
# first get the genos
genos <- read_delim(GENO, delim = " ")

# then prepare a data frame with number of alleles
numa <- dplyr::data_frame(LocIdx = 1:length(CK$loci), NumAlle = sapply(CK$loci, function(x) length(x$freqs)))

# then compute the genotype index of each pair of alleles
wide_genos <- dplyr::left_join(genos, numa) %>%
  dplyr::mutate(GenoIdx = index_ab(alle1, alle2, NumAlle)) %>%
  dplyr::select(type, id, LocIdx, GenoIdx) %>%
  tidyr::spread(data = ., key = LocIdx, value = GenoIdx)


wide_juvies <- wide_genos %>% dplyr::filter(type == "juvie") %>% dplyr::select(-type)
wide_adults <- wide_genos %>% dplyr::filter(type == "pop") %>% dplyr::select(-type)

juvie_mat <- as.matrix(wide_juvies[,-1]) - 1
storage.mode(juvie_mat) <- "integer"
rownames(juvie_mat) <- wide_juvies$id


adult_mat <- as.matrix(wide_adults[,-1]) - 1
storage.mode(adult_mat) <- "integer"
rownames(adult_mat) <- wide_adults$id

## I really should create some missing data in here to make sure it is working!!
set.seed(55)
admask <- sample(1:length(adult_mat), size = floor(length(adult_mat) / 50))
jumask <- sample(1:length(juvie_mat), size = floor(length(juvie_mat) / 50))

adult_mat[admask] <- -1
juvie_mat[jumask] <- -1

#### Now, with that stuff in hand we should be ready to try out our functions
idx <- 1:nrow(juvie_mat)
names(idx) <- idx
system.time({boing <- lapply(idx, function(i) {
    tmp <- comp_ind_pairwise(S = adult_mat, T = juvie_mat, t = i, values = po_logl_flat$probs, nGenos = po_logl_flat$nGenos, Starts = po_logl_flat$base0_locus_starts)
    tmp[rev(top_index(tmp$value, 5)), ]  # just take the top 5 from each
    }) %>%
  dplyr::bind_rows(.id = "offspring") %>%
  dplyr::tbl_df()}
  )

# with roughly 4K adults and 4K offspring, this takes about 45 seconds to a minute on my laptop.  Cool!


# now get the rough list of candidates as those with a logl > 0
candi <- boing %>% dplyr::filter(value > 0) %>%
  dplyr::mutate(offspring = as.integer(offspring))


# and compare these to the true pedigree
ped <- read_delim(PED, delim = " ", col_types = "ccc")

check_em <- candi %>%
  dplyr::mutate(off_name = rownames(juvie_mat)[offspring],
         adult_name = rownames(adult_mat)[ind]) %>%
  dplyr::left_join(ped, by = c("off_name" = "kid")) %>%
  dplyr::mutate(correct = pa == adult_name | ma == adult_name)

# and note that the ones that are wrong have very low logl:
check_em %>% dplyr::filter(correct == FALSE)

# booyah! (of course, there is not genotyping error here, so you expect that it will totally kill.)

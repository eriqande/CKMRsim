library(dplyr)

# quickly make some pedigrees
# name the target individuals 1 and 2.  Missing data are 0's
pedigrees <- list(
  FS = data.frame(
    Kid = c(1, 2, 3, 4),
    Pa = c(3, 3, 0, 0),
    Ma = c(4, 4, 0, 0),
    Sex = c(1, 1, 1, 2),
    Observed = c(1, 1, 0, 0)),
  HS = data.frame(
    Kid = c(1, 2, 3, 4, 5),
    Pa = c(4, 5, 0, 0, 0),
    Ma = c(3, 3, 0, 0, 0),
    Sex = c(1, 1, 2, 1, 1),
    Observed = c(1, 1, 0, 0, 0))
)



# quickly make some mapped marker data which we will put
# onto 5 different chromosomes named 4,9,14,19,and 21
d <- c(4,9,14,19,21)
names(d) <- d
# put n markers on each, where n = 200/the index.  i.e. imagine that they are getting shorter.
# let the length of each chromosome, in morgans be 1000 cM / the index.
set.seed(15)
markers_on_map <- lapply(d, function(x) {

  cm <- 1000 / x
  n <- round(20000/x)

  # simulate positions for them
  pos <- sort(round(runif(n, min = 0, max = cm), digits = 3))

  names(pos) <- paste("chr", x, "-", 1:length(pos), sep ="")

  # simulate allele freqs for them and put those along with the postion in a list
  ret <- lapply(pos, function(x) {
    rr <- list()
    rr$pos <- x
    A <- sample(2:8, 1) # number of alleles
    rg <- rgamma(n = A, shape = 2, scale = 1)
    rr$freqs <- round(rg/sum(rg), digits = 5)
    rr
  })

  ret

})


# let's explore a long format for these guys
long_markers <- lapply(markers_on_map, function(x) {
  lapply(x, function(y) {
    data.frame(Pos = y$pos, Allele = paste("a", 1:length(y$freqs), sep = ""), Freq = y$freqs, stringsAsFactors = FALSE)
  }) %>% bind_rows(.id = "Locus")
}) %>% bind_rows(.id = "Chrom")

long_markers$Chrom <- as.integer(long_markers$Chrom)

# here  we provide an index for each allele at each locus and an index for each locus
long_markers2 <- long_markers %>%
  group_by(Chrom, Locus) %>%
  mutate(AlleIdx = as.integer(1 + n_distinct(Allele) - rank(Freq, ties.method = "first"))) %>%
  group_by(Chrom) %>%
  mutate(LocIdx = as.integer(factor(Pos, levels = sort(unique(Pos))))) %>%
  arrange(Chrom, LocIdx, AlleIdx)


# and now, if we want to efficiently write these to a Morgan file we
# will want to bung the allele freqs together into a string
long_markers2 %>%
  group_by(Chrom, LocIdx) %>%
  summarise(FreqString = paste(Freq, collapse = " ")) %>%
  mutate(setchrom = "set chrom", markers = "markers", allefreqs = "allele freqs") %>%
  select(setchrom, Chrom, markers, LocIdx, allefreqs, FreqString) %>%
  write.table(., file = "/tmp/testit.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# so, when we implement this we will require the user to supply a long format data frame
# that has Chrom (int), LocIdx (int), Pos (dbl),  AlleIdx (int), and Freq (dbl)
# then we don't have to worry about doing the ranking of the alleles, etc---the user will
# have to take care of that.  The AlleIdx will have to correspond to the indexing of the alleles that
# translates directly to the genotypes.  Maybe we should allow for columns of "LocusName"
# and "AlleleName"


# write out a matrix of kappa values for the different relationships we will look at.
kappas <- as.matrix(
  read.table(
    textConnection(
      "Relat  kappa0  kappa1  kappa2
MZ  0 0 1
PO  0 0 1
FS  0.25 0.5 0.25
HS  0.5  0.5  0
GP  0.5  0.5  0
AN  0.5  0.5  0
DFC 0.5625 0.375 0.0625
FC  0.75 0.25 0
HC  0.875 0.125 0
U   0 0 0
      "),
    header = TRUE, row.names = 1
  )
)

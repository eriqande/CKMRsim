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
U   1 0 0
      "),
    header = TRUE, row.names = 1
  )
)



#### Make some microhaplotype data ####
# we just simulate 200 chromosomes using ms and have either 1, 2, 3, 4, or 5 segregating
# sites.  We want 100 of these by the end, so let us do 20 of each number of seg sites
library(stringr)
system("
       source ~/.bashrc;
       echo \"41234 20641 60651\" > seedms;
       rm -f splud;
       for R in 1 2 3 4 5; do
         ms 200 20 -s $R >> splud
       done
")
x <- readLines("splud")
x2 <- x[!str_detect(x, "[a-z/2-9]")]
x3 <- x2[x2!=""]
y <- data.frame(hap01 = x3, stringsAsFactors = FALSE) %>%
  tbl_df
y$LocNum <- rep(1:100, each = 200)

# now me make 0 or 1 be ACGT for each locus.  We make an array for
# what the 0's and 1's mean.  Far out
set.seed(5)
DNA_types <- lapply(1:100, function(x) matrix(unlist(lapply(1:5, function(x) sample(c("A", "C", "G", "T"), 2 ))), byrow=T, ncol = 2))

# now make a column of DNA bases for each haplotype.  We need a function for that.
# this is klugie and not vectorized, but I only need to do it once.
hapseq <- function(x, L) {
  ivec <- as.numeric(str_split(x, "")[[1]]) + 1
  len = length(ivec)
  idx <- cbind(1:len, ivec)
  paste(DNA_types[[L]][idx], collapse = "")
}
# this is klugie but will work....
y$hap <- sapply(1:nrow(y), function(l) hapseq(y$hap01[l], y$LocNum[l]))


# check these:
y %>%
  group_by(LocNum, hap01, hap) %>%
  tally

# yep, it looks good.  So, now, make markers out of them
# we will have some garbage on the chrom name for all of them
y2 <- y %>%
  mutate(Chrom = "ddRAD",
         Locus = paste("ddRAD", LocNum, sep = "_"),
         Pos = LocNum
         ) %>%
  rename(LocIdx = LocNum,
         Allele = hap)



# now compute allele freqs
microhaps <- y2 %>%
  group_by(Chrom, Locus, Pos, LocIdx, Allele) %>%
  tally %>%
  mutate(Freq = n / sum(n)) %>%
  select(-n) %>%
  mutate(AlleIdx = NA) %>%
  reindex_markers()

save(microhaps, file = "data/microhaps.rda")

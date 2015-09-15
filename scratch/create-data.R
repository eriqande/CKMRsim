
# quickly make some pedigrees
pedigrees <- list(
  FS = data.frame(
    Kid = c("i", "j", 1, 2),
    Pa = c("1", "1", 0, 0),
    Ma = c("2", 2, 0, 0),
    Sex = c(0, 0, 1, 2),
    Observed = c(1, 1, 0, 0)),
  HS = data.frame(
    Kid = c("i", "j", 1, 2, 3),
    Pa = c("2", "3", 0, 0, 0),
    Ma = c("1", 1, 0, 0, 0),
    Sex = rep(0, 5),
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
  n <- round(200/x)

  # simulate positions for them
  pos <- sort(runif(n, min = 0, max = cm))

  names(pos) <- paste("chr", x, "-", 1:length(pos), sep ="")

  # simulate allele freqs for them and put those along with the postion in a list
  ret <- lapply(pos, function(x) {
    rr <- list()
    rr$pos <- x
    A <- sample(2:8, 1) # number of alleles
    rg <- rgamma(n = A, shape = 2, scale = 1)
    rr$freqs <- rg/sum(rg)
    rr
  })

  ret

})

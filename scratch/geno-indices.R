
# this is just to make sure that I got my genotype indexing right
library(dplyr)

A <- 15

# a function to compute the index from a and b
indab <- function(a, b, A) {
 2 + (a - 1) * (A + 2) - ( a * (a + 1) / 2) + (b - a)
}

expand.grid(1:A, 1:A)[,c(2,1)] %>%
  setNames(c("a", "b")) %>%
  filter(a <= b) %>%
  mutate(Idx = 1:length(a),
         Indab = indab(a, b, A))

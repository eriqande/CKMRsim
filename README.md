# CKMRsim

This is an R package in development.  Not much else here now.

NOTES:

When you define a pedigree for gene dropping, individuals 1 and 2 are the focal
individauls of the pair and they must have observed parents (because MENDEL output
is hard to parse otherwise...).  Hence PO must look like:
```
  Kid Pa Ma Sex Observed
1   1  2  3   1        1
2   2  4  5   1        1
3   3  0  0   2        0
4   4  0  0   1        0
5   5  0  0   2        0
```
instead of having just three lines.

#!/usr/bin/env Rscript

library(ape)

set.seed(1)
t <- rcoal(100)
write.tree(t, "test.tree")
t <- rcoal(100)
write.tree(t, "test.tree", append=TRUE)

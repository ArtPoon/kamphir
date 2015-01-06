#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

if (length(args) != 2) { stop('Usage: binarize_tree.R [input NWK] [output NWK]')}
input.nwk <- args[1]
output.nwk <- args[2]

if (!file.exists(input.nwk)) {
	stop('input file does not exist')
}
if (!file.exists(output.nwk)) {
	stop('output file does not exist')
}

require(ape, quietly=TRUE)
tree <- read.tree(input.nwk)

# assumes that tip date is last underscore-delimited field
temp <- sapply(tree$tip.label, function(x) strsplit(x, split='_')[[1]])
tip.dates <- as.double(temp[nrow(temp),])

root.tree <- rtt(tree, tip.dates)

write.tree(root.tree, file=output.nwk)

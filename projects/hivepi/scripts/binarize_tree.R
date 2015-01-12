#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

if (length(args) != 2) { stop('Usage: binarize_tree.R [input NWK] [output NWK]')}
input.nwk <- args[1]
output.nwk <- args[2]

if (!file.exists(input.nwk)) {
	stop('input file does not exist')
}

require(ape, quietly=TRUE)
tree <- read.tree(input.nwk)

# assumes that tip date is last underscore-delimited field
temp <- sapply(tree$tip.label, function(x) strsplit(x, split='_')[[1]][[3]])

#tip.dates <- as.double(temp[nrow(temp),])
tip.dates <- (as.double(temp)-1970) * 52.1775

root.tree <- rtt(tree, tip.dates)
root.tree$tip.label <- paste(root.tree$tip.label, as.character(tip.dates), sep='_')

write.tree(root.tree, file=output.nwk)

#~ example: R -f read_master.R --args HIVModel
#~ example: R -f read_master.R --args BDModel2type-2

require(phytools)

clargs <- commandArgs(trailingOnly = TRUE)

nexfile <- scan(paste(clargs[1], ".nexus",sep=''),what=list(character()))
enwk <- nexfile[[1]][7]
# delete reactions and times
nwk <- gsub(",reaction(.*?)\\]","",enwk)
# delete isolated times
nwk <- gsub(",time(.*?)\\]","",nwk)
# get rid of type info
nwk <- gsub("\\[&type=\"","",nwk)
# get rid of quotes
nwk <- gsub("\"","",nwk)

tree <- read.newick(text=nwk)
tree <- collapse.singles(tree)
write.tree(tree, paste(clargs[1], ".nwk",sep='') )

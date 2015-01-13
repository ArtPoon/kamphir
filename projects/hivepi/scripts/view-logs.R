# load kamphir logs
n100 <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n100.RLRootToTip.timetree.log', header=T, sep='\t')
temp <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n100.RLRootToTip.timetree.log.7', header=T, sep='\t')
n100 <- rbind(n100[100:nrow(n100),], temp[1000:nrow(temp),])
temp <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n100.RLRootToTip.timetree.restart.7.log.1', header=T, sep='\t')
n100 <- rbind(n100, temp)

n300 <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n300.RLRootToTip.timetree.log.2', header=T, sep='\t')

n1000 <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n1000.RLRootToTip.timetree.log', header=T, sep='\t')
temp <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n1000.RLRootToTip.timetree.log.3', header=T, sep='\t')
n1000 <- rbind(n1000[100:nrow(n1000),], temp[100:nrow(temp),])
temp <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n1000.RLRootToTip.timetree.restart.log.2', header=T, sep='\t')
n1000 <- rbind(n1000, temp)

n100 <- n100[complete.cases(n100),]
n300 <- n300[complete.cases(n300),]
n1000 <- n1000[complete.cases(n1000),]

vioplot(log10(n100$N), log10(n300$N), log10(n1000$N), col='grey')
points(x=c(1,2,3), y=c(3, log10(3000), 4), cex=2, pch=4, lwd=3, col='yellow')

vioplot(log10(n100$beta), log10(n300$beta), log10(n1000$beta), col='grey')
points(x=c(1,2,3), y=c(-3, log10(0.0003), -4), cex=2, pch=4, lwd=3, col='yellow')

# load BEAST logs
b100 <- read.table('~/git/kamphir/projects/hivepi/beast/SIRTree.n100.log', header=T, sep='\t')
temp <- read.table('~/git/kamphir/projects/hivepi/beast/SIRTree.n100.1.log', header=T, sep='\t')
b100 <- rbind(b100[100:nrow(b100),], temp[100:nrow(temp),])

b300 <- read.table('~/git/kamphir/projects/hivepi/beast/SIRTree.n300.log', header=T, sep='\t')
temp <- read.table('~/git/kamphir/projects/hivepi/beast/SIRTree.n300.1.log', header=T, sep='\t')
b300 <- rbind(b300[100:nrow(b300),], temp[100:nrow(temp),])

b1000 <- read.table('~/git/kamphir/projects/hivepi/beast/SIRTree.n1000.log', header=T, sep='\t')
temp <- read.table('~/git/kamphir/projects/hivepi/beast/SIRTree.n1000.1.log', header=T, sep='\t')
b1000 <- rbind(b1000[100:nrow(b1000),], temp[100:nrow(temp),])

boxplot(b100$S0, b300$S0, b1000$S0, log='y')

vioplot(log10(b100$S0), log10(b300$S0), log10(b1000$S0), col='grey')
points(x=c(1,2,3), y=c(3, log10(3000), 4), cex=2, pch=4, lwd=1.5)



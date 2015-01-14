my.vioplot <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
    horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
    lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
    at, add = FALSE, wex = 1, drawRect = TRUE) 
{
	if(length(col)==1) col <- rep(col,n)
    datas <- list(x, ...)
    n <- length(datas)
    if (missing(at)) 
        at <- 1:n
    upper <- vector(mode = "numeric", length = n)
    lower <- vector(mode = "numeric", length = n)
    q1 <- vector(mode = "numeric", length = n)
    q3 <- vector(mode = "numeric", length = n)
    med <- vector(mode = "numeric", length = n)
    base <- vector(mode = "list", length = n)
    height <- vector(mode = "list", length = n)
    baserange <- c(Inf, -Inf)
    args <- list(display = "none")
    if (!(is.null(h))) 
        args <- c(args, h = h)
    for (i in 1:n) {
        data <- datas[[i]]
        data.min <- min(data)
        data.max <- max(data)
        q1[i] <- quantile(data, 0.25)
        q3[i] <- quantile(data, 0.75)
        med[i] <- median(data)
        iqd <- q3[i] - q1[i]
        upper[i] <- min(q3[i] + range * iqd, data.max)
        lower[i] <- max(q1[i] - range * iqd, data.min)
        est.xlim <- c(min(lower[i], data.min), max(upper[i], 
            data.max))
        smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
            args))
        hscale <- 0.4/max(smout$estimate) * wex
        base[[i]] <- smout$eval.points
        height[[i]] <- smout$estimate * hscale
        t <- range(base[[i]])
        baserange[1] <- min(baserange[1], t[1])
        baserange[2] <- max(baserange[2], t[2])
    }
    if (!add) {
        xlim <- if (n == 1) 
            at + c(-0.5, 0.5)
        else range(at) + min(diff(at))/2 * c(-1, 1)
        if (is.null(ylim)) {
            ylim <- baserange
        }
    }
    if (is.null(names)) {
        label <- 1:n
    }
    else {
        label <- names
    }
    boxwidth <- 0.05 * wex
    if (!add) 
        plot.new()
    if (!horizontal) {
        if (!add) {
            plot.window(xlim = xlim, ylim = ylim)
            axis(2)
            axis(1, at = at, label = label)
        }
        box()
        for (i in 1:n) {
            polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
                c(base[[i]], rev(base[[i]])), col = col[i], border = border, 
                lty = lty, lwd = lwd)
            if (drawRect) {
                lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
                  lty = lty)
                rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
                  q3[i], col = rectCol)
                points(at[i], med[i], pch = pchMed, col = colMed)
            }
        }
    }
    else {
        if (!add) {
            plot.window(xlim = ylim, ylim = xlim)
            axis(1)
            axis(2, at = at, label = label)
        }
        box()
        for (i in 1:n) {
            polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                rev(at[i] + height[[i]])), col = col[i], border = border, 
                lty = lty, lwd = lwd)
            if (drawRect) {
                lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
                  lty = lty)
                rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
                  boxwidth/2, col = rectCol)
                points(med[i], at[i], pch = pchMed, col = colMed)
            }
        }
    }
    invisible(list(upper = upper, lower = lower, median = med, 
        q1 = q1, q3 = q3))
}


# load kamphir logs
n100 <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n100.RLRootToTip.timetree.log', header=T, sep='\t')
temp <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n100.RLRootToTip.timetree.log.7', header=T, sep='\t')
n100 <- rbind(n100[100:nrow(n100),], temp[1000:nrow(temp),])
temp <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n100.RLRootToTip.timetree.restart.7.log.1', header=T, sep='\t')
n100 <- rbind(n100, temp)

n300 <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n300.RLRootToTip.timetree.log.2', header=T, sep='\t')
temp <- read.table('~/git/kamphir/projects/hivepi/data/SIRTree.n300.RLRootToTip.timetree.log.6', header=T, sep='\t')
n300 <- rbind(n300[100:nrow(n300),], temp[100:nrow(temp),])

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
b100$beta <- b100$R0Es * b100$becomeUninfectious / b100$S0Es
b100$phi <- b100$becomeUninfectious * b100$samplingProportion

b300 <- read.table('~/git/kamphir/projects/hivepi/beast/SIRTree.n300.log', header=T, sep='\t')
temp <- read.table('~/git/kamphir/projects/hivepi/beast/SIRTree.n300.1.log', header=T, sep='\t')
b300 <- rbind(b300[100:nrow(b300),], temp[100:nrow(temp),])
b300$beta <- b300$R0Es * b300$becomeUninfectious / b300$S0Es
b300$phi <- b300$becomeUninfectious * b300$samplingProportion

b1000 <- read.table('~/git/kamphir/projects/hivepi/beast/SIRTree.n1000.log', header=T, sep='\t')
temp <- read.table('~/git/kamphir/projects/hivepi/beast/SIRTree.n1000.1.log', header=T, sep='\t')
b1000 <- rbind(b1000[100:nrow(b1000),], temp[100:nrow(temp),])
b1000$beta <- b1000$R0Es * b1000$becomeUninfectious / b1000$S0Es
b1000$phi <- b1000$becomeUninfectious * b1000$samplingProportion

#boxplot(b100$S0, b300$S0, b1000$S0, log='y')


my.vioplot(log10(b1000$S0), log10(n1000$N), log10(b300$S0), log10(n300$N), log10(b100$S0), log10(n100$N),  col=rep(c(rgb(0.8,0.2,0,0.5), rgb(0,0,1,0.5)), times=3), horizontal=T, names=NA)
abline(h=2.5, col='grey40')
abline(h=4.5, col='grey40')
lines(c(3,3), c(4.5,7), lwd=8, col=rgb(0,0.5,0,0.6), lend=1)
lines(rep(log10(3000), 2), c(2.5,4.5), lwd=8, col=rgb(0,0.5,0,0.6), lend=1)
lines(c(4,4), c(0,2.5), lwd=8, col=rgb(0,0.5,0,0.6), lend=1)
text(x=6.5, y=6.1, label='A', cex=3)
text(x=6.5, y=4.1, label='B', cex=3)
text(x=6.5, y=2.1, label='C', cex=3)

#points(x=c(1,2,3), y=c(3, log10(3000), 4), cex=2, pch=4, lwd=1.5)



# generate table contents
df <- cbind(median(n100$N))


######Step 2: Measure consistency between methods
######################i.e. How consistent are the ranks in Genie3, Tigress, and CLR?
pdf('Figure1_consistency.pdf', width=10, height=4 ,pointsize = 10)
par(mfrow=c(1,3));
par(mar=c(4.2, 4.2, 2, 2))
par(cex=1);
plot(order[1:40000], pch=19, col=rgb(0,0,1,0.002), xlab="CLR rank", ylab="Genie3 rank")
plot(order_tig, pch=19, col=rgb(0,0,1,0.002), xlab="Tigress rank", ylab="Genie3 rank")
hist(avg_rank[which(avg_rank<20000)], n=40, xlab="avg rank", ylab="frequency", main="")

dev.off()
# step 1
library(compositions)
library(zCompositions)
library(ALDEx2)
library(CoDaSeq)
library(igraph)
library(car)

# step 2
d=read.table("data/count.table.txt", header=T, check.names=F, row.names=1, sep="\t", comment.char="", stringsAsFactors=FALSE, quote='')
e=read.table("data/subsys.txt", header=T,  sep="\t", comment.char="", stringsAsFactors=FALSE, quote='')

# step 3
d.g <- data.frame(d[,"010B"], d[,"001B"], d[,"009B"], d[,"012B"], d[,"013A"], d[,"014B"], d[,"017B"], d[,"018B"], check.names=FALSE, stringsAsFactors=FALSE)
colnames(d.g) <- c("010B", "001B","009B","012B","013A","014B","017B","018B")
rownames(d.g) <- rownames(d)

# step 4
d.f <- codaSeq.filter(d.g, min.count = 5, samples.by.row=FALSE)

# step 5
d.n0 <- cmultRepl(t(d.f), label=0, method='CZM')

# step 6
d.clr <- codaSeq.clr(d.n0, samples.by.row=TRUE)

# step 7
d.pcx <- prcomp(d.clr)

# step 8
grps=list(c(1:3), c(4:8))
#pdf("biplot.pdf", height=4, width=7)
par(mfrow=c(1,2))
codaSeq.PCAplot(d.pcx, plot.groups=TRUE, grp=grps, grp.col=c("red", "blue"), plot.circles=FALSE, plot.loadings=TRUE)
#dev.off()

# part 1.2
# step 1
conds <- c(rep("A", 3), rep("B", 5))

# step 2
x <- aldex.clr(d.f, conds)

# step 3
x.e <- aldex.effect(x, conds)

#step 4
x.t <- aldex.ttest(x, conds)

# step 5
x.all <- data.frame(x.e,x.t, stringsAsFactors=FALSE)

# step 6
par(mfrow=c(1))
pdf("effect.plot.pdf", height=5, width=5)
aldex.plot(x.all)
dev.off()
# part 1.3
# step 1

x.phi <- codaSeq.phi(x)

# step 2
phi.cutoff <- .01
x.lo.phi <- subset(x.phi, phi <= phi.cutoff)

# generate a graphical object
g <- graph.data.frame(x.lo.phi, directed=FALSE)

# overview of all the proportional relationships
# this can take a long time!!!
V(g)$label.cex <- 1
 plot(g, layout=layout.fruchterman.reingold.grid(g, weight=0.05/E(g)$phi), vertex.size=1, vertex.color="black")

 dg <- decompose.graph(g)

pdf("phi_graph.pdf", height=5, width=5)
plot(dg[[3]])
dev.off()

# # get the clusters from the graph object
 g.clust <- clusters(g)
#
# # data frame containing the names and group memberships of each cluster
 g.df.u <- data.frame(Systematic.name=V(g)$name, cluster=g.clust$membership, cluster.size=g.clust$csize[g.clust$membership], stringsAsFactors=FALSE)
g.df.u <- data.frame(Name=V(g)$name, cluster=g.clust$membership, stringsAsFactors=FALSE)

g.df <- g.df.u[order(g.df.u[,"cluster"]),]

# part 4 stripchart
pdf("stripchart.pdf", height=7, width=7)
par(mfrow=c(1,1))
codaSeq.stripchart(aldex.out=x.all, group.table=e, group.label="SEED1", heir=TRUE, heir.base="SEED4", mar=c(4,22,4,0.5), x.axis="effect", p.cutoff=0.0001, effect.cutoff=2)
dev.off()

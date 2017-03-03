## Example of a lattice plot
## 2017-03-03
## See also: http://stackoverflow.com/questions/24528527/how-to-add-a-title-to-legend-scale-using-levelplot-in-r
require(lattice)

p1 = seq(from=10, to=30, by=2)
p2 = seq(from=20, to=50, by=2)

g1 = expand.grid(p1=p1,p2=p2)
g1$pval = g1$p1 + g1$p2
g1$field = rep("field 1", nrow(g1))

g2 = expand.grid(p1=p1,p2=p2)
g2$pval = g2$p1 + 2*g2$p2
g2$field = rep("field 2", nrow(g2))

df = rbind(g1,g2)
levelplot(pval~p1*p2|field, data=df,
          cuts = 50, xlab="param 1", ylab="param 2",
          aspect="iso",
          main="P values", colorkey = TRUE, region = TRUE)

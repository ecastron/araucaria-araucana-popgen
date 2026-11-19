suppressPackageStartupMessages({library(phangorn); library(ape)})
m <- as.matrix(read.table('figdata/stand_fst_matrix.tsv', header=TRUE, sep='\t', check.names=FALSE))
rownames(m) <- colnames(m)
nn <- neighborNet(as.dist(m))
co <- phangorn:::coords(nn, dim='2D')          # vertex layout phangorn uses internally
cat('tips', length(nn$tip.label), '| vertices', nrow(co), '| edges', nrow(nn$edge), '\n')
write.table(data.frame(v=seq_len(nrow(co)), x=co[,1], y=co[,2]),
            'figdata/nnet_vertices.tsv', sep='\t', quote=FALSE, row.names=FALSE)
write.table(data.frame(v1=nn$edge[,1], v2=nn$edge[,2], len=nn$edge.length),
            'figdata/nnet_edges.tsv', sep='\t', quote=FALSE, row.names=FALSE)
write.table(data.frame(label=nn$tip.label, v=seq_along(nn$tip.label)),
            'figdata/nnet_tips.tsv', sep='\t', quote=FALSE, row.names=FALSE)
cat('written\n')

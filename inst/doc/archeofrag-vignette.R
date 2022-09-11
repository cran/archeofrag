## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.path = "man/figures/",
  comment = "#>"
)

## ----make-empirical-graph,  message=F-----------------------------------------
library(archeofrag)
data(LiangAbu)
abu.frag <- make_frag_object(cr=df.cr, fragments=fragments.info)
abu.g <- make_cr_graph(abu.frag)

## ----frag.relations.by.layers-------------------------------------------------
frag.relations.by.layers(abu.g, "layer")

## ----manipulate-plot-abu, fig.align="center"----------------------------------
par(mar=c(1, 0, 2, 0))
frag.graph.plot(abu.g, layer.attr="layer", main="All layers")

## ----manipulate-get.layer.pair------------------------------------------------
abu.g12 <- frag.get.layers.pair(abu.g, layer.attr="layer",
                                sel.layers=c("1", "2"))

## ----manipulate-plot-abu2, fig.align="center"---------------------------------
par(mar=c(1, 0, 2, 0))
frag.graph.plot(abu.g12, layer.attr="layer", main="Layers 1 and 2")

## ----manipulate-get.layer.pair2-----------------------------------------------
frag.get.layers.pair(abu.g, layer.attr="layer", sel.layers=c("1", "2"),
                     size.mini=2, mixed.components.only=TRUE)

## ----manipulate-frag.get.layers-----------------------------------------------
frag.get.layers(abu.g, layer.attr="layer", sel.layers="1")

## ----weights-g1---------------------------------------------------------------
abu.g12 <- frag.edges.weighting(abu.g12, layer.attr="layer")

## ----cohesion-g12-------------------------------------------------------------
frag.layers.cohesion(abu.g12, layer.attr="layer")

## ----weights-g12morpho--------------------------------------------------------
abu.g12morpho <- frag.edges.weighting(abu.g12,
                                      layer.attr="layer",
                                      morphometry="length")

## ----cohesion-g12morpho-------------------------------------------------------
frag.layers.cohesion(abu.g12morpho, layer.attr="layer")

## ----admixture----------------------------------------------------------------
# topology-based weighting:
frag.layers.admixture(abu.g12, layer.attr="layer")
# topology + morphometry weighting:
frag.layers.admixture(abu.g12morpho, layer.attr="layer")

## ----make-simulatd-graph, message=FALSE---------------------------------------
simul.g <- frag.simul.process(n.components=20, vertices=50)

## ----simulator-all-params, eval=F---------------------------------------------
#  frag.simul.process(initial.layers=1,
#                     n.components=20,
#                     vertices=50,
#                     edges=40,
#                     balance=.4,
#                     components.balance=.4,
#                     disturbance=.1,
#                     aggreg.factor=0,
#                     planar=FALSE,
#                     asymmetric.transport.from="1")

## ----frag.observer.failure, eval=FALSE----------------------------------------
#  frag.observer.failure(abu.g12, likelihood=0.2)

## ----params-------------------------------------------------------------------
params <- frag.get.parameters(abu.g12, layer.attr="layer")
params

## ----simulator-test-----------------------------------------------------------
# for H2:
test.2layers.g <- frag.simul.process(initial.layers=2,
                                    n.components=params$n.components,
                                    vertices=params$vertices,
                                    disturbance=params$disturbance,
                                    aggreg.factor=params$aggreg.factor,
                                    planar=params$planar)
# for H1:
test.1layer.g <- frag.simul.process(initial.layers=1,
                                   n.components=params$n.components,
                                   vertices=params$vertices,
                                   disturbance=params$disturbance,
                                   aggreg.factor=params$aggreg.factor,
                                   planar=params$planar)

## ----simulator-test2, message=FALSE-------------------------------------------
run.test2 <- function(x){
  frag.simul.process(initial.layers=2, # note the different value
                     n.components=params$n.components,
                     vertices=params$vertices,
                     disturbance=params$disturbance,
                     aggreg.factor=params$aggreg.factor,
                     planar=params$planar)
}

## ----simulator-test2-run, message=FALSE---------------------------------------
test2.results <- lapply(1:100, run.test2)

## ----manipulate-param-plot-mar, echo=FALSE------------------------------------
par(mar=c(5, 4, 4, 2))

## ----simulator-test2-edges, message=FALSE, fig.align="center", fig.width=4, fig.height=3----
edges.res <- sapply(test2.results,
                    function(g) frag.get.parameters(g, "layer")$edges)
plot(stats::density(edges.res), main="Edges")
abline(v=params$edges, col="red")

## ----simulator-test2-admix, message=FALSE, fig.align="center", fig.width=4, fig.height=3----
admix.res <- sapply(test2.results,
                    function(g) frag.layers.admixture(g, "layer"))
plot(stats::density(admix.res), main="Admixture")
abline(v=frag.layers.admixture(abu.g12, "layer"), col="red")

## ----simul-compare, message=FALSE---------------------------------------------
compare.res <- frag.simul.compare(abu.g12, layer.attr="layer",
                                  iter=30, summarise=FALSE)
head(compare.res$h1.data)

## ----simul-summarise, message=FALSE-------------------------------------------
frag.simul.summarise(abu.g12, layer.attr="layer",
                     compare.res$h1.data,
                     compare.res$h2.data)

## ----make-similarity----------------------------------------------------------
# make a frag object and generate a similarity graph:
abu.frag <- make_frag_object(sr=df.sr, fragments=fragments.info)
abu.sr <- make_sr_graph(abu.frag)

## ----count-similarity---------------------------------------------------------
# count of similarity relationships in and between layers:
simil.by.layers.df <- frag.relations.by.layers(abu.sr, "layer")
simil.by.layers.df

## ----similarity-perc-tab------------------------------------------------------
# percentage of similarity relationships in and between layers:
round(simil.by.layers.df / sum(simil.by.layers.df, na.rm=T) * 100, 0)

## ----similarity-dendr,  eval=T------------------------------------------------
# turn similarity into distance:
simil.dist <- max(c(simil.by.layers.df), na.rm=T) - simil.by.layers.df
simil.dist <- as.dist(simil.dist)
# hierarchical clustering:
clust.res <- stats::hclust(simil.dist, method="ward.D2")

## ----similarity-dendr-fig, fig.width = 3, fig.height = 3, fig.align="center", fig.cap="Hierarchical clustering of the pottery layers in Liang Abu (distance: based on the number of similarity relationships; clustering method: Ward).", eval=T, message=F----
clust.res$labels <- as.character(factor(clust.res$labels, 
                                        levels=c("0", "1", "2"),
                                        labels=c("layer 0", "layer 1", "layer 2")))
plot(clust.res, hang=-1, axes=F, ann=F)

## ----subgraphs-simul----------------------------------------------------------
# simulate a fragmentation graph:
simul.g <- frag.simul.process(initial.layers=2,
                              n.components=20,
                              vertices=70,
                              balance=.45)
# extract the subgraph of each spatial unit:
simul.g1 <- frag.get.layers(simul.g, layer.attr="layer", sel.layers="1")[[1]]
simul.g2 <- frag.get.layers(simul.g, layer.attr="layer", sel.layers="2")[[1]]

## ----cycles-simul1------------------------------------------------------------
rbind(
  "unit1" = frag.cycles(simul.g1, kmax=5),
  "unit2" = frag.cycles(simul.g2, kmax=5))

## ----pathlengths--------------------------------------------------------------
frag.path.lengths(simul.g1)
frag.path.lengths(simul.g2)
frag.path.lengths(simul.g2, cumulative=T)

## ----diameters----------------------------------------------------------------
frag.diameters(simul.g1)
frag.diameters(simul.g2)


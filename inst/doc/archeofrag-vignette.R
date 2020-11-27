## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.path = "man/figures/",
  comment = "#>"
)

## ----setup, message=F---------------------------------------------------------
library(archeofrag)

## ----make-empirical-graph-----------------------------------------------------
data(LiangAbu)
abu.g <- make_frag_object(df.cr, fragments = fragments.info)
abu.g <- make_cr_graph(abu.g)

## ----make-simulatd-graph, message=FALSE---------------------------------------
simul.g <-frag.simul.process(n.components = 20, vertices = 50)

## ----manipulate-param-plot, echo=FALSE----------------------------------------
pardefault <- par(no.readonly = TRUE)
par(mar=rep(0, 4))

## ----manipulate:plot-abu, fig.align="center"----------------------------------
frag.graph.plot(abu.g, "layer")

## ----manipulate:plot-simul, fig.align="center"--------------------------------
frag.graph.plot(simul.g, "layer")

## ----manipulate:get.layer.pair------------------------------------------------
abu.g12 <- frag.get.layers.pair(abu.g, "layer", c("1", "2"))

## ----manipulate:plot-abu2, fig.align="center"---------------------------------
frag.graph.plot(abu.g12, "layer")

## ----manipulate:frag.get.layers-----------------------------------------------
frag.get.layers(simul.g, "layer", sel.layers = "1")


## ----relations-by-layer-------------------------------------------------------
frag.relations.by.layers(simul.g, "layer")

## ----relations-by-layer2------------------------------------------------------
simul2.g <-frag.simul.process(n.components=20, vertices=50, disturbance=.1)
frag.relations.by.layers(simul2.g, "layer")

## ----get-parameters-----------------------------------------------------------
frag.get.parameters(abu.g12, "layer")

## ----get-parameters-simul-----------------------------------------------------
cbind(
  frag.get.parameters(simul.g, "layer"),
  frag.get.parameters(simul2.g, "layer")
)

## ----weights------------------------------------------------------------------
E(simul.g)$weight
simul.g <- frag.edges.weighting(simul.g, "layer")
E(simul.g)$weight

## ----cohesion-simul1----------------------------------------------------------
frag.layers.cohesion(simul.g, "layer")

## ----cohesion-simul2----------------------------------------------------------
simul2.g <- frag.edges.weighting(simul2.g, "layer")
frag.layers.cohesion(simul2.g, "layer")

## ----admixture----------------------------------------------------------------
frag.layers.admixture(simul.g, "layer")
frag.layers.admixture(simul2.g, "layer")

## ----subgraphs-simul----------------------------------------------------------
simul.l1.g <- induced_subgraph(simul2.g, V(simul2.g)[V(simul2.g)$layer==1])
simul.l2.g <- induced_subgraph(simul2.g, V(simul2.g)[V(simul2.g)$layer==2])

## ----cycles-simul1------------------------------------------------------------
frag.cycles(simul.l1.g, kmax=5)

## ----cycles-simu21------------------------------------------------------------
frag.cycles(simul.l2.g, kmax=5)

## ----pathlengths--------------------------------------------------------------
frag.path.lengths(simul.l1.g)
frag.path.lengths(simul.l2.g, cumulative=T)

## ----diameters----------------------------------------------------------------
frag.diameters(simul.l1.g)
frag.diameters(simul.l2.g)

## ----simulator-all-params-----------------------------------------------------
frag.simul.process(initial.layers = 1,
                   n.components = 20,
                   vertices = 50, edges = 40,
                   disturbance = .1,
                   balance = .4,
                   components.balance = .4,
                   aggreg.factor = 0,
                   planar = T)

## ----params-------------------------------------------------------------------
params <- frag.get.parameters(abu.g12, "layer")

## ----simulator-test-----------------------------------------------------------
test.2layers.g <-frag.simul.process(initial.layers = 2,
                              n.components = params$n.components,
                              vertices = params$vertices,
                              disturbance = params$disturbance,
                              aggreg.factor = params$aggreg.factor,
                              planar = params$planar)

test.1layer.g <-frag.simul.process(initial.layers = 2,
                              n.components = params$n.components,
                              vertices = params$vertices,
                              disturbance = params$disturbance,
                              aggreg.factor = params$aggreg.factor,
                              planar = params$planar)

## ----simulator-test2, message=FALSE-------------------------------------------
run.test2 <- function(x){
  frag.simul.process(initial.layers = 2,
                              n.components = params$n.components,
                              vertices = params$vertices,
                              disturbance = params$disturbance,
                              aggreg.factor = params$aggreg.factor,
                              planar = params$planar)
}

## ----simulator-test2-run, message=FALSE---------------------------------------
test2.results <- lapply(1:40, run.test2)

## ----manipulate:param-plot-mar, echo=FALSE------------------------------------
par(mar=c(5, 4, 4, 2))

## ----simulator-test2-examine, message=FALSE, fig.align="center", fig.width=5----
edges.res <- sapply(test2.results,
                    function(g) frag.get.parameters(g, "layer")$edges)
plot(density(edges.res), main="Edges")
abline(v=params$edges, col="red")

## ----restore-param, echo=FALSE, warning=FALSE, message=FALSE------------------
par(pardefault)


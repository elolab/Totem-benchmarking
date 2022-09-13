set.seed(1234)

suppressMessages(library(Totem))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(cowplot))
suppressMessages(library(dyndimred))
suppressMessages(library(S4Vectors))

dynwrap_object <- readRDS("~/linear_tree/thymus-t-cell-differentiation_mca.rds")

sce <- SingleCellExperiment(assays = list(counts = t(dynwrap_object$counts),
                                          logcounts = t(dynwrap_object$expression)))
sce <- PrepareTotem(sce)

print(dim(sce))

sce <- RunDimRed(object = sce,
                 dim.red.method = "lmds",
                 dim.red.features = NULL,
                 dim.reduction.parameter.list = list(ndim=5))

# extract the generated embedding and set it again into the same slot
own_dim_red <- reducedDim(sce)
reducedDim(sce,type = "lmds") <- own_dim_red # use any name

dim_red <- dimred_mds(dynwrap_object$expression,ndim=2)

sce <- RunClustering(sce,
                         k.range = 3:20,
                         min.cluster.size = 5,
                         N.clusterings=10000)



p1 <- VizCellConnectivity(sce,custom.dim.red = dim_red)

sce <- SelectClusterings(sce,selection.method = 1,
                       selection.N.models = 10,
                       selection.stratified=FALSE,
                       prior.clustering = NULL)


ReturnTrajNames(sce)

VizMST(sce,clustering.names = ReturnTrajNames(sce),custom.dim.red = dim_red)

sce <- RunSmoothing(sce)


ReturnTrajNames(sce)

VizSmoothedTraj(sce,
                traj.names = ReturnTrajNames(sce),
                custom.dim.red = dim_red,plot.pseudotime = FALSE)


ReturnSmoothedTrajNetwork(sce,clustering.name = "4.165")

sce <- ChangeTrajRoot(sce,traj.name="4.165",root.cluster=4)


p2<-VizSmoothedTraj(sce,
                traj.names = "4.165",
                custom.dim.red = dim_red,plot.pseudotime = FALSE)



true_milestones <- unlist(lapply(split(dynwrap_object$milestone_percentages,f = dynwrap_object$milestone_percentages$cell_id),function(x) x$milestone_id[which.max(x$percentage)]))
true_milestones <- true_milestones[dynwrap_object$cell_ids]


# This uses the Slingshot distance method to create the MST
p3 <- VizClustering(sce,clustering = true_milestones,custom.dim.red = dim_red)

# Plot the real trajectory
p4 <- dynplot::plot_dimred(dynwrap_object,dimred = dim_red)

cowplot::plot_grid(p4,p3,p1,p2)




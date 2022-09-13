args <- commandArgs(trailingOnly=TRUE)

cat(args, sep = "\n")

benchmarking_file <- args[1]
output_folder <- args[2]
metric_folder <- args[3]

print(benchmarking_file)
print(output_folder)
print(metric_folder)



file_name <- unlist(strsplit(benchmarking_file,"/"))[c(length(unlist(strsplit(benchmarking_file,"/")))-1,
                                                       length(unlist(strsplit(benchmarking_file,"/"))))]

output_file <- paste0(output_folder,"/",paste(file_name,collapse = "_"))
sce_file <- paste0(metric_folder,"/",paste(file_name,collapse = "_"))
sce_file_ <- gsub(".rds","",sce_file)
sce_file <- paste0(sce_file_,"_sce.rds")

print(sce_file)

if (file.exists(output_file))
{
  print(STOP)
}






library(dynwrap)
library(dyneval)
library(tidyverse)

definition <- definition(
  method = def_method(
    id = "Totem"
  ),
  # leave incomplete. no need to fill.
  parameters = def_parameters(
    dynparam::integer_parameter(
      id = "component",
      default = 1,
      distribution = dynparam::uniform_distribution(1, 10),
      description = "The nth component to use"
    )
  ),
  wrapper = def_wrapper(
    input_required = "expression",
    input_optional = "start_id"
  )
)




run_fun <- function(expression, priors, parameters, seed, verbose) {
  
  set.seed(123456)
  

  library(SingleCellExperiment)
  library(Totem)
  
  sce <- SingleCellExperiment(assays = list(counts=t(expression),
                                            logcounts = t(expression)))
  
  # There's no need to run the algorithm here. skip.
  return(sce)
  
  
}


ti_Totem <- create_ti_method_r(definition, run_fun, package_loaded = c("Totem"))


dynwrap_object <- readRDS(benchmarking_file)

true_milestones <- unlist(lapply(split(dynwrap_object$milestone_percentages,f = dynwrap_object$milestone_percentages$cell_id),function(x) x$milestone_id[which.max(x$percentage)]))
true_milestones <- true_milestones[dynwrap_object$cell_ids]

dataset <- wrap_expression(
  counts = dynwrap_object$counts,
  expression = dynwrap_object$expression
)


library(SingleCellExperiment)
library(Totem)



modelsss <- infer_trajectories(dataset,list(ti_Totem()
),verbose = TRUE)

sce <- readRDS(sce_file)

metrics_list <- list()

for (i in 1:100)
{
  message(i)
  for (selection_method in c(1,2,3,4,5))
  {
    for (stratifiedd in c(TRUE,FALSE))
    {
      sce <- SelectTopModels(sce,selection.method = selection_method,
                             selection.N.models = i,
                             selection.stratified = stratifiedd,
                             prior.clustering = true_milestones)
      
      k_name <- names(sce@metadata$totem$selected.clustering)[i]
      base_clusterings <- sce@metadata$totem$selected.clustering
      sce@metadata$totem$selected.clustering <- list()
      sce@metadata$totem$selected.clustering[[k_name]] <- base_clusterings[[k_name]]
      
      sce <- RunSmoothing(sce)
      
      dynwrap_model <- sce@metadata$totem$dynwrap_trajectory[[1]]
      
      modelsss$model[[1]] <- dynwrap_model
      
      modelsss$model <- map(modelsss$model, add_cell_waypoints)

      metric_ids <- c("correlation","F1_branches","him","featureimp_wcor")
      metrics <- map_dfr(modelsss$model, dyneval::calculate_metrics, dataset = dynwrap_object, metrics = metric_ids)
      metrics$k <- k_name
      metrics$method <- paste0("Totem_",selection_method,"_",stratifiedd)
      
      metrics_list[[paste0("Totem_",selection_method,"_",stratifiedd,"_",k_name)]] <- metrics

      
    }
  }
  
}

res_totem <- metrics_list



metrics<-sce@metadata$totem$all.clustering.scores


i <-  1

ks_to_test <- head(metrics$clustering_name[order(metrics$calinski_harabasz,decreasing = TRUE)],100)
metrics_list <- list()
for (k in ks_to_test)
{
  message(i)
  sce@metadata$totem$selected.clustering <- list()
  sce@metadata$totem$selected.clustering[[k]] <- sce@metadata$totem$all.clustering[[k]]

  sce <- RunSmoothing(sce)

  dynwrap_model <- sce@metadata$totem$dynwrap_trajectory[[1]]

  modelsss$model[[1]] <- dynwrap_model

  modelsss$model <- map(modelsss$model, add_cell_waypoints)
  #

  metric_ids <- c("correlation","F1_branches","him","featureimp_wcor")
  metrics <- map_dfr(modelsss$model, dyneval::calculate_metrics, dataset = dynwrap_object, metrics = metric_ids)
  metrics$k <- k
  metrics$method <- "calinhara"

  metrics_list[[k]] <- metrics
  i <- i +1
}


res_calinhara <- metrics_list









dx <- dist(reducedDim(sce))
asws <- lapply(sce@metadata$totem$all.clustering,function(x) {summary(cluster::silhouette(x,dx))$avg.width})


i <-  1

ks_to_test <- head(names(unlist(asws))[order(unlist(asws),decreasing = TRUE)],100)
metrics_list <- list()
for (k in ks_to_test)
{
  message(i)
  sce@metadata$totem$selected.clustering <- list()
  sce@metadata$totem$selected.clustering[[k]] <- sce@metadata$totem$all.clustering[[k]]
  
  sce <- RunSmoothing(sce)
  
  dynwrap_model <- sce@metadata$totem$dynwrap_trajectory[[1]]
  
  modelsss$model[[1]] <- dynwrap_model
  
  modelsss$model <- map(modelsss$model, add_cell_waypoints)

  metric_ids <- c("correlation","F1_branches","him","featureimp_wcor")
  metrics <- map_dfr(modelsss$model, dyneval::calculate_metrics, dataset = dynwrap_object, metrics = metric_ids)
  metrics$k <- k
  metrics$method <- "asw"
  
  metrics_list[[k]] <- metrics
  i <- i +1
}

res_aws <- metrics_list




i <-  1

ks_to_test <- sample(names(sce@metadata$totem$all.clustering),100,replace = FALSE)
metrics_list <- list()
for (k in ks_to_test)
{
  message(i)
  sce@metadata$totem$selected.clustering <- list()
  sce@metadata$totem$selected.clustering[[k]] <- sce@metadata$totem$all.clustering[[k]]
  
  sce <- RunSmoothing(sce)
  
  dynwrap_model <- sce@metadata$totem$dynwrap_trajectory[[1]]
  
  modelsss$model[[1]] <- dynwrap_model
  
  modelsss$model <- map(modelsss$model, add_cell_waypoints)

  metric_ids <- c("correlation","F1_branches","him","featureimp_wcor")
  metrics <- map_dfr(modelsss$model, dyneval::calculate_metrics, dataset = dynwrap_object, metrics = metric_ids)
  metrics$k <- k
  metrics$method <- "random"
  
  metrics_list[[k]] <- metrics
  i <- i +1
}


res_random <- metrics_list


res <- do.call(rbind,res_totem)
res$metric <- "Totem"

res_ <- do.call(rbind,res_calinhara)
res_$metric <- "calinhara"
res <- rbind(res,res_)

res_ <- do.call(rbind,res_random)
res_$metric <- "random"
res <- rbind(res,res_)

res_ <- do.call(rbind,res_aws)
res_$metric <- "asw"
res <- rbind(res,res_)

res <- as.data.frame(res)
res <- reshape2::melt(res)
res <- res[!grepl("time",res$variable),]


saveRDS(res,file = output_file)


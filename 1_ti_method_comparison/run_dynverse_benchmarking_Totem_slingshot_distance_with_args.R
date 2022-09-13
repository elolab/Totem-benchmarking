args <- commandArgs(trailingOnly=TRUE)

cat(args, sep = "\n")

benchmarking_file <- args[1]
output_folder <- args[2]

print(benchmarking_file)
print(output_folder)



file_name <- unlist(strsplit(benchmarking_file,"/"))[c(length(unlist(strsplit(benchmarking_file,"/")))-1,
                                                       length(unlist(strsplit(benchmarking_file,"/"))))]

output_file <- paste0(output_folder,"/",paste(file_name,collapse = "_"))

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
        dynparam::character_parameter(
            id = "algorithm",
            values = c("clara"),
            default = "clara",
            description = "clustering algorithm"
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


dataset <- readRDS(benchmarking_file)

true_milestones <- unlist(lapply(split(dataset$milestone_percentages,f = dataset$milestone_percentages$cell_id),function(x) x$milestone_id[which.max(x$percentage)]))
true_milestones <- true_milestones[dataset$cell_ids]


data <- wrap_expression(
    counts = dataset$counts,
    expression = dataset$expression
)

modelsss <- infer_trajectories(data,list(ti_Totem()
),verbose = TRUE)

library(SingleCellExperiment)
library(Totem)

sce <- SingleCellExperiment(assays = list(counts=t(dataset$expression),
                                          logcounts = t(dataset$expression)))

sce <- RunTotem(sce,dim.red.method = "lmds",dim.reduction.parameter.list = list(ndim=5),k.range = 3:20,
                min.cluster.size = 5,selection.method = 1,
                selection.N.models = 1,selection.stratified = FALSE,N.clusterings = 10000)

output_file_ <- gsub(".rds","",output_file)
output_file_ <- paste0(output_file_,"_sce",".rds")
saveRDS(sce,file = output_file_)


modelsss$model[[1]] <- sce

# run different selection methods. We use type==1 in the actual evaluation.

for (type1 in c(1,2,3,4,5))
{
    for (type2 in c("single"))
    {
        modelss <- modelsss
        sce <- modelss$model[[1]]
        
        if (type1 >= 4 && type2=="single") {
          sce <- SelectClusterings(sce,selection.method = type1,prior.clustering = true_milestones)
        } else if (type1 < 4 && type2=="single") {
          sce <- SelectClusterings(sce,selection.method = type1)
        }
        
        sce <- RunSmoothing(sce)
        
        if (type2=="single")
        {
            modelss$model[[1]] <- list()
            modelss$model[[1]][[1]] <- sce@metadata$totem$dynwrap_trajectory[[1]]
            names(modelss$model[[1]]) <- names(sce@metadata$totem$dynwrap_trajectory)
            
        }
        
        
        metrics_list <- list()
        overalls <- c()
        print(names(modelss))
        print(names(modelss$model[[1]]))
        

        for (model_name in names(modelss$model[[1]]))
        {
            models <- modelss
            models$model[[1]] <- modelss$model[[1]][[model_name]] 
            
            dataset <- add_cell_waypoints(dataset)
            models$model <- map(models$model, add_cell_waypoints)
            
            message("2. add waypoints... DONE")
            
            metric_ids <- dyneval::metrics %>% pull(metric_id)
            metrics <- map_dfr(models$model, dyneval::calculate_metrics, dataset = dataset, metrics = metric_ids)
            
            
            
            message("3. calculate metrics... DONE")
            
            metrics$method <- c("Totem")
            
            metrics$overall <- (metrics$correlation)*(metrics$featureimp_wcor)*(metrics$F1_branches)*(metrics$him)
            metrics$overall <- (metrics$overall)^(1/4)
            
            metrics_list[[model_name]] <- metrics
            overalls <- c(overalls,metrics$overall)
        }
        
        metrics <- metrics_list[[which.max(overalls)]]
        
        dir.create(output_folder,recursive=TRUE)
        
        output_file_ <- gsub(".rds","",output_file)
        output_file_ <- paste0(output_file_,"_type_",type1,"_",type2,".rds")
        
        saveRDS(metrics,file = output_file_)
        
        
        
    }
}

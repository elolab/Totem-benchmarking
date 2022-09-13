args <- commandArgs(trailingOnly=TRUE)

cat(args, sep = "\n")

benchmarking_file <- args[1]
output_folder <- args[2]

print(benchmarking_file)
print(output_folder)

library(dynwrap)
library(dyneval)
library(tidyverse)

definition <- definition(
    method = def_method(
        id = "Tinga"
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
    
    
    library(tidyverse)
    library(dyno)
    library(dynwrap)
    
    data <- wrap_expression(
        expression = expression,
        counts = expression
    )
    
    traj <- infer_trajectory(
        data, 
        method = gng_param2(), 
        seed = 123456
    )
    
    
    
    dynwrap::wrap_data(cell_ids = rownames(expression)) %>%
        dynwrap::add_trajectory(milestone_ids = traj$milestone_ids,
                                milestone_network = traj$milestone_network,
                                divergence_regions = traj$divergence_regions,
                                milestone_percentages = traj$milestone_percentages,
                                progressions = traj$progressions)
}



ti_tinga <- create_ti_method_r(definition, run_fun, package_loaded = c("TInGa"))

dataset <- readRDS(benchmarking_file)

models <- infer_trajectories(dataset,list(ti_tinga()
),verbose = TRUE)

message("1. trajectory inference... DONE")


dataset <- add_cell_waypoints(dataset)
models$model <- map(models$model, add_cell_waypoints)

message("2. add waypoints... DONE")

metric_ids <- dyneval::metrics %>% pull(metric_id)
metrics <- map_dfr(models$model, dyneval::calculate_metrics, dataset = dataset, metrics = metric_ids)

message("3. calculate metrics... DONE")
metrics$method <- c("scShaper")


file_name <- unlist(strsplit(benchmarking_file,"/"))[c(length(unlist(strsplit(benchmarking_file,"/")))-1,
                                          length(unlist(strsplit(benchmarking_file,"/"))))]

output_file <- paste0(output_folder,"/",paste(file_name,collapse = "_"))

dir.create(output_folder,recursive=TRUE)

saveRDS(metrics,file = output_file)











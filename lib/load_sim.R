load_sim <- function(folder, setting, method) {
  # folder <- "Simulation BARTBART"
  # method <- "BART"
  # setting <- "BART"
  
  files <- paste0(folder, "/", list.files(folder))
  
  indv_files <- files %>% str_subset("/indv_seed")
  avg_files <- files %>% str_subset("/avg_seed")
  tree_files <- files %>% str_subset("/tree_subgroup")
  subgroup_files <- files %>% str_subset("/subgroup")
  
  indv_results <- map_df(indv_files, readRDS) %>%
    mutate(Method = method, Setting = setting)
  avg_results <- map_df(avg_files, readRDS) %>%
    mutate(Method = method, Setting = setting)
  tree_results <- map_df(tree_files, readRDS) %>%
    mutate(Method = method, Setting = setting)
  subgroup_results <- map_df(subgroup_files, readRDS) %>%
    mutate(Method = method, Setting = setting)
  
  return(
    list(
      indv_results     = indv_results,
      avg_results      = avg_results,
      tree_results     = tree_results,
      subgroup_results = subgroup_results
    )
  )
  
}
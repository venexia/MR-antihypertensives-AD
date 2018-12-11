 # Required libraries: TwoSampleMR, tidyverse, furrr
# For group MR a backend must be specified from the future package. 
# For sequential use plan(sequential) and for multicore plan(multiprocess).

group_MR <- function (df, outcome.data = NULL, outcome = NULL, ...) {
  
  group_var <- quos(...)

  # Add unique ID to dataframe for grouping variables
  
  group_df <- df %>%
    as_tibble() %>%
    mutate(index = group_indices(df, !!!group_var))
  
  # Perform MR function for each group

  if (outcome.data=="mrbase") {
  
  output <- group_df %>% 
    group_by(index) %>%
    nest() %>% 
    mutate(output = future_map(data, ~ suppressMessages(MR(., outcome = outcome, wait = TRUE)), .progress = TRUE)) %>% 
    dplyr::select(index, output) %>%
    unnest(output) 
  
  }
  
  if (outcome.data=="file") {
    
    output <- group_df %>% 
      group_by(index) %>%
      nest() %>% 
      mutate(output = future_map(data, ~ suppressMessages(MR_file(., outcome = outcome, wait = TRUE)), .progress = TRUE)) %>% 
      dplyr::select(index, output) %>%
      unnest(output) 
    
  }
  
  output <- group_df %>% 
    dplyr::select(index = index, !!!group_var) %>%
    distinct() %>% 
    full_join(output, by = "index") %>%
    dplyr::select(index, !!!group_var, everything())
     
  return(output)
  
}


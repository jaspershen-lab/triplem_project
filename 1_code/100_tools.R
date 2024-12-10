lm_adjust <-
  function(expression_data,
           sample_info,
           threads = 5) {
    library(future)
    library(furrr)
    # plan(strategy = multisession(workers = threads))
    new_expression_data <-
      rownames(expression_data) %>%
      purrr::map(function(name) {
        # cat(name, " ")
        x = as.numeric(expression_data[name,])
        temp_data =
          data.frame(x = x,
                     sample_info)
        
        
        temp_data$Gender[temp_data$Gender == 'Female'] = 0
        temp_data$Gender[temp_data$Gender == 'Male'] = 1
        temp_data$Gender = as.numeric(temp_data$Gender)
        
        temp_data$IRIS[temp_data$IRIS == 'IR'] = 1
        temp_data$IRIS[temp_data$IRIS == 'IS'] = 2
        temp_data$IRIS[temp_data$IRIS == 'Unknown'] = 0
        temp_data$IRIS = as.numeric(temp_data$IRIS)
        
        adjusted_x <-
          lm(x ~ Gender + BMI + IRIS, data = temp_data) %>%
          residuals()
        adjusted_x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(new_expression_data) <-
      colnames(expression_data)
    
    rownames(new_expression_data) <-
      rownames(expression_data)
    new_expression_data
  }





body_site_color = c(
  "gut" = "#edd064",
  "skin" = "#f2ccac",
  "oral" = "#a1d5b9",
  "nasal" = "#a17db4"
)
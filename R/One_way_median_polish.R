###########################################################################
# One way median polish
##########################################################################
#' Dependencies from other functions
# source("aux_FMP.R")
#' Parameters for the  One_way_median_polish
#' @param  n.  It represents the total number of functional curves.
#' @param  p.  It represents the grid size.
#' @param  Y.  A  matrix with dimension n by p. The functional data.
#' @param  n_prefectures. The number of prefectures, states or departments.
#' @param  year. Vector with the years considered in the population.
#' @param  age. Vector with the ages considered in each year.
#'

One_way_median_polish <- function(Y, n_prefectures=51, year=1959:2020, age=0:100)
{
    #Number of years considered in the population
    n_year = length(year)
    #Number of ages considered in each year
    n_age = length(age)

    #######
    #Get the indexes for the row partitions by prefectures
    #######
    part_list = list()
    for(ik in 1:n_prefectures) {
      part_list[[ik]] = (n_year*ik-(n_year-1)):(n_year*ik)
    }
    row_partition_index=part_list


    if(length(class(Y)) == 2 && class(Y)[1] == "matrix" && class(Y)[2] == "array") 
    {
      sample_size <- n_year
      grid_size <- ncol(Y)
    } 
    else if(all(!inherits(Y, "matrix"), !inherits(Y, "array")))
    {
      stop("Y should be functional data")
    }

    Functional_grand_effect <- rep(0, length.out = grid_size)
    Functional_row_effect <- matrix(0, nrow = length(row_partition_index), ncol = grid_size)
    new_Functional_row_effect <- Functional_row_effect + 1
    Y_iter<-Y
    Iteration <- 1
    while (sum(new_Functional_row_effect) != 0) {
      cat("Iteration=", Iteration,"\n")
      cat("Sum =", sum(new_Functional_row_effect),"\n")

      partition_Y <- lapply(1:length(row_partition_index), function(k) {
        Y_iter[row_partition_index[[k]], ]})

      #Compute the median in each row with the curves in the row. Depth based selection

      Median_at_each_row<-lapply(partition_Y,get_median)

      #Get the functional median of the row functional medians
      Medians_medians_row<-get_median(list.rbind( Median_at_each_row))
      Functional_grand_effect=Functional_grand_effect+Medians_medians_row

      #Get the functional row effect
      new_Functional_row_effect<-t(sapply(Median_at_each_row, function(k) {k -Medians_medians_row}))
      # sum(new_Functional_row_effect)
      Functional_row_effect=Functional_row_effect+new_Functional_row_effect

      # Substract the row-median in each row from the rest of the curves in the row

      partition_Y<-lapply(partition_Y,Substraction)
      Y_iter<-list.rbind(partition_Y)
      Iteration <- Iteration + 1
    }
  return(list(grand_effect = Functional_grand_effect,
              row_effect = Functional_row_effect))
}

#' Return from One_way_median_polish
#'
#' @return Return a list containing grand effect, and a row effect.
#' @return grand_effect, a vector of dimension p.
#' @return row_effect, a matrix of dimension length(row_partition_index) by p.

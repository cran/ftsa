
###########################################################################
# Two way median polish 
##########################################################################
#' Dependencies from other functions
# source("aux_FMP.R")
#' Parameters for the  Two_way_median_polish
#' @param  n.  It represents the total number of functional curves.
#' @param  p.  It represents the grid size.
#' @param  Y.  A  matrix with dimension n by 2p. The functional data.
#' @param  n_prefectures. The number of prefectures, states or departments. 
#' @param  year. Vector with the years considered in each population. 
#' @param  age. Vector with the ages considered in each year. 
#' @param  n_populations. Number of populations. 
Two_way_median_polish <- function(Y, year=1959:2020, age=0:100, n_prefectures=51, n_populations=2)
{
  #Number of years considered in each population
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
  #######
  #Get the indexes for the column partitions by populations.
  #######
  part_list_c = list()
  for(ik in 1:n_populations) {
    part_list_c[[ik]] = (n_age*ik-(n_age-1)):(n_age*ik)
  }
  column_partition_index=part_list_c
  
  if (length(class(Y)) == 2 && class(Y)[1] == "matrix" && class(Y)[2] == "array")
  {
      sample_size <- 2*nrow(Y)
      grid_size <- ncol(Y)/2
  }
  else if(!is(Y, "matrix") && !is(Y, "array"))
  {
      stop("Y should be functional data")
  } 
  
  if (length(row_partition_index) < 2|length(column_partition_index) < 2)
  {
      stop("partition_index is not a qualified partition")
  }  
  Functional_grand_effect <- rep(0, length.out = grid_size)
  Functional_row_effect <- matrix(0, nrow = length(row_partition_index), ncol = grid_size)
  Functional_column_effect <- matrix(0, nrow = length(column_partition_index), ncol = grid_size)
  new_Functional_row_effect <- Functional_row_effect + 1
  Y_iter<-Y
  Iteration <- 1
  while (sum(new_Functional_row_effect) != 0) {
    cat("Iteration=", Iteration,"\n")
    cat("Sum =", sum(new_Functional_row_effect),"\n")
    
    # Partition of the data according to the row partition index
    partition_row_Y<-partition_row_Y_f(Y_iter,column_partition_index,row_partition_index)
    
    #Compute the median in each row with the curves in the row. Depth based selection
    Median_at_each_row<-lapply(partition_row_Y,get_median)
    
    
    #Get the functional median of the row functional medians
    Medians_medians_row<-get_median(list.rbind( Median_at_each_row))
    Functional_grand_effect<-Functional_grand_effect+Medians_medians_row
    
    #Get the functional row effect
    new_Functional_row_effect<-t(sapply(Median_at_each_row, function(k) {k -Medians_medians_row}))
    #sum(new_Functional_row_effect)
    Functional_row_effect<-Functional_row_effect+new_Functional_row_effect
    
    # Subtract the row-median in each row from the rest of the curves in the row
    partition_row_Y<-lapply(partition_row_Y,Substraction)
    
    #Match again Y_iter
    Y_iter<-f(partition_row_Y,n_year)
    
    # Partition of the data according to the column partition index
    partition_column_Y<-partition_column_Y_f(Y_iter,n_age)
    #get the median at each column
    Median_at_each_column<-lapply(partition_column_Y,get_median)
    
    #Get the functional median of the row functional medians
    Medians_medians_col<-get_median(list.rbind( Median_at_each_column))
    
    #Add the FGE
    Functional_grand_effect=Functional_grand_effect+Medians_medians_col
    
    #Get the functional column effect
    new_Functional_column_effect<-t(sapply(Median_at_each_column, function(k) {k -Medians_medians_col}))
    Functional_column_effect<- Functional_column_effect + new_Functional_column_effect
    
    #data
    partition_column_Y<-lapply(partition_column_Y,Substraction)
    
    Y_iter<-list.cbind(partition_column_Y)
    sum(new_Functional_row_effect)
    Iteration <- Iteration + 1    
  }
  return(list(grand_effect = Functional_grand_effect,
              row_effect = Functional_row_effect,
              col_effect = Functional_column_effect))
}

#' Return from Two_way_median_polish
#' 
#' @return Return a list containing grand effect,row effect and column effect.
#' @return grand_effect, a vector of dimension p
#' @return row_effect, a matrix of dimension length(row_partition_index) by p.
#' @return col_effect, a matrix of dimension length(column_partition_index) by p


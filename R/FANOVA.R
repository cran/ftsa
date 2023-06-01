###############################################################################
# Two-way Functional Anova
###############################################################################
#' Description for the function FANOVA
#' This is a general function for the Functional ANOVA based on means. It provides
#' the functional grand effect, the functional row and column effects.
#' Dependencies from other functions

#' Parameters for the function FANOVA
#' @param  n.  It represents the total number of functional curves.
#' @param  p.  It represents the grid size.
#' @param  data_pop1. It's a p by n matrix
#' @param  data_pop2. It's a p by n matrix
#' @param  year. Vector with the years considered in each population. 
#' @param  age. Vector with the ages considered in each year. 
#' @param  n_prefectures. The number of prefectures, states or departments. 
#' @param  n_populations. Number of populations. 

FANOVA<-function(data_pop1,data_pop2,year=1959:2020,age= 0:100,n_prefectures=51,n_populations=2)
{
  #Number of years considered in each population
  n_year = length(year)
  #Number of ages considered in each year
  n_age = length(age)
 
  #######################################################
  #Get the indexes for the row partitions by prefectures
  #######################################################
  
  part_list = list()
  for(ik in 1:n_prefectures)
  {
      part_list[[ik]] = (n_year*ik-(n_year-1)):(n_year*ik)
      rm(ik)
  }
  row_par = part_list
  
  #functional grand effect
  FGE_means <- mu_hat(data_pop1, data_pop2, n_year, n_prefectures, n_age, n_populations, row_par)
  
  #functional row effect
  FRE_mean <- FRE_means(data_pop1, data_pop2, n_year, n_prefectures, n_age, n_populations, row_par)
  
  #Functional column effect
  FCE_mean <- FCE_means(data_pop1, data_pop2, n_year, n_prefectures, n_age, n_populations, row_par)
  return(list(FGE_mean=FGE_means,FRE_mean=t(FRE_mean),FCE_mean=t(FCE_mean)))
}

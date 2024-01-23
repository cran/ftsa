###########################################################################
# One way median polish  Residuals
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

One_way_Residuals <- function(Y, n_prefectures=51, year=1959:2020, age=0:100)
{
    #Number of years considered in the population
    n_year = length(year)
    #Number of ages considered in each year
    n_age = length(age)

    Pop1 <- One_way_median_polish(Y,n_prefectures=n_prefectures,year=year,age=age)
    FGE = Pop1$grand_effect
    FRE = Pop1$row_effect

    residuals_b1<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
    for(i in 1:nrow(residuals_b1))
    {
        residuals_b1[i,] = Y[i,] - FGE
    }

    # Remove the row effect

    residuals_b1r <- matrix(0, nrow = n_prefectures*n_year, ncol = n_age)
    for(j in 1:n_prefectures)
    {
        residuals_b1r[(n_year*j-(n_year-1)):(n_year*j),]=residuals_b1[(n_year*j-(n_year-1)):(n_year*j),]-t(replicate(n_year,FRE[j,]))
    }
    return(residuals1 = residuals_b1r)
}

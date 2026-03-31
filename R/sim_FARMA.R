sim_FARMA <-
function(n = 500, d = 0.1, no_grid = 101, seed_number, process = c("FAR", "FARMA"))
{
    process = match.arg(process)
    set.seed(123 + seed_number)
  
    # simulate stochastic process realizations 
  
    step_1 = BrownMat(N = 2 * n + 100, refinement = no_grid)

    # calculating beta values
    
    beta_val = vector("numeric", n)
    for(i in 1:n)
    {
        beta_val[i] = exp(lgamma(i + d) - lgamma(i + 1) - lgamma(d))
    }
    
    # Step 2
    
    eta_star = matrix(NA, no_grid, (n+100))
    for(t_val in (n+1):(2*n+100))
    {
        inner_sum = matrix(NA, no_grid, n)
        for(ik in 1:n)
        {
            inner_sum[,ik] = (beta_val[ik] * step_1[,t_val - ik])
        }
        eta_star[,t_val - n] = step_1[,t_val] + rowSums(inner_sum)
    }

    # Step 3
    
    if(process == "FAR")
    {
	      sim_X = funARMat(refinement = no_grid, eta_star_val = eta_star)
	  }
	
	  else if(process == "FARMA")
	  {
		    eta_star_MA = funMAMat(refinement = no_grid, eta_star_val = eta_star)	
		    sim_X = funARMat(refinement = no_grid, eta_star_val = eta_star_MA)
	  }
    else
    {
        warning("DGP should either be FAR or FARMA process.")
    }
        
    # Step 4 (last n observations)
    
    sim_X_record = sim_X[,(ncol(sim_X)-(n-1)):ncol(sim_X)]
    return(sim_X_record)
}

stop_time_sim_data <-
function(sample_size, omega, seed_number)
{
    set.seed(123 + seed_number)
  
    # set a set of grid points

    grid = seq(0, 1, by = 0.01)
    n_grid = length(grid)

    # compute mean function

    mt = err = vector("numeric", n_grid)
    for(i in 1:length(grid))
    {
        mt[i] = 10 * grid[i] * (1 - grid[i])
    }

    sim_mat = matrix(NA, n_grid, sample_size)
    err = as.numeric(rwiener(end = 1, frequency = n_grid))
    sim_mat[,1] = mt + omega * err
    for(ik in 2:ceiling(sample_size/2))
    {
        sim_mat[,ik] = 0.2 * sim_mat[,ik-1] + omega * as.numeric(rwiener(end = 1, frequency = n_grid))
    }

    # introducing a stopping time

    for(ik in (ceiling(sample_size/2) + 1):sample_size)
    {
        sim_mat[,ik] = (0.2 + 0.7) * sim_mat[,ik-1] + omega * as.numeric(rwiener(end = 1, frequency = n_grid))
    }

    # transforming a non-stationary series to stationary

    sim_mat_standard = matrix(NA, n_grid, (sample_size - 1))
    for(ik in 1:(sample_size - 1))
    {
        sim_mat_standard[,ik] = abs(sim_mat[,ik+1] - sim_mat[,ik])/(abs(sim_mat[,ik]) + 0.1)
    }

    sample_size = (sample_size - 1)
    colnames(sim_mat_standard) = 1:sample_size
    data_comb_fts = fts(grid, sim_mat_standard, xname = "Grid point", yname = "Simulated data")
    return(data_comb_fts)
}

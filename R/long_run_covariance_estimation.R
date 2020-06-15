  long_run_covariance_estimation <- function (dat, C0 = 3, H = 3)
  {
      T = ncol(dat)
      no_grid = nrow(dat)
      center_dat = t(scale(t(dat), center = TRUE, scale = FALSE))
      cov_l <- function(porder, band, nval, kern_type) {
          cov_sum = gamma_l(0, nval)
          if (kern_type == "FlatTop") {
              for (ik in 1:(nval - 1)) {
                  cov_sum = cov_sum + FlatTop(ik/band) * abs(ik)^(porder) *
                    (gamma_l(ik, nval) + t(gamma_l(ik, nval)))
              }
          }
          else {
              for (ik in 1:(nval - 1)) {
                  cov_sum = cov_sum + kweights(ik/band, kernel = kern_type) *
                    abs(ik)^(porder) * (gamma_l(ik, nval) + t(gamma_l(ik,
                    nval)))
              }
          }
          return(cov_sum)
      }
      gamma_l <- function(lag, T) {
          gamma_lag_sum = 0
          if (lag >= 0) {
              for (ij in 1:(T - lag)) {
                  gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,
                    ij]) %*% t(as.matrix(center_dat[, (ij + lag)])))
              }
          }
          else {
              for (ij in 1:(T + lag)) {
                  gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,
                    ij - lag]) %*% t(as.matrix(center_dat[, ij])))
              }
          }
          return(gamma_lag_sum/T)
      }
      gamma_hat = list()
      for (ik in 1:T) {
          gamma_hat[[ik]] = gamma_l(lag = ik - 1, T = T)
      }
      rho_hat <- function(gamma, ik) {
          a = sum(diag(gamma[[1]]))
          denom = a * a
          b = gamma[[ik + 1]] * gamma[[ik + 1]]
          num = sum(b)
          c = num/denom
          calc = sqrt(c)
          return(calc)
      }
      rho = vector(, T - 1)
      for (ik in 1:(T - 1)) {
          rho[ik] = rho_hat(gamma_hat, ik)
      }
      h_hat <- function(rho, N, H) {
          l <- -1
          c <- 1
          end <- N - H
          const <- C0 * sqrt(log10(N)/N)
          while (c <= end) {
              if (rho[c] <= const) {
                  c2 <- 1
                  while (c2 <= H - 1) {
                    if (rho[c + c2] <= const) {
                      c2 <- c2 + 1
                      l <- c - 1
                    }
                    else {
                      c2 <- H
                      l <- -1
                      c <- c + 1
                    }
                  }
                  if (l != -1) {
                    c <- end + 1
                  }
              }
              else {
                  c <- c + 1
              }
          }
          return(l)
      }
      h_adaptive_ini <- h_hat(rho = rho, N = T, H = H)
      h_opt_fun <- function(band_ini, T, kern_type, kern_type_ini) {
          kern_name_kweights_ini = switch(kern_type_ini, BT = "Bartlett",
              PR = "Parzen", TH = "Tukey-Hanning", QS = "Quadratic Spectral",
              FT = "FlatTop")
          kern_name_kweights = switch(kern_type, BT = "Bartlett",
              PR = "Parzen", TH = "Tukey-Hanning", QS = "Quadratic Spectral")
          w_weights = switch(kern_type, BT = 1, PR = 6, TH = pi^2/4,
              QS = 18 * pi * pi/125)
          q = switch(kern_type, BT = 1, PR = 2, TH = 2, QS = 2)
          C_0 = cov_l(porder = 0, band = band_ini, nval = T, kern_type = kern_name_kweights_ini)
          C_2 = w_weights * cov_l(porder = q, band = band_ini,
              nval = T, kern_type = kern_name_kweights_ini)
          first_part = (2 * q * sum(C_2^2))^(1/(1 + 2 * q))
          kernel_square_int = switch(kern_type, BT = 2/3, PR = 0.539285,
              TH = 3/4, QS = 1, FT = 4/3)
          second_part = ((sum(C_0^2) + sum(diag(C_0))^2) * kernel_square_int)^(-1/(1 +
              2 * q))
          c_0 = first_part * second_part
          hat_h_opt = c_0 * (T^(1/(1 + 2 * q)))
          C_0_est = cov_l(porder = 0, band = hat_h_opt, nval = T,
              kern_type = kern_name_kweights)
          return(list(hat_h_opt = hat_h_opt, C_0_est = C_0_est))
      }
      h_opt_FT_BT = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "BT",
          kern_type_ini = "FT")
      return(h_opt_FT_BT$C_0_est)
  }

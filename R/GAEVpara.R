GAEVpara = function(data,kdf1,kdf2)
{
  Tperiod = nrow(data)
  ngrid = ncol(data)

  loc_est_para = matrix(0,nrow = Tperiod,ncol=kdf1)
  logscale_est_para = matrix(0, nrow=Tperiod,ncol=kdf2)
  logscale_est = loc_est = matrix(0,nrow = Tperiod,ncol=ngrid)
  shape_est = matrix(0,nrow = Tperiod,ncol=1)

  for (t in 1:Tperiod){

    yeardata = data.frame(temp=data[t,],day=1:ngrid, kdf1, kdf2)
    form1 <- as.formula(glue("temp ~ s(day, k = {kdf1}, bs = 'cr')"))
    form2 <- as.formula(glue("~ s(day, k = {kdf2}, bs = 'cr')"))
    fmla_gev = list(form1, form2, ~ 1)

    # fmla_gev = list(temp ~ s(day,k=kdf1.global,bs= 'cr'), ~ s(day,k=kdf2.global,bs= 'cr'), ~ 1)
    m_gev = evgam(formula = fmla_gev, data = yeardata, family='gev')
    loc_est[t,] = m_gev$location$fitted
    loc_est_para[t,] =  m_gev$location$coefficients

    logscale_est[t,] = m_gev$logscale$fitted
    logscale_est_para[t,] = m_gev$logscale$coefficients
    shape_est[t,] = m_gev$shape$fitted[1]
  }
  para.est = list(loc = list(est = loc_est, para = loc_est_para, basis = m_gev$location$X),
                  logscale = list(est = logscale_est, para = logscale_est_para, basis = m_gev$logscale$X),
                  shape = shape_est)

  return (para.est)
}

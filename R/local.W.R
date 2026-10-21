local.W <-
function(data, m, int = c(-0.5, 2.5), taper = c("none", "Velasco", "HC"), 
                    diff_param = 1, l = 1)
{
    taper<-taper[1]
    if((taper%in%c("none","Velasco","HC")) == FALSE) 
    {
        stop('taper must be either "none", "Velasco", or "HC".')
    }
    T <- length(data)
    if(taper == "Velasco")
    {
        p <- 3
        ht <- cos_bell((1:T)/T)
        data <- ht*data
        peri <- per(data)[-1]
        d.hat <- optimize(f=R.lw.tapered, interval=int, peri=peri,  m=m, T=T, p=p)$minimum
        se <- sqrt(p*PHI_lw(ht=ht,T=T,p=p)/(4*m))
    }
    if(taper == "HC")
    {
        data <- diff(data, differences = diff_param)
        T <- length(data)
        ht <- cos_bell_cmplx((1:T)/T)
        data <- ht*data
        peri <- per(data)[-1]
        d.hat <- optimize(f = R.lw.hc, interval = int, peri = peri, m = m, T = T)$minimum + 1
        se <- sqrt(1.5/(4*m))
    }
    if(taper == "none")
    {
        peri <- per(data)[-1]
        d.hat <- optimize(f=R.lw, interval=int, peri=peri,  m=m, T=T, l=1)$minimum
        se <- 1/(2*sqrt(m))
    }
    return(list("d" = d.hat, "s.e." = se))
}

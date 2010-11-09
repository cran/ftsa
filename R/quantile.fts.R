quantile.fts = function(x, probs = c(0.25, 0.75), ...)
{
  if (class(x)[1] == "fts"|class(x)[1] == "fds"|class(x)[1] =="sfts"){
      y = x$y
      p = dim(y)[1]
      q = matrix(, p, length(probs))
      for(i in 1:p){
          q[i,] = quantile(y[i,], probs = probs)
      }
      colnames(q) = c(0.25,0.75)
      rownames(q) = x$x
      if (class(x)[1] == "fds"){
          warning("Object is not a functional time series.")
      }
      return(q)
   }
   else {
        stop("Not a functional object.")
   }
}

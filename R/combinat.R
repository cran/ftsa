combinat <-
function(n,p){
  if (n<p){combinat=0}
  else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
}

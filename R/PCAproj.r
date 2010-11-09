PCAproj=function (x, k = 2, method = c("mad", "sd", "qn"), CalcMethod = c("eachobs",
    "lincomb", "sphere"), nmax = 1000, update = TRUE, scores = TRUE,
    maxit = 5, maxhalf = 5, scale = NULL, center = l1median,
    control)
{
    if (!missing(control))
        ParseControlStructure(control, c("k", "method", "CalcMethod",
            "nmax", "update", "scores", "maxit", "maxhalf"))
    CalcMethod = CalcMethod[1]
    method = method[1]
    if (!any(CalcMethod == c("eachobs", "lincomb", "sphere")))
        stop(paste("Unknown calcmethod:", CalcMethod))
    n = nrow(x)
    p = ncol(x)
    if (k > min(n, p))
        stop("k too large")
    if (p > n) {
        svdx = svd(t(x))
        x = svdx$v %*% diag(svdx$d)
        pold = p
        p = n
    }
    else pold = p
    DataObj = ScaleAdv(x, scale = scale, center = center)
    if (pold > n) {
        DataObj$center <- as.vector(svdx$u %*% DataObj$center)
        DataObj$scale <- ScaleAdv(x %*% t(svdx$u), center = NULL,
                                  scale = scale)$scale
    }
    y = DataObj$x
    if (scores)
        scoresize = n * k
    else scoresize = 0
    if (CalcMethod == "lincomb") {
        update = FALSE
        if (nmax > n) {
            aux = matrix(runif((nmax - n) * n), nrow = nmax -
                n)
            y = rbind(y, t(t(aux %*% as.matrix(x)) - DataObj$center))
        }
    }
    else if (CalcMethod == "sphere") {
             update = FALSE
             if (nmax > n)
                 y = rbind(y, rmvnorm(nmax - n, rep(0, p), diag(p)))
    }
    nn = nrow(y)
    if (update)
        ret = .C("rpcnup", PACKAGE = "pcaPP", as.double(y), as.integer(c(nn,
            p, k, ParseDevString(method), scores, maxit, maxhalf,
            n)), scores = double(scoresize), loadings = double(p *
            k), lambda = double(k))
    else ret = .C("rpcn", PACKAGE = "pcaPP", as.double(y), as.integer(c(nn,
        p, k, ParseDevString(method), scores, n)), scores = double(scoresize),
        loadings = double(p * k), lambda = double(k))
    veig = matrix(ret$loadings, ncol = k)
    if (pold > n)
        veig = svdx$u %*% veig
    if (scores)
        DataPostProc(DataObj, ret$lambda, veig, matrix(ret$scores,
                     ncol = k), match.call(), scores)
    else DataPostProc(DataObj, ret$lambda, veig, NULL, match.call(),
                      scores)
}

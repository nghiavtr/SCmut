######################### 
scfdr <- function(sc.ref, sc.alt,  mut.sites, germstat=NULL, germAA=NULL, 
     seed=NULL, nrep=100, qlim= 0.999, nr =25,  vaf.cut=0.2, 
     n0=0.1, m=2, smooth = 0.9, verb = TRUE,  ...) 
## sc.ref, sc.alt = reference and alternative counts statistics from the pipeline
##                  include only single-cell data
## mut.sites      = names of candidate mutation sites 
##                  names must match with the rownames for sc.ref
## germstat = name of germline count variable in the full, then use germfun() to extract germline AA sites
## germAA = germline AA sites (a boolean vector, corresponding to each row of sc.ref) input directly from users
##
## seed = seed for set.seed() 
## nrep = number of random replicates of the null/non-mutated datasets
## qlim = quantile limit of the total reads to avoid outliers
## 
## nr = number of grid points: 2D-area is split to nr x nr
## n0 = proportion of points near the origin where the null is enforced (so fdr=1)
## vaf.cut = threshold of the vaf: If vaf<cut then fdr=1
## m = moderating constant to compute VAFs, also avoiding 0/0
## smooth = smoothing parameter
##
## output: 
## (x,y,fdr.xy) triplet of fdr2d  on the nrxnr 2D-region
## ifdr = individual site-cell fdr estimate: a matrix of mut.sites x ncell
##  
{
    vcat = function(...) if (verb) 
        cat(...)
   if ((smooth < 0.01) | (smooth > 0.99)) {
        smooth = min(max(smooth, 0.01), 0.99)
        warning("smooth must be between 0.01 and 0.99 - reset to ", 
            smooth)
    }

  ## germ-line AA sites: used in the filter
  if (is.null(germAA))  germAA = germfun(germstat)

  ## observed-mutation sites
  x0.obs = c(sc.ref[mut.sites,])
     x.obs = x0.obs[!is.na(x0.obs)]
   y0.obs = c(sc.alt[mut.sites,])
     y.obs = y0.obs[!is.na(x0.obs)]
   nread.obs = x.obs+y.obs
   vaf.obs = y.obs/nread.obs

  ## null statistics: random sample from presumed non-mutation sites, 
  ##   must be reasonable germline AA 
  ## note the seed is used here
   filt = germAA & !(rownames(sc.ref) %in% mut.sites)
   n = length(mut.sites)   ## number of mutation sites to test
   if (!is.null(seed)) set.seed(seed)
   ran = sample(rownames(sc.ref)[filt],n*nrep)
   x0 = c(sc.ref[ran,])
      ranX = x0[!is.na(x0)]
   y0 = c(sc.alt[ran,])
      ranY = y0[!is.na(x0)]

## statistics for fdr
   nread.null = ranX+ ranY   ## nreads
   pick=nread.null>0
   ranX=ranX[pick]
   ranY=ranY[pick]
   nread.null=nread.null[pick]
   
   vaf.null = (ranY)/nread.null    ## vaf here!!
   x = c(nread.obs, nread.null); xnull= nread.null
   y = c(vaf.obs, vaf.null); ynull = vaf.null

## break points
   xbreaks = MidBreaks(x[x<quantile(x,qlim, na.rm=TRUE) & x>quantile(x,1-qlim, na.rm=TRUE)], nr)
   ybreaks = MidBreaks(y[y<quantile(y,qlim, na.rm=TRUE) & y>quantile(y,1-qlim, na.rm=TRUE)], nr)
   xmid = xbreaks[-1] - diff(xbreaks)/2
   ymid = ybreaks[-1] - diff(ybreaks)/2

## tabulated values
    cut.stat.x = cut(x, xbreaks, include.lowest = TRUE)
    cut.stat.y = cut(y, ybreaks, include.lowest = TRUE)
    cut.null.x = cut(xnull, xbreaks, include.lowest = TRUE)
    cut.null.y = cut(ynull, ybreaks, include.lowest = TRUE)
    tab.stat = table(cut.stat.x, cut.stat.y)
    tab.null = table(cut.null.x, cut.null.y)
 
### need to equalize near the (0,0) point to get fdr=1
    N0 = n0*nr; N0 = ceiling(N0)
    tab.stat[1:N0,1:N0]= tab.null[1:N0,1:N0]= 
        (tab.stat[1:N0,1:N0]+ tab.null[1:N0,1:N0])/2   

## smoothing fdr2d
    df.smooth = (1 - smooth) * (nr - 1)^2
    b = smooth2d.yn(tab.null, tab.stat, df = df.smooth)
    ## note: bfit ~y/n ~ nrep*f0/(nrep*f0 + f)
    ## ==> nrep*f0 = bfit*nrep*f0+ bfit*f)
    ## ==> f0/f = bfit/(1-bfit)/nrep
    f0fz = b$fit/(1 - b$fit)/nrep 
       fdr = ifelse(f0fz>1,1,f0fz)
    dimnames(fdr) = dimnames(tab.stat)

##  check code with SGI BC data, run with sc-runfdr.r
check = FALSE
if (check){
 contour(xmid, ymid,f0fz, levels=seq(0.1,1,len=10), xlab='scTotal reads', ylab='scVAF')
   points(nread.obs, vaf.obs, pch=8, col='red')
 tum = cellType$index=='Tumor'
   tum.mat = matrix(rep(tum,length(mut.sites)),nrow=length(mut.sites), byrow=TRUE)
   tum.indic = c(tum.mat)[!is.na(x0.obs)]
   points(nread.obs[tum.indic], vaf.obs[tum.indic], pch=16, col='red')
   points(nread.obs[!tum.indic], vaf.obs[!tum.indic], pch=16, col='green')
}

## indiv fdr values, if needed: cut the fc to reduce computation
    pick = vaf.obs> vaf.cut
    ifdr = rep(1,length(vaf.obs))
    ifdr[pick] = approx2d(xmid, ymid, fdr, nread.obs[pick], vaf.obs[pick])

## put fdr values to the original cell positions in result matrix
    ifdr.vec = c(sc.ref[mut.sites,])
      ifdr.vec[!is.na(ifdr.vec)] = ifdr  ## otherwise keep NA
      ifdr.mat = matrix(ifdr.vec, nrow=length(mut.sites))
        rownames(ifdr.mat) = mut.sites
        colnames(ifdr.mat) = colnames(sc.ref)

if (check){
  signif = c(ifdr.mat<0.2)[!is.na(x0.obs)]
  points(nread.obs[signif], vaf.obs[signif], pch=0, cex=1.5)
}

## output
    out = list(x= xmid, y=ymid, fdr.xy=fdr, ifdr=ifdr.mat)
    return(out)
}

## Confidence of a real zero of alt-allele, 
##    accounting for coverage and errors
##
## conf0 = -10*log10(upper-confidence bound of prob)
## based on observing x successes out of n trials
## Example:
#> conf0(0,5); conf0(0,10); conf0(0,30)
#[1] 3.460869 [1] 5.869316 [1] 10.22103
#
conf0 = function(x,n,alpha=0.05){
fn = function(th,x,n, alpha){
    pleft=  pbinom(x,n,th) 
    return(pleft-alpha)
  }
 if (is.na(x)) return(0)
 if (x >= n)  return(0)
 run = uniroot(fn,c(x/n,1),x=x,n=n, alpha=alpha)
 return(-10*log10(run$root))
}


## statistics from the count matrix
## germstat= germ counts in DNA sequencing, Nx2 matrix of (lab[1]='alt', lab[2]='total')_reads
## gnlim = coverage limit to get confident AA call
## conf.fun = evaluate conf0 score, quite slow, so default=FALSE
##
germfun = function(germstat,lab = c('alt','total'), gnlim= 10, conf.fun =FALSE, conf.lim=6){
  gx = germstat[,lab[1]]    ## alternative allele counts
  gn = germstat[,lab[2]]    ## germline coverage
 
 if (!conf.fun){
   germAA =  (gn> gnlim) & (gx==0)  ## germline AA
     germAA = ifelse(is.na(germAA),FALSE, germAA)
  }
 if (conf.fun){  ## rather slow!  
   gconf0 =  mapply(conf0, gx, gn)
   germAA = gconf0 > conf.lim
     germAA = ifelse(is.na(germAA),FALSE, germAA)
  }
 return(germAA)
}


# from OCplus
MidBreaks <- function (x, nr) 
{
    mids = seq(min(x), max(x), length = nr - 1)
    d = mids[2] - mids[1]
    breaks = c(mids - d/2, mids[nr - 1] + d/2)
    breaks
}

## adapted from OCplus simul2d.direct, but with a modification on n==0
## smoothing y/n= proportion of nulls
## NOTE: if n==0 set y/n=1  (all null, so fdr=1)
##
smooth2d.yn<- 
function (y, n, df = 100, err = 0.01, edge.count = 3, mid.start = 6) 
{
    nx = nrow(y)
    ny = ncol(y)
    if (nx < mid.start | edge.count > ny) 
        stop("Grid size too small")
    ## n = ifelse(n <= 1, 1, n)  ## wrong default here!
    edge = 1:edge.count
    mid = mid.start:(nx - mid.start + 1)
    y[mid, edge] = ifelse(n[mid, edge] == 1, 1, y[mid, edge])
    p = y/n
      p = ifelse(n<1, 1, p) ## if n==0 set p=1  (results= null)
    sv2 = df2var(nx, df = df)
    fit = smooth2d.basic(p, sv2 = sv2, err = err)
    list(fit = fit, df = df, sv2 = sv2)
}

## from OCplus
approx2d<- 
function (x, y, z, xout, yout, ...) 
{
    n = length(x)
    m = length(y)
    no = length(xout)
    zout = rep(NA, no)
    for (i in 2:n) {
        for (j in 2:m) {
            ndx = x[i - 1] <= xout & xout <= x[i] & y[j - 1] <= 
                yout & yout <= y[j]
            lower = z[i - 1, j - 1] + (xout[ndx] - x[i - 1]) * 
                (z[i, j - 1] - z[i - 1, j - 1])/(x[i] - x[i - 
                1])
            upper = z[i - 1, j] + (xout[ndx] - x[i - 1]) * (z[i, 
                j] - z[i - 1, j])/(x[i] - x[i - 1])
            zout[ndx] = lower + (yout[ndx] - y[j - 1]) * (upper - 
                lower)/(y[j] - y[j - 1])
        }
    }
    undef = which(is.na(zout))
    if (length(undef > 0)) {
        gg = as.matrix(expand.grid(x, y))
        for (i in undef) {
            d2 = (xout[i] - gg[, 1])^2 + (yout[i] - gg[, 2])
            nn = which.min(d2)
            zout[i] = z[nn]
        }
    }
    zout
}

smooth2d.basic <- function (y, sv2 = 0.01, se2 = 1, err = 0.01) 
{
    W = 1/se2
    nx = nrow(y)
    ny = ncol(y)
    Y = c(y)
    R = rmatrix(nx, ny)
    d = R$d/sv2 + W
    xx = R$x
    ybar = mean(Y)
    B = W * (Y - ybar)
    v = gausei(d, xx, B, sv2, err = err, maxiter = 100)
    est = ybar + v
    matrix(est, ncol = ny)
}

gausei <- function (d, xx, B, sv2, err = 0.01, maxiter = 10) 
{
    ng = length(d)
    xx[xx < 1] = ng + 1
    v = c(B/(d + 1e-06), 0)
    vold = 3 * v
    iter = 1
    while ((sd(vold - v, na.rm = TRUE)/sd(v, na.rm = TRUE) > 
        err) & (iter < maxiter)) {
        vold = v
        vmat = matrix(v[xx], ncol = 4)
        v[1:ng] = (B + c(vmat %*% rep(1, 4))/sv2)/d
        iter = iter + 1
    }
    v[1:ng]
}

rmatrix <- function (nx, ny) 
{
    xymat = matrix(0, nrow = nx, ncol = ny)
    grid = c((row(xymat) - 1) * 1000 + col(xymat) - 1)
    neighbor = c(grid + 1000, grid - 1000, grid + 1, grid - 1)
    loc = match(neighbor, grid, 0)
    loc = matrix(loc, ncol = 4)
    num = c((loc > 0) %*% rep(1, 4))
    list(x = loc, d = num)
}

df2var <- function (nx, df = nx, lower = 1e-04, upper = 1000) 
{
    ff = function(x) dgfree(nx, x) - df
    out = uniroot(ff, lower = lower, upper = upper)
    return(out$root)
}

dgfree <- function (nx, sv2 = 1, ny = nx) 
{
    imp = matrix(0, nrow = nx, ncol = ny)
    imp[nx/2, ny/2] = 1
    imp.resp = smooth2d.basic(imp, sv2 = sv2)
    df = nx * ny * imp.resp[nx/2, ny/2]
    return(df)
}



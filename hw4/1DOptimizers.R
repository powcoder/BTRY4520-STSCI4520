NewtonRaphson = function(fn,dfn,x0,tol=1e-8,maxit=100){
  # fn = function for which a root is needed. 
  # df = derivative of the function for which a root is needed. 
  # x0 = starting guess for the root
  # tol = tolerance, stop when |fn(x)| < tol
  # maxit = maximum number of iterations   

  f0 = fn(x0); df0 = dfn(x0);    # Initialization
  
  tol.met = FALSE    # No tolerance met
  iter = 0           # No iterations

  while(!tol.met){
    x0 = x0 - f0/df0
    f0 = fn(x0); df0 = dfn(x0)
    
    iter = iter + 1 # Update iterations and tolerance
    if( abs(f0) < tol | iter > maxit ){
      tol.met=TRUE
    }
#    print(paste(c(iter,x0,f0)))
  }
  return(list(sol=x0,iter=iter))
}




GoldenSection = function(fn,xl,xr,tol=1e-8,maxit=100){
  gr = (1 + sqrt(5))/2  # Pre-calculate golden ratio

  xm = xl + (xr-xl)/(1+gr)  
  fl = fn(xl); fr = fn(xr); fm = fn(xm)

  tol.met = FALSE    # No tolerance met
  iter = 0           # No iterations

  while(!tol.met){
    iter = iter + 1
    if( (xr-xm) > (xm-xl) ){  # Working on the right-hand side
      y = xm + (xr-xm)/(1+gr); fy = fn(y);
      if( fy > fm){ xl = xm; fl = fm; xm = y; fm = fy }
      else{ xr = y; fr = fy }
    }
    else{
      y = xm - (xm-xl)/(1+gr); fy = fn(y);  
      if( fy > fm){ xr = xm; fr = fm; xm = y; fm = fy }
      else{ xl = y; fl = fy }
    }
#    print(c(iter,xm))
    if( (xr-xm) < tol | iter > maxit ){ tol.met=TRUE }
  }
  return(list(xm=xm,iter=iter))
}

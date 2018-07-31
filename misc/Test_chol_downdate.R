chol_update1 = function(L,x,sign = 1){
  n = nrow(L)
  x = as.matrix(x)
  for(i in 1:ncol(x)) {
    for(k in 1:n){
      L_kk = L[k,k]
      x_k = x[k,i]
      r = sqrt(L_kk^2 + sign * x_k^2)
      c = r/L_kk
      s = x_k / L_kk
      L[k,k] = r
      if(k < n) {
        index = (k+1):n
        L[index,k] = (L[index,k] + sign *  s*x[index,i])/c
        x[index,i] = c*x[index,i] - s*L[index,k]
      }
    }
  }
  return(L)
}

x  =rstdnorm_mat(1000,2000)
X = tcrossprod(x) + diag(1,nrow(x))
L = t(chol(X))

y = x[,1:5,drop=FALSE]
Xresid = X - y %*% t(y)
L2 = t(chol(Xresid))
chol_update(L,y,sign = -1) - L2
microbenchmark(
t(chol(Xresid)),
chol_update1(L,y,sign = -1),
chol_update(L,y,sign = -1),
times=10)

chol_dropRows1 = function(L,start_row, end_row) {
  n = nrow(L)
  if(start_row > 1) {
    top_rows = 1:(start_row-1)
  } else{
    top_rows = c()
  }
  if(end_row < n) {
    bottom_rows = (end_row+1):n
  } else{
    bottom_rows = c()
  }
  
  L11 = L[top_rows,top_rows]
  L31 = L[bottom_rows,top_rows]
  S33 = L[bottom_rows,bottom_rows]
  S23 = L[bottom_rows,-c(top_rows,bottom_rows)]
  
  n_new = length(c(top_rows,bottom_rows))
  Lnew = matrix(0,n_new,n_new)
  Lnew[top_rows,top_rows] = L11
  if(length(bottom_rows)>0) {
    Lnew_bottomRows = length(top_rows) + 1:length(bottom_rows)
    Lnew[Lnew_bottomRows,top_rows] = L31
    L33 = chol_update1(S33,S23,1)
    Lnew[Lnew_bottomRows,Lnew_bottomRows] = L33
  }
  
  return(Lnew)
}

Ls = t(chol(X[-c(6:10),-c(6:10)]))
chol_dropRows1(L,6,10) - Ls
chol_dropRows(L,6,5) - Ls

microbenchmark(chol(X[-c(600+(1:10)),-c(600+(1:10))]),chol_dropRows(L,600,11),times=10)

               
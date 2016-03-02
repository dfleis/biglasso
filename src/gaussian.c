#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

// c++ functions for big.matrix
int get_row_bm(SEXP xP);
int get_col_bm(SEXP xP);
double crossprod_bm(SEXP xP, double *y_, int *row_idx_, double center_, 
                    double scale_, int n_row, int j);
double get_elem_bm(SEXP xP, double center_, double scale_, int i, int j);
double crossprod_resid(SEXP xP, double *y_, double sumY_, int *row_idx_, 
                       double center_, double scale_, int n_row, int j);
double sum(double *x, int n);
int checkConvergence(double *beta, double *beta_old, double eps, int l, int J);
double lasso(double z, double l1, double l2, double v);

// Memory handling, output formatting (Gaussian)
SEXP cleanupG(double *a, double *r, int *e1, int *e2, double *z, 
              SEXP beta, SEXP loss, SEXP iter) {
  Free(a);
  Free(r);
  Free(e1);
  Free(e2);
  Free(z);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, loss);
  SET_VECTOR_ELT(res, 2, iter);
  UNPROTECT(4);
  return(res);
}

// Gaussian loss
double gLoss(double *r, int n) {
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
}

// Coordinate descent for gaussian models
SEXP cdfit_gaussian(SEXP X_, SEXP y_, SEXP row_idx_, SEXP center_, SEXP scale_, 
                    SEXP lambda, SEXP eps_, SEXP max_iter_, 
                    SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_) {
  // Declarations
  int n = length(row_idx_);
  int p = get_col_bm(X_);
  int L = length(lambda);
  
  SEXP res, beta, loss, iter;
  PROTECT(beta = allocVector(REALSXP, L*p));
  double *b = REAL(beta);
  for (int j=0; j<(L*p); j++) b[j] = 0;
  PROTECT(loss = allocVector(REALSXP, L));
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  double *a = Calloc(p, double); // Beta from previous iteration
  for (int j=0; j<p; j++) a[j]=0;

  double *y = REAL(y_);
  int *row_idx = INTEGER(row_idx_);
  double *center = REAL(center_);
  double *scale = REAL(scale_);
  
// const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  double *lam = REAL(lambda);
  double eps = REAL(eps_)[0];
  int max_iter = INTEGER(max_iter_)[0];
//  double gamma = REAL(gamma_)[0];
  double *m = REAL(multiplier);
  double alpha = REAL(alpha_)[0];
  int dfmax = INTEGER(dfmax_)[0];
  int user = INTEGER(user_)[0];
  double *r = Calloc(n, double);
  for (int i=0; i<n; i++) r[i] = y[i];
  double *z = Calloc(p, double);
  
  double sumResid = sum(r, n);
  
  for (int j=0; j<p; j++) {
    z[j] = crossprod_bm(X_, r, row_idx, center[j], scale[j], n, j)/n;
  }
  
  int *e1 = Calloc(p, int);
  for (int j=0; j<p; j++) e1[j] = 0;
  int *e2 = Calloc(p, int);
  for (int j=0; j<p; j++) e2[j] = 0;
  double cutoff, l1, l2;
  int converged, lstart;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user) {
    lstart = 0;
  } else {
    REAL(loss)[0] = gLoss(r,n);
    lstart = 1;
  }
  
  // Path
  for (int l=lstart;l<L;l++) {
    if (l != 0) {
      // Assign a
      
      for (int j=0;j<p;j++) {
        a[j] = b[(l-1)*p+j];
      }

      // Check dfmax
      int nv = 0;
      for (int j=0; j<p; j++) {
	      if (a[j] != 0) nv++;
      }
      if (nv > dfmax) {
      	for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
      	res = cleanupG(a, r, e1, e2, z, beta, loss, iter);
      	return(res);
      }

      // Determine eligible set
      cutoff = 2*lam[l] - lam[l-1];
      
      for (int j=0; j<p; j++) {
        if (fabs(z[j]) > (cutoff * alpha * m[j])) e2[j] = 1;
      } 

    } else {
        // Determine eligible set
        double lmax = 0;
        for (int j=0; j<p; j++) if (fabs(z[j]) > lmax) lmax = fabs(z[j]);
        cutoff = 2*lam[l] - lmax;
        for (int j=0; j<p; j++) if (fabs(z[j]) > (cutoff * alpha * m[j])) e2[j] = 1;
    }

    while (INTEGER(iter)[l] < max_iter) {
      while (INTEGER(iter)[l] < max_iter) {
      	while (INTEGER(iter)[l] < max_iter) {
      	  // Solve over the active set
      	  INTEGER(iter)[l]++;
      	  for (int j=0; j<p; j++) {
      	    if (e1[j]) {
      	      //z[j] = crossprod_bm(X_, r, row_idx, center[j], scale[j], n, j)/n + a[j];
      	      z[j] = crossprod_resid(X_, r, sumResid, row_idx, center[j], scale[j], n, j) / n + a[j];
      	      // Update beta_j
      	      l1 = lam[l] * m[j] * alpha;
      	      l2 = lam[l] * m[j] * (1-alpha);
      	      
      	      b[l*p+j] = lasso(z[j], l1, l2, 1);
      	      
      	      // Update r
      	      double shift = b[l*p+j] - a[j];
      	      if (shift !=0) {
      	        for (int i=0;i<n;i++) r[i] -= shift * get_elem_bm(X_, center[j], scale[j], row_idx[i], j);
      	        sumResid = sum(r, n); //update sum of residual
      	      }
      	    }
      	  }
      	  // Check for convergence
      	  converged = checkConvergence(b, a, eps, l, p);
      	  
      	  for (int j=0; j<p; j++) {
      	    a[j] = b[l*p+j];   
      	  }
      	  if (converged) break;
      	}

      	// Scan for violations in strong set
      	int violations = 0;

        for (int j=0; j<p; j++) {
          if (e1[j]==0 & e2[j]==1) {
            //z[j] = crossprod_bm(X_, r, row_idx, center[j], scale[j], n, j)/n;
            z[j] = crossprod_resid(X_, r, sumResid, row_idx, center[j], scale[j], n, j) / n;
            
            // Update beta_j
            l1 = lam[l] * m[j] * alpha;
            l2 = lam[l] * m[j] * (1-alpha);
            
            b[l*p+j] = lasso(z[j], l1, l2, 1);
            
            // If something enters the eligible set, update eligible set & residuals
            if (b[l*p+j] !=0) {
              e1[j] = e2[j] = 1;
              for (int i=0; i<n; i++) r[i] -= b[l*p+j] * get_elem_bm(X_, center[j], scale[j], row_idx[i], j);
              sumResid = sum(r, n); //update sum of residual
              
              a[j] = b[l*p+j];
              violations++;
            }
          }
        }

    	  if (violations==0) break;
      }

      // Scan for violations in rest
      int violations = 0;

      for (int j=0; j<p; j++) {
        if (e2[j]==0) {
          //z[j] = crossprod_bm(X_, r, row_idx, center[j], scale[j], n, j)/n;
          z[j] = crossprod_resid(X_, r, sumResid, row_idx, center[j], scale[j], n, j) / n;
          
          // Update beta_j
          l1 = lam[l] * m[j] * alpha;
          l2 = lam[l] * m[j] * (1-alpha);
          b[l*p+j] = lasso(z[j], l1, l2, 1);
          // If something enters the eligible set, update eligible set & residuals
          if (b[l*p+j] !=0) {
            e1[j] = e2[j] = 1;
            for (int i=0; i<n; i++) r[i] -= b[l*p+j] * get_elem_bm(X_, center[j], scale[j], row_idx[i], j);
            sumResid = sum(r, n);
            
            a[j] = b[l*p+j];
            violations++;
          }
        }
      }
      
      if (violations==0) {
      	REAL(loss)[l] = gLoss(r, n);
      	break;
      }
    }
  }
  res = cleanupG(a, r, e1, e2, z, beta, loss, iter);
  return(res);
}

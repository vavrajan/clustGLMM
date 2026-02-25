/*
 * Mathematical functions needed
 */

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "structures.h"

double logit_inv(double x){
  double expx;
  
  if(x < 0){
    expx = exp(x);
    return(expx/(1.0+expx));
  }else{
    expx = exp(-x);
    return(1.0/(1.0+expx));
  }
}

double safe_log1pexp(double x) {
  if (x > 700) {
    // e^700 overflows double; log(1+e^x) is approx x
    return x; 
  } else if (x > 0) {
    // log(1+e^x) = log(e^x * (e^-x + 1)) = x + log(1+e^-x)
    return x + log1p(exp(-x));
  } else if (x < -37) {
    // e^-37 is very small; log(1+e^x) is approx e^x
    return exp(x);
  } else {
    // Small to medium negative values: use log1p for precision
    return log1p(exp(x));
  }
}
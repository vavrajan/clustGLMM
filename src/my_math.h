#ifndef MY_MATH_H
#define MY_MATH_H

double logit_inv(double x);
double safe_log1pexp(double x);
double ord_qk(int k, double eta, double* c, int* Kord);

#endif
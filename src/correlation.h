#ifndef CORRELATION_H
#define CORRELATION_H

double pearson (const double data1[], const double data2[], const size_t n);

double spearman (const double data1[], const double data2[], const size_t n, double work[]);

double spearman_permutation_test (const double data1[], const double data2[], const size_t n,
		double work1[], double work2[], unsigned int *seed, const double rho);

#endif /* correlation.h */

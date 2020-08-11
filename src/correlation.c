/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2020 Thiago L. A. Miller <tmiller@mochsl.org.br
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"

#include <stdlib.h>
#include <math.h>
#include "correlation.h"

#define PERMUTATION_SIZE 1001

/*
 * FROM sort/sortvec_source.c - GSL
 *
 * Implement Heap sort -- direct and indirect sorting
 * Based on descriptions in Sedgewick "Algorithms in C"
 *
 * Copyright (C) 1999  Thomas Walter
 *
 * 18 February 2000: Modified for GSL by Brian Gough
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 3, or (at your option) any
 * later version.
 *
 * This source is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 */

static inline void
my_downheap2 (double *data1, double *data2, const size_t N, size_t k)
{
	double v1 = data1[k];
	double v2 = data2[k];

	while (k <= N / 2)
		{
			size_t j = 2 * k;

			if (j < N && data1[j] < data1[(j + 1)])
				{
					j++;
				}

			if (!(v1 < data1[j]))  /* avoid infinite loop if nan */
				{
					break;
				}

			data1[k] = data1[j];
			data2[k] = data2[j];

			k = j;
		}

	data1[k] = v1;
	data2[k] = v2;
}

static inline void
sort2 (double *data1, double *data2, const size_t n)
{
	size_t N;
	size_t k;

	if (n == 0)
		{
			return;   /* No data to sort */
		}

	/* We have n_data elements, last element is at 'n_data-1', first at
	'0' Set N to the last element number. */

	N = n - 1;

	k = N / 2;
	k++;           /* Compensate the first use of 'k--' */
	do
		{
			k--;
			my_downheap2 (data1, data2, N, k);
		}
	while (k > 0);

	while (N > 0)
		{
			/* first swap the elements */
			double tmp;

			tmp = data1[0];
			data1[0] = data1[N];
			data1[N] = tmp;

			tmp = data2[0];
			data2[0] = data2[N];
			data2[N] = tmp;

			/* then process the heap */
			N--;

			my_downheap2 (data1, data2, N, 0);
		}
}

/*
 * FROM statistics/covar_source.c - GSL
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Jim Davies, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

double
pearson (const double data1[], const double data2[],
		const size_t n)
{
	size_t i;
	long double sum_xsq = 0.0;
	long double sum_ysq = 0.0;
	long double sum_cross = 0.0;
	long double ratio;
	long double delta_x, delta_y;
	long double mean_x, mean_y;
	long double r;

	/*
	* Compute:
	* sum_xsq = Sum [ (x_i - mu_x)^2 ],
	* sum_ysq = Sum [ (y_i - mu_y)^2 ] and
	* sum_cross = Sum [ (x_i - mu_x) * (y_i - mu_y) ]
	* using the above relation from Welford's paper
	*/

	mean_x = data1[0];
	mean_y = data2[0];

	for (i = 1; i < n; ++i)
		{
			ratio = i / (i + 1.0);
			delta_x = data1[i] - mean_x;
			delta_y = data2[i] - mean_y;
			sum_xsq += delta_x * delta_x * ratio;
			sum_ysq += delta_y * delta_y * ratio;
			sum_cross += delta_x * delta_y * ratio;
			mean_x += delta_x / (i + 1.0);
			mean_y += delta_y / (i + 1.0);
		}

	r = sum_cross / (sqrt (sum_xsq) * sqrt (sum_ysq));

	return r;
}

static void
compute_rank (double *v, const size_t n)
{
	size_t i = 0;

	while (i < n - 1)
		{
			double vi = v[i];

			if (vi == v[i + 1])
				{
					size_t j = i + 2;
					size_t k;
					double rank = 0.0;

					/* we have detected a tie, find number of equal elements */
					while (j < n && vi == v[j])
						++j;

					/* compute rank */
					for (k = i; k < j; ++k)
						rank += k + 1.0;

					/* divide by number of ties */
					rank /= (double) (j - i);

					for (k = i; k < j; ++k)
						v[k] = rank;

					i = j;
				}
			else
				{
					/* no tie - set rank to natural ordered position */
					v[i] = i + 1.0;
					++i;
				}
		}

	if (i == n - 1)
		v[n - 1] = (double) n;
}

double
spearman (const double data1[], const double data2[],
		const size_t n, double work[])
{
	size_t i;
	double *ranks1 =  &work[0];
	double *ranks2 =  &work[n];

	double r;

	for (i = 0; i < n; ++i)
		{
			ranks1[i] = data1[i];
			ranks2[i] = data2[i];
		}

	/* sort data1 and update data2 at same time; compute rank of data1 */
	sort2 (ranks1, ranks2, n);
	compute_rank (ranks1, n);

	/* now sort data2, updating ranks1 appropriately; compute rank of data2 */
	sort2 (ranks2, ranks1, n);
	compute_rank (ranks2, n);

	/* compute correlation of rank vectors in double precision */
	r = pearson (ranks1, ranks2, n);

	return r;
}

// Fisher–Yates shuffle Algorithm
static void
shuffle (double *base, const size_t n, unsigned int *seed)
{
	size_t i = 0;
	size_t j = 0;

	double temp = 0.0;

	for (i = n - 1; i > 0; i--)
		{
			j = rand_r (seed) % (i + 1);

			temp = base[j];
			base[j] = base[i];
			base[i] = temp;
		}
}

double
spearman_permutation_test (const double data1[], const double data2[],
		const size_t n, double work1[], double work2[], unsigned int *seed,
		const double rho)
{
	const size_t wn = 2 * n;
	size_t i = 0;
	size_t acm = 0;
	double p_value = 0.0;

	for (i = 0; i < n; i++)
		{
			work1[i] = data1[i];
			work1[n + i] = data2[i];
		}

	for (i = 0; i < PERMUTATION_SIZE; i++)
		{
			shuffle (work1, wn, seed);

			double rho2 = spearman (work1, &work1[n], n, work2);

			if (fabs (rho2) < fabs (rho))
				acm++;
		}

	p_value = (double) (PERMUTATION_SIZE - acm) / PERMUTATION_SIZE;
	return p_value;
}

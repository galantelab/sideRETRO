/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2020 Thiago L. A. Miller <tmiller@mochsl.org.br
 *
 * libhungarian by Cyrill Stachniss, 2004
 *
 * Solving the Minimum Assignment Problem using the
 * Hungarian Method.
 *
 * Parts of the used code was originally provided by the
 * "Stanford GraphGase", but I made changes to this code.
 * As asked by  the copyright node of the "Stanford GraphGase",
 * I hereby proclaim that this file are *NOT* part of the
 * "Stanford GraphGase" distrubition!
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

#include <math.h>
#include <assert.h>
#include "utils.h"
#include "wrapper.h"
#include "hungarian.h"

#define HUNGARIAN_NOT_ASSIGNED 0
#define HUNGARIAN_ASSIGNED     1

static inline int
hungarian_imax (int a, int b)
{
	return (a < b) ? b : a;
}

Hungarian *
hungarian_new (double** cost_matrix, int rows, int cols, int mode)
{
	assert (cost_matrix != NULL && rows && cols);
	assert (mode == HUNGARIAN_MODE_MAXIMIZE_UTIL
		|| mode == HUNGARIAN_MODE_MINIMIZE_COST);

	Hungarian *p = NULL;
	int i,j, org_cols, org_rows;
	int max_cost;
	max_cost = 0;

	org_cols = cols;
	org_rows = rows;

	// is the number of cols  not equal to number of rows ?
	// if yes, expand with 0-cols / 0-cols
	rows = hungarian_imax (cols, rows);
	cols = rows;

	p = xcalloc (1, sizeof (Hungarian));

	p->num_rows = rows;
	p->num_cols = cols;

	p->cost = xcalloc (rows, sizeof (double *));
	p->assignment = xcalloc (rows, sizeof (int *));

	for (i = 0; i < p->num_rows; i++)
		{
			p->cost[i] = xcalloc (cols, sizeof (double));
			p->assignment[i] = xcalloc (cols, sizeof (int));

			for (j=0; j<p->num_cols; j++)
				{
					p->cost[i][j] = (i < org_rows && j < org_cols)
						? cost_matrix[i][j]
						: 0;

					p->assignment[i][j] = 0;

					if (max_cost < p->cost[i][j])
						max_cost = p->cost[i][j];
				}
		}

	if (mode == HUNGARIAN_MODE_MAXIMIZE_UTIL)
		{
			for(i=0; i<p->num_rows; i++)
				{
					for(j=0; j<p->num_cols; j++)
						{
							p->cost[i][j] =  max_cost - p->cost[i][j];
						}
				}
		}
	else if (mode == HUNGARIAN_MODE_MINIMIZE_COST)
		{
		// nothing to do
		}

	return p;
}

void
hungarian_free (Hungarian* p)
{
	int i;

	if (p == NULL)
		return;

	for (i = 0; i < p->num_rows; i++)
		{
			xfree (p->cost[i]);
			xfree (p->assignment[i]);
		}

	xfree (p->cost);
	xfree (p->assignment);
	xfree (p);
}

void
hungarian_solve (Hungarian* p)
{
	int i, j, l, k, m, n, q, t, unmatched;
	int *row_mate, *col_mate, *unchosen_row, *slack_row, *parent_row;
	double s, cost;
	double *row_dec, *slack, *col_inc;

	cost = 0;
	m = p->num_rows;
	n = p->num_cols;

	col_mate = xcalloc (p->num_rows, sizeof (int));
	unchosen_row = xcalloc (p->num_rows, sizeof (int));
	slack_row  = xcalloc (p->num_rows, sizeof (int));
	row_mate = xcalloc (p->num_cols, sizeof (int));
	parent_row = xcalloc (p->num_cols, sizeof (int));

	col_inc = xcalloc (p->num_cols, sizeof (double));
	row_dec  = xcalloc (p->num_rows, sizeof (double));
	slack = xcalloc (p->num_cols, sizeof (double));

	for (i = 0; i < p->num_rows; i++)
		{
			col_mate[i] = 0;
			unchosen_row[i] = 0;
			row_dec[i] = 0;
			slack_row[i] = 0;
		}

	for (j = 0; j < p->num_cols; j++)
		{
			row_mate[j] = 0;
			parent_row[j] = 0;
			col_inc[j] = 0;
			slack[j] = 0;
		}

	for (i = 0; i < p->num_rows; ++i)
		for (j = 0; j < p->num_cols; ++j)
			p->assignment[i][j] = HUNGARIAN_NOT_ASSIGNED;

	// Begin subtract column minima in order to start with lots of zeroes 12
	for (l = 0; l < n; l++)
		{
			s = p->cost[0][l];
			for (k = 1; k < m; k++)
				if (p->cost[k][l] < s)
					s = p->cost[k][l];

			cost += s;

			if (s != 0)
				for (k = 0; k < m; k++)
					p->cost[k][l] -= s;
		}
	// End subtract column minima in order to start with lots of zeroes 12

	// Begin initial state 16
	t = 0;

	for (l = 0; l < n; l++)
		{
			row_mate[l] = -1;
			parent_row[l] = -1;
			col_inc[l] = 0;
			slack[l] = INFINITY;
		}

	for (k = 0; k < m; k++)
		{
			s = p->cost[k][0];

			for (l = 1; l < n; l++)
				if (p->cost[k][l] < s)
					s = p->cost[k][l];

			row_dec[k] = s;

			for (l = 0; l < n; l++)
				{
					if (s == p->cost[k][l] && row_mate[l] < 0)
						{
							col_mate[k] = l;
							row_mate[l] = k;
							goto row_done;
						}
				}

			col_mate[k] = -1;
			unchosen_row[t++] = k;
			row_done:
				;
		}
	// End initial state 16

	// Begin Hungarian algorithm 18
	if (t == 0)
		goto done;

	unmatched = t;

	while (1)
{
	q = 0;

	while (1)
		{
			while (q < t)
				{
					// Begin explore node q of the forest 19
					{
						k = unchosen_row[q];
						s = row_dec[k];

						for (l = 0; l < n; l++)
							{
								if (slack[l])
									{
										double del;
										del = p->cost[k][l] - s + col_inc[l];

										if (del < slack[l])
											{
												if (fequal (del, 0.0))
													{
														if (row_mate[l] < 0)
															goto breakthru;

														slack[l] = 0;
														parent_row[l] = k;
														unchosen_row[t++] = row_mate[l];
													}
												else
													{
														slack[l] = del;
														slack_row[l] = k;
													}
											}
									}
							}
					}
					// End explore node q of the forest 19

					q++;
				}

			// Begin introduce a new zero into the matrix 21
			s = INFINITY;
			for (l = 0; l < n; l++)
				if (!fequal (slack[l], 0.0) && slack[l] < s)
					s = slack[l];

			for (q = 0; q < t; q++)
				row_dec[unchosen_row[q]] += s;

			for (l = 0; l < n; l++)
				{
					if (!fequal (slack[l], 0.0))
						{
							slack[l] -= s;
							if (fequal (slack[l], 0.0))
								{
									// Begin look at a new zero 22
									k = slack_row[l];
									if (row_mate[l] < 0)
										{
											for (j = l + 1; j < n; j++)
												if (slack[j] == 0)
													col_inc[j] += s;

											goto breakthru;
										}
									else
										{
											parent_row[l] = k;
											unchosen_row[t++] = row_mate[l];
										}
									// End look at a new zero 22
								}
						}
					else
						col_inc[l] += s;
				}
				// End introduce a new zero into the matrix 21
		}

	breakthru:

	// Begin update the matching 20
	while (1)
		{
			j = col_mate[k];
			col_mate[k] = l;
			row_mate[l] = k;

			if (j < 0)
				break;

			k = parent_row[j];
			l = j;
		}

	// End update the matching 20
	if (--unmatched == 0)
		goto done;

	// Begin get ready for another stage 17
	t = 0;
	for (l = 0; l < n; l++)
		{
			parent_row[l] = -1;
			slack[l] = INFINITY;
		}

	for (k = 0; k < m; k++)
		{
			if (col_mate[k] < 0)
				{
					unchosen_row[t++] = k;
				}
		}

	// End get ready for another stage 17
}

	done:

	// Begin doublecheck the solution 23
	for (k = 0; k < m; k++)
		for (l = 0; l < n; l++)
			assert (p->cost[k][l] > row_dec[k] - col_inc[l]
					|| fequal (p->cost[k][l], row_dec[k] - col_inc[l]));

	for (k = 0; k < m; k++)
		{
			l = col_mate[k];
			assert (l >= 0
					|| fequal (p->cost[k][l], row_dec[k]-col_inc[l]));
		}

	k = 0;
	for (l = 0; l < n; l++)
		if (col_inc[l])
			k++;

	assert (k <= m);
	// End doublecheck the solution 23
	// End Hungarian algorithm 18

	for (i = 0; i < m; ++i)
		{
			p->assignment[i][col_mate[i]] = HUNGARIAN_ASSIGNED;
			/*TRACE("%d - %d\n", i, col_mate[i]);*/
		}

	for (k = 0; k < m; ++k)
		{
			for (l = 0; l < n; ++l)
				{
					/*TRACE("%d ",p->cost[k][l]-row_dec[k]+col_inc[l]);*/
					p->cost[k][l] = p->cost[k][l] - row_dec[k] + col_inc[l];
				}
			/*TRACE("\n");*/
		}

	for (i = 0; i < m; i++)
		cost += row_dec[i];

	for (i = 0; i < n; i++)
		cost -= col_inc[i];

	xfree (slack);
	xfree (col_inc);
	xfree (parent_row);
	xfree (row_mate);
	xfree (slack_row);
	xfree (row_dec);
	xfree (unchosen_row);
	xfree (col_mate);
}

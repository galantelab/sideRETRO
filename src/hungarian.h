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

#ifndef HUNGARIAN_H
#define HUNGARIAN_H

#define HUNGARIAN_MODE_MINIMIZE_COST 0
#define HUNGARIAN_MODE_MAXIMIZE_UTIL 1

struct _Hungarian
{
  int      num_rows;
  int      num_cols;
  double **cost;
  int    **assignment;
};

typedef struct _Hungarian Hungarian;

/** This method initialize the hungarian_problem structure and init
 *  the  cost matrices (missing lines or columns are filled with 0).
 *  It returns the size of the quadratic(!) assignment matrix. **/
Hungarian * hungarian_new (double **cost_matrix, int rows, int cols, int mode);

/** Free the memory allocated by init. **/
void hungarian_free (Hungarian *p);

/** This method computes the optimal assignment. **/
void hungarian_solve (Hungarian *p);

#define hungarian_assignment(p) ((p)->assignment)

#endif /* hungarian.h */

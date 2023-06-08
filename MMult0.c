/* Create macros so that the matrices are stored in column-major order */
#include <math.h>
#define A(i, j) a[(i)*n + (j)]

/* Routine for computing C = A * B + C */

int MY_MMult(double *a, int n, double Tol, int *P)
{
  int i, j, k, imax;
  double maxA, absA;

  for (i = 0; i <= n; i++)
    P[i] = i; // Unit permutation matrix, P[N] initialized with N

  for (i = 0; i < n; i++)
  {
    maxA = 0.0;
    imax = i;

    for (k = i; k < n; k++)
      if ((absA = fabs(A(k, i))) > maxA)
      {
        maxA = absA;
        imax = k;
      }

    if (maxA < Tol)
      return 0; // failure, matrix is degenerate

    if (imax != i)
    {
      // pivoting P
      j = P[i];
      P[i] = P[imax];
      P[imax] = j;

      // Pivoting rows of A
      for (j = 0; j < n; j++)
      {
        double temp = A(i, j);
        A(i, j) = A(imax, j);
        A(imax, j) = temp;
      }

      // Counting pivots starting from N (for determinant)
      P[n]++;
    }

    for (j = i + 1; j < n; j++)
    {
      A(j, i) /= A(i, i);

      for (k = i + 1; k < n; k++)
        A(j, k) -= A(j, i) * A(i, k);
    }
  }

  return 1; // decomposition done
}

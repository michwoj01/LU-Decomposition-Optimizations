#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"
#include <sys/time.h>
#include <time.h>
#include <math.h>

#define A(i, j) a[(i)*n + (j)]
#define B(i, j) b[(i)*n + (j)]
#define abs(x) ((x) < 0.0 ? -(x) : (x))
static double gtod_ref_time_sec = 0.0;

int LUPDecompose(double *a, int n, double Tol, int *P)
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

double check(int n, double *a)
{
  int i, j;
  double temp = 0.0;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      temp = temp + A(i, j);
    }
  }
  return temp;
}

int MY_MMult(double *, int, double, int *);

void random_matrix(int n, double *a)
{
  double drand48();
  int i, j;

  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
      A(i, j) = 2.0 * drand48() - 1.0;
}

void copy_matrix(int n, double *a, double *b)
{
  int i, j;

  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
      B(i, j) = A(i, j);
}

double dclock()
{
  double the_time, norm_sec;
  struct timeval tv;

  gettimeofday(&tv, NULL);

  if (gtod_ref_time_sec == 0.0)
    gtod_ref_time_sec = (double)tv.tv_sec;

  norm_sec = (double)tv.tv_sec - gtod_ref_time_sec;

  the_time = norm_sec + tv.tv_usec * 1.0e-6;

  return the_time;
}

int main()
{
  int p, rep;

  int *b;

  double
      dtime,
      dtime_best,
      gflops,
      original,
      after;

  double
      *a,
      *aold,
      *aref;

  printf("MY_MMult = [\n");

  for (p = PFIRST; p <= PLAST; p += PINC)
  {

    gflops = 2.0 * p * p * p * 1.0e-09;

    /* Allocate space for the matrices */
    /* Note: I create an extra column in A to make sure that
       prefetching beyond the matrix does not cause a segfault */
    a = (double *)malloc(p * p * sizeof(double));
    aold = (double *)malloc(p * p * sizeof(double));
    aref = (double *)malloc(p * p * sizeof(double));
    b = (int *)malloc((p + 1) * sizeof(int));

    /* Generate random matrices Aold  */
    random_matrix(p, aold);

    copy_matrix(p, aold, aref);

    if (LUPDecompose(aref, p, 0.0, b) != 1)
    {
      printf("Error with ref decomposision");
      exit(-1);
    }

    original = check(p, aref);

    /* Run the reference implementation so the answers can be compared */

    /* Time the "optimized" implementation */
    for (rep = 0; rep < NREPEATS; rep++)
    {
      copy_matrix(p, aold, a);

      /* Time your implementation */
      dtime = dclock();

      if (MY_MMult(a, p, 0.0, b) != 1)
      {
        printf("ERROR with new decomposing");
        exit(-1);
      }

      dtime = dclock() - dtime;

      if (rep == 0)
      {
        dtime_best = dtime;
        after = check(p, a);
      }
      else
      {
        dtime_best = (dtime < dtime_best ? dtime : dtime_best);
        after = check(p, a);
      }
    }

    printf("%d %le %le \n", p, gflops / dtime_best, after - original);
    fflush(stdout);

    free(a);
    free(aref);
  }

  printf("];\n");

  exit(0);
}

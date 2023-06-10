#define A(i, j) a[(i) * (n) + (j)]
#define IDX(i, j) ((j) + (i) * (n))
#define myabs(x) ((x) < 0 ? (-x) : (x))
#define MAX(a, b) ((a > b) ? a : b)
#include <immintrin.h>

int MY_MMult(double *a, int n, double Tol, int *P)
{
    register unsigned i, j, k, imax;
    register double maxA, absA, multiplier, divider, temp;
    register __m256d devider;
    register __m256d tmp0, tmp1, tmp4, tmp5;

    for (i = 0; i <= n; i++)
        P[i] = i;

    for (i = 0; i < n; i++)
    {
        maxA = 0.0;
        imax = i;

        for (j = i; j < n; j++)
            if ((absA = myabs(A(j, i))) > maxA)
            {
                maxA = absA;
                imax = j;
            }

        if (maxA < Tol)
            return 0;

        if (imax != i)
        {
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            for (j = 0; j < n; j++)
            {
                temp = A(i, j);
                A(i, j) = A(imax, j);
                A(imax, j) = temp;
            }
            P[n]++;
        }
        multiplier = A(i, i);

        for (j = i + 1; j < n; j++)
        {
            A(j, i) /= multiplier;
            divider = A(j, i);
            devider[0] = divider;
            devider[1] = divider;
            devider[2] = divider;
            devider[3] = divider;

            for (k = i + 1; k < n;)
            {
                if (k < MAX(n - 8, 0))
                {
                    tmp0 = _mm256_loadu_pd(a + IDX(j, k));
                    tmp4 = _mm256_loadu_pd(a + IDX(j, k + 4));

                    tmp1 = _mm256_loadu_pd(a + IDX(i, k));
                    tmp5 = _mm256_loadu_pd(a + IDX(i, k + 4));

                    tmp1 = _mm256_mul_pd(tmp1, devider);
                    tmp5 = _mm256_mul_pd(tmp5, devider);

                    tmp0 = _mm256_sub_pd(tmp0, tmp1);
                    tmp4 = _mm256_sub_pd(tmp4, tmp5);

                    _mm256_storeu_pd(a + IDX(j, k), tmp0);
                    _mm256_storeu_pd(a + IDX(j, k + 4), tmp4);
                    
                    k += 8;
                }
                else
                {
                    A(j, k) -= A(i, k) * divider;
                    k++;
                }
            }
        }
    }

    return 1;
}

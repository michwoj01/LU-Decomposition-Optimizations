#define A(i, j) a[(i)*n + (j)]
#define abs(x) ((x) < 0 ? (-x) : (x))
#define MAX(a, b) ((a > b) ? a : b)

int MY_MMult(double *a, int n, double Tol, int *P)
{
    register int i, j, k, imax;
    register double maxA, absA, multiplier, divider, temp;

    for (i = 0; i <= n; i++)
        P[i] = i;

    for (i = 0; i < n; i++)
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < n; k++)
            if ((absA = abs(A(k, i))) > maxA)
            {
                maxA = absA;
                imax = k;
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

            for (k = i + 1; k < n;)
            {
                if (k < MAX(n - 8, 0))
                {
                    A(j, k) -= A(i, k) * divider;
                    A(j, k + 1) -= A(i, k + 1) * divider;
                    A(j, k + 2) -= A(i, k + 2) * divider;
                    A(j, k + 3) -= A(i, k + 3) * divider;
                    A(j, k + 4) -= A(i, k + 4) * divider;
                    A(j, k + 5) -= A(i, k + 5) * divider;
                    A(j, k + 6) -= A(i, k + 6) * divider;
                    A(j, k + 7) -= A(i, k + 7) * divider;
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

#define A(i, j) a[(i)*n + (j)]
#define abs(x) ((x) < 0 ? (-x) : (x))

int MY_MMult(double *a, int n, double Tol, int *P)
{
    register unsigned i, j, k, imax;
    register double maxA, absA;

    for (i = 0; i <= n; i++)
        P[i] = i;

    for (i = 0; i < n; i++)
    {
        maxA = 0.0;
        imax = i;

        for (j = i; j < n; j++)
            if ((absA = abs(A(j, i))) > maxA)
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
                double temp = A(i, j);
                A(i, j) = A(imax, j);
                A(imax, j) = temp;
            }
            P[n]++;
        }

        for (j = i + 1; j < n; j++)
        {
            A(j, i) /= A(i, i);

            for (k = i + 1; k < n; k++)
                A(j, k) -= A(j, i) * A(i, k);
        }
    }

    return 1;
}

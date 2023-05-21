#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <sstream>
#include <iostream>
#include <iomanip>

double func(double x)
{
    return sin(x + 2) / (0.4 + cos(x));
}

int main(int argc, char** argv)
{
    int rank, size, i;
    long N = 200000000;
    double a_, h, I = 0, Ifinal = 0, eps, time1, time2, x, a = -1, b = 1;

    h = (b - a) / N;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        time1 = MPI_Wtime();

        printf("N = %d, size = %d\n", N, size);
    }

    a_ = a + (b - a) * rank / size;

    for (i = 0; i <= (N / size) - 1; i++)
    {
        x = a_ + i * h;
        I += func(x);
    }

    MPI_Reduce(&I, &Ifinal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        time2 = MPI_Wtime();
        Ifinal = (Ifinal + (func(a) + func(b)) / 2) * h;
        eps = fabs(1.2246281058 - Ifinal);

        printf("\nIntegral = %.10f\n\nError = %.10f\n\nTime = %f\n\n", Ifinal, eps, time2 - time1);
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}

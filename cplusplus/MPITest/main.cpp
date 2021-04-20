#include <iostream>
#include <cmath>
#include "mpi.h"

int main(int argc, char** argv)
{
	int N = 100000000;
	MPI_Init(&argc, &argv);

	int rank;
	int numprocess;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocess);

	char process_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(process_name, &name_len);

	double sum = 0;
	for (int i = rank; i <= N; i += numprocess)
	{
		sum += std::sqrt(1 - (double(i) / N) * (double(i) / N)) / N * 4;
	}
	std::cout << "rank: " << rank << ", name: " << process_name << ", part ressult: " << sum << std::endl;

	double* each_sum = nullptr;
	if (rank == 0)
	{
		each_sum = new double[numprocess];
	}

	MPI_Gather(&sum, 1, MPI_DOUBLE, each_sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		double tot_sum = 0;
		for (int i = 0; i != numprocess; ++i)
		{
			tot_sum += each_sum[i];
		}
		std::cout << "pi: " << tot_sum << std::endl;
	}

	MPI_Finalize();

}

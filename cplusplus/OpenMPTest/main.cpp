#include <iostream>
#include <cmath>
#include <omp.h>
#include <itp/timer>
#include <fmt/format.h>

using namespace std;

int main(int argc, char** argv)
{
	long long sum = 0;
	int N = 20000;
	int thread_num = stoi(argv[1]);
	itp::Timer t;

	t.start();

	for (int i = 0; i < 100; i++)
	{
#pragma omp parallel for reduction(+:sum) num_threads(thread_num)
		for (int j = 0; j < N; j++)
		{
			for (int j = 0; j < N; j++)
			{
				sum++;
			}
		}

	}

	t.stop();

	fmt::print("time: {}, sum is {}\n", t.span(), sum);
}

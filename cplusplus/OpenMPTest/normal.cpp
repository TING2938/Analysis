#include <iostream>
#include <cmath>
#include <itp/timer>
#include <fmt/format.h>

using namespace std;

int main(int argc, char** argv)
{
	long double sum = 0;
	int N = 20000;
	itp::Timer t;

	t.start();

	for (int i = 0; i < 100; i++)
	{
		for (long j = 0; j < N; j++)
		{
			for (long k = 0; k < N; k++)
			{
				sum += j * k - j * j * k;
				sum *= 0.000005;
			}
		}

	}

	t.stop();

	fmt::print("time: {}, sum is {}\n", t.span(), sum);
}

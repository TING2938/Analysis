#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fmt/format.h>
#include <itp/timer>

using namespace std;

int main()
{
	itp::Timer t;
	long long sum = 0;
	int N1 = 1000000;
	int N2 = 100;

	t.start();
	for (int i = 0; i != N1; ++i)
	{
		for (int j = 0; j != N2; ++j)
			sum += (i + j);
	}
	t.stop();
	fmt::print("sum1 is {}, time1 is {}\n", sum, t.span());

	sum = 0;
	t.start();
	for (int i = 0; i != N2; ++i)
	{
		for (int j = 0; j != N1; ++j)
			sum += (i + j);
	}
	t.stop();
	fmt::print("sum2 is {}, time2 is {}\n", sum, t.span());
}

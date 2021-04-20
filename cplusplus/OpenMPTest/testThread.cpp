#include <thread>
#include <vector>
#include <itp/timer>
#include <fmt/format.h>

using namespace std;

void func(long double& sum, int nthread, int rank)
{
	constexpr int N = 20000;

	for (int i = rank; i < 100; i += nthread)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				sum += j * k - j * j * k;
				sum *= 0.000005;
			}
		}
	}
}

int main(int argc, char** argv)
{
	vector<thread> td;

	long double sum = 0;
	int N = 20000;
	int thread_num = stoi(argv[1]);
	itp::Timer t;
	
	t.start();

	for (int i = 0; i < thread_num; i++)
	{
		td.emplace_back(func, std::ref(sum), thread_num, i);
	}

	for (auto&& i : td)
	{
		if (i.joinable())
			i.join();
	}

	t.stop();


	fmt::print("time: {}, sum is {}\n", t.span(), sum);
}

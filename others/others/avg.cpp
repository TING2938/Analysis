#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
using std::string;
using std::vector;

class GetOpt {
public:
	GetOpt(int gc, char** gv) : argc(gc), argv(gv) {}

	template<typename T>
	void operator()(string str, T &x) {
		auto p = find(argv + 1, argv + argc, str);
		if (p != argv + argc) {
			std::stringstream ss{ *(p + 1) };
			ss >> x;
		}
	}

	template <typename T>
	void operator()(string str, vector<T> &x) {
		auto b = find(argv + 1, argv + argc, str);
		if (b != argv + argc) {
			x.clear();
			auto e = std::find_if(b + 1, argv + argc, [](char *chr) { return chr[0] == '-'; });
			T tmp = 0;
			for (size_t i = 1; i != e - b; ++i) {
				std::stringstream ss{ *(b + i) };
				ss >> tmp;
				x.push_back(tmp);
			}
		}
	}

private:
	int argc;
	char **argv;
};
vector<vector<double>> readXVG(string fileName, vector<size_t> &col);

int main(int argc, char *argv[]) {
	if (argc == 1 || argv[1] == string("-h")) {
		std::cout << "Usage: avg <fileName> <-b <begin>> <-e <end>> <-c <column>>" << std::endl;
		return 0;
	}

	GetOpt getopt(argc, argv);
	vector<size_t> col{ 2 };
	getopt("-c", col);
	sort(col.begin(), col.end());

	auto data = readXVG(argv[1], col);

	auto t1 = data[0][0];
	auto t2 = *(data[0].end() - 1);
	auto l = data[0].size();
	auto begTime = t1, endTime = t2;
	getopt("-b", begTime);
	getopt("-e", endTime);

	if (begTime < t1 || begTime > t2)
		begTime = t1;

	if (begTime > endTime || endTime < t1 || endTime > t2)
		endTime = t2;

	auto b = static_cast<size_t>((begTime - t1) * (l - 1) / (t2 - t1));
	auto e = static_cast<size_t>((endTime - t1) * (l - 1) / (t2 - t1));

	vector<double> sum;
	double tmp = 0;
	for (size_t i = 0; i != data[0].size(); ++i) {
		for (size_t j = 1; j != col.size() + 1; ++j)
			tmp += data[j][i];
		sum.push_back(tmp);
		tmp = 0;
	}
	data.push_back(sum);

	vector<double> avg, min, max, Std;

	for (size_t j = 0; j != col.size() + 1; ++j) {
		auto begit = data[j + 1].begin() + b;
		auto endit = data[j + 1].begin() + e + 1;
		avg.push_back(accumulate(begit, endit, 0.0) / (e - b + 1));
		min.push_back(*min_element(begit, endit));
		max.push_back(*max_element(begit, endit));
		double x2 = 0;
		for (auto &&i = begit; i != endit; ++i)
			x2 += *i * *i;
		Std.push_back(sqrt((x2 - (e - b + 1) * avg[j] * avg[j]) / (e - b)));
	}

	printf("\nStatistics from %.2f to %.2f, over %d points\n", data[0][b], data[0][e], e-b+1);
	printf("\n%8s%13s%13s%13s%13s\n", "Column", "Average", "Std", "Max", "Min");
	std::cout << string(13 * 4 + 8, '-') << std::endl;
	size_t i;
	for (i = 0; i != col.size(); ++i)
		printf("%8zd%13.3f%13.3f%13.3f%13.3f\n", col[i], avg[i], Std[i], max[i], min[i]);
	if (col.size() != 1) {
		std::cout << string(13 * 4 + 8, '-') << std::endl;
		printf("%8s%13.3f%13.3f%13.3f%13.3f\n", "Totol", avg[i], Std[i], max[i], min[i]);
	}
	printf("\n");
	//system("pause");

}

vector<vector<double>> readXVG(string fileName, vector<size_t>& col)
{
	std::ifstream file(fileName);
	vector<vector<double>> data(col.size() + 1);
	string line;
	std::stringstream ss;
	double tmp = 0;
	while (getline(file, line)) {
		if (line.size() > 1 && line[0] != '#' && line[0] != '@') {
			ss << line;
			ss >> tmp;
			data[0].push_back(tmp);
			size_t beg = 1, end = col[0], j = 0;
			while (true)
			{
				for (size_t i = 0; i != end - beg; ++i)
					ss >> tmp;
				data[j + 1].push_back(tmp);
				if (++j == col.size())
					break;
				beg = col[j - 1];
				end = col[j];
			}
			ss.clear();
			ss.str("");
		}
	}
	return data;
}

#include <string>
#include <iostream>
#include <itp/getopt>
#include <vector>
#include <sstream>

using namespace std;

int main(int argc, char** argv)
{
	vector<int> vec = {3, 5, 7, 2, 7};
	stringstream ss;
	ss << "[";
	for (auto&& v : vec)
	{
		ss << v << ", ";
	}
	ss << "\b\b]";
	cout << ss.str() << endl;


}
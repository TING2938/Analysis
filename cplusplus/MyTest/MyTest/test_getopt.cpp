#include <itp/getopt>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

int main(int argc, char** argv)
{
	int aa = 12;
	double bb = 3.24;
	std::vector<std::string> vec;
	
	std::string str = "-h";
	itp::Getopt p(str);
	p(aa, "-a", true, "aa");
	p(bb, "-b", true, "bb");
	p.getArray(vec, "-s", false, "vec");
	p.finish();

}

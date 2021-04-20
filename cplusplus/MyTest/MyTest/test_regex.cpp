#include <regex>
#include <itp/core>

int main()
{
	std::regex r;
	r.assign(R"((\d{6})(\d{4})(\d{2})(\d{2})(\d{3}[X\d]))");
	std::string str = "42070319960228293X";
	std::smatch sm;

	auto ad = std::regex_replace(str, r, "$1/$2/$3/$4/$5");
	fmt::print("match: {}\n", ad);


}

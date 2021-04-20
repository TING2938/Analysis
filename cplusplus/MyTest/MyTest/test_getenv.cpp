#include <cstdlib>
#include <iostream>
#include <string>
#include <iomanip>

int main()
{
	int str = std::atoi(std::getenv("TEST2"));
	std::cout << str << std::endl;

	if (!str.empty())
		std::cout << str << std::endl;
	else
		std::cout << "Not Found!" << std::endl;
		
}
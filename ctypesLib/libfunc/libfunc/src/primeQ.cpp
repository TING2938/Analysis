#include "header.h"

unsigned long long primeQ(unsigned long long num) {

	for (size_t i = 2; i <= sqrt(num); ++i) {
		if (num % i == 0) {
			return i;
		} 
	}
	return num;
}

unsigned long long speed(unsigned long long num) {
	auto i = num - 1;
	while (i > 1) {
		if (num % i == 0) 
			return i;
		--i;
	}
	return num;
}
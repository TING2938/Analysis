#include "header.h"

double** SUM(size_t a, size_t b) {
	auto s = CreatMatrix<double>(a, b);
	for (size_t i = 0; i != a; ++i)
		for (size_t j = 0; j != b; ++j)
			s[i][j] = i + j;
	return s;
}
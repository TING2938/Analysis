#include <gtest/gtest.h>
#include <itp/core>

template<typename T>
bool equalQ(const itp::Vector<T>& v1, const itp::Vector<T>& v2)
{
	if (v1.size() != v2.size())
		return false;
	for (size_t i = 0; i != v1.size(); ++i) {
		if (v1[i] != v2[i])
			return false;
	}
	return true;
}

template<typename T>
bool equalQ(const itp::Matrix<T>& m1, const itp::Matrix<T>& m2)
{
	if (m1.size(0) != m2.size(0))
		return false;
	if (m1.size(1) != m2.size(1))
		return false;
	for (size_t i = 0; i != m1.size(0) * m1.size(1); ++i) {
		if (m1(i) != m2(i))
			return false;
	}
	return true;
}


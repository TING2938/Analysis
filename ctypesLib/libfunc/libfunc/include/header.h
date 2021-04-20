#ifndef __HEADER_H__
#define __HEADER_H__
#ifdef _WINDLL
#define DLL_API __declspec(dllexport)
 // __declspec(dllimport)
#else
#define DLL_API 
#endif

#include <iostream>
using namespace std;

template <typename T> T* CreatVector(size_t a) {
	return new T[a]();
}

template <typename T> T** CreatMatrix(size_t a, size_t b) {
	T **s = new T*[a];
	for (size_t i = 0; i != a; ++i)
		s[i] = new T[b]();
	return s;
}

template <typename T> T*** CreatBox(size_t a, size_t b, size_t c) {
	T ***s = new T**[a];
	for (size_t i = 0; i != a; ++i)
		s[i] = CreatMatrix<T>(b, c);
	return s;
}

extern "C"
{
	DLL_API unsigned long long primeQ(unsigned long long num);

	DLL_API double** SUM(size_t a, size_t b);
	DLL_API unsigned long long speed(unsigned long long num);
}


#endif // !__HEADER_H__
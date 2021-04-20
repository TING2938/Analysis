#include <iostream>
using namespace std;

class A
{
public:
	A()
	{
		cout << "Generate class A." << endl;
	}
	~A()
	{
		cout << "Destroy class A" << endl;
	}

	int aaa;
	int bbb;
};

int main()
{
	bool bFirst = true;

	for (int i = 0; i != 5; ++i)
	{

		static A a;

		if (bFirst)
		{
			a.aaa = 34;
			bFirst = false;
		}

		cout << "i is " << i << ", aaa is " << a.aaa << endl;
		
	}
	



}
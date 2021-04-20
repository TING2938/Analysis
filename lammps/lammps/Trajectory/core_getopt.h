#ifndef __GETOPT__
#define __GETOPT__

#include <string>

template <typename ARRAY, typename T>
inline bool contains(const ARRAY& vec, const T& str)
{
	for (auto&& i : vec) {
		if (str == i)
			return true;
	}
	return false;
}

class Getopt
{
public:

	Getopt(int gc, char** gv) : argc(gc), argv(gv) {}

	std::string toString()
	{
		std::string str(argv[0]);
		for (size_t i = 1; i != argc; ++i) {
			str += " ";
			str += argv[i];
		}
		return str;
	}

	template<typename T, typename... S>
	void operator()(T& x, S&& ... str)
	{
		auto p = findArgs(str...);
		if (p != argv + argc) {
			std::stringstream ss{ *(p + 1) };
			ss >> x;
		}
	}

	template<typename... S>
	void operator()(bool& x, S&& ... str)
	{
		auto p = findArgs(str...);
		if (p != argv + argc)
			x = true;
	}

	template <typename T, typename... S>
	void getArray(T& x, S&& ... str)
	{
		auto b = findArgs(str...);
		if (b != argv + argc) {
			x.clear();
			auto e = std::find_if(b + 1, argv + argc,
				[](char* chr) { return chr[0] == '-'; });
			typename T::value_type tmp = 0;
			for (size_t i = 1; i != e - b; ++i) {
				std::stringstream ss{ *(b + i) };
				ss >> tmp;
				x.append(tmp);
			}
		}
	}

protected:

	/*! \brief
	 *  return argument position, or (argv+argc) if not found.
	 */
	template <typename... S>
	char** findArgs(S&& ... str)
	{
		const char* opt[]{ str... };
		return std::find_if(argv + 1, argv + argc, [&opt](std::string&& i) {
			return contains(opt, i);
			});
	}

protected:
	int argc;
	char** argv;
}; // class Getopt;


#endif // !__GETOPT__
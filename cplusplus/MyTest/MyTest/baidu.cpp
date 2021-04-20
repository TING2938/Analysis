#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <fmt/format.h>

using namespace std;

int lengthOfLongestSubstring(string s)
{
	if (s.size() == 1)
		return 1;
	int lng = 0;
	vector<int> veclng;
	unordered_set<int> tmp;

	for (int i = 0; i != s.size(); i++)
	{
		if (!tmp.count(s[i]))
		{
			tmp.insert(s[i]);
			lng++;
		}
		else
		{
			veclng.push_back(lng);
			lng = 0;
			tmp.clear();
			i = s.find_first_not_of(s[i], i) - 1;
			if (i == string::npos - 1)
				break;
			tmp.insert(s[i]);
			lng++;
		}
	}
	if (veclng.empty())
		return lng;
	return max(*max_element(veclng.begin(), veclng.end()), lng);
}

int main()
{
	auto ss = lengthOfLongestSubstring("dvdf");
	fmt::print("{}\n", ss);
}

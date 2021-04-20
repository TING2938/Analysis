//给定一个非负整数数组和一个整数 m，你需要将这个数组分成 m 个非空的连续子数组。设计一个算法使得这 m 个子数组各自和的最大值最小。
//
//注意 :
//数组长度 n 满足以下条件:
//
//1 ≤ n ≤ 1000
//1 ≤ m ≤ min(50, n)
//示例 :
//
//	输入 :
//	nums = [7, 2, 5, 10, 8]
//	m = 2
//
//	输出 :
//	18
//
//	解释 :
//	一共有四种方法将nums分割为2个子数组。
//	其中最好的方式是将其分为[7, 2, 5] 和[10, 8]，
//	因为此时这两个子数组各自的和的最大值为18，在所有情况中最小。
//
//	来源：力扣（LeetCode）
//	链接：https ://leetcode-cn.com/problems/split-array-largest-sum
//著作权归领扣网络所有。商业转载请联系官方授权，非商业转载请注明出处。

#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>

using namespace std;

class Solution
{
public:

	int N;
	map<int, int> mp;

	inline int calc_hash(int m, int beg)
	{
		return ((m & 0x3FF) << 30) | ((beg & 0x3FF) << 20);
	}

	int splitArray(vector<int>& nums, int m)
	{
		N = nums.size();

		return split(m, nums, 0);
	}

	int split(int m, vector<int>& nums, int beg)
	{

		int hash = calc_hash(m, beg);
		auto it = mp.find(hash);
		if (it != mp.end())
		{
			return it->second;
		}

		int sum = 0;
		int ret = INT_MAX;

		if (m == 1)
		{
			ret = accumulate(nums.begin() + beg, nums.end(), 0);
			mp[hash] = ret;
			return ret;
		}

		if (m == nums.size() - beg)
		{
			ret = *max_element(nums.begin() + beg, nums.end());
			mp[hash] = ret;
			return ret;
		}

		for (int end = beg + 1; end != N - 1; ++end)
		{
			sum = accumulate(nums.begin() + beg, nums.begin() + end, 0);
			if (sum >= ret)
			{
				break;
			}
			sum = max(split(m - 1, nums, end), sum);
			ret = min(ret, sum);
		}
		mp[hash] = ret;
		return ret;
	}

};

int main()
{
	vector<int> nums = { 2, 3, 1, 2, 4, 3 };
	int m = 5;
	Solution sol;
	cout << sol.splitArray(nums, m) << endl;
}

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

using namespace std;

class Solution
{
public:
	int splitArray(vector<int>& nums, int m)
	{
		int N = nums.size();
		if (m == N)
		{
			return *max_element(nums.begin(), nums.end());
		}

		if (m == 1)
		{
			return accumulate(nums.begin(), nums.end(), 0);
		}

		vector<int> kn(m - 1);

		vector<int> tmpsum(m);
		int maxnum;


		maxnum = *max_element(nums.begin(), nums.begin() + m - 1);
		maxnum = max(maxnum, accumulate(nums.begin() + m - 1, nums.end(), 0));

		for (kn[0] = 1; kn[0] <= N - m + 1; kn[0]++)
		{
			tmpsum[0] = accumulate(nums.begin(), nums.begin() + kn[0], 0);
			if (tmpsum[0] >= maxnum)
			{
				break;
			}
			loopForKn(1, kn, maxnum, tmpsum, nums, N, m);
		}

		return maxnum;

	}

	void loopForKn(int k, vector<int>& kn, int& maxnum, vector<int>& tmpsum, vector<int>& nums, int N, int m)
	{
		if (k == m - 1)
		{
			tmpsum[m - 1] = accumulate(nums.begin() + kn[m - 2], nums.end(), 0);
			maxnum = *max_element(tmpsum.begin(), tmpsum.end());
			return;
		}

		for (kn[k] = kn[k - 1] + 1; kn[k] <= N - m + k + 1; kn[k]++)
		{
			tmpsum[k] = accumulate(nums.begin() + kn[k - 1], nums.begin() + kn[k], 0);
			if (tmpsum[k] >= maxnum)
			{
				break;
			}
			loopForKn(k + 1, kn, maxnum, tmpsum, nums, N, m);
		}
	}
};

int main()
{
	vector<int> nums = { 5334, 6299, 4199, 9663, 8945, 3566, 9509 };
	int m = 5;
	Solution sol;
	cout << sol.splitArray(nums, m) << endl;
}

//给出两个 非空 的链表用来表示两个非负的整数。其中，它们各自的位数是按照 逆序 的方式存储的，并且它们的每个节点只能存储 一位 数字。
//
//如果，我们将这两个数相加起来，则会返回一个新的链表来表示它们的和。
//
//您可以假设除了数字 0 之外，这两个数都不会以 0 开头。
//
//示例：
//
//输入：(2 -> 4 -> 3) + (5 -> 6 -> 4)
//输出：7 -> 0 -> 8
//原因：342 + 465 = 807
//
//来源：力扣（LeetCode）
//链接：https ://leetcode-cn.com/problems/add-two-numbers
//著作权归领扣网络所有。商业转载请联系官方授权，非商业转载请注明出处。

#include <iostream>

// Definition for singly-linked list.
struct ListNode
{
	int val;
	ListNode* next;
	ListNode(int x) : val(x), next(nullptr)
	{
	}
};

class Solution
{
public:
	ListNode* addTwoNumbers(ListNode* l1, ListNode* l2)
	{
		auto p1 = l1;
		auto p2 = l2;

		int sum, res;

		sum = p1->val + p2->val;
		res = (sum - sum % 10) / 10;
		auto p3 = new ListNode(sum % 10);
		auto l3 = p3;

		if (p1)	p1 = p1->next;
		if (p2) p2 = p2->next;

		while (p1 || p2 || res)
		{
			sum = res;
			if (p1)
			{
				sum += p1->val;
				p1 = p1->next;
			}
			if (p2)
			{
				sum += p2->val;
				p2 = p2->next;
			}

			res = (sum - sum % 10) / 10;

			p3->next = new ListNode(sum % 10);
			p3 = p3->next;
		}
		return l3;
	}
};

int main()
{
	ListNode* l1 = new ListNode(3);
	l1->next = new ListNode(7);

	ListNode* l2 = new ListNode(4);
	l2->next = new ListNode(5);

	Solution sol;
	auto l3 = sol.addTwoNumbers(l1, l2);

	ListNode* tmp = l3;
	while (tmp)
	{
		std::cout << tmp->val << " ";
		tmp = tmp->next;
	}
}

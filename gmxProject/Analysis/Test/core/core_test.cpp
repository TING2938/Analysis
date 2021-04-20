#include "test_config.h"

TEST(ITP_CORE, VECTOR) {

	veci vec1;
	veci vec2{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	veci vec3 = vec2 + 1;

	for (int i = 0; i != vec2.size(); ++i) {
		ASSERT_EQ(vec2[i], i);
		ASSERT_EQ(vec2(i), i);
	}
	ASSERT_EQ(vec2.size(), 10);
	ASSERT_EQ(vec2.capacity(), 10);
	ASSERT_EQ(vec2.max(), 9);
	ASSERT_EQ(vec2.min(), 0);
	ASSERT_EQ(vec2.argmax(), 9);
	ASSERT_EQ(vec2.argmin(), 0);
	ASSERT_DOUBLE_EQ(vec2.mean(), 4.5);
	ASSERT_EQ(vec2.sum(), 45);
	ASSERT_DOUBLE_EQ(vec2.stdev(), 3.0276503540974917);
	ASSERT_EQ(vec2.front(), 0);
	ASSERT_EQ(vec2.back(), 9);
	ASSERT_NE(vec3.size(), 11);

	ASSERT_TRUE(vec2.contains(2));
	ASSERT_FALSE(vec2.contains(12));
	ASSERT_FALSE(vec2.isRef());
	ASSERT_TRUE(equalQ(vec2.flip(), {9, 8, 7, 6, 5, 4, 3, 2, 1, 0}));
	ASSERT_TRUE(equalQ(vec2.part(Rg(1, END)).cumprod(), {1, 2, 6, 24, 120, 720, 5040, 40320, 362880}));
	ASSERT_TRUE(equalQ(vec2.cumsum(), {0, 1, 3, 6, 10, 15, 21, 28, 36, 45}));
	ASSERT_TRUE(equalQ(vec2.diff(), {1, 1, 1, 1, 1, 1, 1, 1, 1}));
	ASSERT_TRUE(equalQ(vec2.part(Rg(2, 4)), {2, 3}));

	veci vec4(vec2.data(), 3);
	ASSERT_TRUE(equalQ(vec4 + 2, {2, 3, 4}));
	ASSERT_TRUE(equalQ(vec4 - 2, {-2, -1, 0}));
	ASSERT_TRUE(equalQ(vec4 * 2, {0, 2, 4}));
	ASSERT_TRUE(equalQ(vec4 / 2, {0, 0, 1}));
	ASSERT_TRUE(equalQ(vec4 += 2, {2, 3, 4}));
	ASSERT_TRUE(equalQ(vec4 -= 2, {0, 1, 2}));
	ASSERT_TRUE(equalQ(vec4 *= 2, {0, 2, 4}));
	ASSERT_TRUE(equalQ(vec4 /= 2, {0, 1, 2}));

	veci vec5(vec2.data() + 2, 3);
	ASSERT_TRUE(equalQ(vec4.setDifference(vec5), {0, 1}));
	ASSERT_TRUE(equalQ(vec4.setUnion(vec5), {0, 1, 2, 3, 4}));
	ASSERT_TRUE(equalQ(vec4.setIntersection(vec5), {2}));
	ASSERT_TRUE(equalQ(vec4 > 1, {false, false, true}));
	ASSERT_TRUE(equalQ(vec4 < 1, {true, false, false}));
	ASSERT_TRUE(equalQ(vec4 >= 1, {false, true, true}));
	ASSERT_TRUE(equalQ(vec4 <= 1, {true, true, false}));
	ASSERT_TRUE(equalQ(vec4 == 1, {false, true, false}));
	ASSERT_TRUE(equalQ(vec4 != 1, {true, false, true}));

	veci vec6{3, 1, 2};
	ASSERT_TRUE(equalQ(vec4 > vec6, {false, false, false}));
	ASSERT_TRUE(equalQ(vec4 >= vec6, {false, true, true}));
	ASSERT_TRUE(equalQ(vec4 < vec6, {true, false, false}));
	ASSERT_TRUE(equalQ(vec4 <= vec6, {true, true, true}));
	ASSERT_TRUE(equalQ(vec4 == vec6, {false, true, true}));
	ASSERT_TRUE(equalQ(vec4 != vec6, {true, false, false}));
	
	ASSERT_TRUE(equalQ(vec4 + vec6, {3, 2, 4}));
	ASSERT_TRUE(equalQ(vec4 - vec6, {-3, 0, 0}));
	ASSERT_TRUE(equalQ(vec4 * vec6, {0, 1, 4}));
	ASSERT_TRUE(equalQ(vec4 / vec6, {0, 1, 1}));
	ASSERT_TRUE(equalQ(vec4 += vec6, {3, 2, 4}));
	ASSERT_TRUE(equalQ(vec4 -= vec6, {0, 1, 2}));
	ASSERT_TRUE(equalQ(vec4 *= vec6, {0, 1, 4}));
	ASSERT_TRUE(equalQ(vec4 /= vec6, {0, 1, 2}));

	vec6.append(4);
	vec6.append(6);
	ASSERT_EQ(vec6.size(), 5);
	ASSERT_EQ(vec6.capacity(), 6);
	ASSERT_EQ(vec6[3], 4);
	ASSERT_EQ(vec6[4], 6);
	vec6.append(4);
	ASSERT_TRUE(equalQ(itp::find(vec6 == 4), {3, 5}));
	ASSERT_TRUE(equalQ(itp::find(vec6 < 4), {0, 1, 2}));
	ASSERT_TRUE(equalQ(vec6.toMatrix(2, 2), {{3, 1},{2, 4}}));
	
	vecb vec7{true, true, false};
	vecb vec8{false, true, false};
	ASSERT_TRUE(equalQ(vec7 && vec8, {false, true, false}));
	ASSERT_TRUE(equalQ(vec7 || vec8, {true, true, false}));
	ASSERT_TRUE(equalQ(!vec8, {true, false, true}));
	
}


TEST(ITP_CORE, MATRIX)
{
	mati m1;
	mati m2(2, 3, {1, 2, 3,
				   4, 5, 6});
	for (size_t i = 0; i != m2.nrow(); ++i)
		for (size_t j = 0; j != m2.ncol(); ++j) {
			ASSERT_EQ(m2[i][j], m2(i*m2.ncol() + j));
		}
	for (size_t i = 0; i != m2.ncol() * m2.nrow(); ++i) {
		ASSERT_EQ(m2(i), i + 1);
	}
	ASSERT_EQ(m2.nrow(), 2);
	ASSERT_EQ(m2.ncol(), 3);
	ASSERT_EQ(m2.capacity(), 6);
	ASSERT_EQ(m2.size(0), 2);
	ASSERT_EQ(m2.size(1), 3);
	ASSERT_FALSE(m2.isRef());
	ASSERT_TRUE(equalQ(m2.max(), {4, 5, 6}));
	ASSERT_TRUE(equalQ(m2.min(), {1, 2, 3}));
	ASSERT_TRUE(equalQ(m2.mean(), {2.5, 3.5, 4.5}));
	ASSERT_TRUE(equalQ(m2.sum(), {5, 7, 9}));
    ASSERT_TRUE(equalQ(m2.stdev(), {2.1213203435596424, 2.1213203435596424, 2.1213203435596424}));
	ASSERT_TRUE(equalQ(m2.max(ROW), {3, 6}));
	ASSERT_TRUE(equalQ(m2.min(ROW), {1, 4}));
	ASSERT_TRUE(equalQ(m2.mean(ROW), {2, 5}));
	ASSERT_TRUE(equalQ(m2.sum(ROW), {6, 15}));
	ASSERT_TRUE(equalQ(m2.stdev(ROW), {1, 1}));
	ASSERT_TRUE(equalQ(m2.part(Rg(1, END), Rg(1, END)), {{5, 6}}));
	ASSERT_TRUE(equalQ(m2.toVector(), {1, 2, 3, 4, 5, 6}));
	ASSERT_TRUE(equalQ(m2.transpose(), {{1, 4}, {2, 5}, {3, 6}}));
	
	mati m3{{1, 2}, {3, 4}};
	ASSERT_TRUE(equalQ(m3 + 1, {{2, 3}, {4, 5}}));
	ASSERT_TRUE(equalQ(m3 - 1, {{0, 1}, {2, 3}}));
	ASSERT_TRUE(equalQ(m3 * 1, {{1, 2}, {3, 4}}));
	ASSERT_TRUE(equalQ(m3 / 1, {{1, 2}, {3, 4}}));
	ASSERT_TRUE(equalQ(m3 += 1, {{2, 3}, {4, 5}}));
	ASSERT_TRUE(equalQ(m3 -= 1, {{1, 2}, {3, 4}}));
	ASSERT_TRUE(equalQ(m3 *= 1, {{1, 2}, {3, 4}}));
	ASSERT_TRUE(equalQ(m3 /= 1, {{1, 2}, {3, 4}}));
	
	ASSERT_TRUE(equalQ(m3 > 2, {{false, false}, {true, true}}));
	ASSERT_TRUE(equalQ(m3 < 2, {{true, false}, {false, false}}));
	ASSERT_TRUE(equalQ(m3 >= 2, {{false, true}, {true, true}}));
	ASSERT_TRUE(equalQ(m3 <= 2, {{true, true}, {false, false}}));
	ASSERT_TRUE(equalQ(m3 == 2, {{false, true}, {false, false}}));
	ASSERT_TRUE(equalQ(m3 != 2, {{true, false}, {true, true}}));

	mati m4{{2, 2}, {1, 5}};
	ASSERT_TRUE(equalQ(m3 + m4, {{3, 4}, {4, 9}}));
	ASSERT_TRUE(equalQ(m3 - m4, {{-1, 0}, {2, -1}}));
	ASSERT_TRUE(equalQ(m3 * m4, {{2, 4}, {3, 20}}));
	ASSERT_TRUE(equalQ(m3 / m4, {{0, 1}, {3, 0}}));
	ASSERT_TRUE(equalQ(m3 += m4, {{3, 4}, {4, 9}}));
	ASSERT_TRUE(equalQ(m3 -= m4, {{1, 2}, {3, 4}}));
	ASSERT_TRUE(equalQ(m3 *= m4, {{2, 4}, {3, 20}}));
	ASSERT_TRUE(equalQ(m3 /= m4, {{1, 2}, {3, 4}}));

	ASSERT_TRUE(equalQ(m3 > m4, {{false, false}, {true, false}})); 
	ASSERT_TRUE(equalQ(m3 < m4, {{true, false}, {false, true}}));
	ASSERT_TRUE(equalQ(m3 >= m4, {{false, true}, {true, false}}));
	ASSERT_TRUE(equalQ(m3 <= m4, {{true, true}, {false, true}}));
	ASSERT_TRUE(equalQ(m3 == m4, {{false, true}, {false, false}}));
	ASSERT_TRUE(equalQ(m3 != m4, {{true, false}, {true, true}}));

	ASSERT_TRUE(equalQ(itp::find(m4 == 2), {0, 1}));
	ASSERT_TRUE(equalQ(itp::find(m4 < 2), {2}));
	ASSERT_TRUE(equalQ(itp::find(m4 > 1), {0, 1, 3}));
	ASSERT_TRUE(equalQ(itp::find(m4 >= 2), {0, 1, 3}));

	ASSERT_TRUE(equalQ(m3.dot(m4), {{4, 12}, {10, 26}}));


}



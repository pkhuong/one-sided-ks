#include "one-sided-ks.h"

#include <cfloat>
#include <climits>
#include <cmath>
#include <tuple>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {
using ::testing::DoubleNear;
using ::testing::Lt;

TEST(OneSidedKs, ConstantsOk)
{
	EXPECT_EQ(one_sided_ks_check_constants(), 0);
}

TEST(OneSidedKs, ThresholdMonotonicGolden)
{
	// The paper says eps = 0.05, m = 6 -> b = 4
	for (size_t i = 6; i < 100; ++i) {
		// The threshold should be f(x) / x,
		// f(x) = sqrt((x + 1)(2 log x + log 4))
		const double fx
		    = std::sqrt((i + 1) * (2 * std::log(i) + std::log(4)));

		EXPECT_THAT(one_sided_ks_threshold(i, 6, std::log(0.05)),
		    DoubleNear(fx / i, 1e-15));
	}
}

// If we increase min count, the threshold should decrease
TEST(OneSidedKs, ThresholdMonotonicMinCount)
{
	EXPECT_THAT(one_sided_ks_threshold(1000, 100, -1),
	    Lt(one_sided_ks_threshold(1000, 10, -1)));
}

// If we allow more false positives, the threshold should decrease.
TEST(OneSidedKs, ThresholdMonotonicEps)
{
	EXPECT_THAT(one_sided_ks_threshold(10000, 1000, -3),
	    Lt(one_sided_ks_threshold(10000, 1000, -4)));
}

// As we get more data points the threshold should also decrease.
TEST(OneSidedKs, ThresholdMonotonicCount)
{
	EXPECT_THAT(one_sided_ks_threshold(100000, 1000, -4),
	    Lt(one_sided_ks_threshold(10000, 1000, -4)));
}

TEST(OneSidedKs, MinCountGolden)
{
	// Paper says 6.
	EXPECT_TRUE(one_sided_ks_min_count_valid(6, std::log(0.05)));

	// Going higher should still be fine.
	EXPECT_TRUE(one_sided_ks_min_count_valid(7, std::log(0.05)));
	EXPECT_TRUE(
	    one_sided_ks_min_count_valid(SIZE_MAX - 1, std::log(0.05)));

	// Lower should not be good
	EXPECT_FALSE(one_sided_ks_min_count_valid(5, std::log(0.05)));
	EXPECT_FALSE(one_sided_ks_min_count_valid(1, std::log(0.05)));
	EXPECT_FALSE(one_sided_ks_min_count_valid(0, std::log(0.05)));
}

TEST(OneSidedKs, FindMinCountGolden)
{
	// Even if we want eps ~= 1, we should always return a valid
	// value.
	//
	// 2 is the bare minimum to avoid singularities.  3 is the
	// actual value here (and for any eps bounded away from 1, I
	// think).
	EXPECT_EQ(one_sided_ks_find_min_count(-DBL_MIN), 3);

	// Paper says 0.05 - 6
	EXPECT_EQ(one_sided_ks_find_min_count(std::log(0.05)), 6);

	// Don't overflow.
	EXPECT_EQ(one_sided_ks_find_min_count(-HUGE_VAL), SIZE_MAX);
}

TEST(OneSidedKs, ExpectedIterEdgeCase)
{
	// Numerical trickery makes this computation extra conservative, but
	// at least it makes sense.
	EXPECT_THAT(one_sided_ks_expected_iter(6, std::log(0.05), 1.0),
	    DoubleNear(100.0, 0.1));
	EXPECT_THAT(one_sided_ks_expected_iter(1000, -1, 0.0),
	    DoubleNear(DBL_MAX, 1.0));
	EXPECT_THAT(one_sided_ks_expected_iter(1000, -1, DBL_MIN),
	    DoubleNear(DBL_MAX, 1.0));
}
} // namespace

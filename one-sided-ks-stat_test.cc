#include "one-sided-ks.h"

#include <cfloat>
#include <climits>
#include <cmath>
#include <iostream>
#include <random>
#include <tuple>
#include <vector>

#include "external/csm/csm.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {
double inv_occurrence(const std::vector<size_t> &x)
{
	size_t sum = 0;
	for (const size_t count : x) {
		sum += count;
	}

	return 1.0 / std::max<size_t>(1, sum);
}

double max_cdf_delta(
    const std::vector<size_t> &x, const std::vector<size_t> &y)
{
	double max_delta = 0.0;
	const double x_scale = inv_occurrence(x);
	const double y_scale = inv_occurrence(y);

	size_t sum_x = 0;
	size_t sum_y = 0;
	for (size_t i = 0; i < std::max(x.size(), y.size()); ++i) {
		if (i < x.size()) {
			sum_x += x[i];
		}

		if (i < y.size()) {
			sum_y += y[i];
		}

		const double delta
		    = std::abs(x_scale * sum_x - y_scale * sum_y);
		max_delta = std::max(delta, max_delta);
	}

	return max_delta;
}

double max_uniform_cdf_delta(const std::vector<size_t> &x)
{
	double max_delta = 0.0;
	const double x_scale = inv_occurrence(x);

	size_t sum_x = 0;
	double sum_y = 0;
	for (size_t i = 0; i < x.size(); ++i) {
		sum_x += x[i];

		sum_y += 1.0 / x.size();

		const double delta = std::abs(x_scale * sum_x - sum_y);
		max_delta = std::max(delta, max_delta);
	}

	return max_delta;
}

// Returns whether the evidence lets us reject the equality hypothesis.
bool uniform_eq_test(
    size_t range, size_t repeat, size_t min_count, double log_eps)
{
	std::random_device dev;
	std::mt19937 rng(dev());

	std::uniform_int_distribution<size_t> dist(0, range - 1);
	std::vector<size_t> x(range, 0);
	std::vector<size_t> y(range, 0);

	for (size_t i = 0; i < repeat; ++i) {
		++x[dist(rng)];
		++y[dist(rng)];

		const double delta = max_cdf_delta(x, y);
		const double threshold = one_sided_ks_pair_threshold_fast(
		    i + 1, min_count, log_eps);
		if (delta > threshold) {
			return true;
		}
	}

	return false;
}

// Same thing, but test directly against the uniform distribution.
bool uniform_distribution_eq_test(
    size_t range, size_t repeat, size_t min_count, double log_eps)
{
	std::random_device dev;
	std::mt19937 rng(dev());

	std::uniform_int_distribution<size_t> dist(0, range - 1);
	std::vector<size_t> x(range, 0);

	for (size_t i = 0; i < repeat; ++i) {
		++x[dist(rng)];

		const double delta = max_uniform_cdf_delta(x);
		const double threshold
		    = one_sided_ks_distribution_threshold_fast(
			i + 1, min_count, log_eps);
		if (delta > threshold) {
			return true;
		}
	}

	return false;
}

// Compare identical uniform distributions for 500K iterations.  We
// should have a false positive rate less than the eps of 0.01.
TEST(OneSidedKs, UniformPair)
{
	size_t total = 0;
	size_t failures = 0;

#ifndef NDEBUG
	std::cout << "This test suite needs a few minutes in optimized mode. "
		     "Debug mode may take for approximately ever."
		  << std::endl;
#endif
	for (size_t i = 0; i < 10000; ++i) {
		++total;
		if (uniform_eq_test(
			10, 500000, 100, std::log(0.01) + one_sided_ks_eq)) {
			++failures;
		}

		if (csm(total, 0.01, failures, std::log(1e-4), nullptr)
		    != 0) {
			std::cout << "Actual rate " << 1.0 * failures / total
				  << ": " << failures << " / " << total
				  << "\n";
			EXPECT_LE(1.0 * failures / total, 0.01)
			    << failures << " / " << total;
			return;
		}
	}

	EXPECT_TRUE(false) << "Too many iterations " << total << "("
			   << failures << ")";
}

TEST(OneSidedKs, UniformDistribution)
{
	size_t total = 0;
	size_t failures = 0;

#ifndef NDEBUG
	std::cout << "This test suite needs a few minutes in optimized mode. "
		     "Debug mode may take for approximately ever."
		  << std::endl;
#endif
	for (size_t i = 0; i < 10000; ++i) {
		++total;
		if (uniform_distribution_eq_test(
			10, 500000, 100, std::log(0.01) + one_sided_ks_eq)) {
			++failures;
		}

		if (csm(total, 0.01, failures, std::log(1e-4), nullptr)
		    != 0) {
			std::cout << "Actual rate " << 1.0 * failures / total
				  << ": " << failures << " / " << total
				  << "\n";
			EXPECT_LE(1.0 * failures / total, 0.01)
			    << failures << " / " << total;
			return;
		}
	}

	EXPECT_TRUE(false) << "Too many iterations " << total << "("
			   << failures << ")";
}

constexpr double kDiscrepancyRate = 0.025;

// Like the EQ test, but differ in kDiscrepancyRate of cases.
std::pair<bool, size_t> uniform_neq_test(
    size_t range, size_t repeat, size_t min_count, double log_eps)
{
	std::random_device dev;
	std::mt19937 rng(dev());

	// First, generate all the sets.
	std::bernoulli_distribution error(kDiscrepancyRate);
	std::uniform_int_distribution<size_t> dist(0, range - 1);
	std::vector<size_t> x(range, 0);
	std::vector<size_t> y(range, 0);

	for (size_t i = 0; i < repeat; ++i) {
		++x[dist(rng)];
		if (error(rng)) {
			++y[range - 1];
		} else {
			++y[dist(rng)];
		}

		const double delta = max_cdf_delta(x, y);
		const double threshold = one_sided_ks_pair_threshold_fast(
		    i + 1, min_count, log_eps);
		if (delta > threshold) {
			return { true, i + 1 };
		}
	}

	return { false, SIZE_MAX };
}

std::pair<bool, size_t> uniform_distribution_neq_test(
    size_t range, size_t repeat, size_t min_count, double log_eps)
{
	std::random_device dev;
	std::mt19937 rng(dev());

	// First, generate all the sets.
	std::bernoulli_distribution error(kDiscrepancyRate);
	std::uniform_int_distribution<size_t> dist(0, range - 1);
	std::vector<size_t> x(range, 0);

	for (size_t i = 0; i < repeat; ++i) {
		if (error(rng)) {
			++x[range - 1];
		} else {
			++x[dist(rng)];
		}

		const double delta = max_uniform_cdf_delta(x);
		const double threshold
		    = one_sided_ks_distribution_threshold_fast(
			i + 1, min_count, log_eps);
		if (delta > threshold) {
			return { true, i + 1 };
		}
	}

	return { false, SIZE_MAX };
}

// Compare a slightly non-uniform and an uniform distribution for 100K
// iterations.  We should have a ridiculous false negative rate.
TEST(OneSidedKs, NonUniformPair)
{
	size_t total = 0;
	size_t successes = 0;
	double total_iter = 0;

	for (size_t i = 0; i < 20000; ++i) {
		++total;

		bool different;
		size_t num_iter;
		std::tie(different, num_iter) = uniform_neq_test(
		    10, 100000, 100, std::log(0.01) + one_sided_ks_eq);
		if (different) {
			total_iter += num_iter;
			++successes;
		}

		if (csm(total, 0.999, successes, std::log(1e-4), nullptr)
		    != 0) {
			std::cout
			    << "Actual rate " << 1.0 * successes / total
			    << " " << successes << " / " << total << " - "
			    << total_iter / std::max<size_t>(1, successes)
			    << "\n";
			EXPECT_GE(1.0 * successes / total, 0.99)
			    << successes << " / " << total;
			return;
		}
	}

	EXPECT_EQ(total, successes)
	    << "Should not have any failure after " << total << " tests ("
	    << total_iter / std::max<size_t>(1, successes) << " iter/test)";
}

TEST(OneSidedKs, NonUniformDistribution)
{
	size_t total = 0;
	size_t successes = 0;
	double total_iter = 0;

	for (size_t i = 0; i < 20000; ++i) {
		++total;

		bool different;
		size_t num_iter;
		std::tie(different, num_iter) = uniform_distribution_neq_test(
		    10, 100000, 100, std::log(0.01) + one_sided_ks_eq);
		if (different) {
			total_iter += num_iter;
			++successes;
		}

		if (csm(total, 0.999, successes, std::log(1e-4), nullptr)
		    != 0) {
			std::cout
			    << "Actual rate " << 1.0 * successes / total
			    << " " << successes << " / " << total << " - "
			    << total_iter / std::max<size_t>(1, successes)
			    << "\n";
			EXPECT_GE(1.0 * successes / total, 0.99)
			    << successes << " / " << total;
			return;
		}
	}

	EXPECT_EQ(total, successes)
	    << "Should not have any failure after " << total << " tests ("
	    << total_iter / std::max<size_t>(1, successes) << " iter/test)";
}

// Compare the same distributions. We should usually do better than
// the expected iteration count.
TEST(OneSidedKs, NonUniformPairExpectedIter)
{
	const double expected_iter = one_sided_ks_expected_iter(
	    100, std::log(0.01) + one_sided_ks_eq, kDiscrepancyRate);

	size_t total = 0;
	size_t successes = 0;
	double total_iter = 0;

	// Given the long-tailed nature of the number of iterations,
	// the median should be lower than the expected value.
	for (size_t i = 0; i < 10000; ++i) {
		++total;

		bool different;
		size_t num_iter;
		std::tie(different, num_iter) = uniform_neq_test(
		    10, 100000, 100, std::log(0.01) + one_sided_ks_eq);
		total_iter += num_iter;
		if (different && num_iter < expected_iter) {
			++successes;
		}

		if (csm(total, 0.5, successes, std::log(1e-4), nullptr)
		    != 0) {
			std::cout << "Average iter " << total_iter / total
				  << " expected " << expected_iter
				  << " hit ratio: " << 1.0 * successes / total
				  << "\n";
			EXPECT_GE(1.0 * successes / total, 0.5)
			    << successes << " / " << total;
			return;
		}
	}

	EXPECT_TRUE(false) << "Too many iterations " << total << "("
			   << successes << " in " << total_iter << ")";
}

} // namespace

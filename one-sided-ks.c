#include "one-sided-ks.h"

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

/* Pairwise <= test is the base case. */
const double one_sided_ks_pair_le = 0;

/* -log 2 rounded away from 0. */
const double one_sided_ks_pair_eq = -0.6931471805599454;

/* -log(2 sqrt(2)) rounded away from 0. */
const double one_sided_ks_fixed_le = -1.039720770839918;

/* -log(4 sqrt(2)) rounded away from 0. */
const double one_sided_ks_fixed_eq = -1.7328679513998635;

/* -log(4 sqrt(2)) rounded away from 0. */
const double one_sided_ks_class = -1.7328679513998635;

static inline uint64_t float_bits(double x)
{
	uint64_t bits;
	uint64_t mask;

	memcpy(&bits, &x, sizeof(bits));
	/* extract the sign bit. */
	mask = (int64_t)bits >> 63;
	/*
	 * If negative, flip the significand bits to convert from
	 * sign-magnitude to 2's complement.
	 */
	return bits ^ (mask >> 1);
}

static inline double bits_float(uint64_t bits)
{
	double ret;
	uint64_t mask;

	mask = (int64_t)bits >> 63;
	/* Undo the bit-flipping above. */
	bits ^= (mask >> 1);
	memcpy(&ret, &bits, sizeof(ret));
	return ret;
}

static inline double next_k(double x, uint64_t delta)
{
	return bits_float(float_bits(x) + delta);
}

__attribute__((__unused__)) static inline double next(double x)
{
	return next_k(x, 1);
}

static inline double prev_k(double x, uint64_t delta)
{
	return bits_float(float_bits(x) - delta);
}

__attribute__((__unused__)) static inline double prev(double x)
{
	return prev_k(x, 1);
}

/* Assume libm is off by < 4 ULPs. */
static const uint64_t libm_error_limit = 4;

static inline double log_up(double x)
{
	return next_k(log(x), libm_error_limit);
}

static inline double log_down(double x)
{
	return prev_k(log(x), libm_error_limit);
}

static inline double sqrt_up(double x)
{
	/* sqrt is supposed to be rounded correctly. */
	return next(sqrt(x));
}

static inline double sqrt_down(double x)
{
	/* sqrt is supposed to be rounded correctly. */
	return prev(sqrt(x));
}

int one_sided_ks_check_constants(void)
{
	int ret = 0;

	/*
	 * Use memcpy instead of float_bits to directly compare bit
	 * patterns in sign-magnitude instead of float_bits's
	 * conversion to 2's complement.
	 */

	size_t index = 0;
#define CHECK(NAME, EXPECTED)                                                \
	do {                                                                 \
		assert(index < CHAR_BIT * sizeof(int) - 1);                  \
		uint64_t actual;                                             \
		memcpy(&actual, &NAME, sizeof(actual));                      \
		if (actual != (uint64_t)EXPECTED) {                          \
			ret |= 1 << index;                                   \
		}                                                            \
		++index;                                                     \
	} while (0)

	/* pair_le is the default: no adjustment. */
	CHECK(one_sided_ks_pair_le, 0);
	CHECK(one_sided_ks_pair_eq, -4618953502541334032ULL);
	CHECK(one_sided_ks_fixed_le, -4616010731606004876ULL);
	CHECK(one_sided_ks_fixed_eq, -4612889074221922196ULL);
	CHECK(one_sided_ks_class, -4612889074221922196ULL);
#undef CHECK

	return ret;
}

/*
 * f(x) / x, where f(x) = ((x + 1)(2 log x + log b))^1/2, rounded up.
 */
static double threshold_up(double x, double log_b_up)
{
	/* Exact up to 2^53. */
	const double xp1 = x + 1;
	/*
	 * compute f(x)^2 = (x + 1)(2 log x + log b).
	 *
	 * x = 1 is exact, and so is the multiplication by 2.
	 */
	const double f_x2 = next(xp1 * next(2 * log_up(x) + log_b_up));

	return next(sqrt_up(f_x2) / x);
}

static double threshold_down(double x, double log_b_down)
{
	/* Exact up to 2^53. */
	const double xp1 = x + 1;
	/*
	 * compute f(x)^2 = (x + 1)(2 log x + log b).
	 *
	 * x = 1 is exact, and so is the multiplication by 2.
	 */
	const double f_x2 = prev(xp1 * prev(2 * log_down(x) + log_b_down));

	return prev(sqrt_down(f_x2) / x);
}

/*
 * b = 1/[eps (min_count - 1)]
 *
 * log(b) = -log(eps) - log(min_count - 1).
 */
static double log_b_up(size_t min_count, double log_eps)
{
	return next(-log_down(min_count - 1) - log_eps);
}

static double log_b_down(size_t min_count, double log_eps)
{
	return prev(-log_up(min_count - 1) - log_eps);
}

double one_sided_ks_threshold(size_t n, size_t min_count, double log_eps)
{
	assert(log_eps < 0);

	if (log_eps >= 0) {
		return -HUGE_VAL;
	}

	if (!one_sided_ks_min_count_valid) {
		min_count = one_sided_ks_find_min_count(log_eps);
	}

	if (n < min_count) {
		return HUGE_VAL;
	}

	return threshold_up(n, log_b_up(min_count, log_eps));
}

/*
 * We need eps exp(min_count - 1) >= min_count + 1.
 *
 * In log space, that's log(eps) + min_count - 1 >= log (min_count + 1)
 */
bool one_sided_ks_min_count_valid(size_t min_count, double log_eps)
{
	assert(log_eps < 0);
	if (log_eps >= 0) {
		return true;
	}

	if (min_count <= 2) {
		return false;
	}

	return prev(log_eps + (min_count - 1)) >= log_up(min_count + 1.0);
}

size_t one_sided_ks_find_min_count(double log_eps)
{
	assert(log_eps < 0);
	if (log_eps >= 0) {
		return 0;
	}

	/*
	 * Galloping search for a min_count that satisfies
	 * one_sided_ks_min_count_valid.
	 */
	size_t i;
	for (i = 1; i < 64; ++i) {
		if (one_sided_ks_min_count_valid(1ULL << i, log_eps)) {
			break;
		}
	}

	/* Nothing to bsearch: m is at least 2. */
	if (i == 1) {
		return 2;
	}

	/* That's practically as good as infinity. */
	if (i >= 64) {
		return SIZE_MAX;
	}

	/*
	 * Binary search for the min count.
	 *
	 * Invariant: low is invalid, and high is valid.
	 */
	size_t low = 1ULL << (i - 1);
	size_t high = 1ULL << i;
	while (low + 1 < high) {
		const size_t pivot = low + (high - low) / 2;
		if (one_sided_ks_min_count_valid(pivot, log_eps)) {
			high = pivot;
		} else {
			low = pivot;
		}
	}

	return high;
}

/*
 * threshold is a monotonically decreasing function of x > 0
 *
 * If rounding up, find the min x s.t. threshold(x, log_b) <= target.
 *
 * If rounding down, find the max x s.t. threshold(x, log_b) >= target
 */
static double invert_threshold(size_t min_count, double target, bool up,
    double threshold(double x, double log_b), double log_b)
{
	if (threshold(min_count, log_b) <= target) {
		return min_count;
	}

	if (threshold(DBL_MAX, log_b) >= target) {
		return DBL_MAX;
	}

	/*
	 * Invariant: threshold(low, log_b) > target
	 *            threshold(high, log_b) < target;
	 */
	double low = min_count;
	double high = DBL_MAX;
	for (size_t i = 0; i < 64; ++i) {
		const uint64_t low_bits = float_bits(low);
		const uint64_t high_bits = float_bits(high);
		const double pivot
		    = bits_float(low_bits + (high_bits - low_bits) / 2);
		const double fx = threshold(pivot, log_b);
		if (fx == target) {
			return pivot;
		}

		if (fx < target) {
			high = pivot;
		} else {
			low = pivot;
		}
	}

	return (up ? high : low);
}

/*
 * Over approximate g(target) = inv[f(x)/x].
 */
static double invert_threshold_up(
    double target, size_t min_count, double log_eps)
{
	return invert_threshold(min_count, target, true, threshold_up,
	    log_b_up(min_count, log_eps));
}

/*
 * Under approximate g(target).
 */
static double invert_threshold_down(
    double target, size_t min_count, double log_eps)
{
	return invert_threshold(min_count, target, false, threshold_down,
	    log_b_down(min_count, log_eps));
}

/*
 * E[N] <= g(d - m / g(d)) <= g_up(d - m / g_down(d))
 */
double one_sided_ks_expected_iter(
    size_t min_count, double log_eps, double delta)
{
	if (min_count == 0 || delta <= 0) {
		return DBL_MAX;
	}

	if (!one_sided_ks_min_count_valid(min_count, log_eps)) {
		return -1;
	}

	assert(log_eps < 0);
	if (log_eps >= 0) {
		log_eps = prev(-0.0);
	}

	/*
	 * The formula doesn't hold if the expected difference exceeds
	 * our first threshold.  Return a conservative count by
	 * clamping `delta` to much less than our first threshold.
	 */
	const double first_threshold
	    = one_sided_ks_threshold(min_count, min_count, log_eps);
	if (delta > first_threshold / 2) {
		delta = prev(first_threshold / 2);
	}

	const double g_delta
	    = invert_threshold_down(delta, min_count, log_eps);
	const double inner = delta - next(1.0 * min_count / g_delta);
	return invert_threshold_up(prev(inner), min_count, log_eps);
}

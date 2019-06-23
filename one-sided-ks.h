#ifndef ONE_SIDED_KS_H
#define ONE_SIDED_KS_H
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Add these constants to `log_eps` dependong on the type of comparison.*/

/* For pairwise <= test (same as not adding anything). */
extern const double one_sided_ks_pair_le;

/* For pairwise equality test. */
extern const double one_sided_ks_pair_eq;

/* <= and equality test with a specific distribution. */
extern const double one_sided_ks_fixed_le;
extern const double one_sided_ks_fixed_eq;

/* Equality test with a family of distributions. */
extern const double one_sided_ks_class;

/*
 * Returns 0 if the constants were correctly compiled, non-zero otherwise.
 *
 * If non-zero the return value is a bitmask with ones for each
 * constant with an incorrect value, in order:
 *
 * bit 0: pair_le
 * bit 1: pair_eq
 * bit 2: fixed_le
 * bit 3: fixed_eq
 * bit 4: class
 *
 */
int one_sided_ks_check_constants(void);

/*
 * Given a sample size of `n` pairs of datapoints, of which the first
 * `min_count` were accumulated without performing any significance
 * test, returns a `threshold` such that, if the supremum of the
 * difference between the two empirical CDFs exceeds `threshold`, we
 * can conclude that the first set of datapoints is not sampled from a
 * distribution that is always less than or equal to the second.
 *
 * The threshold is such that, if we apply this process to an infinite
 * stream of data, the probability of false positive is at most
 * `exp(log_eps)`.
 *
 * log_eps must be negative.
 */
double one_sided_ks_threshold(size_t n, size_t min_count, double log_eps);

/*
 * Determines whether `min_count` is high enough to achieve a log
 * error rate of at most `log_eps`.
 */
bool one_sided_ks_min_count_valid(size_t min_count, double log_eps);

/* Computes the min `min_count` for `log_eps`, which must be negative. */
size_t one_sided_ks_find_min_count(double log_eps);

/*
 * Upper-bounds the expected number of iterations to reject the
 * hypothesis if the actual distance from then null hypothesis is
 * `delta`.
 *
 * Extremely conservative when delta is large.
 *
 * min_count must be valid for log_eps, and log_eps must be negative.
 */
double one_sided_ks_expected_iter(
    size_t min_count, double log_eps, double delta);
#ifdef __cplusplus
} /* extern "C" */
#endif
#endif /* !ONE_SIDED_KS_H */

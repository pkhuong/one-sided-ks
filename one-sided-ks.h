#ifndef ONE_SIDED_KS_H
#define ONE_SIDED_KS_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
/*
 * The functions declared in this header compute one-sided confidence
 * sequences on the Kolmogorov-Smirnov statistic.  Add the following
 * constants to `log_eps` to implement other comparisons.
 *
 * For example, to perform a two-sided test at p < 0.01, pass in
 * `log_eps = ln(0.01) + one_sided_ks_eq`.
 */

/* For the one-sided (<=) test, the default case. */
extern const double one_sided_ks_le;

/* For the two-sided (equality) test. */
extern const double one_sided_ks_eq;

/*
 * For the two-sided one-sample test against a family of distribution,
 * where the distance statistic is the supremum for the distribution
 * that minimises the statistic.
 *
 * See Darling and Robbins's Nonparametric sequential tests with
 * power one (https://www.pnas.org/content/pnas/61/3/804.full.pdf)
 * for technical conditions.
 *
 * Add this constant to `log_eps` and use
 * `one_sided_ks_distribution_threshold`.
 */
extern const double one_sided_ks_class;

/*
 * Returns 0 if the constants were definitely compiled correctly,
 * non-zero otherwise.
 *
 * When non-zero, the return value is a bitmask with ones for each
 * constant with an incorrect value, in order:
 *
 * bit 0: le
 * bit 1: eq
 * bit 2: class
 * bit 3: internal constant (sqrt 1/2).
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
 * `min_count` should be valid for `log_eps`, and `log_eps` must be
 * negative.  When `n < min_count`, this function immediately returns
 * +infty.
 */
double one_sided_ks_pair_threshold(
    uint64_t n, uint64_t min_count, double log_eps);

/*
 * Same as `one_sided_ks_pair_threshold`, without any safety check.
 */
double one_sided_ks_pair_threshold_fast(
    uint64_t n, uint64_t min_count, double log_eps);

/*
 * Given a sample size of `n` datapoints, of which the first
 * `min_count` were accumulated without performing any significance
 * test, returns a `threshold` such that the supremum of the
 * difference between the empirical and theoretical CDFs should only
 * exceed `threshold` with low probability.
 *
 * The threshold is such that, if we apply this process to an infinite
 * stream of data, the probability of exceeding the threshold even
 * once is at most `exp(log_eps)`.
 *
 * `min_count` should be valid for `log_eps`, and `log_eps` must be
 * negative.  When `n < min_count`, this function immediately returns
 * +infty.
 */
double one_sided_ks_distribution_threshold(
    uint64_t n, uint64_t min_count, double log_eps);

/*
 * Same as `one_sided_ks_distribution_threshold`, without any safety check.
 */
double one_sided_ks_distribution_threshold_fast(
    uint64_t n, uint64_t min_count, double log_eps);

/*
 * Determines whether `min_count` is high enough to achieve a log
 * error rate of at most `log_eps`, which must be negative.
 *
 * Returns non-zero if valid, 0 if invalid.
 */
int one_sided_ks_min_count_valid(uint64_t min_count, double log_eps);

/* Computes the min `min_count` for `log_eps`, which must be negative. */
uint64_t one_sided_ks_find_min_count(double log_eps);

/*
 * Upper-bounds the expected number of iterations to reject the
 * hypothesis if the actual distance from then null hypothesis is
 * `delta`.
 *
 * Extremely conservative when delta is large, and for comparisons
 * against a distribution.
 *
 * `min_count` must be valid for `log_eps`, and `log_eps` must be negative.
 */
double one_sided_ks_expected_iter(
    uint64_t min_count, double log_eps, double delta);

#ifdef __cplusplus
} /* extern "C" */
#endif
#endif /* !ONE_SIDED_KS_H */

cc_library(
    name = "one-sided-ks",
    srcs = ["one-sided-ks.c"],
    hdrs = ["one-sided-ks.h"],
    visibility = ["//visibility:public"],
    deps = [],
)

cc_test(
    name = "one-sided-ks_test",
    srcs = ["one-sided-ks_test.cc"],
    deps = [
        ":one-sided-ks",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "one-sided-ks-stat_test",
    srcs = ["one-sided-ks-stat_test.cc"],
    size = "enormous",  # we need a lot of data points
    shard_count = 3,
    deps = [
        ":one-sided-ks",
        "@com_google_googletest//:gtest_main",
        "@csm//:csm",
    ],
)
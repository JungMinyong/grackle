Overview
========

The contents of this directory are used to run the test suite associated with the core Grackle library.
These tests are completely unrelated to the ``pytest`` testing suite.

Zero UV-background rates
=======================

With ``GRACKLE_BUILD_TESTS=ON``, build the ``testUVBackgroundZeros`` target
and run ``ctest --test-dir <build> -R '^UVBackgroundZeros$' --output-on-failure``.
This test exercises the production reader using synthetic UV tables without
external data files. It checks all 14 interpolated quantities with zero,
one-zero and positive endpoint pairs, table endpoints, and background-off
redshifts. Zero intervals must stay exactly zero, one-zero intervals must
remain nonnegative, and positive pairs retain logarithmic interpolation.

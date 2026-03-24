"""Run selected FDTD1D unit tests.

This script executes the following tests together:
- test_fdtd_solves_basic_propagation
- test_fdtd_PMC_boundary_conditions

Usage:
    python run_fdtd_tests.py
"""

import pytest

if __name__ == "__main__":
    # Run only the two specified tests in test_fdtd1d.py
    pytest.main([
        "-q",
        "test_fdtd1d.py",
        "-k",
        "test_fdtd_solves_basic_propagation or test_fdtd_PMC_boundary_conditions",
    ])

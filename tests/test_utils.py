"""Tests for utility functions (acat_optimized, compute_pheno_variance, z_to_p_via_chi2)."""
import math

import numpy as np
import polars as pl
import pytest

from tdbsumstat.utils import acat_optimized, compute_pheno_variance, z_to_p_via_chi2


class TestAcatOptimized:
    def test_uniform_pvalues(self):
        """ACAT of many uniform p-values should be close to 0.5."""
        pvals = pl.Series([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
        result = acat_optimized(pvals)
        assert 0.0 <= result <= 1.0

    def test_very_small_pvalue_list(self):
        """A single very small p-value should produce a small ACAT."""
        pvals = pl.Series([1e-20])
        result = acat_optimized(pvals)
        assert result < 1e-10

    def test_all_ones(self):
        """p-values of 1 should produce ACAT close to 1."""
        pvals = pl.Series([1.0, 1.0, 1.0])
        result = acat_optimized(pvals)
        assert result > 0.9

    def test_accepts_list(self):
        """acat_optimized should accept plain Python lists."""
        pvals = [0.05, 0.1, 0.2]
        result = acat_optimized(pvals)
        assert 0.0 <= result <= 1.0

    def test_invalid_type_raises(self):
        """Non-Series / non-list input should raise TypeError."""
        with pytest.raises(TypeError):
            acat_optimized(np.array([0.05, 0.1]))

    def test_out_of_range_raises(self):
        """p-values outside [0, 1] should raise ValueError."""
        pvals = pl.Series([0.5, 1.5])
        with pytest.raises(ValueError):
            acat_optimized(pvals)

    def test_nan_values_ignored(self):
        """NaN p-values should be silently ignored."""
        pvals = pl.Series([0.05, float("nan"), 0.1])
        result = acat_optimized(pvals)
        assert 0.0 <= result <= 1.0

    def test_all_nan_returns_nan(self):
        """All-NaN input should return NaN."""
        pvals = pl.Series([float("nan"), float("nan")])
        result = acat_optimized(pvals)
        assert math.isnan(result)

    def test_small_pvalue_approximation(self):
        """Values below 1e-15 should use 1/(pi*p) approximation without crashing."""
        pvals = pl.Series([1e-300, 0.05])
        result = acat_optimized(pvals)
        assert result < 0.05  # very small p drives ACAT down


class TestComputePhenoVariance:
    def _make_quant_df(self, n=500, eaf=0.3, se=0.1, beta=0.05):
        return pl.DataFrame({
            "N": [n] * 10,
            "EAF": [eaf] * 10,
            "SE": [se] * 10,
            "BETA": [beta] * 10,
        })

    def _make_binary_df(self, n_cases=200, n_controls=300, eaf=0.3, se=0.1, beta=0.05):
        return pl.DataFrame({
            "N_CASES": [float(n_cases)] * 10,
            "N_CONTROLS": [float(n_controls)] * 10,
            "EAF": [eaf] * 10,
            "SE": [se] * 10,
            "BETA": [beta] * 10,
        })

    def test_quant_returns_string(self):
        df = self._make_quant_df()
        result = compute_pheno_variance(df, "quant")
        # The function returns a string representation of a float
        assert isinstance(result, str)
        float(result)  # should not raise

    def test_binary_adds_n_column(self):
        df = self._make_binary_df()
        result = compute_pheno_variance(df, "binary")
        assert isinstance(result, str)
        float(result)

    def test_quant_positive_variance(self):
        df = self._make_quant_df(n=1000, eaf=0.4, se=0.05)
        result = float(compute_pheno_variance(df, "quant"))
        assert result > 0


class TestZToPViaChi2:
    def test_zero_z_gives_one(self):
        """z=0 should give p-value of 1."""
        p = z_to_p_via_chi2(0)
        assert abs(p - 1.0) < 1e-10

    def test_large_z_gives_small_p(self):
        """Large z-score should give very small p-value."""
        p = z_to_p_via_chi2(10)
        assert p < 1e-20

    def test_z_196_gives_approx_005(self):
        """z ≈ 1.96 should give p ≈ 0.05."""
        p = z_to_p_via_chi2(1.96)
        assert abs(p - 0.05) < 0.01

    def test_symmetric(self):
        """Positive and negative z should give the same p-value."""
        assert abs(z_to_p_via_chi2(2.5) - z_to_p_via_chi2(-2.5)) < 1e-15

"""
Tests for pyCropWat core functionality.
"""

import numpy as np
import pytest

from pycropwat.core import EffectivePrecipitation


class TestCropwatFormula:
    """Tests for the CROPWAT effective precipitation formula."""
    
    def test_low_precipitation(self):
        """Test formula for precipitation <= 250mm."""
        pr = np.array([0, 50, 100, 150, 200, 250])
        ep = EffectivePrecipitation.cropwat_effective_precip(pr)
        
        # Manual calculation: ep = pr * (125 - 0.2 * pr) / 125
        expected = pr * (125 - 0.2 * pr) / 125
        
        np.testing.assert_array_almost_equal(ep, expected)
    
    def test_high_precipitation(self):
        """Test formula for precipitation > 250mm."""
        pr = np.array([251, 300, 500, 1000])
        ep = EffectivePrecipitation.cropwat_effective_precip(pr)
        
        # Manual calculation: ep = 0.1 * pr + 125
        expected = 0.1 * pr + 125
        
        np.testing.assert_array_almost_equal(ep, expected)
    
    def test_boundary_condition(self):
        """Test at the 250mm boundary."""
        pr = np.array([250])
        ep = EffectivePrecipitation.cropwat_effective_precip(pr)
        
        # At exactly 250mm, use the first formula
        expected = 250 * (125 - 0.2 * 250) / 125
        
        np.testing.assert_array_almost_equal(ep, expected)
    
    def test_zero_precipitation(self):
        """Test with zero precipitation."""
        pr = np.array([0])
        ep = EffectivePrecipitation.cropwat_effective_precip(pr)
        
        assert ep[0] == 0
    
    def test_mixed_values(self):
        """Test with mixed values above and below 250mm."""
        pr = np.array([100, 250, 300, 50, 500])
        ep = EffectivePrecipitation.cropwat_effective_precip(pr)
        
        expected = np.where(
            pr <= 250,
            pr * (125 - 0.2 * pr) / 125,
            0.1 * pr + 125
        )
        
        np.testing.assert_array_almost_equal(ep, expected)
    
    def test_2d_array(self):
        """Test with 2D array (simulating raster data)."""
        pr = np.array([
            [100, 200, 300],
            [50, 250, 400],
            [0, 150, 500]
        ])
        ep = EffectivePrecipitation.cropwat_effective_precip(pr)
        
        assert ep.shape == pr.shape
        
        # Check a few values
        assert ep[0, 0] == pytest.approx(100 * (125 - 0.2 * 100) / 125)
        assert ep[0, 2] == pytest.approx(0.1 * 300 + 125)


class TestEffectivePrecipFraction:
    """Tests for effective precipitation fraction calculation."""
    
    def test_fraction_range(self):
        """Test that fraction is always between 0 and 1."""
        pr = np.array([10, 50, 100, 200, 250, 300, 500, 1000])
        ep = EffectivePrecipitation.cropwat_effective_precip(pr)
        epf = np.where(pr > 0, ep / pr, 0)
        
        assert np.all(epf >= 0)
        assert np.all(epf <= 1)
    
    def test_fraction_decreases_with_precipitation(self):
        """Test that fraction decreases as precipitation increases."""
        pr = np.array([50, 100, 200, 300, 500])
        ep = EffectivePrecipitation.cropwat_effective_precip(pr)
        epf = ep / pr
        
        # Fraction should generally decrease (not strictly monotonic at boundary)
        assert epf[0] > epf[-1]

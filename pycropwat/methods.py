"""
Effective precipitation calculation methods.

This module provides various methods for calculating effective precipitation
from total precipitation data.
"""

import numpy as np
from typing import Literal

# Type alias for method names
PeffMethod = Literal[
    "cropwat",
    "fao_aglw", 
    "fixed_percentage",
    "dependable_rainfall"
]


def cropwat_effective_precip(pr: np.ndarray) -> np.ndarray:
    """
    Calculate effective precipitation using the USDA SCS/CROPWAT method.
    
    This is the default method used in FAO CROPWAT software.
    
    Formula:
        - If P <= 250 mm: Peff = P × (125 - 0.2 × P) / 125
        - If P > 250 mm: Peff = 0.1 × P + 125
    
    Parameters
    ----------
    pr : np.ndarray
        Precipitation in mm.
        
    Returns
    -------
    np.ndarray
        Effective precipitation in mm.
        
    References
    ----------
    Smith, M. (1992). CROPWAT: A computer program for irrigation planning
    and management. FAO Irrigation and Drainage Paper No. 46.
    
    Muratoglu, A., et al. (2023). Performance analyses of effective rainfall
    estimation methods. Water Research, 238, 120011.
    """
    ep = np.where(
        pr <= 250,
        pr * (125 - 0.2 * pr) / 125,
        0.1 * pr + 125
    )
    return ep.astype(np.float32)


def fao_aglw_effective_precip(pr: np.ndarray) -> np.ndarray:
    """
    Calculate effective precipitation using the FAO/AGLW formula.
    
    This method is used by FAO's Land and Water Division (AGLW).
    
    Formula:
        - If P <= 250 mm: Peff = 0.6 × P - 10 (but not less than 0)
        - If P > 250 mm: Peff = 0.8 × P - 25
    
    Parameters
    ----------
    pr : np.ndarray
        Precipitation in mm.
        
    Returns
    -------
    np.ndarray
        Effective precipitation in mm.
        
    References
    ----------
    FAO. (1986). Yield response to water. FAO Irrigation and Drainage
    Paper No. 33.
    """
    ep = np.where(
        pr <= 250,
        np.maximum(0.6 * pr - 10, 0),
        0.8 * pr - 25
    )
    return ep.astype(np.float32)


def fixed_percentage_effective_precip(
    pr: np.ndarray,
    percentage: float = 0.7
) -> np.ndarray:
    """
    Calculate effective precipitation using a fixed percentage method.
    
    This simple method assumes a constant fraction of precipitation
    is effective. Common values range from 70-80%.
    
    Formula:
        Peff = P × percentage
    
    Parameters
    ----------
    pr : np.ndarray
        Precipitation in mm.
    percentage : float, optional
        Fraction of precipitation that is effective (0-1).
        Default is 0.7 (70%).
        
    Returns
    -------
    np.ndarray
        Effective precipitation in mm.
    """
    if not 0 <= percentage <= 1:
        raise ValueError(f"Percentage must be between 0 and 1, got {percentage}")
    
    ep = pr * percentage
    return ep.astype(np.float32)


def dependable_rainfall_effective_precip(
    pr: np.ndarray,
    probability: float = 0.75
) -> np.ndarray:
    """
    Calculate effective precipitation using the FAO Dependable Rainfall method.
    
    This method estimates the amount of rainfall that can be depended upon
    at a given probability level. Based on the FAO method which considers
    that effective rainfall at 75% probability is approximately:
    
    Formula (approximation for 75% probability):
        - If P < 100 mm: Peff = 0.6 × P - 10 (but not less than 0)
        - If P >= 100 mm: Peff = 0.8 × P - 25
    
    For other probability levels, a scaling factor is applied.
    
    Parameters
    ----------
    pr : np.ndarray
        Monthly precipitation in mm.
    probability : float, optional
        Probability level (0.5-0.9). Default is 0.75 (75%).
        Higher probability = more conservative estimate.
        
    Returns
    -------
    np.ndarray
        Effective precipitation in mm.
        
    References
    ----------
    FAO. (1992). CROPWAT - A computer program for irrigation planning
    and management. FAO Irrigation and Drainage Paper No. 46.
    
    Notes
    -----
    The scaling factors are approximations based on typical rainfall
    distributions. For more accurate results, site-specific analysis
    of historical rainfall data is recommended.
    """
    if not 0.5 <= probability <= 0.9:
        raise ValueError(f"Probability must be between 0.5 and 0.9, got {probability}")
    
    # Base calculation at 75% probability
    ep_base = np.where(
        pr < 100,
        np.maximum(0.6 * pr - 10, 0),
        0.8 * pr - 25
    )
    
    # Apply probability scaling
    # At 50% probability, multiply by ~1.2
    # At 75% probability, multiply by 1.0 (base case)
    # At 90% probability, multiply by ~0.8
    prob_scale = 1.0 + (0.75 - probability) * 0.8
    
    ep = ep_base * prob_scale
    return np.maximum(ep, 0).astype(np.float32)


def get_method_function(method: PeffMethod):
    """
    Get the effective precipitation function for a given method name.
    
    Parameters
    ----------
    method : str
        Method name: 'cropwat', 'fao_aglw', 'fixed_percentage', 
        or 'dependable_rainfall'.
        
    Returns
    -------
    callable
        The effective precipitation calculation function.
        
    Raises
    ------
    ValueError
        If method name is not recognized.
    """
    methods = {
        'cropwat': cropwat_effective_precip,
        'fao_aglw': fao_aglw_effective_precip,
        'fixed_percentage': fixed_percentage_effective_precip,
        'dependable_rainfall': dependable_rainfall_effective_precip,
    }
    
    if method not in methods:
        raise ValueError(
            f"Unknown method '{method}'. Available methods: {list(methods.keys())}"
        )
    
    return methods[method]


def list_available_methods() -> dict:
    """
    List all available effective precipitation methods with descriptions.
    
    Returns
    -------
    dict
        Dictionary mapping method names to descriptions.
    """
    return {
        'cropwat': 'USDA SCS method as implemented in FAO CROPWAT (default)',
        'fao_aglw': 'FAO/AGLW formula from FAO Irrigation Paper No. 33',
        'fixed_percentage': 'Simple fixed percentage method (default 70%)',
        'dependable_rainfall': 'FAO Dependable Rainfall at specified probability',
    }

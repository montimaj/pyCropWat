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
    "dependable_rainfall",
    "farmwest",
    "usda_scs"
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


def farmwest_effective_precip(pr: np.ndarray) -> np.ndarray:
    """
    Calculate effective precipitation using the FarmWest method.
    
    This is a simple empirical formula used by Washington State University's
    FarmWest program for irrigation scheduling in the Pacific Northwest.
    
    Formula:
        Peff = (P - 5) × 0.75 (but not less than 0)
    
    The method assumes the first 5 mm is lost to interception/evaporation,
    and 75% of the remaining precipitation is effective.
    
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
    FarmWest. Effective Precipitation. Washington State University.
    https://farmwest.com/climate/calculator-information/et/effective-precipitation/
    """
    ep = np.maximum((pr - 5) * 0.75, 0)
    return ep.astype(np.float32)


def usda_scs_effective_precip(
    pr: np.ndarray,
    eto: np.ndarray,
    awc: np.ndarray,
    rooting_depth: float = 1.0
) -> np.ndarray:
    """
    Calculate effective precipitation using the USDA-SCS method with AWC.
    
    This method accounts for soil water holding capacity and evaporative
    demand to estimate effective precipitation. It is based on the USDA
    Soil Conservation Service method that considers soil storage factors.
    
    Formula:
        1. Calculate soil storage depth: d = AWC × 0.5 × rooting_depth (inches)
        2. Calculate storage factor: sf = 0.531747 + 0.295164×d - 0.057697×d² + 0.003804×d³
        3. Calculate effective precipitation:
           Peff = sf × (P^0.82416 × 0.70917 - 0.11556) × 10^(ETo × 0.02426)
        4. Peff is clamped to be between 0 and min(P, ETo)
    
    Note: Internal calculations are done in inches, output is converted to mm.
    
    Parameters
    ----------
    pr : np.ndarray
        Total precipitation in mm.
    eto : np.ndarray
        Reference evapotranspiration in mm.
    awc : np.ndarray
        Available Water Capacity in inches/inch (volumetric fraction).
        For SSURGO data, this is typically in units of cm/cm.
    rooting_depth : float, optional
        Crop rooting depth in meters. Default is 1.0 m.
        
    Returns
    -------
    np.ndarray
        Effective precipitation in mm.
        
    References
    ----------
    USDA SCS. (1993). Chapter 2 Irrigation Water Requirements. In Part 623
    National Engineering Handbook. USDA Soil Conservation Service.
    https://www.wcc.nrcs.usda.gov/ftpref/wntsc/waterMgt/irrigation/NEH15/ch2.pdf
    
    Notes
    -----
    - AWC data for U.S.: projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite
    - AWC data for global: projects/sat-io/open-datasets/FAO/HWSD_V2_SMU (band 'AWC')
    - ETo data for U.S.: projects/openet/assets/reference_et/conus/gridmet/monthly/v1
    - ETo data for global: projects/climate-engine-pro/assets/ce-ag-era5-v2/daily
    """
    # Convert mm to inches for calculation
    pr_inches = pr / 25.4
    eto_inches = eto / 25.4
    
    # Convert rooting depth to inches (1 meter = 39.37 inches)
    rz_inches = rooting_depth * 39.37
    
    # Calculate soil storage depth (d term for eq. 2-85)
    # d = AWC × 0.5 × rooting_depth_inches
    d = awc * 0.5 * rz_inches
    
    # Calculate storage factor (sf) using polynomial equation (eq. 2-85)
    sf = 0.531747 + 0.295164 * d - 0.057697 * np.power(d, 2) + 0.003804 * np.power(d, 3)
    
    # Calculate base effective precipitation term
    # (P^0.82416 × 0.70917 - 0.11556)
    pr_term = np.power(np.maximum(pr_inches, 0.001), 0.82416) * 0.70917 - 0.11556
    pr_term = np.maximum(pr_term, 0)  # Ensure non-negative
    
    # Calculate ETo adjustment term: 10^(ETo × 0.02426)
    eto_term = np.power(10, eto_inches * 0.02426)
    
    # Calculate effective precipitation (eq. 2-84) in inches
    ep_inches = sf * pr_term * eto_term
    
    # Clamp EP: must be <= P and <= ETo, and >= 0
    ep_inches = np.minimum(ep_inches, pr_inches)
    ep_inches = np.minimum(ep_inches, eto_inches)
    ep_inches = np.maximum(ep_inches, 0)
    
    # Convert back to mm
    ep = ep_inches * 25.4
    
    return ep.astype(np.float32)


def get_method_function(method: PeffMethod):
    """
    Get the effective precipitation function for a given method name.
    
    Parameters
    ----------
    method : str
        Method name: 'cropwat', 'fao_aglw', 'fixed_percentage', 
        'dependable_rainfall', 'farmwest', or 'usda_scs'.
        
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
        'farmwest': farmwest_effective_precip,
        'usda_scs': usda_scs_effective_precip,
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
        'farmwest': 'FarmWest method: Peff = (P - 5) × 0.75',
        'usda_scs': 'USDA-SCS method with AWC and ETo (requires GEE assets)',
    }

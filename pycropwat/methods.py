"""
Effective precipitation calculation methods.

This module provides various methods for calculating effective precipitation
from total precipitation data. Effective precipitation represents the portion
of total rainfall that is available for crop use after accounting for losses
due to surface runoff, deep percolation, and evaporation.

Available Methods
-----------------
* ensemble: Ensemble of all methods except TAGEM-SuET (default). Returns the mean of CROPWAT, FAO/AGLW,
Fixed Percentage (70%), Dependable Rainfall (75%), FarmWest, and USDA-SCS methods.
Requires AWC and ETo data.

* cropwat: CROPWAT method - The method used in FAO CROPWAT software.
    
* fao_aglw: FAO Land and Water Division formula from FAO Irrigation Paper No. 33. A two-part linear approximation.
    
* fixed_percentage: Simple fixed percentage method. Assumes a constant fraction (default 70%) of precipitation is 
effective.
    
* dependable_rainfall: FAO Dependable Rainfall method. Estimates rainfall that can be depended
    upon at a given probability level (default 75%).
    
* farmwest: FarmWest method (https://farmwest.com/climate/calculator-information/et/effective-precipitation/). 
Simple empirical formula assuming 5mm interception loss and 75% effectiveness.
    
* usda_scs: USDA-SCS soil moisture depletion method. Accounts for soil water holding capacity (AWC) and 
evaporative demand (ETo). Requires additional GEE assets for AWC and ETo data.

* suet: TAGEM-SuET (Turkish Irrigation Management and Plant Water Consumption System).
Calculates effective precipitation based on the difference between P and ETo. Requires ETo data.

* pcml: Physics-Constrained Machine Learning method (Hasan et al., 2025). Pre-computed Peff
available as a GEE asset for the Western U.S. (17 states) from Jan 2000 to Sep 2024.
Asset: projects/ee-peff-westus-unmasked/assets/effective_precip_monthly_unmasked
Note: Only Western U.S. vectors overlapping the 17-state extent can be used with PCML.

Note
----
For methods requiring ETo (usda_scs, suet, ensemble), crop evapotranspiration (ETc) can be
used instead of grass reference ET (ETo) for more accurate crop-specific estimates.
ETc = ETo × Kc, where Kc is the crop coefficient.

Examples
--------
```python
from pycropwat.methods import (
    cropwat_effective_precip,
    list_available_methods
)
import numpy as np

# Calculate effective precipitation
precip = np.array([50, 100, 200, 300])
eff_precip = cropwat_effective_precip(precip)

# List all available methods
methods = list_available_methods()
for name, description in methods.items():
    print(f"{name}: {description}")
```

References
----------
Smith, M. (1992). CROPWAT: A computer program for irrigation planning
    and management. FAO Irrigation and Drainage Paper No. 46.
    
FAO. (1986). Yield response to water. FAO Irrigation and Drainage
    Paper No. 33.
    
USDA SCS. (1993). Chapter 2 Irrigation Water Requirements. In Part 623
    National Engineering Handbook.

Hasan, M. F., Smith, R. G., Majumdar, S., Huntington, J. L., Alves Meira Neto, A.,
    & Minor, B. A. (2025). Satellite data and physics-constrained machine learning
    for estimating effective precipitation in the Western United States and
    application for monitoring groundwater irrigation. Agricultural Water Management,
    319, 109821. https://doi.org/10.1016/j.agwat.2025.109821
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
    "usda_scs",
    "suet",
    "pcml",
    "ensemble"
]


def cropwat_effective_precip(pr: np.ndarray) -> np.ndarray:
    r"""
    Calculate effective precipitation using the CROPWAT method.
    
    This is the method used in FAO CROPWAT software.
    
    Formula
    -------
    $$
    P_{eff} = \begin{cases}
    P \times \frac{125 - 0.2P}{125} & \text{if } P \leq 250 \text{ mm} \\
    0.1P + 125 & \text{if } P > 250 \text{ mm}
    \end{cases}
    $$
    
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
    r"""
    Calculate effective precipitation using the FAO/AGLW Dependable Rainfall formula.
    
    This method is the FAO Dependable Rainfall method based on 80% probability
    exceedance, used by FAO's Land and Water Division (AGLW). Also known as
    the dependable rainfall method.
    
    Formula
    -------
    $$
    P_{eff} = \begin{cases}
    \max(0.6P - 10, 0) & \text{if } P \leq 70 \text{ mm} \\
    0.8P - 24 & \text{if } P > 70 \text{ mm}
    \end{cases}
    $$
    
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
        pr <= 70,
        np.maximum(0.6 * pr - 10, 0),
        0.8 * pr - 24
    )
    return ep.astype(np.float32)


def fixed_percentage_effective_precip(
    pr: np.ndarray,
    percentage: float = 0.7
) -> np.ndarray:
    r"""
    Calculate effective precipitation using a fixed percentage method.
    
    This simple method assumes a constant fraction of precipitation
    is effective. Common values range from 70-80%.
    
    Formula
    -------
    $$P_{eff} = P \times f$$
    
    where $f$ is the effectiveness fraction (default 0.7).
    
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
    r"""
    Calculate effective precipitation using the FAO Dependable Rainfall method.
    
    This method is based on the FAO/AGLW formula (80% probability exceedance)
    but allows adjustment for different probability levels. The formula estimates
    the amount of rainfall that can be depended upon at a given probability level.
    
    Formula (at 80% probability exceedance)
    ----------------------------------------
    $$
    P_{eff} = \begin{cases}
    \max(0.6P - 10, 0) & \text{if } P \leq 70 \text{ mm} \\
    0.8P - 24 & \text{if } P > 70 \text{ mm}
    \end{cases}
    $$
    
    For other probability levels, a scaling factor is applied.
    
    Parameters
    ----------
    
    pr : np.ndarray
        Monthly precipitation in mm.
    
    probability : float, optional
        Probability level (0.5-0.9). Default is 0.75 (75%).
        Higher probability = lower Peff (plan for more irrigation).
        Lower probability = higher Peff (riskier, less irrigation needed).
        
    Returns
    -------
    np.ndarray
        Effective precipitation in mm.
        
    References
    ----------
    FAO. (1986). Yield response to water. FAO Irrigation and Drainage
    Paper No. 33.
    
    Notes
    -----
    The scaling factors are approximations based on typical rainfall
    distributions. For more accurate results, site-specific analysis
    of historical rainfall data is recommended.
    """
    if not 0.5 <= probability <= 0.9:
        raise ValueError(f"Probability must be between 0.5 and 0.9, got {probability}")
    
    # Base calculation at 80% probability exceedance
    ep_base = np.where(
        pr <= 70,
        np.maximum(0.6 * pr - 10, 0),
        0.8 * pr - 24
    )
    
    # Apply probability scaling
    # At 50% probability, multiply by 1.3
    # At 80% probability, multiply by 1.0 (base case)
    # At 90% probability, multiply by 0.9
    prob_scale = 1.80 - probability
    
    ep = ep_base * prob_scale
    return np.maximum(ep, 0).astype(np.float32)


def farmwest_effective_precip(pr: np.ndarray) -> np.ndarray:
    r"""
    Calculate effective precipitation using the FarmWest method.
    
    This is a simple empirical formula used by the
    FarmWest program for irrigation scheduling in the Pacific Northwest.
    
    Formula
    -------
    $$P_{eff} = \max((P - 5) \times 0.75, 0)$$
    
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
    FarmWest. Effective Precipitation.
    https://farmwest.com/climate/calculator-information/et/effective-precipitation/
    """
    ep = np.maximum((pr - 5) * 0.75, 0)
    return ep.astype(np.float32)


def usda_scs_effective_precip(
    pr: np.ndarray,
    eto: np.ndarray,
    awc: np.ndarray,
    rooting_depth: float = 1.0,
    mad_factor: float = 0.5
) -> np.ndarray:
    r"""
    Calculate effective precipitation using the USDA-SCS method with AWC.
    
    This method accounts for soil water holding capacity and evaporative
    demand to estimate effective precipitation. It is based on the USDA
    Soil Conservation Service method that considers soil storage factors.
    
    Formula
    -------
    1. Calculate soil storage depth: $d = AWC \times MAD \times D_r$ (rooting depth in inches)
    2. Calculate storage factor: $SF = 0.531747 + 0.295164 \cdot d - 0.057697 \cdot d^2 + 0.003804 \cdot d^3$
    3. Calculate effective precipitation:
       $P_{eff} = SF \times (P^{0.82416} \times 0.70917 - 0.11556) \times 10^{ET_o \times 0.02426}$
    4. $P_{eff}$ is clamped to be between 0 and $\min(P, ET_o)$
    
    Note: Internal calculations are done in inches, output is converted to mm.
    
    Parameters
    ----------
    
    pr : np.ndarray
        Total precipitation in mm.
    
    eto : np.ndarray
        Reference evapotranspiration in mm.
    
    awc : np.ndarray
        Available Water Capacity. For SSURGO data (U.S.), this is in inches
        (total for 0-152cm profile). For FAO HWSD data (global), this is in mm/m.
    
    rooting_depth : float, optional
        Crop rooting depth in meters. Default is 1.0 m.
        
    mad_factor : float, optional
        Management Allowed Depletion factor (0-1). Controls what fraction of
        the soil water storage is considered available. Default is 0.5.
        
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
    - AWC data for global: projects/sat-io/open-datasets/FAO/HWSD_V2_SMU (band: 'AWC')
    - ETo data for U.S.: projects/openet/assets/reference_et/conus/gridmet/monthly/v1 (band: 'eto')
    - ETo data for global: projects/climate-engine-pro/assets/ce-ag-era5-v2/daily (band: 'ReferenceET_PenmanMonteith_FAO56')
    """
    # Validate mad_factor
    if not 0 <= mad_factor <= 1:
        raise ValueError(f"mad_factor must be between 0 and 1, got {mad_factor}")
    
    # Convert mm to inches for calculation
    pr_inches = pr / 25.4
    eto_inches = eto / 25.4
    
    # Convert rooting depth to inches (1 meter = 39.37 inches)
    rz_inches = rooting_depth * 39.37
    
    # Calculate soil storage depth (d term for eq. 2-85)
    # d = AWC × MAD × rooting_depth_inches
    d = awc * mad_factor * rz_inches
    
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


def suet_effective_precip(
    pr: np.ndarray,
    eto: np.ndarray
) -> np.ndarray:
    r"""
    Calculate effective precipitation using the TAGEM-SuET method.
    
    TAGEM-SuET (Türkiye'de Sulanan Bitkilerin Bitki Su Tüketimleri) is the
    Turkish Irrigation Management and Plant Water Consumption System.
    This method calculates effective precipitation based on the difference
    between precipitation and reference evapotranspiration. When precipitation
    exceeds ETo, the excess becomes effective precipitation, with a non-linear
    reduction for large excesses.
    
    For more accurate crop-specific estimates, use crop evapotranspiration (ETc)
    instead of grass reference ET (ETo) when available.
    
    Formula
    -------
    $$
    P_{eff} = \begin{cases}
    0 & \text{if } P \leq ET_o \\
    P - ET_o & \text{if } P > ET_o \text{ and } (P - ET_o) < 75 \\
    75 + 0.0011(P - ET_o - 75)^2 + 0.44(P - ET_o - 75) & \text{otherwise}
    \end{cases}
    $$
    
    Parameters
    ----------
    
    pr : np.ndarray
        Total precipitation in mm.
    
    eto : np.ndarray
        Reference evapotranspiration (ETo) in mm. For more accurate crop-specific
        estimates, use crop evapotranspiration (ETc) instead when available.
        
    Returns
    -------
    np.ndarray
        Effective precipitation in mm.
    """
    # Calculate P - ETo
    p_minus_eto = pr - eto
    
    # Case 1: P <= ETo -> Peff = 0
    # Case 2: P > ETo and (P - ETo) < 75 -> Peff = P - ETo
    # Case 3: P > ETo and (P - ETo) >= 75 -> Peff = 75 + 0.0011*(P-ETo-75)^2 + 0.44*(P-ETo-75)
    
    excess = p_minus_eto - 75  # (P - ETo - 75) term
    
    ep = np.where(
        pr <= eto,
        0,
        np.where(
            p_minus_eto < 75,
            p_minus_eto,
            75 + 0.0011 * np.power(excess, 2) + 0.44 * excess
        )
    )
    
    return np.maximum(ep, 0).astype(np.float32)


def pcml_effective_precip(peff: np.ndarray) -> np.ndarray:
    r"""
    Physics-Constrained Machine Learning (PCML) effective precipitation.
    
    This method uses pre-computed effective precipitation from a machine learning
    model trained with physics constraints. The PCML Peff data is available as a
    GEE Image asset for the Western United States (17 states).
    
    Unlike other methods in this module, PCML does not compute Peff from precipitation.
    Instead, it retrieves pre-computed Peff values from the GEE asset. This function
    serves as a pass-through that validates and returns the input PCML Peff data.
    
    GEE Asset
    ---------
    ``projects/ee-peff-westus-unmasked/assets/effective_precip_monthly_unmasked``
    
    Coverage
    --------
    - **Region**: Western United States (17 states: AZ, CA, CO, ID, KS, MT, NE, NV,
      NM, ND, OK, OR, SD, TX, UT, WA, WY)
    - **Temporal**: January 2000 - September 2024 (monthly)
    - **Resolution**: ~2 km (native scale retrieved dynamically from GEE asset)
    - **Band Format**: ``bYYYY_M`` where M is month without leading zero
      (e.g., ``b2015_9`` for September 2015, ``b2016_10`` for October 2016)
    
    Parameters
    ----------
    
    peff : np.ndarray
        Pre-computed effective precipitation from PCML GEE asset in mm.
        This is retrieved from the GEE asset, not calculated.
        
    Returns
    -------
    np.ndarray
        Effective precipitation in mm (pass-through).
        
    References
    ----------
    Hasan, M. F., Smith, R. G., Majumdar, S., Huntington, J. L., Alves Meira Neto, A.,
    & Minor, B. A. (2025). Satellite data and physics-constrained machine learning
    for estimating effective precipitation in the Western United States and
    application for monitoring groundwater irrigation. Agricultural Water Management,
    319, 109821. https://doi.org/10.1016/j.agwat.2025.109821
    
    Examples
    --------
    Using PCML via CLI (no asset/band/geometry required - uses defaults):
    
    ```bash
    pycropwat process --method pcml \\
        --start-year 2000 --end-year 2024 \\
        --output ./WesternUS_PCML --workers 8
    ```
    
    Or with a custom geometry to subset the region:
    
    ```bash
    pycropwat process --method pcml \\
        --geometry pacific_northwest.geojson \\
        --start-year 2000 --end-year 2024 \\
        --output ./PacificNW_PCML --workers 8
    ```
    
    Notes
    -----
    - PCML is trained on Western U.S. data and should only be used in that region.
    - **Geometry requirement**: Only Western U.S. vectors that overlap with the extent of the
      17 states (AZ, CA, CO, ID, KS, MT, NE, NV, NM, ND, OK, OR, SD, TX, UT, WA, WY) can be used.
    - The asset contains monthly effective precipitation values, not raw precipitation.
    - When using ``--method pcml``, the default PCML asset is automatically used.
    - Bands are automatically selected based on year/month (e.g., b2015_9 for Sep 2015).
    - **Only annual (water year, Oct-Sep)** effective precipitation fractions are available for PCML (not monthly),
      loaded directly from a separate GEE asset (band format: ``bYYYY``, e.g., ``b2020``).
    """
    # Pass-through function - PCML Peff is pre-computed in the GEE asset
    return np.asarray(peff).astype(np.float32)


def ensemble_effective_precip(
    pr: np.ndarray,
    eto: np.ndarray,
    awc: np.ndarray,
    rooting_depth: float = 1.0,
    fixed_percentage: float = 0.7,
    dependable_probability: float = 0.75
) -> np.ndarray:
    r"""
    Calculate effective precipitation using ensemble of all methods except TAGEM-SuET.
    
    This method computes the mean of 6 effective precipitation methods:
    CROPWAT, FAO/AGLW, Fixed Percentage, Dependable Rainfall, FarmWest, and USDA-SCS.
    The TAGEM-SuET method is excluded as it tends to underperform in arid/semi-arid
    climates (Muratoglu et al., 2023).
    
    Formula
    -------
    $$P_{eff}^{ensemble} = \frac{1}{6} \sum_{i=1}^{6} P_{eff}^{(i)}$$
    
    where the six methods are: CROPWAT, FAO/AGLW, Fixed %, Dependable Rain, FarmWest, USDA-SCS.
    
    Parameters
    ----------
    
    pr : np.ndarray
        Total precipitation in mm.
    
    eto : np.ndarray
        Reference evapotranspiration (ETo) in mm. For more accurate crop-specific
        estimates, use crop evapotranspiration (ETc) instead when available.
    
    awc : np.ndarray
        Available Water Capacity in inches/inch (volumetric fraction).
        SSURGO asset provides this directly.
    
    rooting_depth : float, optional
        Crop rooting depth in meters. Default is 1.0 m.
        
    fixed_percentage : float, optional
        Fraction for Fixed Percentage method (0-1). Default is 0.7.
        
    dependable_probability : float, optional
        Probability for Dependable Rainfall method (0.5-0.9). Default is 0.75.
        
    Returns
    -------
    np.ndarray
        Ensemble mean effective precipitation in mm.
        
    References
    ----------
    Muratoglu, A., Bilgen, G. K., Angin, I., & Kodal, S. (2023). Performance
    analyses of effective rainfall estimation methods for accurate quantification
    of agricultural water footprint. Water Research, 238, 120011.
    """
    # Calculate Peff using each method
    peff_cropwat = cropwat_effective_precip(pr)
    peff_fao = fao_aglw_effective_precip(pr)
    peff_fixed = fixed_percentage_effective_precip(pr, fixed_percentage)
    peff_dependable = dependable_rainfall_effective_precip(pr, dependable_probability)
    peff_farmwest = farmwest_effective_precip(pr)
    peff_usda = usda_scs_effective_precip(pr, eto, awc, rooting_depth)
    
    # Stack and compute mean
    peff_stack = np.stack([
        peff_cropwat,
        peff_fao,
        peff_fixed,
        peff_dependable,
        peff_farmwest,
        peff_usda
    ], axis=0)
    
    ep = np.nanmean(peff_stack, axis=0)
    
    return ep.astype(np.float32)


def get_method_function(method: PeffMethod):
    """
    Get the effective precipitation function for a given method name.
    
    Parameters
    ----------
    
    method : str
        Method name: 'cropwat', 'fao_aglw', 'fixed_percentage', 
        'dependable_rainfall', 'farmwest', 'usda_scs', 'suet', or 'ensemble'.
        
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
        'suet': suet_effective_precip,
        'pcml': pcml_effective_precip,
        'ensemble': ensemble_effective_precip,
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
        'ensemble': 'Ensemble mean of 6 methods - default (excludes SuET and PCML, requires AWC and ETo)',
        'cropwat': 'CROPWAT method from FAO',
        'fao_aglw': 'FAO/AGLW Dependable Rainfall (80% exceedance)',
        'fixed_percentage': 'Simple fixed percentage method (default 70%)',
        'dependable_rainfall': 'FAO Dependable Rainfall at specified probability',
        'farmwest': r'FarmWest method: $P_{eff} = (P - 5) \times 0.75$',
        'usda_scs': 'USDA-SCS method with AWC and ETo (requires GEE assets)',
        'suet': 'TAGEM-SuET method based on P - ETo (requires ETo asset)',
        'pcml': 'Physics-Constrained ML (Western U.S. 17 states only, Jan 2000 - Sep 2024, GEE asset; only Western U.S. vectors can be used)',
    }

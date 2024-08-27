"""fib MC2020 Chapter 15.10."""

import math
from typing import Literal

import numpy as np


def alpha_T_s() -> float:
    """Returns the coefficient of thermal expansion for reinforcing steel.

    fib Model Code 2020, section 15.5.4

    Returns:
        float: Coefficient of thermal expansion alpha_T_s in K-1.
    """
    return 10e-6


def Es() -> float:
    """Returns the elastic modulus for steel reinforcement.

    fib Model Code 2020, section 15.10.1

    Returns:
        float: Elastic modulus in MPa
    """
    return 200000


def rho_s() -> float:
    """Returns the density of the steel.

    fib Model Code 2020, section 15.4.

    Returns:
        float: the steel density in kg/m3.
    """
    return 7850


def determine_ductility_class(ftk: float, fyk: float, euk: float) -> str:
    """Determine the ductility class of reinforcement steel based on the
    characteristic tensile strength (ftk), yield strength (fyk),
        and strain at maximum force (euk).

    fib Model Code 2020, section 15.5

    Args:
        ftk (float): Characteristic tensile strength in MPa.
        fyk (float): Characteristic yield strength in MPa.
        euk (float): Characteristic strain at maximum stress in percentage (%).

    Returns:
        str: Ductility class ('A', 'B', 'C', or 'D').

    Raises:
        ValueError: If input values are not within valid
            ranges or do not match any ductility class.
    """
    if fyk > 700:
        raise ValueError(
            'Characteristic yield strength fyk must be ≤ 700 MPa'
            + ' for ductility classification.'
        )

    ratio_ft_fy = ftk / fyk

    if ratio_ft_fy >= 1.05 and euk >= 2.5:
        if ratio_ft_fy >= 1.25 and ratio_ft_fy <= 1.45 and euk >= 8:
            return 'D'
        if ratio_ft_fy >= 1.15 and ratio_ft_fy <= 1.35 and euk >= 7.5:
            return 'C'
        if ratio_ft_fy >= 1.08 and euk >= 5:
            return 'B'
        if ratio_ft_fy >= 1.05 and euk >= 2.5:
            return 'A'

    raise ValueError(
        'The provided values do not match any defined ductility class.'
    )


def p_duct_cold_work(eps_u: float, ft: float, fy: float) -> float:
    """Calculate the equivalent steel ductility
        parameter p for cold-worked steel.

    fib Model Code 2020, eq. (15.5-1a)

    Args:
        eps_uu (float): Strain at maximum force (εu) in percentage (%).
        ft (float): Tensile strength (ft) in MPa.
        fy (float): Yield strength (fy) in MPa.

    Returns:
        float: Equivalent steel ductility parameter p for cold-worked steel.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if eps_u <= 0:
        raise ValueError('Strain at maximum force εu must be positive.')
    if ft <= 0 or fy <= 0:
        raise ValueError(
            'Tensile strength ft and yield strength fy must be positive.'
        )
    if fy > 700:
        raise ValueError(
            'Yield strength fy must be ≤ 700 MPa for ductility classification.'
        )

    return eps_u**0.75 * (ft / fy - 1) ** 0.8


def p_duct_hot_rolled(
    eps_u: float, eps_y: float, ft: float, fy: float
) -> float:
    """Calculate the equivalent steel ductility parameter p for hot-rolled steel.

    fib Model Code 2020, eq. (15.5-1b)

    Args:
        eps_u (float): Strain at maximum force (εu) in percentage (%).
        eps_y (float): Yield strain (εy) in percentage (%).
        ft (float): Tensile strength (ft) in MPa.
        fy (float): Yield strength (fy) in MPa.

    Returns:
        float: Equivalent steel ductility parameter p for hot-rolled steel.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if eps_u <= 0 or eps_y <= 0:
        raise ValueError(
            'Strain at maximum force εu and yield strain εy must be positive.'
        )
    if ft <= 0 or fy <= 0:
        raise ValueError(
            'Tensile strength ft and yield strength fy must be positive.'
        )
    if fy > 700:
        raise ValueError(
            'Yield strength fy must be ≤ 700 MPa for ductility classification.'
        )

    return ((eps_u - eps_y) + 3 * eps_y) ** 0.75 * (ft / fy - 1) ** 0.8


def Es_theta(
    T: float,
    steel_type: Literal['hot-rolled', 'cold-worked'],
    Es: float = 200000,
) -> float:
    """Computes the modulus of elasticity at temperature T of
        reinforcing steel.

    fib Model Code 2020, Table (15.10-1)

    Args:
        T (float): Temperature in degrees Celsius (°C).
        steel_type (str): Type of steel ('hot-rolled' or 'cold-worked').
        Es (float). Modulus of elasticity of the steel in MPa.
            Defaults to 200000 MPa.

    Returns:
        float: the modulus of elasticity of reinforcing steel in MPa.

    Raises:
        ValueError: If the temperature is not
        within the valid range or steel_type is invalid.
    """
    if Es < 0:
        raise ValueError('Es cannot be less than 0')
    data = {
        'hot-rolled': {
            20: 1.00,
            100: 1.00,
            200: 0.81,
            300: 0.61,
            400: 0.42,
            500: 0.36,
            600: 0.18,
            700: 0.07,
            800: 0.05,
            900: 0.04,
            1000: 0.02,
            1100: 0.01,
            1200: 0,
        },
        'cold-worked': {
            20: 1.00,
            100: 0.96,
            200: 0.92,
            300: 0.81,
            400: 0.63,
            500: 0.44,
            600: 0.26,
            700: 0.08,
            800: 0.06,
            900: 0.05,
            1000: 0.03,
            1100: 0.02,
            1200: 0,
        },
    }

    if steel_type not in data:
        raise ValueError("Steel type must be 'hot-rolled' or 'cold-worked'.")

    temperatures = np.array(sorted(data[steel_type].keys()))
    values = np.array([data[steel_type][t] for t in temperatures])

    if temperatures[0] > T or temperatures[-1] < T:
        raise ValueError(
            'Temperature must be within the '
            + f'range {temperatures[0]} to {temperatures[-1]} °C.'
        )

    return float(np.interp(T, temperatures, values)) * Es


def fsp_theta(
    T: float, steel_type: Literal['hot-rolled', 'cold-worked'], fyk: float
) -> float:
    """Computes the proportional limit at temperature T of
        reinforcing steel.

    fib Model Code 2020, Table (15.10-1)

    Args:
        T (float): Temperature in degrees Celsius (°C).
        steel_type (str): Type of steel ('hot-rolled' or 'cold-worked').
        fyk (float). Characteristic yield stength in MPa.

    Returns:
        float: the proportional limit at temperature
            of reinforcing steel in MPa.

    Raises:
        ValueError: If the temperature is not
        within the valid range or steel_type is invalid.
    """
    if fyk < 0:
        raise ValueError('fyk cannot be less than 0')
    data = {
        'hot-rolled': {
            20: 1.00,
            100: 1.00,
            200: 0.81,
            300: 0.61,
            400: 0.42,
            500: 0.36,
            600: 0.18,
            700: 0.07,
            800: 0.05,
            900: 0.04,
            1000: 0.02,
            1100: 0.01,
            1200: 0,
        },
        'cold-worked': {
            20: 1.00,
            100: 0.96,
            200: 0.92,
            300: 0.81,
            400: 0.63,
            500: 0.44,
            600: 0.26,
            700: 0.08,
            800: 0.06,
            900: 0.05,
            1000: 0.03,
            1100: 0.02,
            1200: 0,
        },
    }

    if steel_type not in data:
        raise ValueError("Steel type must be 'hot-rolled' or 'cold-worked'.")

    temperatures = np.array(sorted(data[steel_type].keys()))
    values = np.array([data[steel_type][t] for t in temperatures])

    if temperatures[0] > T or temperatures[-1] < T:
        raise ValueError(
            'Temperature must be within the '
            + f'range {temperatures[0]} to {temperatures[-1]} °C.'
        )

    return float(np.interp(T, temperatures, values)) * fyk


def fsy_theta(
    T: float, steel_type: Literal['hot-rolled', 'cold-worked'], fyk: float
) -> float:
    """Computes the maximum stress at temperature T of
        reinforcing steel.

    fib Model Code 2020, Table (15.10-1)

    Args:
        T (float): Temperature in degrees Celsius (°C).
        steel_type (str): Type of steel ('hot-rolled' or 'cold-worked').
        fyk (float). Characteristic yield stength in MPa.

    Returns:
        float: the maximum stress at temperature
            of reinforcing steel in MPa.

    Raises:
        ValueError: If the temperature is not
        within the valid range or steel_type is invalid.
    """
    if fyk < 0:
        raise ValueError('fyk cannot be less than 0')
    data = {
        'hot-rolled': {
            20: 1.00,
            100: 1.00,
            200: 1.00,
            300: 1.00,
            400: 1.00,
            500: 0.78,
            600: 0.47,
            700: 0.23,
            800: 0.11,
            900: 0.06,
            1000: 0.04,
            1100: 0.02,
            1200: 0,
        },
        'cold-worked': {
            20: 1.00,
            100: 1.00,
            200: 1.00,
            300: 1.00,
            400: 0.94,
            500: 0.67,
            600: 0.40,
            700: 0.12,
            800: 0.11,
            900: 0.08,
            1000: 0.05,
            1100: 0.03,
            1200: 0,
        },
    }

    if steel_type not in data:
        raise ValueError("Steel type must be 'hot-rolled' or 'cold-worked'.")

    temperatures = np.array(sorted(data[steel_type].keys()))
    values = np.array([data[steel_type][t] for t in temperatures])

    if temperatures[0] > T or temperatures[-1] < T:
        raise ValueError(
            'Temperature must be within the '
            + f'range {temperatures[0]} to {temperatures[-1]} °C.'
        )

    return float(np.interp(T, temperatures, values)) * fyk


def sigma_theta(
    eps_theta: float, Es_theta: float, fsp_theta: float, fsy_theta: float
) -> float:
    """Computes the stress under certain deformation and temperature.

    fib Model Code 2020, Figg (15.10-2)

    Args:
        eps_theta (float): strain to evaluate the stress.
        Es_theta (float): elastic modulus at temperature in MPa.
        fsp_theta (float): proportional limit at temperature in MPa.
        fsy_theta (float): maximum stress at temperature in MPa.
    """
    if eps_theta < 0 or Es_theta < 0 or fsp_theta < 0 or fsy_theta < 0:
        raise ValueError('All parameters should be larger or equal to zero.')

    eps_sp_theta = fsp_theta / Es_theta
    eps_sy_theta = 0.02
    eps_st_theta = 0.15
    eps_su_theta = 0.20

    # Comptue a,b and c parameters
    c = (fsy_theta - fsp_theta) ** 2 / (
        (eps_sy_theta - eps_sp_theta) * Es_theta - 2 * (fsy_theta - fsp_theta)
    )
    a = math.sqrt(
        (eps_sy_theta - eps_sp_theta)
        * (eps_sy_theta - eps_sp_theta + c / Es_theta)
    )

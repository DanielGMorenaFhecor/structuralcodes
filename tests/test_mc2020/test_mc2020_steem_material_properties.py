"""Tests for fib MC2020 Chapter 15.10."""

import pytest

from structuralcodes.codes.mc2020 import _steel_material_properties


def test_alpha_T_s():
    """Test alpha_T_s."""
    result = _steel_material_properties.alpha_T_s()
    assert result == 10e-6


def test_Es():
    """Test Es."""
    result = _steel_material_properties.Es()
    assert result == 200000


def test_rho_s():
    """Test rho_s."""
    result = _steel_material_properties.rho_s()
    assert result == 7850


@pytest.mark.parametrize(
    'ftk, fyk, euk, expected_class',
    [
        (525, 500, 2.5, 'A'),
        (540, 500, 5.0, 'B'),
        (600, 500, 7.5, 'C'),
        (625, 500, 8.0, 'D'),
    ],
)
def test_determine_ductility_class(ftk, fyk, euk, expected_class):
    """Test determine_ductility_class."""
    result = _steel_material_properties.determine_ductility_class(
        ftk, fyk, euk
    )
    assert result == expected_class


@pytest.mark.parametrize(
    'ftk, fyk, euk',
    [
        (600, 800, 8.0),
        (500, 500, 2.0),
        (500, 500, 8.0),
    ],
)
def test_determine_ductility_class_invalid(ftk, fyk, euk):
    """Test determination of ductility class with invalid input values."""
    with pytest.raises(ValueError):
        _steel_material_properties.determine_ductility_class(ftk, fyk, euk)


@pytest.mark.parametrize(
    'eu, ft, fy, expected',
    [
        (5.0, 600, 500, 0.9226),
    ],
)
def test_p_duct_cold_work(eu, ft, fy, expected):
    """Test p_duct_cold_work."""
    result = _steel_material_properties.p_duct_cold_work(eu, ft, fy)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'eu, ey, ft, fy, expected',
    [
        (10.0, 2.0, 600, 500, 1.9972),
    ],
)
def test_p_duct_hot_rolled(eu, ey, ft, fy, expected):
    """Test p_duct_hot_rolled."""
    result = _steel_material_properties.p_duct_hot_rolled(eu, ey, ft, fy)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'temperature, steel_type, expected',
    [
        (20, 'hot-rolled', 200000),
        (500, 'cold-worked', 88000),
        (450, 'hot-rolled', 78000),
        (1100, 'hot-rolled', 2000),
    ],
)
def test_Es_theta(temperature, steel_type, expected):
    """Test Es_theta."""
    result = _steel_material_properties.Es_theta(temperature, steel_type)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'temperature, steel_type, fyk, expected',
    [
        (20, 'hot-rolled', 500, 500),
        (400, 'cold-worked', 400, 252),
        (350, 'hot-rolled', 500, 257.5),
        (1000, 'cold-worked', 450, 13.5),
    ],
)
def test_fsp_theta(temperature, steel_type, fyk, expected):
    """Test fsp_theta."""
    result = _steel_material_properties.fsp_theta(temperature, steel_type, fyk)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'temperature, steel_type, fyk, expected',
    [
        (20, 'hot-rolled', 500, 500),
        (500, 'hot-rolled', 400, 312),
        (750, 'cold-worked', 500, 57.5),
        (900, 'cold-worked', 450, 36),
    ],
)
def test_fsy_theta(temperature, steel_type, fyk, expected):
    """Test fsy_theta."""
    result = _steel_material_properties.fsy_theta(temperature, steel_type, fyk)
    assert result == pytest.approx(expected, rel=1e-2)

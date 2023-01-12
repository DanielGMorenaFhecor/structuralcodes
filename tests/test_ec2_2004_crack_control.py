"""Tests for EUROCODE 2-1-1:2004 Chapter 7.3 Crack Control"""
import math

import pytest
from structuralcodes.codes.ec2_2004 import _crack_control


@pytest.mark.parametrize(
    'test_exposure_class, test_load_combination, expected',
    [
        ('X0', 'f', 0.2),
        ('x0', 'F', 0.2),
        ('X0', 'qp', 0.4),
        ('x0', 'QP', 0.4),
        ('XC2', 'f', 0.2),
        ('xc2', 'F', 0.2),
        ('XC3', 'f', 0.2),
        ('xc3', 'F', 0.2),
        ('XC4', 'f', 0.2),
        ('xc4', 'F', 0.2),
        ('XC2', 'qp', 0.3),
        ('xc2', 'QP', 0.3),
        ('XC3', 'qp', 0.3),
        ('xc3', 'QP', 0.3),
        ('XC4', 'qp', 0.3),
        ('xc4', 'QP', 0.3),
        ('XD1', 'qp', 0.3),
        ('xd1', 'QP', 0.3),
        ('XD2', 'qp', 0.3),
        ('xd2', 'QP', 0.3),
        ('XS1', 'qp', 0.3),
        ('xs1', 'QP', 0.3),
        ('XS2', 'qp', 0.3),
        ('xs2', 'QP', 0.3),
        ('XS3', 'qp', 0.3),
        ('xs3', 'QP', 0.3),
    ],
)
def test_w_max_returns_expected_values(
    test_exposure_class, test_load_combination, expected
):
    """Test that the w_max function returns expected values"""
    w_max = _crack_control.w_max(test_exposure_class, test_load_combination)
    assert w_max == expected


@pytest.mark.parametrize(
    'test_exposure_class, test_load_combination',
    [('dummy1', 'f'), ('dummy2', 'qp'), ('XD1', 'dummy3'), ('XS1', 'dummy4')],
)
def test_w_max_not_valid_input_raises_valueerror(
    test_exposure_class, test_load_combination
):
    """Test that not valid input returns ValueError"""
    with pytest.raises(ValueError):
        _crack_control.w_max(test_exposure_class, test_load_combination)


@pytest.mark.parametrize(
    'h, expected',
    [
        (200, 1),
        (300, 1),
        (800, 0.65),
        (1000, 0.65),
        (400, 0.93),
        (500, 0.86),
        (600, 0.79),
        (700, 0.72),
    ],
)
def test_k_crack_min_steel_area_returns_expected_values(h, expected):
    """Test the k_crack_min_steel_area function"""
    k = _crack_control.k_crack_min_steel_area(h)
    assert math.isclose(k, expected)


def test_k_crack_min_steel_area_raises_valueerror():
    """Test that not valid input returns ValueError exeption"""
    with pytest.raises(ValueError):
        h = -100
        _crack_control.k_crack_min_steel_area(h)


def test_kc_crack_min_steel_area_pure_tension_returns_expected_values():
    """Test the kc_crack_min_steel_area_pure_tension function"""
    assert 1 == _crack_control.kc_crack_min_steel_area_pure_tension()


@pytest.mark.parametrize(
    'h, b, fct_eff, n_ed, expected',
    [
        (600, 400, 3, 20, 0.3925926),
        (600, 400, 3, -20, 0.4166667),
        (400, 200, 4, 3, 0.397500),
        (200, 50, 5, -80, 1),
        (200, 50, 5, 80, 0),
    ],
)
def test_kc_crack_min_steel_area_rectangular_returns_expected_values(
    h, b, fct_eff, n_ed, expected
):
    """Test the kc_crack_min_steel_area_rectangular"""
    kc = _crack_control.kc_crack_min_steel_area_rectangular(
        h,
        b,
        fct_eff,
        n_ed,
    )
    assert math.isclose(kc, expected, rel_tol=0.000001)


def test_kc_crack_min_steel_area_rectangular_raises_valueerror():
    """Test the kc_crack_min_steel_area_rectangular raises Value
    Error for not correct input values for b and h"""
    with pytest.raises(ValueError):
        _crack_control.kc_crack_min_steel_area_rectangular(
            h=-100, b=100, fct_eff=100, n_ed=10
        )
        _crack_control.kc_crack_min_steel_area_rectangular(
            h=100, b=-100, fct_eff=100, n_ed=10
        )


@pytest.mark.parametrize(
    'f_cr, a_ct, fct_eff, expected',
    [
        (30, 10000, 5, 0.54),
        (20, 5000, 3, 1.2),
        (55, 7500, 4, 1.65),
        (55, 50000, 4, 0.5),
    ],
)
def test_kc_crack_min_steel_area_flanges(f_cr, a_ct, fct_eff, expected):
    """Test the kc_crack_min_steel_area_flanges function"""
    kc = _crack_control.kc_crack_min_steel_area_flanges(f_cr, a_ct, fct_eff)
    assert math.isclose(kc, expected, rel_tol=0.000001)


@pytest.mark.parametrize(
    'a_ct, s_steel, fct_eff, k, kc, expected',
    [
        (10000, 500, 3, 1, 1, 60),
        (80000, 500, 5, 0.65, 0.5, 260),
        (80000, 400, 4, 0.9, 0.75, 540),
    ],
)
def test_crack_min_steel_area_returns_expected_values(
    a_ct, s_steel, fct_eff, k, kc, expected
):
    """Test the crack_min_steel_area returns expected values"""
    as_min = _crack_control.crack_min_steel_area(a_ct, s_steel, fct_eff, k, kc)
    assert math.isclose(as_min, expected, rel_tol=10e-6)


@pytest.mark.parametrize(
    'a_ct, s_steel, fct_eff, k, kc',
    [
        (-10000, 100, 3, 0.7, 0.67),
        (10000, -100, 3, 0.7, 0.65),
        (10000, 100, 3, 0.5, 0.65),
        (10000, 100, 3, 1.1, 0.65),
        (10000, 100, 3, 0.7, -0.1),
        (10000, 100, 3, 0.7, 1.1),
    ],
)
def test_crack_min_steel_area_raises_valueerror(a_ct, s_steel, fct_eff, k, kc):
    """Test the crack_min_steel_area raises value error"""
    with pytest.raises(ValueError):
        _crack_control.crack_min_steel_area(a_ct, s_steel, fct_eff, k, kc)


@pytest.mark.parametrize(
    (
        'a_ct, s_steel, fct_eff, k, kc, ap, d_steel, d_press, e,    '
        ' incr_stress, expected'
    ),
    [
        (80000, 400, 4, 0.9, 0.75, 500, 10, 10, 0.5, 10, 531.161),
        (50000, 500, 3, 0.7, 0.4, 700, 10, 30, 0.8, 20, 69.541),
        (50000, 500, 4, 1, 1, 1000, 0, 20, 0.8, 20, 364.223),
    ],
)
def test_crack_min_steel_area_with_press_tendons_returns_expected_values(
    a_ct,
    s_steel,
    fct_eff,
    k,
    kc,
    ap,
    d_steel,
    d_press,
    e,
    incr_stress,
    expected,
):
    """Test the crack_min_steel_area returns expected values"""
    as_min = _crack_control.crack_min_steel_area_with_prestresed_tendons(
        a_ct, s_steel, fct_eff, k, kc, ap, d_steel, d_press, e, incr_stress
    )
    assert math.isclose(as_min, expected, rel_tol=10e-6)


@pytest.mark.parametrize(
    'a_ct, s_steel, fct_eff, k, kc, ap, d_steel, d_press, e,     incr_stress',
    [
        (-80000, 400, 4, 0.9, 0.75, 500, 10, 10, 0.5, 10),
        (80000, -400, 4, 0.9, 0.75, 500, 10, 10, 0.5, 10),
        (80000, 400, 4, 0.5, 0.75, 500, 10, 10, 0.5, 10),
        (80000, 400, 4, 1.1, 0.75, 500, 10, 10, 0.5, 10),
        (80000, 400, 4, 0.9, -0.1, 500, 10, 10, 0.5, 10),
        (80000, 400, 4, 0.9, 1.1, 500, 10, 10, 0.5, 10),
        (80000, 400, 4, 0.9, 0.75, -500, 10, 10, 0.5, 10),
        (80000, 400, 4, 0.9, 0.75, 500, -10, 10, 0.5, 10),
        (80000, 400, 4, 0.9, 0.75, 500, 10, 0, 0.5, 10),
        (80000, 400, 4, 0.9, 0.75, 500, 10, 10, 0.1, 10),
        (80000, 400, 4, 0.9, 0.75, 500, 10, 10, 0.9, 10),
    ],
)
def test_crack_min_steel_area_with_press_tendons_raise_valueerror(
    a_ct, s_steel, fct_eff, k, kc, ap, d_steel, d_press, e, incr_stress
):
    """Test the crack_min_steel_area raise ValueError for non valid values"""
    with pytest.raises(ValueError):
        _crack_control.crack_min_steel_area_with_prestresed_tendons(
            a_ct, s_steel, fct_eff, k, kc, ap, d_steel, d_press, e, incr_stress
        )


@pytest.mark.parametrize(
    'wk, s_steel, fct_eff, h_cr, h, d, incr_stress, kc, exp_phi,  exp_sep',
    [
        (0.3, 240, 2.9, 200, 400, 360, 40, 0.4, 25, 250),
        (0.2, 260, 2.9, 200, 400, 360, 40, None, 8.75, 125),
        (0.35, 360, 2.9, 200, 400, 360, 40, 0.4, 11, 125),
        (0.35, 360, 2.9, 200, 400, 360, 40, None, 6.875, 125),
    ],
)
def test_crack_min_steel_without_direct_calculation_returns_expected_values(
    wk,
    s_steel,
    fct_eff,
    h_cr,
    h,
    d,
    incr_stress,
    kc,
    exp_phi,
    exp_sep,
):
    """Test the crack_min_steel_area raise ValueError for non valid values"""
    phi, sep = _crack_control.crack_min_steel_without_direct_calculation(
        wk, s_steel, fct_eff, h_cr, h, d, incr_stress, kc
    )
    assert math.isclose(phi, exp_phi, rel_tol=10e-6)
    assert math.isclose(sep, exp_sep, rel_tol=10e-6)


@pytest.mark.parametrize(
    'wk, s_steel, fct_eff, h_cr, h, d, incr_stress, kc',
    [
        (-0.1, 200, 3, 250, 300, 280, 0, 0.7),
        (0.2, 200, -3, 250, 300, 280, 0, 0.7),
        (0.2, 200, 3, 250, 300, 280, 0, 1.1),
        (0.2, 200, 3, -250, 300, 280, 0, 0.7),
        (0.2, 200, 3, -250, -300, 280, 0, 0.7),
        (0.2, 200, 3, -250, -300, -280, 0, 0.7),
        (0.2, 360, 2.9, 200, 400, 360, 0, 0.4),
        (0.5, 200, 2.9, 200, 400, 360, 0, 0.4),
    ],
)
def test_crack_min_steel_without_direct_calculation_raise_valueerror(
    wk, s_steel, fct_eff, h_cr, h, d, incr_stress, kc
):
    """Test the crack_min_steel_area raise ValueError for non valid values"""
    with pytest.raises(ValueError):
        _crack_control.crack_min_steel_without_direct_calculation(
            wk, s_steel, fct_eff, h_cr, h, d, incr_stress, kc
        )


@pytest.mark.parametrize(
    'e, d_press, d_steel, expected',
    [
        (0.8, 20, 0, 0.894427),
        (0.6, 25, 10, 0.489898),
        (0.5, 10, 10, 0.707107),
    ],
)
def test_adjusted_bond_length_return_expected_values(
    e, d_press, d_steel, expected
):
    """Test the adjusted_bond_length_function returns expected values"""
    assert math.isclose(
        _crack_control.adjusted_bond_strength(e, d_press, d_steel),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    'e, d_press, d_steel',
    [
        (0.1, 20, 0),
        (-2, 25, 10),
        (1.15, 10, 10),
        (0.6, -10, 10),
        (0.6, 10, -10),
    ],
)
def test_adjusted_bond_length_raise_valuerror(e, d_press, d_steel):
    """Test the adjusted_bond_length_function raises exceptions"""
    with pytest.raises(ValueError):
        _crack_control.adjusted_bond_strength(e, d_press, d_steel)


@pytest.mark.parametrize(
    'h, d, x, expected',
    [
        (400, 200, 100, 100),
        (400, 200, 150, 83.333333),
        (550, 150, 150, 133.33333),
    ],
)
def test_hc_eff_concrete_tension_returns_expected_values(h, d, x, expected):
    """Test the hc_eff_concrete_tension returns expected results"""
    assert math.isclose(
        _crack_control.hc_eff_concrete_tension(h, d, x),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    'h, d, x',
    [
        (-50, 200, 100),
        (50, -200, 100),
        (50, 200, -100),
        (400, 450, 100),
        (400, 200, 450),
    ],
)
def test_hc_eff_concrete_tension_raise_exceptions(h, d, x):
    """Test hc_eff_concrete tension raises expected exceptions"""
    with pytest.raises(ValueError):
        _crack_control.hc_eff_concrete_tension(h, d, x)


@pytest.mark.parametrize(
    'es, ecm, expected',
    [
        (10e9, 10e5, 1e4),
    ],
)
def test_alpha_e_returns_expected_values(es, ecm, expected):
    """Test alpha_e returns expected values"""
    assert math.isclose(
        _crack_control.alpha_e(es, ecm),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    'es, ecm',
    [
        (-10e9, 10e5),
        (100e9, -10e-5),
    ],
)
def test_alpha_e_raise_exceptions(es, ecm):
    """Test alpha_e raises exceptions"""
    with pytest.raises(ValueError):
        _crack_control.alpha_e(es, ecm)


@pytest.mark.parametrize(
    'a_s, e1, a_p, ac_eff, expected',
    [
        (200, 0.8, 125, 600, 0.46666667),
        (125, 1.5, 125, 1200, 0.33854),
    ],
)
def test_rho_p_eff_returns_expected_values(a_s, e1, a_p, ac_eff, expected):
    """Test rho_p_eff returns expeceted values"""
    assert math.isclose(
        _crack_control.rho_p_eff(a_s, e1, a_p, ac_eff), expected, rel_tol=10e-5
    )


@pytest.mark.parametrize(
    'a_s, e1, a_p, ac_eff',
    [
        (-200, 0.8, 125, 600),
        (200, -0.8, 125, 600),
        (200, 0.8, -125, 600),
        (200, 0.8, 125, -600),
    ],
)
def test_rho_p_eff_raise_value_error(a_s, e1, a_p, ac_eff):
    """Test rho_p_eff raise exceptions"""
    with pytest.raises(ValueError):
        _crack_control.rho_p_eff(a_s, e1, a_p, ac_eff)


@pytest.mark.parametrize(
    'load_type, expected',
    [
        ('short', 0.6),
        ('long', 0.4),
    ],
)
def test_kt_load_duration_returns_expected_values(load_type, expected):
    """Test kt_load_duration returns expected values"""
    assert _crack_control.kt_load_duration(load_type) == expected


def test_kt_load_duration_raise_value_errors():
    """Test kt_load_duration raise value errors"""
    with pytest.raises(TypeError):
        _crack_control.kt_load_duration(load_type=123)

    with pytest.raises(ValueError):
        _crack_control.kt_load_duration(load_type='asdf')

import structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage
import structuralcodes.codes.ec2_2004._concrete_material_properties
import structuralcodes.codes.ec2_2004._reinforcement_material_properties
import structuralcodes.codes.ec2_2004._section_7_3_crack_control
import structuralcodes.codes.ec2_2004.shear
import structuralcodes.codes.ec2_2023._annexB_time_dependent
import structuralcodes.codes.ec2_2023._section5_materials
import structuralcodes.codes.ec2_2023._section9_sls
import xlwings as xw


def main():
    """Example function to demonstrate xlwings UDFs."""
    wb = xw.Book.caller()
    sheet = wb.sheets[0]
    # Toggle the cell value between two states
    if sheet['A1'].value == 'It works!':
        sheet['A1'].value = 'Still works!'
    else:
        sheet['A1'].value = 'It works!'


@xw.func
def SC_ec2_2023_annexB_time_dependent_alpha_c(fcm_28):
    """Returns the coefficient to obtain the tangent modulus of elasticity from
    the secant modulus of elasticity.

    EN 1992-1-1, Table 5.1

    Args:
    fcm_28 (float): The characteristic compressive strength at 28 days
        in MPa.

    Returns:
    float: Coefficient alphac so that Ec=alphac*Ecm,28.
    """
    return structuralcodes.codes.ec2_2023._annexB_time_dependent.alpha_c(
        fcm_28
    )


@xw.func
def SC_ec2_2023_section5_materials_A_phi_correction_exp(hn, atm_conditions):
    """Computes the correction exponent for the modification for the phi_50y_t0
    with respect the fck value.

    EN 1992-1-1:2023, Table 5.2.

    Args:
    hn (float): The notional size in mm.
    atm_conditions (str): 'dry' or 'humid'.

    Returns:
    float: The correction exponent value.

    Raises:
    ValueError: If hn is less than 100 or greater than 1000.
    ValueError: If atm_conditions is not 'dry' or 'humid'.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.A_phi_correction_exp(
        hn, atm_conditions
    )


@xw.func
def SC_ec2_2023_section5_materials_Ecm(fcm, kE):
    """Computes the secant modulus between sigma_c=0 and sigma_c=0.4*fcm.

    EN 1992-1-1:2023, Eq. (5.1).

    Args:
    fcm (float): The mean compressive strength in MPa.

    Keyword Args:
    kE (float): Coefficient relating the aggregates used in concrete.
        Default value is 9500, but it can vary from 5000 to 13000.

    Returns:
    float: The secant concrete modulus in MPa.

    Raises:
    ValueError: If fcm is less than 0.
    ValueError: If kE is not between 5000 and 13000.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.Ecm(fcm, kE)


@xw.func
def SC_ec2_2023_section5_materials_Es():
    """Returns the value of the modulus of elasticity for weldable reinforcing
    steel.

    EN 1992-1-1:2023, 5.2.4-3.

    Returns:
    float: Modulus of elasticity in MPa.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.Es()


@xw.func
def SC_ec2_2023_section5_materials_alpha_c_th():
    """Returns the linear coefficient of thermal expansion in 1/Cº for
    concrete.

    EN 1992-1-1:2023, 5.1.6-6.

    Returns:
    float: The linear coefficient of thermal expansion in 1/Cº for
    concrete.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.alpha_c_th()


@xw.func
def SC_ec2_2023_section5_materials_alpha_s_th():
    """Returns the linear coefficient of thermal expansion in 1/Cº for weldable
    reinforced steel.

    EN 1992-1-1:2023, 5.2.4-5

    Returns:
    float: The linear coefficient of thermal expansion in 1/Cº for weldable
    reinforcement steel.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.alpha_s_th()


@xw.func
def SC_ec2_2023_section5_materials_eps_c1(fcm):
    """Computes the strain at maximum compressive strength of concrete (fcm)
    for the Sargin constitutive law.

    EN 1992-1-1:2023, Eq. (5.9).

    Args:
    fcm (float): The mean strength of concrete in MPa.

    Returns:
    float: The strain at maximum compressive strength of concrete.

    Raises:
    ValueError: If fcm is less than 12+8MPa.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.eps_c1(fcm)


@xw.func
def SC_ec2_2023_section5_materials_eps_c2():
    """The strain at maximum compressive stress of concrete for the
    parabolic-rectangular law.

    EN 1992-1-1:2023, Eq. 8.4

    Returns:
    float: The strain at maximum compressive stress, absolute value, no
    unit.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.eps_c2()


@xw.func
def SC_ec2_2023_section5_materials_eps_cs_50y(
    fck_28, atm_conditions, hn, strength_dev_class
):
    """Computes the nominal total shrinkage in ‰ for concrete after a duration
    of drying of 50 years.

    EN 1992-1-1:2023, Table 5.3.

    Args:
    fck_28 (float): Characteristic strength at 28 days in MPa
    atm_conditions (str): 'dry' or 'humid'.
    hn (float): The notional size in mm.
    strength_dev_class (str): Strength development class 'CS', 'CN', 'CR',
        'slow', 'normal' or 'rapid'.

    Returns:
    float: The nominal shrinkage value in percent.

    Raises:
    ValueError: If fck_28 is less than 20 MPa or larger than 80 MPa.
    ValueError: If atm_conditions is not 'dry' or 'humid'.
    ValueError: If hn is less than 100 or larger than 1000.
    ValueError: If strength_dev_class is not CS', 'CN', 'CR',
        'slow', 'normal' or 'rapid'.
    ValueError: If combination of fck_28 and hn is out of scope.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.eps_cs_50y(
        fck_28, atm_conditions, hn, strength_dev_class
    )


@xw.func
def SC_ec2_2023_section5_materials_eps_cu1(fcm):
    """Computes the strain at concrete failure of concrete.

    EN 1992-1-1:2023, Eq. (5.10).

    Args:
    fcm (float): The mean strength of concrete in MPa.

    Returns:
    float: The maximum strength at failure of concrete.

    Raises:
    ValueError: If fcm is less than 12+8MPa.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.eps_cu1(fcm)


@xw.func
def SC_ec2_2023_section5_materials_eps_cu2():
    """The ultimate strain of the parabolic-rectangular law.

    EN 1992-1-1:2023, Eq. 8.4

    Returns:
    float: The ultimate strain, absolute value, no unit.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.eps_cu2()


@xw.func
def SC_ec2_2023_section5_materials_eps_ud(eps_uk, gamma_s):
    """Design value for the ultimate limit strain welding reinforcing steel.

    EN 1992-1-1:2023, 5.2.4-2.

    Args:
    eps_uk (float): Characteristic ultimate limit strain.
    gamma_s (float): Safety coefficient.

    Returns:
    float: Design ultimate strain limit.

    Raises:
    ValueError: If eps_uk is less than 0.
    ValueError: If gamma_s is less than 1.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.eps_ud(
        eps_uk, gamma_s
    )


@xw.func
def SC_ec2_2023_section5_materials_eta_cc(fck, fck_ref):
    """Computes the factor to measure the difference between the undistributed
    compressive strength of a cylinder and the effective compressive strength
    in a structural member.

    EN 1992-1-1:2023, Eq. (5.4).

    Args:
    fck (float): The characterisitic compressive strength in MPa.

    Keyword Args:
    fck_ref (float, optional): The reference compressive strength MPa.

    Returns:
    float: The value of the factor eta_cc.

    Raises:
    ValueError: If fck is less than 12 MPa.
    ValueError: If fkc_ref is less or equal to 0.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.eta_cc(
        fck, fck_ref
    )


@xw.func
def SC_ec2_2023_section5_materials_fcd(fck, eta_cc, k_tc, gamma_c):
    """Computes the value of the design compressive strength of concrete.

    EN 1992-1-1:2023, Eq. (5.3).

    Args:
    fck (float): Characteristic compressive strength in MPa.
    eta_cc (float): Factor for measuring the difference between the
        undistributed compressive strength of a cylinder and the effective
        compressive strength in the real structural member.
    k_tc (float): Factor for taking into consideration high sustained
        loads and of time of loading.
    gamma_c (float): Partial factor of concrete.

    Returns:
    float: The design compressive strength of concrete in MPa.

    Raises:
    ValueError: If fck is less than 12 MPa.
    ValueError: If _etc_cc is not between 0 and 1.
    ValueError: If gamma_c is less than 1.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.fcd(
        fck, eta_cc, k_tc, gamma_c
    )


@xw.func
def SC_ec2_2023_section5_materials_fcm(fck, delta_f):
    """Determines the mean strength of concrete from its characteristic value.

    EN 1992-1-1:2023, Table 5.1.

    Args:
    fck (float): Is the characteristic compressive strength in MPa.

    Keyword Args:
    delta_s (float): The increment in MPa to compute the mean compressive
        strength.

    Returns:
    float: The mean compressive strength in MPa
    """
    return structuralcodes.codes.ec2_2023._section5_materials.fcm(fck, delta_f)


@xw.func
def SC_ec2_2023_section5_materials_fctd(fctk_5, k_tt, gamma_c):
    """Computes the value of the design tensile strength of concrete.

    EN 1992-1-1:2023, Eq. (5.5).

    Args:
    fctk_5 (float): The 5% mean concrete tensile strength fractile in MPa.
    k_tt (float): The factor for considering the effect of high sustained
        loads and of time of loading on concrete tensile strength.
    gamma_c (float): Partial factor of concrete.

    Returns:
    float: The design tensile strength of concrete in MPa.

    Raises:
    ValueError: If fctk_5 is less than 0.
    ValueError: If gamma_c is less than 1.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.fctd(
        fctk_5, k_tt, gamma_c
    )


@xw.func
def SC_ec2_2023_section5_materials_fctk_5(fctm):
    """Compute the 5% mean concrete tensile strength fractile.

    EN 1992-1-1:2023, Table 5.1.

    Args:
    fctm (float): The mean concrete tensile strength in MPa.

    Returns:
    float: The 5% mean concrete tensile strength fractile in MPa.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.fctk_5(fctm)


@xw.func
def SC_ec2_2023_section5_materials_fctk_95(fctm):
    """Compute the 95% mean concrete tensile strength fractile.

    EN 1992-1-1:2023, Table 5.1.

    Args:
    fctm (float): The mean concrete tensile strength in MPa.

    Returns:
    float: The 5% mean concrete tensile strength fractile in MPa.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.fctk_95(fctm)


@xw.func
def SC_ec2_2023_section5_materials_fctm(fck):
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength.

    EN 1992-1-1:2023, Table 5.1.

    Args:
    fck (float): The characteristic compressive strength in MPa.

    Returns:
    float: The mean tensile strength in MPa.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.fctm(fck)


@xw.func
def SC_ec2_2023_section5_materials_fpd(fp01k, gamma_P):
    """Computes the design value for the prestressing steel stress.

    EN 1992-1-1:2023, 5.3.3.

    Args:
    fp01k (float): The 0.1% proof stress in MPa.
    gamma_P (float): The safety coefficient.

    Returns:
    float: The design value for the design prestressing steel stress in
    MPa.

    Raises:
    ValueError: If fp01k is less than 0.
    ValueError: If gamma_P is less than 1.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.fpd(
        fp01k, gamma_P
    )


@xw.func
def SC_ec2_2023_section5_materials_fyd(fyk, gamma_s):
    """Design value for the yielding stress for welding reinforcing steel.

    EN 1992-1-1:2023, Eq (5.11).

    Args:
    fyk (float): Characteristic yield stress for the steel in MPa.
    gamma_s (float): Safety coefficient.

    Returns:
    float: Design yielding stress for steel in MPa.

    Raises:
    ValueError: If fyk is less than 0.
    ValueError: If gamma_s is less than 1.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.fyd(fyk, gamma_s)


@xw.func
def SC_ec2_2023_section5_materials_hn(Ac, u):
    """Computes the notional size of a given concrete cross-section.

    EN 1992-1-1:2023, Table 5.2.

    Args:
    Ac (float): The concrete cross-sectional area in (mm2).
    u (float): The perimeter exposed to drying (mm).

    Returns:
    float: The notional size (mm).

    Raises:
    ValueError: If Ac is less than 0.
    ValueError: If u is less than 0.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.hn(Ac, u)


@xw.func
def SC_ec2_2023_section5_materials_k_sargin(Ecm, fcm, eps_c1):
    """Computes the coefficient k for Sargin constitutive law in compression.

    EN 1992-1-1:2003, eq. (5.7)

    Args:
    Ecm (float): the secant modulus between sigma_c=0 and sigma_c=0.4*fcm
        in MPa
    fcm (float): the mean compressive strength of concrete in MPa
    eps_c1 (float): the strain of concrete at stress fcm

    Returns:
    float: the coefficient k for Sargin constitutive law

    Raises:
    ValueError: if Ecm is less or equal to 0
    ValueError: if fcm is less than 12+8MPa
    ValueError: if eps_c1 is less or equal to 0
    """
    return structuralcodes.codes.ec2_2023._section5_materials.k_sargin(
        Ecm, fcm, eps_c1
    )


@xw.func
def SC_ec2_2023_section5_materials_k_tc(t_ref, t0, strength_dev_class):
    """Computes the factor for considering the effect of high sustained loads
    and of time of loading on concrete compressive strength.

    EN 1992-1-1:2023, Eq. (5.3).

    Args:
    t_ref (float): The reference time in days.
    t0 (float): Age at loading in days.
    strength_dev_class (str): Strength development class 'CS', 'CN', 'CR',
        'slow', 'normal', 'rapid'.

    Returns:
    float: The factor value.

    Raises:
    ValueError: If t_ref is less than 0.
    ValueError: If t0 is less than 0.
    ValueError: If strength_dev_class is not 'CS', 'CN', 'CR', 'slow',
        'normal', 'rapid'.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.k_tc(
        t_ref, t0, strength_dev_class
    )


@xw.func
def SC_ec2_2023_section5_materials_k_tt(t_ref, strength_dev_class):
    """Computes the factor for considering the effect of high sustained loads
    and of time of loading on concrete tensile strength.

    EN 1992-1-1:2023, Eq. (5.5).

    Args:
    t_ref (float): The reference time in days.
    strength_dev_class (str): Strength development class 'CS', 'CN', 'CR',
        'slow', 'normal' or 'rapid'.

    Returns:
    float: The factor value.

    Raises:
    ValueError: If t_ref is less than 0.
    ValueError: If strength_dev_class is not 'CS', 'CN', 'CR', 'slow',
        'normal' or 'rapid'.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.k_tt(
        t_ref, strength_dev_class
    )


@xw.func
def SC_ec2_2023_section5_materials_n_parabolic_rectangular():
    """The exponent in the parabolic-rectangular law.

    EN 1992-1-1:2023, Eq. 8.4
    Returns:
    float: The exponent n, absolute value, no unit.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.n_parabolic_rectangular()


@xw.func
def SC_ec2_2023_section5_materials_p_steel_stress_params(
    prestress_class, element
):
    """Computes the stress-diagram parameters fp01k and fpk.

    EN 1992-1-1:2023, 5.3.3.

    Args:
    prestress_class (str): Possible values: Y1560, Y1670,
        Y1770, Y1860, Y1770, Y1860, Y1960, Y2060,
        Y1030, Y1050, Y1100 and Y1230
    element (str): Element type, 'W' for Wires, 'S' for Strands, and 'B'
        for Bars.

    Returns:
    Tuple(float, float): With the value of fp01k and fpk in MPa.

    Raises:
    ValueError: If combination of prestress_class and element is not a
        possible value from the range.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.p_steel_stress_params(
        prestress_class, element
    )


@xw.func
def SC_ec2_2023_section5_materials_phi_50y_t0(
    t0, atm_conditions, hn, strength_dev_class
):
    """Computes the creep coefficient of plain concrete at 50 years of loading.
    Interpolation is linear between values.

    EN 1992-1-1:2023, Table 5.2.

    Args:
    t0 (float): Age at loading [days].
    atm_conditions (str): 'dry' or 'humid'.
    hn (float): The notional size in mm.
    strength_dev_class (str): Strength development class 'CS', 'CN', 'CR',
        'slow', 'normal' or 'rapid'.

    Returns:
    float: The creep coefficient.

    Raises:
    ValueError: If t0 is less than 1.
    ValueError: If atm_conditions is not 'dry' or 'humid'.
    ValueError: If hn is less than 100 or larger than 1000.
    ValueError: If strength_dev_class is not 'CS', 'CN',
        'CR', 'slow', 'normal' or 'rapid'.
    ValueError: If combination of t0 and hn is out of scope.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.phi_50y_t0(
        t0, atm_conditions, hn, strength_dev_class
    )


@xw.func
def SC_ec2_2023_section5_materials_phi_correction_factor(fck, A_exponent):
    """Computes the correction factor for the computation of the phi_50y_t0.

    EN 1992-1-1:2023, Table 5.2.

    Args:
    fck (float): Characteristic strength of concrete in MPa.
    A_exponent (float): The A correction exponent value.

    Returns:
    float: The correction factor value.

    Raises:
    ValueError: If fck is not between 12 and 100 MPa.
    ValueError: If A_exponent is not between 0.64 and 0.82.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.phi_correction_factor(
        fck, A_exponent
    )


@xw.func
def SC_ec2_2023_section5_materials_reinforcement_duct_props(
    fyk, ductility_class
):
    """Return a dict with the minimum characteristic ductility properties for
    reinforcement ductility class.

    EUROCODE 2 1992-1-1:2023, Tab. 5.5.

    Args:
    fyk (float): The characteristic yield strength.
    ductility_class (Literal['A', 'B', 'C']): The reinforcement ductility
        class designation.

    Returns:
    Dict[str, float]: A dict with the characteristik strain value at the
    ultimate stress level (epsuk), and the characteristic ultimate stress
    (ftk).
    """
    return structuralcodes.codes.ec2_2023._section5_materials.reinforcement_duct_props(
        fyk, ductility_class
    )


@xw.func
def SC_ec2_2023_section5_materials_sigma_p(eps, fpy, fpu, eps_u, Ep):
    """Computes the stress for prestressing steel as a function of the strain.

    EN 1992-1-1:2023, 5.3.3.

    Args:
    eps (float): Strain value.
    fpy (float): Yielding stress of the steel in MPa. Use fd for design
        stress values, and fp01k for nominal stress values.
    fpu (float): The maximum stress at eps_u in MPa. Use fpd for design
        stress values, fpk for nominal stress values, and fpu == fpy for
        horizontal post-elastic branch without strain limit.

    Keyword Args:
    eps_u (float): Ultimate strain. Use eps_uk = 0.035 for nominal ultimate
        strain, and eps_ud for design ultimate strain.
    Ep (float): Modulus of elasticity of prestressing steel in MPa.

    Raises:
    ValueError: If eps is less than 0 or larger than eps_u.
    ValueError: If fpy is less or equal to 0.
    ValueError: If fpu is less than fpy.
    ValueError: If eps_u is lower or equal to 0.
    ValueError: If _Ep is less or equal to 0.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.sigma_p(
        eps, fpy, fpu, eps_u, Ep
    )


@xw.func
def SC_ec2_2023_section5_materials_sigma_s(eps, fy, k, eps_u, Es):
    """Compute the stress for welded reinforcing steel in MPa for a given
    strain.

    EN 1992-1-1:2023, 5.2.4.

    Args:
    eps (float): The strain value.
    fy (float): The yielding stress in MPa. Use fyd for the design
        strength, and fyk for the characteristic strength.
    k (float): Curve parameter. Ratio between the ultimate stress and the
        yielding stress. k = 1 for horizontal post-elastic branch without
        strain limit.
    eps_u (float): Ultimate strain at failure. Use eps_ud for the design
        ultimate strain, and eps_uk for the characteristic ultimate strain.

    Keyword Args:
    Es (float): The modulus of elasticity for reinforcing steel.

    Returns:
    float: The nominal stress in MPa.

    Raises:
    ValueError: If eps is less than 0 or larger than eps_uk.
    ValueError: If Es is less or equal to 0.
    ValueError: If fy is less or equal to 0.
    ValueError: If k is less than 1.
    ValueError: If eps_u is less or equal to 0 and k > 1.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.sigma_s(
        eps, fy, k, eps_u, Es
    )


@xw.func
def SC_ec2_2023_section5_materials_weight_c(concrete_type):
    """Returns the mean unit weight of concrete in kN/m3.

    EN 1992-1-1:2023, 5.1.6-5.

    Args:
    concrete_type (str): 'nc' for normal concrete, or 'npc' for normal
        plain concrete.

    Returns:
    float: Mean unit weight in kN/m3.

    Raises:
    ValueError: If concrete_type is not 'nc' or 'npc'.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.weight_c(
        concrete_type
    )


@xw.func
def SC_ec2_2023_section5_materials_weight_s():
    """Returns the mean unit weight of reinforced steel for the purposes of
    design in kN/m3.

    EN 1992-1-1:2023.2.4-4.

    Returns:
    float: The mean unit weight in kN/m3.
    """
    return structuralcodes.codes.ec2_2023._section5_materials.weight_s()


@xw.func
def SC_ec2_2023_section9_sls_As_min_y(NEd, b, h, fct_eff, fyk):
    """Returns the minimum reinforcement to avoid yielding of steel. Box or T
    sections are to be divided into rectangles.

    EN 1992-1-1:2023, Eq. (9.4)

    Eq. (9.2) and (9.3) are particular cases of the general equation

    Eq. (9.2) is valid for pure bending, hence NEd=0

    Eq. (9.3) is valid for pure tension. The general expression has an upper
    limit that equals the values of Eq. (9.3)

    Args:
    NEd (float): SLS axial force applied on the section or rectangle
        (compressions are negative) in kN.
    b (float): The width of the section or rectangle in meters.
    h (float): The height of the section or rectange in meters.
    fct_eff (float): Effective tension strength of concrete (can normally
        be taken as the mean tensile strength) in MPa.
    fyk (float): Characteristic yield strength of steel in MPa.

    Returns:
    tuple(float, float): The minimum tensile reinforcement to avoid
    yielding of steel on the most tensioned fibre of the rectangle
    (As_min_y1) in cm2, and the minimum tensile reinforcement to avoid
    yielding of steel on the most tensioned fibre of the rectangle
    (As_min_y2) in cm2.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.As_min_y(
        NEd, b, h, fct_eff, fyk
    )


@xw.func
def SC_ec2_2023_section9_sls_Ec_eff(fcm, phi, kE):
    """Returns de effective modulus of elasticity from fcm and phi.

    EN 1992-1-1:2023, Eq. (9.1).

    Args:
    fcm (float): The mean compressive strength in MPa.
    phi (float): The creep coefficient.

    Keyword Args:
    kE (float): Constant to account for the type of aggregate.

    Returns:
    float: The effective modulus of elastiticy in MPa.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.Ec_eff(fcm, phi, kE)


@xw.func
def SC_ec2_2023_section9_sls_alpha_c(fcm_28):
    """Returns the coefficient to obtain the tangent modulus of elasticity from
    the secant modulus of elasticity.

    EN 1992-1-1, Table 5.1

    Args:
    fcm_28 (float): The characteristic compressive strength at 28 days
        in MPa.

    Returns:
    float: Coefficient alphac so that Ec=alphac*Ecm,28.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.alpha_c(fcm_28)


@xw.func
def SC_ec2_2023_section9_sls_delta_simpl(
    delta_loads, delta_shr, fck1, phi1, b1, h, d, As1, Mk
):
    """Simplified calculation of the deflection for rectangular sections.

    EN1992-1-1:2023, Eq. (9.23).

    Args:
    delta_loads (float): Linear elastic deflection due to loads.
    delta_shr (float): Linear elastic deflection due to shrinkage.
    fck1 (float): Characteristic concrete strength in MPa.
    phi1 (float): Weighted mean value of the creep coefficient.
    b1 (float): Width of rectangular cross-section in m.
    h (float): Height of rectanguar cross-section in m.
    d (float): Effective height of cross-section in m.
    As1 (float): Tension reinforcement at centre span for continuous in cm2
        beams or at the embedment for a cantilever.
    Mk (float): Characteristic moment at centre span for continuous
        beams or at the embedment for a cantilever.

    Returns:
    float: The deflection of the beam in units consistent with delta_loads
    and delta_shr.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.delta_simpl(
        delta_loads, delta_shr, fck1, phi1, b1, h, d, As1, Mk
    )


@xw.func
def SC_ec2_2023_section9_sls_epssm_epscm(
    sigma_s, kt, fct_eff, rho_eff, alphae, Es
):
    """Returns the mean strain difference between steel and concrete along 2
    transfer lengths.

    EN 1992-1-1:2023, Eq. (9.11).

    Args:
    sigma_s (float): The stress in steel at the section of the crack.
    kt (float): An integration factor to account for the variation in
        strain in steel and concrete it is to be taken as 0.6 for short
        term loading or instantaneous loading and equal to 0.4 for long
        term or repeated loading.
    fct_eff (float): The effective cracking stress, which can be taken
        equal to the mean tensile strength of concrete, fctm.
    rho_eff (float): The effective reinforcement ratio in the tension zone.
    alphae (float): The equivalence factor equal to Es/Ecm.
    Es (float): The modulus of elasticity of steel, normally taken as 200
        GPa.

    Returns:
    float: The mean strain difference bewteen steel and concrete along 2
    transfer lengths.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.epssm_epscm(
        sigma_s, kt, fct_eff, rho_eff, alphae, Es
    )


@xw.func
def SC_ec2_2023_section9_sls_fcm(fck, delta_f):
    """Determines the mean strength of concrete from its characteristic value.

    EN 1992-1-1:2023, Table 5.1.

    Args:
    fck (float): Is the characteristic compressive strength in MPa.

    Keyword Args:
    delta_s (float): The increment in MPa to compute the mean compressive
        strength.

    Returns:
    float: The mean compressive strength in MPa
    """
    return structuralcodes.codes.ec2_2023._section9_sls.fcm(fck, delta_f)


@xw.func
def SC_ec2_2023_section9_sls_fctm(fck):
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength.

    EN 1992-1-1:2023, Table 5.1.

    Args:
    fck (float): The characteristic compressive strength in MPa.

    Returns:
    float: The mean tensile strength in MPa.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.fctm(fck)


@xw.func
def SC_ec2_2023_section9_sls_k_1_r(h, x, ay):
    """Returns k1/r factor to account for increase in crack width due to
    curvature of the section in bending.

    EN 1992-1-1:2023, Eq. (9.9).

    Args:
    h (float): Height of the section in consistent units (e.g. meters).
    x (float): Distance from most compressed fibre to neutra axis in
        consistent units (e.g. meters).
    ay (float): Cover to centre of tensioned reinforcement closest to most
        tensioned face in consistent units (e.g. meters).

    Returns:
    float: Factor k1/r (non-dimensional) which accounts for the increase in
    crack width due to curvature of the section in bending.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.k_1_r(h, x, ay)


@xw.func
def SC_ec2_2023_section9_sls_kfl(h, xg, hceff):
    """Returns factor kfl which accounts for the distribution of stresses
    before cracking.

    EN 1992-1-1:2023, Eq. (9.17).

    Args:
    h (float): Height of the cross section.
    xg (float): Distance from the compressed fibre to the centroid of the
        uncracked section.
    hceff (float): Height of the effective tension area.

    Returns:
    float: Returns factor kfl which accounts for the distribution of
    stresses before cracking.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.kfl(h, xg, hceff)


@xw.func
def SC_ec2_2023_section9_sls_kh(b, h):
    """Returns factor kh, which reduces the tensile strength of concrete to
    account for imposed restrained deformations due to shrinkage.

    EN 1992-1-1:2023, Eq. (9.5).

    Args:
    b (float): Width of the rectangle in meters.
    h (float): Height of the rectangle in meters.

    Returns:
    float: Factor kh which reduces the tensile strength of concrete to
    account for imposed restrained deformations due to shrinkage.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.kh(b, h)


@xw.func
def SC_ec2_2023_section9_sls_srm_cal(c, kfl_, kb, phi, rho_eff, kw, h, x):
    """Returns the mean crack spacing.

    EN 1992-1-1:2023, Eq. (9.15).

    Args:
    c (float): Concrete cover of reinforcement to bar surface. Larger value
        of lateral and vertical cover should be applied.
    kfl (float): Factor accounting for distribution of stresses prior to
        cracking.
    kb (float): Factor accounting for bond conditions.
    phi (float): Bar diameter.
    rho_eff(float): Effective reinforcement ratio in the tension zone.
    kw (float): Factor converting the mean crack spacing into a
        characteristic crack spacing, with a reocmmended value of 1.3
        (NDP).
    h (float): Height of the cross section.
    x (float): Depth of the neutral axis measured form the most compressed
        fibre.

    Returns:
    float: The mean crack spacing in units consistent with c and phi.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.srm_cal(
        c, kfl_, kb, phi, rho_eff, kw, h, x
    )


@xw.func
def SC_ec2_2023_section9_sls_wk_cal(
    kw, h, xg, hc_eff, c, kb, phi, rho_eff, x, sigma_s, kt, fct_eff, alphae, Es
):
    """Returns the characteristic crack width, wk,cal, as well as auxiliary
    variables, 1/r, srm,cal and epssm-epscm.

    EN1992-1-1:2023 Eq. (9.8), complemented with Eq. (9.11), Eq. (9.15), Eq.
    (9.17).

    Args:
    kw (float): Factor that converts the mean crack spacing to a
        characteristic value.
    h (float): Height of cross section.
    xg (float): Depth of centroid of section measured from compressed
        fibre.
    hc_eff (float): Height of the effective tensioned concrete area.
    c (float): Concrete cover of reinforcement to bar surface. Larger
        value of lateral and vertical cover should be applied.
    kb (float): Factor account for bond conditions of bar.
    phi (float): Diameter of tensioned bars (for different bar diameters,
        equivalent diameter according to Eq. (9.19).
    rho_eff (float): Effective tension reinforcement ratio.
    x (float): Depth of the neutral axis of the cracked section measured
        from compressed fibre.
    sigma_s (float): Tension in most tensioned bar according to fully
        cracked analysis.
    kt (float): Factor accounting for tension stiffening.
    fct_eff (float): Effective tensile strength of concrete.
    alphae (float): Modular ratio Es/Ecm.
    Es (float): Modulus of elasticity of steel bars (normally Es=200 MPa).

    Returns:
    Tuple[float, float, float, float]: The characteristic crack width,
    wk,cal, in consistent units, as well as auxiliary variables, 1/r,
    srm,cal and epssm-epscm.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.wk_cal(
        kw,
        h,
        xg,
        hc_eff,
        c,
        kb,
        phi,
        rho_eff,
        x,
        sigma_s,
        kt,
        fct_eff,
        alphae,
        Es,
    )


@xw.func
def SC_ec2_2023_section9_sls_wk_cal2(kw, k_1_r, srm_cal, epssm_epscm):
    """Returns the calculated characteristic crack width.

    EN 1992-1-1:2023, Eq. (9.8).

    Args:
    kw (float): Factor that converts the mean crack spacing to a
        characteristic value.
    k_1_r (float): Factor accounting for the effect of curvature on crack
        width - can be determined using the function k_1_r.
    srm_cal (float): Mean crack spacing - can be determined using the
        function srm_cal.
    epssm_epscm (float): Mean diference of strain between steel and
        concrete - can be determined using the function epssm_epscm.

    Returns:
    float: The calculated characteristic crack width in in units consistent
    with srm_cal.
    """
    return structuralcodes.codes.ec2_2023._section9_sls.wk_cal2(
        kw, k_1_r, srm_cal, epssm_epscm
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage__get_alpha_ds_dict(cement_class):
    """Return a dictionary with values for aplha_ds1 and alpha_ds2."""
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage._get_alpha_ds_dict(
        cement_class
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_alpha_1(fcm):
    """Return a factor taking into account the effect of concrete strength.

    EN 1992-1-1:2004, Eq. (B.8c).

    Args:
    fcm (float): The mean concrete strength in MPa.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.alpha_1(
            fcm
        )
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_alpha_2(fcm):
    """Return a factor taking into account the effect of concrete strength.

    EN 1992-1-1:2004, Eq. (B.8c).

    Args:
    fcm (float): The mean concrete strength in MPa.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.alpha_2(
            fcm
        )
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_alpha_3(fcm):
    """Return a factor taking into account the effect of concrete strength.

    EN 1992-1-1:2004, Eq. (B.8c).

    Args:
    fcm (float): The mean concrete strength in MPa.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.alpha_3(
            fcm
        )
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_alpha_cement(cement_class):
    """Return an exponent that depends on the cement type.

    Args:
    cement_class (str): The cement class, either 'S', 'N' or 'R'.

    Returns:
    float: The exponent that depends on the cement type.

    Raises:
    ValueError: If an invalid cement class is provided.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.alpha_cement(
        cement_class
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_alpha_ds1(cement_class):
    """Return a coefficient depending on the cement class.

    EN 1992-1-1:2004, Sec. B.2.

    Args:
    cement_class (str): The cement class, either 'S', 'N' or 'R'.

    Returns:
    float: The exponent that depends on the cement type.

    Raises:
    ValueError: If an invalid cement class is provided.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.alpha_ds1(
            cement_class
        )
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_alpha_ds2(cement_class):
    """Return a coefficient depending on the cement class.

    EN 1992-1-1:2004, Sec. B.2.

    Args:
    cement_class (str): The cement class, either 'S', 'N' or 'R'.

    Returns:
    float: The exponent that depends on the cement type.

    Raises:
    ValueError: If an invalid cement class is provided.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.alpha_ds2(
            cement_class
        )
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_beta_H(h_0, fcm, RH, alpha_3):
    """Calculate the effect of relative humidity and the effective thickness
    of the structural element.

    EN 1992-1-1:2004, Eq. (B.8a and b).

    Args:
    h_0 (float): The effective cross sectional thickness in mm, Equation
        (B.6).
    fcm (float): The mean concrete strength in MPa.
    RH (float): The relative humidity in percent.
    alpha_3 (float): A factor describing the effect of concrete strength
        defined in Eq. B.8c.

    Returns:
    float: The effect of humidity and the effective thickness of the
    element.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.beta_H(
        h_0, fcm, RH, alpha_3
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_beta_RH(RH, RH_0):
    """Calculate the factor describing the effect of relative humidity.

    EN 1992-1-1:2004, Eq. (B.12).

    Args:
    RH (float): The relative humidity in percent.

    Keyword Args:
    RH_0 (float): The reference relative humidity, default: 100%.

    Returns:
    float: The factor taking into account the relative humidity.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.beta_RH(
            RH, RH_0
        )
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_beta_as(t):
    """Calculate the factor describing the development of autogenous
    shrinkage.

    EN 1992-1-1:2004, Eq. (3.13).

    Args:
    t (npt.ArrayLike): The age of the concrete in days.

    Returns:
    npt.ArrayLike: The factor describing the development of autogenous
    shrinkage.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.beta_as(t)
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_beta_c(t0, t, beta_H):
    """Calculate the factor that describes the creep development as a function
    of time after loading.

    EN 1992-1-1:2004, Eq. (B.7).

    Args:
    t0 (float): The concrete age in days a the time of loading.
    t (ArrayLike): The concrete age in days at the evaluated time.
    beta_H (float): Parameter defined in (B.8).

    Returns:
    float: Parameter defined by Equation (B.7), beta_c.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.beta_c(
        t0, t, beta_H
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_beta_ds(t, t_s, h_0):
    """Calculate the coefficient taking into account the time of drying.

    EN 1992-1-1:2004, Eq. (3.10).

    Args:
    t (npt.ArrayLike): The age of the concrete in days.
    t_s (float): The age of the concrete in days at the start of drying
        (normally this is the point in time when curing measures end).
    h_0 (float): The effective cross section thichkness in mm.

    Returns:
    npt.ArrayLike: The coefficient taking into account the time of drying.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.beta_ds(
            t, t_s, h_0
        )
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_beta_fcm(fcm):
    """Calculate the effect of the concrete strength on the standardized creep
    number.

    EN 1992-1-1:2004, Eq. (B.4).

    Args:
    fcm (float): The mean concrete strength in MPa.

    Returns:
    float: The effect of concrete strength.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.beta_fcm(
            fcm
        )
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_beta_t0(t0):
    """Calculate the effect of age at loading on the notional creep
    coefficient.

    EN 1992-1-1:2004, Eq. (B.5).

    Args:
    t0 (float): The age at loading in days.

    Returns:
    float: The effect of age at loading.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.beta_t0(
            t0
        )
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_eps_ca(beta_as, eps_ca_inf):
    """Calculate the autogenous shrinkage.

    EN 1992-1-1:2004, Eq. (3.11).

    Args:
    beta_as (npt.ArrayLike): A factor describing the autogenous shrinkage
        development.
    eps_ca_inf (float): The final autogenous shrinkage.

    Returns:
    npt.ArrayLike: The autogenous shrinkage.

    Note:
    In EC2 (2004), the shrinkage strain is calculated as a positive number.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.eps_ca(
        beta_as, eps_ca_inf
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_eps_ca_inf(fck):
    """Calculate the final autogenous shrinkage.

    EN 1992-1-1:2004, Eq. 3.12.

    Args:
    fck (float): The characteristic compressive strength in MPa.

    Returns:
    float: The final autogenous shrinkage.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.eps_ca_inf(
        fck
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_eps_cd(beta_ds, k_h, eps_cd_0):
    """Calculate the drying shrinkage.

    EN 1992-1-1:2004, Eq. (3.9).

    Args:
    beta_ds (npt.ArrayLike): A coefficient taking into account the time of
        drying defined in Eq. (3.10).
    k_h (float): A coefficient depending on the effective thickness of the
        section defined in Tab. 3.3.
    eps_cd_0 (float): The nominal value of drying shrinkage defined in Eq.
        (B.11).

    Returns:
    npt.ArrayLike: The drying shrinkage.

    Note:
    In EC2 (2004), the shrinkage strain is calculated as a positive number.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.eps_cd(
        beta_ds, k_h, eps_cd_0
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_eps_cd_0(
    alpha_ds1, alpha_ds2, fcm, beta_RH, fcm_0
):
    """Calculate the nominal value of drying shrinkage.

    EN 1992-1-1:2004, Eq. (B.11).

    Args:
    alpha_ds1 (float): A coefficient depending on the cement type, defined
        in EC2 (2004), Sec. B.2.
    alpha_ds2 (float): A coefficient depending on the cement type, defined
        in EC2 (2004), Sec. B.2.
    fcm (float): The mean compressive strength in MPa.
    beta_RH (float): A factor describing the effect of relative humidity,
        defined in Eq. (B.12).

    Keyword Args:
    fcm_0 (float): A reference strength in MPa, default 10 MPa.

    Returns:
    float: The nominal value of drying shrinkage.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.eps_cd_0(
            alpha_ds1, alpha_ds2, fcm, beta_RH, fcm_0
        )
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_eps_cs(eps_cd, eps_ca):
    """Calculate the total shrinkage strain.

    EN 1992-1-1:2004, Eq. (3.8).

    Args:
    eps_cd (npt.ArrayLike): The drying shrinkage defined in Eq. (3.9).
    eps_ca (npt.ArrayLike): The autogenous shrinkage defined in Eq. (3.11).

    Returns:
    npt.ArrayLike: The total shrinkage strain.

    Note:
    In EC2 (2004), the shrinkage strain is calculated as a positive number.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.eps_cs(
        eps_cd, eps_ca
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_h_0(Ac, u):
    """Calculate the effective thickness of the cross section.

    EN 1992-1-1:2004, Eq. (B.6).

    Args:
    Ac (float): The cross section area.
    u (float): The part of the circumference of the cross section subject
        to drying.

    Returns:
    float: The effective thickness.

    Note:
    The unit of the return will be consistent with the input. E.g. if Ac is
    mm ** 2 and u is mm, the return is mm.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.h_0(
        Ac, u
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_k_h(h_0):
    """Calculate the coefficient depending on the notional size.

    EN 1992-1-1:2004, Tab. 3.3.

    Args:
    h_0 (float): The notional size of the cross-section in mm.

    Returns:
    float: The coefficient depending on the notional size.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.k_h(
        h_0
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_phi(phi_0, beta_c):
    """Calculate the creep number.

    EN 1992-1-1:2004, Eq. (B.1).

    Args:
    phi_0 (float): The standardized creep number defined in Eq. B.2.
    beta_c (npt.ArrayLike): A factor taking into account the creep
        development as a function of time after loading defined in Eq.
        (B.7).

    Returns:
    float: The creep number.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.phi(
        phi_0, beta_c
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_phi_0(phi_RH, beta_fcm, beta_t0):
    """Calculate the standardized creep number.

    EN 1992-1-1:2004, Eq. (B.2).

    Args:
    phi_RH (float): The effect of relative humidity defined in Eq. B.3.
    beta_fcm (float): The effect of the concrete strength defined in Eq.
        B.4.
    beta_t0 (float): The effect of the age at loading defined in Eq. B.5.

    Returns:
    float: The standardized creep number.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.phi_0(
        phi_RH, beta_fcm, beta_t0
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_phi_RH(
    h_0, fcm, RH, alpha_1, alpha_2
):
    """Calculate the effect of relative humidity on the standardized creep
    number.

    EN 1992-1-1:2004, Eq. (B.3).

    Args:
    h_0 (float): The effective cross sectional thickness in mm, Equation
        (B.6).
    fcm (float): The mean concrete strength in MPa.
    RH (float): The relative humidity in percent.
    alpha_1 (float): A factor describing the effect of concrete strength
        defined in Eq. (B.8c).
    alpha_2 (float): A factor describing the effect of concrete strength
        defined in Eq. (B.8c).

    Returns:
    float: The calculation parameter (B.3).
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.phi_RH(
        h_0, fcm, RH, alpha_1, alpha_2
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_t0_adj(t0, alpha_cement):
    """Calculate the adjusted age of the concrete.

    EN 1992-1-1:2004, Eq. (B.9).

    Args:
    t0 (float): The concrete age in days at the time of loading.
    alpha_cement (float): Exponent derived from the sement type.

    Returns:
    float: The adjusted age of the concrete.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.t0_adj(
        t0, alpha_cement
    )


@xw.func
def SC_ec2_2004_concrete_creep_and_shrinkage_t_T(T, dt):
    """Calculate the maturity of the concrete.

    EN 1992-1-1:2004, Eq. (B.10).

    Args:
    T (npt.ArrayLike): The curing temperature history in degrees Celcius.
    dt (npt.ArrayLike): The number of days with temperature T.

    Returns:
    float: The maturity of the concrete.

    Note:
    The two arrays T and dt should have the same length. Each item in dt
    represents a time interval, and the corresponding item in T represents
    the average temperature in that time interval.
    """
    return structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage.t_T(
        T, dt
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_Ecm(fcm):
    """The secant modulus of concrete.

    EN 1992-1-1:2004, Table 3.1.

    Args:
    fcm (float): The mean compressive strength of concrete in MPa.

    Returns:
    float: The secant modulus of concrete in MPa.
    """
    return structuralcodes.codes.ec2_2004._concrete_material_properties.Ecm(
        fcm
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_eps_c1(fcm):
    """The strain at maximum compressive stress of concrete (fcm) for the
    Sargin constitutive law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
    fcm (float): The mean compressive strength of concrete in MPa.

    Returns:
    float: The strain at maximum compressive stress, absolute value, no
    unit.
    """
    return structuralcodes.codes.ec2_2004._concrete_material_properties.eps_c1(
        fcm
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_eps_c2(fck):
    """The strain at maximum compressive stress of concrete for the
    parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
    fck (float): The characteristic compressive strength of concrete in
        MPa.

    Returns:
    float: The strain at maximum compressive stress, absolute value, no
    unit.
    """
    return structuralcodes.codes.ec2_2004._concrete_material_properties.eps_c2(
        fck
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_eps_c3(fck):
    """The strain at maximum compressive stress of the bi-linear law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
    fck (float): The characteristic compressive strength of concrete in
        MPa.

    Returns:
    float: The strain at maximum compressive stress, absolute value, no
    unit.
    """
    return structuralcodes.codes.ec2_2004._concrete_material_properties.eps_c3(
        fck
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_eps_cu1(fck):
    """The ultimate strain for the Sargin constitutive law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
    fck (float): The characteristic compressive strength of concrete in
        MPa.

    Returns:
    float: The ultimate strain, absolute value, no unit.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_material_properties.eps_cu1(
            fck
        )
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_eps_cu2(fck):
    """The ultimate strain of the parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
    fck (float): The characteristic compressive strength of concrete in
        MPa.

    Returns:
    float: The ultimate strain, absolute value, no unit.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_material_properties.eps_cu2(
            fck
        )
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_eps_cu3(fck):
    """The ultimate strain of the bi-linear law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
    fck (float): The characteristic compressive strength of concrete in
        MPa.

    Returns:
    float: The ultimate strain, absolute value, no unit.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_material_properties.eps_cu3(
            fck
        )
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_fcd(fck, alpha_cc, gamma_c):
    """The design compressive strength of concrete.

    EN 1992-1-1:2004, Eq. (3.15)

    Args:
    fck (float): The characteristic compressive strength in MPa.
    alpha_cc (float): A factor for considering long-term effects on the
        strength, and effects that arise from the way the load is applied.
    gamma_c (float): The partial factor of concrete.

    Returns:
    float: The design compressive strength of concrete in MPa
    """
    return structuralcodes.codes.ec2_2004._concrete_material_properties.fcd(
        fck, alpha_cc, gamma_c
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_fcm(fck, delta_f):
    """The mean compressive strength of concrete.

    EN 1992-1-1:2004, Table 3.1.

    Args:
    fck (float): The characteristic compressive strength of concrete in
        MPa.

    Keyword Args:
    delta_f (float): The difference between the mean and the
        characteristic strength.

    Returns:
    float: The mean compressive strength in MPa.
    """
    return structuralcodes.codes.ec2_2004._concrete_material_properties.fcm(
        fck, delta_f
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_fctk_5(fctm):
    """The 5% fractile of the tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
    fctm (float): The mean tensile strength of concrete in MPa.

    Returns:
    float: The 5% fractile of the tensile strength in MPa.
    """
    return structuralcodes.codes.ec2_2004._concrete_material_properties.fctk_5(
        fctm
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_fctk_95(fctm):
    """The 95% fractile of the tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
    fctm (float): The mean tensile strength of concrete in MPa.

    Returns:
    float: The 95% fractile of the tensile strength in MPa.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_material_properties.fctk_95(
            fctm
        )
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_fctm(fck):
    """The mean tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
    fck (float): The characteristic compressive strength of concrete in
        MPa.

    Returns:
    float: The mean tensile strength in MPa.
    """
    return structuralcodes.codes.ec2_2004._concrete_material_properties.fctm(
        fck
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_k_sargin(Ecm, fcm, eps_c1):
    """Computation of k parameter for Sargin constitutive Law.

    EN 1992-1-1:2004, Eq. (3.14)

    Args:
    Ecm (float): the mean elastic modulus of concrete in MPa.
    fcm (float): the mean compressive strength in MPa.
    eps_c1 (float): the strain corresponding to peak stress.
    """
    return (
        structuralcodes.codes.ec2_2004._concrete_material_properties.k_sargin(
            Ecm, fcm, eps_c1
        )
    )


@xw.func
def SC_ec2_2004_concrete_material_properties_n_parabolic_rectangular(fck):
    """The exponent in the parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
    fck (float): The characteristic compressive strength of concrete in
        MPa.

    Returns:
    float: The exponent n, absolute value, no unit.
    """
    return structuralcodes.codes.ec2_2004._concrete_material_properties.n_parabolic_rectangular(
        fck
    )


@xw.func
def SC_ec2_2004_reinforcement_material_properties_epsud(epsuk, gamma_eps):
    """Calculate the design value of the reinforcement ultimate strain.

    EUROCDE 2 1992-1-1:2004, Fig 3.8

    Args:
    epsuk (float): The characteristic ultimate strain

    Keyword Args:
    gamma_eps (float): The partial factor specified in NA.
        Default value 0.9.

    Returns:
    float: The design ultimate strain

    Raises:
    ValueError: if epsuk is less than 0
    ValueError: if gamma_eps is greater than 1
    """
    return structuralcodes.codes.ec2_2004._reinforcement_material_properties.epsud(
        epsuk, gamma_eps
    )


@xw.func
def SC_ec2_2004_reinforcement_material_properties_fyd(fyk, gamma_s):
    """Calculate the design value of the reinforcement yield strength.

    EUROCODE 2 1992-1-1:2004, Fig. 3.8

    Args:
    fyk (float): The characteristic yield strength in MPa.
    gamma_s (float): The partial factor.

    Returns:
    float: The design yield strength in MPa.

    Raises:
    ValueError: if fyk is less than 0
    ValueError: if gamma_s is less than 1
    """
    return (
        structuralcodes.codes.ec2_2004._reinforcement_material_properties.fyd(
            fyk, gamma_s
        )
    )


@xw.func
def SC_ec2_2004_reinforcement_material_properties_reinforcement_duct_props(
    fyk, ductility_class
):
    """Return a dict with the minimum characteristic ductility properties for
    reinforcement ductility class.

    EUROCODE 2 1992-1-1:2004, Tab. C.1

    Args:
    fyk (float): The characteristic yield strength.
    ductility_class (Literal['A', 'B', 'C']): The reinforcement ductility
        class designation.

    Returns:
    Dict[str, float]: A dict with the characteristik strain value at the
    ultimate stress level (epsuk), and the characteristic ultimate stress
    (ftk).

    Raises:
    ValueError: when the ductility_class does not define a valid ductility
        class
    """
    return structuralcodes.codes.ec2_2004._reinforcement_material_properties.reinforcement_duct_props(
        fyk, ductility_class
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_As_min(
    A_ct, sigma_s, fct_eff, k, kc
):
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas.

    EUROCODE 2 1992-1-1:2004, Eq. (7.1).

    Args:
    A_ct (float): is the area of concrete within the tensile zone in mm2.
        The tensile zone is that part of the section which is calculated
        to be in tension just before the formation of the first crack.
    sigma_s (float): is the absolute value of the maximum stress in MPa
        permitted in the reinforcement immediately after the formation
        of the crack. This may be taken as theyield strength of the
        reinforcement, fyk. A lower value may, however, be needed to
        satisfy the crack width limits according to the maximum
        bar size of spacing (see 7.3.3 (2)).
    fct_eff (float): is the mean value of the tensile strength in MPa of
        the concrete effective at the time when the cracks may first be
        expected to occur: fct,eff=fct or lower (fct(t)), is cracking
        is expected earlier than 28 days.
    k (float): is the coefficient which allow for the effect of
        non-uniform self-equilibrating stresses, which lead to a
        reduction of restraint forces.
        k=1 for webs w<=300mm or flanges widths less than 300mm
        k=0.65 for webs w>=800mm or flanges with widths greater than 800mm
        Intermediate values may be interpolated.
    kc (float): is a coefficient which takes account of the stress
        distribution within the section immediately prior to cracking and
        the change of the lever arm.

    Returns:
    float: The minimum area of reinforcing steel within the tensile zone in
    mm2.

    Raises:
    ValueError: if k value is not between 0.65 and 1 or kc is not
        larger than 0 and lower than 1.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.As_min(
        A_ct, sigma_s, fct_eff, k, kc
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_As_min_2(
    wk, sigma_s, fct_eff, h_cr, h, d, delta_s, kc
):
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas.

    EUROCODE 2 1992-1-1:2004, Table (7.2N), Table (7.3N).

    Args:
    wk (float): the characteristic crack width value in mm.
    sigma_s (float): the steel stress value in MPa under the relevant
        combination of actions.
    fct_eff (float): is the mean value of the tensile strength in MPa of
        the concrete effective at the time when the cracks may first be
        expected to occur: fct,eff=fct or lower (fct(t)), is cracking is
        expected earlier than 28 days.
    h_cr (float): is the depth of the tensile zone immediately prior to
        cracking, considering the characteristic values of prestress and
        axial forces under the quasi-permanent combination of actions.
    h (float): the overall depth of the section in mm.
    d (float): is the effective depth to the centroid of the outer layer of
        the reinforcement.

    Keyword Args:
    delta_s (float, optional): value of prestressed stress in MPa if
        applicable.
    kc (float, optional): is a coefficient which takes account of the
        stress distribution within the section immediately prior to
        cracking and the change of the lever arm in a bending section. None
        for pure tensile uniform axial section.

    Returns:
    tuple(float, float): With the value of the maximum bar diameters in mm
    in the first position and the maximum bar spacing in mm in the second
    position.

    Raises:
    ValueError: If wk, fct_eff, h_cr, h or d are less than 0.
    ValueError: If kc is not between 0 and 1.
    ValueError: If combination of wk and stress values are out of scope.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.As_min_2(
        wk, sigma_s, fct_eff, h_cr, h, d, delta_s, kc
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_As_min_p(
    A_ct, sigma_s, fct_eff, k, kc, Ap, phi_s, phi_p, xi, delta_s
):
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas in addition with bonded tendons.

    EUROCODE 2 1992-1-1:2004, Eq. (7.1).

    Args:
    A_ct (float): is the area of concrete within the tensile zone in mm2.
        The tensile zone is that part of the section which is calculated to
        be in tension just before the formation of the first crack.
    sigma_s (float): is the absolute value of the maximum stress in MPa
        permitted in the reinforcement immediately after the formation of
        the crack. This may be taken as the yield strength of the
        reinforcement, fyk. A lower value may, however, be needed to
        satisfy the crack width limits according to the maximum bar size of
        spacing (see 7.3.3 (2)).
    fct_eff (float): is the mean value of the tensile strength in MPa of
        the concrete effective at the time when the cracks may first be
        expected to occur: fct,eff=fct or lower (fct(t)), is cracking is
        expected earlier than 28 days.
    k (float): is the coefficient which allow for the effect of non-
        uniform self-equilibrating stresses, which lead to a reduction of
        restraint forces. k=1 for webs w<=300mm or flanges widths less than
        300mm. k=0.65 for webs w>=800mm or flanges with widths greater than
        800mm. Intermediate values may be interpolated.
    kc (float): is a coefficient which takes account of the stress
        distribution within the section immediately prior to cracking and
        the change of the lever arm.
    Ap (float): is the area in mm2 of pre or post-tensioned tendons within
        ac_eff.
    phi_s (float): largest bar diameter in mm of reinforcing steel. Equal
        to 0 if only prestressing is used in control cracking.
    phi_p (float): equivalent diameter in mm of tendon according to 6.8.2.
    xi (float): ratio of bond strength of prestressing and reinforcing
        steel, according to Table 6.2 in 6.8.2.
    delta_s (float): stress variation in MPa in prestressing tendons from
        the state of zero strain of the concrete at the same level.

    Returns:
    float: The minimm area of reinforcing steel within the tensile zone in
    mm2.

    Raises:
    ValueError: If k value is not between 0.65 and 1 or kc is not larger
    than 0 and lower than 1. If diameters phi_s or phi_p are lower than 0.
    If ratio of bond xi strength is less than 0.15 or larger than 0.8. If
    stress variation incr_stress is less than 0.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.As_min_p(
        A_ct, sigma_s, fct_eff, k, kc, Ap, phi_s, phi_p, xi, delta_s
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_alpha_e(Es, Ecm):
    """Compute the ratio between the steel and mean concrete elastic modules.

    EUROCODE 2 1992-1-1:2004, Section 7.3.4-2

    Args:
    Es (float): Steel elastic modulus in MPa.
    Ecm (float): Concrete mean elastic modulus in MPa.

    Returns:
    float: Ratio between modules.

    Raises:
    ValueError: If any of es or ecm is lower than 0.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.alpha_e(
        Es, Ecm
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_eps_sm_eps_cm(
    sigma_s, alpha_e, rho_p_eff, kt, fct_eff, Es
):
    """Returns the strain difference (epsilon_sm - epsilon_cm) needed to
    compute the crack width. esm is the mean strain in the reinforcement under
    the relevant combination of loads of imposed deformations and taking into
    account the effects of tension stiffening. Only the additional tensile
    strain beyond the state of zero strain of the concrete is considered.
    epsilon_cm is the mean strain in the concrete between the cracks.

    EUROCODE 2 1992-1-1:2004, Eq. (7.9).

    Args:
    sigma_s (float): Is the stress in MPa in the tension reinforcement
        assuming a cracked section. For pretensioned members, s_steel may
        be replaced by increment of s_steel stress variation in
        prestressing tendons from the state of zero strain of the concrete
        at the same level.
    alpha_e (float): Is the ratio Es/Ecm.
    rho_p_eff (float): Effective bond ratio between areas given by Eq.
        (7.10).
    kt (float): Is a factor dependent on the load duration.
    fct_eff (float): Is the mean value of the tensile strength in MPa of
        the concrete effective at the time when the cracks may first be
        expected to occur: fct_eff=fctm or fctm(t) if crack is expected
        earlier than 28 days.
    Es (float): Steel elastic modulus in MPa.

    Returns:
    float: The strain difference between concrete and steel.

    Raises:
    ValueError: If any sigma_s, alpha_e, rho_p_eff, fct_eff or Es is less
        than 0.
    ValueError: if kt is not 0.6 and not 0.4.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.eps_sm_eps_cm(
        sigma_s, alpha_e, rho_p_eff, kt, fct_eff, Es
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_hc_eff(h, d, x):
    """Returns the effective height of concrete in tension surrounding the
    reinforcement or prestressing tendons.

    EUROCODE 2 1992-1-1:2004, Section (7.3.2-3).

    Args:
    h (float): total depth of the element in mm.
    d (float): distance in mm to the level of the steel centroid.
    x (float): distance in mm to the zero tensile stress line.

    Returns:
    float: The effective height in mm.

    Raises:
    ValueError: If any of h, d or x is lower than zero.
    ValueError: If d is greater than h.
    ValueError: If x is greater than h.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.hc_eff(
        h, d, x
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_k(h):
    """Is the coefficient which allow for the effect of
    non-uniform self-equilibrating stresses, which lead to a
    reduction of restraint forces.
    k=1 for webs w<=300mm or flanges widths less than 300mm
    k=0.65 for webs w>=800mm or flanges with widths greater than 800mm.

    EUROCODE 2 1992-1-1:2004, Eq. (7.1).

    Args:
    h (float): flange length or flange width in mm

    Returns:
    float: k coefficient value.

    Raises:
    ValueError: if h is less than 0
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.k(h)


@xw.func
def SC_ec2_2004_section_7_3_crack_control_k1(bond_type):
    """Get the k1 coefficient which takes account of the bond properties of the
    bounded reinforcement.

    EUROCODE 2 1992-1-1:2004, Eq. (7.11-k1).

    Args:
    bond_type (str): The bond property of the reinforcement. High bond bars
        (bond), or bars with an effectively plain surface (plain).

    Returns:
    float: Value of the k1 coefficient.

    Raises:
    ValueError: If bond_type is neither 'bond' nor 'plain'.
    TypeError: If bond_type is not an str.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.k1(
        bond_type
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_k2(eps_r):
    """Computes a coeff. which takes into account the distribution of strain.

    EUROCODE 2 1992-1-1:2004, Eq. (7.13).

    Args:
    eps_r (float): ratio epsilon_2/epsilon_1 where epsilon_1 is the greater
        and epsilon_2 is the lesser strain at the boundaries of the section
        considered, assessed on the basis of a cracked section. epsilon_r=0
        for bending and epsilon_r=1 for pure tension.

    Returns:
    float: The k2 coefficient value.

    Raises:
    ValueError: If eps_r is not between 0 and 1.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.k2(eps_r)


@xw.func
def SC_ec2_2004_section_7_3_crack_control_k3():
    """Returns the k3 coefficient for computing sr_max.

    Returns:
    float: Value for the coefficient.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.k3()


@xw.func
def SC_ec2_2004_section_7_3_crack_control_k4():
    """Returns the k4 coefficient for computing sr_max.

    Returns:
    float: Value for the coefficient.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.k4()


@xw.func
def SC_ec2_2004_section_7_3_crack_control_kc_flanges_area(f_cr, A_ct, fct_eff):
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm for bending+axial combination
    in rectangular sections for flanges of box sections and T-sections.

    EUROCODE 2 1992-1-1:2004, Eq. (7.3).

    Args:
    f_cr: is the absolute value in kN of the tensile force within the
        flange immediately prior to cracking due to cracking moment
        calculated with fct,eff.
    A_ct (float): is the area of concrete within the tensile zone in mm2.
        The tensile zone is that part of the section which is calculated
        to be in tension just before the formation of the first crack.
    fct_eff (float): is the mean value of the tensile strength in MPa of
        the concrete effective at the time when the cracks may first be
        expected to occur: fct,eff=fct or lower (fct(t)), is cracking
        is expected earlier than 28 days.

    Returns:
    float: Value of the kc coefficient.

    Raises:
    ValueError: If A_ct is less than 0mm2.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.kc_flanges_area(
        f_cr, A_ct, fct_eff
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_kc_rect_area(h, b, fct_eff, N_ed):
    """Computes the coefficient which takes account of the stress distribution
    within the section immediately prior to cracking and the change of the
    lever arm for bending+axial combination in rectangular sections and webs of
    box sections and T-sections.

    EUROCODE 2 1992-1-1:2004, Eq. (7.2).

    Args:
    h (float): heigth of the element in mm
    b (float): width of the element in mm
    fct_eff (float): is the mean value of the tensile strength in MPa of
        the concrete effective at the time when the cracks may first be
        expected to occur: fct,eff=fct or lower (fct(t)), is cracking is
        expected earlier than 28 days.
    N_ed (str): axial force at the serviceability limit state acting on the
        part of the cross-section under consideration (compressive force
        positive). n_ed should be determined considering the characteristic
        values of prestress and axial forces under the relevant combination
        of actions.

    Returns:
    float: Value of the kc coefficient.

    Raises:
    ValueError: If h or b are less than 0.
    """
    return (
        structuralcodes.codes.ec2_2004._section_7_3_crack_control.kc_rect_area(
            h, b, fct_eff, N_ed
        )
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_kc_tension():
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm in pure tension.

    EUROCODE 2 1992-1-1:2004, Eq. (7.1).

    Returns:
    float: Value of the kc coefficient in pure tension.
    """
    return (
        structuralcodes.codes.ec2_2004._section_7_3_crack_control.kc_tension()
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_kt(load_type):
    """Returns the kt factor dependent on the load duration for the crack width
    calculation.

    Args:
    load_type (str): The load type, 'short' for term loading, 'long' for
        long term loading.

    Returns:
    float: With the kt factor.

    Raises:
    ValueError: If load_type is not 'short' and not 'long'.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.kt(
        load_type
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_phi_eq(n1, n2, phi1, phi2):
    """Computes the equivalent diameter. For a section with n1 bars of diameter
    phi1 and n2 bars of diameter phi2.

    EUROCODE 2 1992-1-1:2004, Sect. (7.12).

    Args:
    n1 (int): Number of bars with diameter phi1.
    n2 (int): Number of bars with diameter phi2.
    phi1 (float): Diameter of n1 bars in mm.
    phi2 (float): Diamater of n2 bars in mm.

    Returns:
    float: The equivalent diameter in mm.

    Raises:
    ValueError: If any of n1 or n2 is less than 0.
    ValueError: If any of phi1 or phi2 is less than 0.
    TypeError: If any of n1 or n2 is not an integer.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.phi_eq(
        n1, n2, phi1, phi2
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_rho_p_eff(As_, xi1, Ap, Ac_eff):
    """Effective bond ratio between areas.

    EUROCODE 2 1992-1-1:2004, Eq. (7.10).

    Args:
    As (float): Steel area in mm2.
    xi1 (float): The adjusted ratio of bond according to expression (7.5).
    Ap (float): The area in mm2 of post-tensioned tendons in ac_eff.
    Ac_eff (float): Effective area of concrete in tension surrounding the
        reinforcement or prestressing tendons of depth hc_eff.

    Returns:
    float: With the retio between areas.

    Raises:
    ValueError: If any of As, xi1, Ap or Ac_eff is less than 0.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.rho_p_eff(
        As_, xi1, Ap, Ac_eff
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_sr_max_close(
    c, phi, rho_p_eff, k1, k2, k3, k4
):
    """Computes the maximum crack spacing in cases where bonded reinforcement
    is fixed at reasonably close centres within the tension zone
    (w_spacing<=5(c+phi/2)).

    EUROCODE 2 1992-1-1:2004, Eq. (7.11).

    Args:
    c (float): Is the cover in mm of the longitudinal reinforcement.
    phi (float): Is the bar diameter in mm. Where mixed bar diameters used,
        then it should be replaced for an equivalent bar diameter.
    rho_p_eff (float): Effective bond ratio between areas given by Eq.
        (7.10).
    k1 (float): Coefficient that takes into account the bound properties
        of the bonded reinforcement.
    k2 (float): Coefficient that takes into account the distribution of of
        the strain.

    Keyword Args:
    k3 (float, optional): Coefficient from the National Annex. If not
        specified then k3=3.4.
    k4 (float): Coefficient from the National Annex. If not specified then
        k4=0.425.

    Returns:
    float: The maximum crack spaing in mm.

    Raises:
    ValueError: If one or more of c, phi, rho_p_eff, k3 or k4
        is lower than zero.
    ValueError: If k1 is not 0.8 or 1.6.
    ValueError: If k2 is not between 0.5 and 1.0.
    """
    return (
        structuralcodes.codes.ec2_2004._section_7_3_crack_control.sr_max_close(
            c, phi, rho_p_eff, k1, k2, k3, k4
        )
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_sr_max_far(h, x):
    """Computes the maximum crack spacing in cases where bonded reinforcement
    exceeds (w_spacing>5(c+phi/2)) or where there is no bonded reinforcement at
    all.

    EUROCODE 2 1992-1-1:2004, Eq. (7.14).

    Args:
    h (float): Total depth of the beam in mm.
    x (float): Distance to non tension area of the element mm.

    Returns:
    float: Maximum crack spacing in mm.

    Raises:
    ValueError: If one of h or x is less than zero.
    ValueError: If x is greater than h.
    """
    return (
        structuralcodes.codes.ec2_2004._section_7_3_crack_control.sr_max_far(
            h, x
        )
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_sr_max_theta(
    sr_max_y, sr_max_z, theta
):
    """Computes the crack spacing sr_max when there is an angle between the
    angle of principal stress and the direction of the reinforcement, for
    members in two orthogonal directions, that is significant (> 15 degrees).

    EUROCODE 2 1992-1-1:2004, Eq. (7.15).

    Args:
    sr_max_y (float): Crack spacing in mm in the y-direction.
    sr_max_z (float): Crack spacing in mm in the z-direction.
    theta (float): Angle in radians between the reinforcement in the
        y-direction and the direction of the principal tensile stress.

    Returns:
    float: The crack spacing in mm.

    Raises:
    ValueError: If sr_max_y or sr_max_z is negative.
    ValueError: If theta is not between 0 and pi/2.
    """
    return (
        structuralcodes.codes.ec2_2004._section_7_3_crack_control.sr_max_theta(
            sr_max_y, sr_max_z, theta
        )
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_w_max(
    exposure_class, load_combination
):
    """Computes the recommended value of the maximum crack width.

    EUROCODE 2 1992-1-1:2004, Table (7.1N).

    Args:
    exposure_class (str): The exposure class. Possible values: X0, XC1,
        XC2, XC3, XC4, XD1, XD2, XS1, XS2, XS3.
    load_combination (str): The characteristic of the load combination.
        Frequent (f), or quasi-permanent (qp).

    Returns:
    float: The maximum recommended value for the crack width wmax in mm.

    Raises:
    ValueError: if not valid exposure_class or load_combination values.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.w_max(
        exposure_class, load_combination
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_w_spacing(c, phi):
    """Computes the distance threshold from which the maximum crack spacing is
    constant.

    EUROCODE 2 1992-1-1:2004, Sect. (7.3.4-3).

    Args:
    c (float): Cover of the longitudinal reinforcement in mm.
    phi (float): Is the bar diameter in mm. Where mixed bar diameters used,
        then it should be replaced for an equivalent bar diameter.

    Returns:
    float: Threshold distance in mm.

    Raises:
    ValueError: If any of c or phi is less than 0.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.w_spacing(
        c, phi
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_wk(sr_max, eps_sm_eps_cm):
    """Computes the crack width.

    EUROCODE 2 1992-1-1:2004, Eq. (7.8).

    Args:
    sr_max (float): The maximum crack length spacing in mm.
    eps_sm_eps_cm (float): the difference between the mean strain in the
        reinforcement under relevant combination of loads, including the
        effect of imposed deformations and taking into account tension
        stiffening and the mean strain in the concrete between cracks.

    Returns:
    float: Crack width in mm.

    Raises:
    ValueError: If any of sr_max or esm_ecm is less than zero.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.wk(
        sr_max, eps_sm_eps_cm
    )


@xw.func
def SC_ec2_2004_section_7_3_crack_control_xi1(xi, phi_p, phi_s):
    """Computes the adjusted ratio of bond strength taking into account
    the different diameters of prestressing and reinforcing steel.

    EUROCODE 2 1992-1-1:2004, Eq. (7.5).

    Args:
    xi (float): ratio of bond strength of prestressing and reinforcing
        steel, according to Table 6.2 in 6.8.2.
    phi_p (float): largest bar diameter in mm of reinforcing steel.
        Equal to 0 if only prestressing is used in control cracking.
    phi_s (float): equivalent diameter in mm of tendon acoording
        to 6.8.2.

    Returns:
    float: With the value of the ratio.

    Raises:
    ValueError: If diameters phi_s or phi_p are lower than 0. If ratio of
        bond strength xi is less than 0.15 or larger than 0.8.
    """
    return structuralcodes.codes.ec2_2004._section_7_3_crack_control.xi1(
        xi, phi_p, phi_s
    )


@xw.func
def SC_ec2_2004shear_Asw_max(fcd, fck, bw, s, fywd, NEd, Ac, alpha):
    """Calculate the maximum cross-sectional area of the shear reinforcement,
    based on the assumption 1/tan(theta) == 1.

    EN 1992-1-1 (2005). Eq. (6.13)

    Args:
    fcd (float): The design strength of the concrete in MPa.
    fck (float): The characteristic compressive strength in MPa.
    bw (float): The smallest width of the cross-section in tension in mm.
    s (float): The centre-to-centre distance of the shear reinforcement in
        mm.
    fwyd (float): The design strength of the shear reinforcement steel in
        MPa.
    NEd (float): The normal force in the cross-section due to loading or
        prestress (NEd > 0 for compression) in N.
    Ac (float): The cross-sectional area of the concrete in mm2.

    Keyword Args:
    alpha (float): The angle of the shear reinforcement with respect to the
        neutral axis in degrees. Default value = 90 degrees.

    Returns:
    float: The maximum allowable cross-sectional area of the shear
    reinforcement in mm2.

    Raises:
    ValueError: When sigma_cp > fcd.
    """
    return structuralcodes.codes.ec2_2004.shear.Asw_max(
        fcd, fck, bw, s, fywd, NEd, Ac, alpha
    )


@xw.func
def SC_ec2_2004shear_VEdmax_unreinf(bw, d, fck, fcd):
    """Calculate the maximum allowable shear force for cross-sections without
    shear reinforcement.

    En 1992-1-1 (2005), Eq. (6.5).

    Args:
    bw (float): The smallest width of the cross-section in tension in mm.
    d (float): The effective depth of the cross-section in mm.
    fck (float): The characteristic compressive strength in MPa.
    fcd (float): The design compressive strength in MPa.

    Returns:
    float: The maximum allowable shear force in the cross-section in N.
    When a reduced shear force may be considered for the calculations, the
    unreduced shear force has to comply to this value.
    """
    return structuralcodes.codes.ec2_2004.shear.VEdmax_unreinf(bw, d, fck, fcd)


@xw.func
def SC_ec2_2004shear_VRdc(fck, d, Asl, bw, NEd, Ac, fcd, k1, gamma_c, CRdc):
    """Compute the design strength of the shear resistance.

    EN 1992-1-1 (2005), Eq. (6.2)

    Args:
    fck (float): The characteristic compressive strength in MPa.
    d (float): The effective depth of the cross-section in mm.
    Asl (float): The cross-sectional area of the tensile reinforcement,
        anchored atleast (lbd + d) beyond the considered cross-section, in
        mm2.
    bw (float): The smallest width of the cross-section in tension in mm.
    NEd (float): The normal force in the cross-section due to loading or
        prestress (NEd > 0 for compression) in N.
    Ac (float): The cross-sectional area of the concrete in mm2.
    fcd (float): The design compressive strength in MPa.

    Keyword Args:
    k1 (float): Factor used to include the effect of the normal stress
        into the shear resistance of the concrete. Default value = 0.15,
        value might differ between National Annexes.
    gamma_c (float): Partial factor for concrete. Default value = 1.5,
        value might differ between National Annexes.
    CRdc (Optional[float]): Scaling factor for the shear resistance.
        Default value is 0.18 / gamma_c.

    Returns:
    float: The concrete shear resistance in MPa.
    """
    return structuralcodes.codes.ec2_2004.shear.VRdc(
        fck, d, Asl, bw, NEd, Ac, fcd, k1, gamma_c, CRdc
    )


@xw.func
def SC_ec2_2004shear_VRdc_prin_stress(Iy, bw, S, fctd, NEd, Ac, L_x, L_pt2):
    """Calculate the shear resistance in uncracked, prestressed elements
    without shear reinforcement, value is determined via Mohr's circle.

    The maximal value of the principle tensile stress does no necessarily lay
    at the centre of gravity. If this is the ase the minimum value of the shear
    resistance and corresponding stress needs to be found at the relevant
    location.

    EN 1992-1-1 (2005), Eq. (6.4).

    Args:
    Iy (float): The second moment of area of the considered cross-section
        in mm4.
    bw (float): The width of the cross-section at the centre of gravity.
    S (float): The first moment of area of the considered cross-section of
        the part above the centre of gravity, and with respect to the
        centre of gravity in mm3.
    fctd (float): Design value of the tensile strength of the concrete.
    NEd (float): The normal force in the cross-section due to loading or
        prestress (NEd > 0 for compression) in N.
    Ac (float): The cross-sectional area of the concrete in mm2.

    Keyword Args:
    L_x (float): Distance from the considered cross-section until the
        starting point of the transference length of the prestress steel.
        This value should be provided when the prestressing steel is
        prestreched. Default value is None.
    L_pt2 (float): Maximum value of the transference length of the
        prestress steel, according to Eq. (8.18). This value should be
        provided when the prestressing steel is prestreched. Default value
        is None.

    Returns:
    float: The maximum allowable shear force in N for an uncracked,
    prestressed element without shear reinforcement, determined from
    maximum allowable principle stress.
    """
    return structuralcodes.codes.ec2_2004.shear.VRdc_prin_stress(
        Iy, bw, S, fctd, NEd, Ac, L_x, L_pt2
    )


@xw.func
def SC_ec2_2004shear_VRdmax(bw, z, fck, theta, NEd, Ac, fcd, alpha, limit_fyd):
    """Calculate the maximum shear strength of the compression strut.

    EN 1992-1-1 (2005). Eq. (6.9)

    Args:
    bw (float): The smallest width of the cross-section in tension in mm.
    z (float): The inner lever arm of internal forces in mm.
    fck (float): The characteristic compressive strength in MPa.
    theta (float): The angle of the compression strut in degrees.
    NEd (float): The normal force in the cross-section due to loading or
        prestress (NEd > 0 for compression) in N.
    Ac (float): The cross-sectional area of the concrete in mm2.
    fcd (float): The design compressive strength in MPa.

    Keyword Args:
    alpha (float): The angle of the shear reinforcement with respect to the
        neutral axis in degrees. Default value = 90 degrees.
    limit_fyd (bool): Flag to indicate if the design yield stress is
        limited to 0.8 * fyk or not. This controls whether the stress
        reduction factor of concrete is given by Eq. (6.6) (False) or
        (6.10) (True).

    Returns:
    float: The shear strength of the shear reinforcement in N.

    Raises:
    ValueError: When theta < 21.8 degrees or theta > 45 degrees.
    ValueError: When sigma_cp > fcd.
    """
    return structuralcodes.codes.ec2_2004.shear.VRdmax(
        bw, z, fck, theta, NEd, Ac, fcd, alpha, limit_fyd
    )


@xw.func
def SC_ec2_2004shear_VRds(Asw, s, z, theta, fyk, alpha, gamma_s):
    """Calculate the shear resistance of vertical shear reinforcement.

    EN 1992-1-1 (2005). Eq. (6.8)

    Args:
    Asw (float): the cross-sectional area of the shear reinforcement in
        mm2.
    s (float): The centre-to-centre distance of the shear reinforcement in
        mm.
    z (float): The inner lever arm of internal forces in mm.
    theta (float): The angle of the compression strut in degrees.
    fyk (float): The characteristic strength of the reinforcement steel in
        MPa.

    Keyword Args:
    alpha (float): The angle of the shear reinforcement with respect to the
        neutral axis in degrees. Default value = 90 degrees.
    gamma_s (float): Partial factor of the reinforcement steel. Default
        value = 1.15. Value might differ between National Annexes.

    Returns:
    float: The shear resistance of the shear reinforcement in N.

    Raises:
    ValueError: When theta < 21.8 degrees or theta > 45 degrees.
    """
    return structuralcodes.codes.ec2_2004.shear.VRds(
        Asw, s, z, theta, fyk, alpha, gamma_s
    )


@xw.func
def SC_ec2_2004shear__alpha_l(L_x, L_pt2):
    """Compute the relative anchorage length for prestreched prestressing
    steel.

    Defined in EN 1992-1-1 (2005), Eq. (6.4).

    Args:
    L_x (float): Distance from the considered cross-section until the
        starting point of the transference length of the prestress steel.
    L_pt2 (float): Maximum value of the transference length of the
        prestress steel, according to Eq. (8.18).

    Returns:
    float: Fraction (relative anchorage length) for determining the amount
    of prestress that may be used when determining the shear resistance
    using Mohr's circle.
    """
    return structuralcodes.codes.ec2_2004.shear._alpha_l(L_x, L_pt2)


@xw.func
def SC_ec2_2004shear__k(d):
    """Compute a correction factor.

    Defined in EN 1992-1-1 (2005), Eq. (6.2).

    Args:
    d (float): The effective depth of the cross-section in mm.

    Returns:
    float: Correction factor to account for the cross-sectional size on the
    shear resistance.
    """
    return structuralcodes.codes.ec2_2004.shear._k(d)


@xw.func
def SC_ec2_2004shear__rho_L(Asl, bw, d):
    """Compute the longitudinal reinforcement ratio.

    Defined in EN 1992-1-1 (2005), Eq. (6.2).

    Args:
    Asl (float): The cross-sectional area of the tensile reinforcement,
        anchored at least (lbd + d) beyond the considered cross-section, in
        mm2.
    bw (float): The smallest width of the cross-section in tension in mm.
    d (float): The effective depth of the cross-section in mm.

    Returns:
    float: The maximum allowable reinforcement ratio of the longitudinal
    reinforcement, unitless.
    """
    return structuralcodes.codes.ec2_2004.shear._rho_L(Asl, bw, d)


@xw.func
def SC_ec2_2004shear__sigma_cp(NEd, Ac, fcd):
    """Calculate the average prestress stress in the cross-section.

    Defined in EN 1992-1-1 (2005), Eq. (6.2).

    Args:
    NEd (float): The normal force in the cross-section due to loading or
        prestress (NEd > 0 for compression) in N.
    Ac (float): The cross-sectional area of the concrete in mm2.
    fcd (float): The design compressive strength in MPa.

    Returns:
    float: The maximum allowable average prestress in the cross-section in
    MPa.
    """
    return structuralcodes.codes.ec2_2004.shear._sigma_cp(NEd, Ac, fcd)


@xw.func
def SC_ec2_2004shear__theta(theta, cot_min, cot_max):
    """Check if the provided angle theta is within the bounds provided by the
    code.

    EN 1992-1-1 (2005). Eq. (6.7N)

    Args:
    theta (float): The chosen angle of the compression strut in degrees.

    Keyword Args:
    cot_min (float): The minimum value for cot(theta). Default value is
        1.0. Different value might be provided in the National Annexes.
    cot_max (float): The maximum value for cot(theta). Default value is
        2.5. Different value might be provided in the National Annexes.

    Raises:
    ValueError if the chosen angle is not within the given bounds.
    """
    return structuralcodes.codes.ec2_2004.shear._theta(theta, cot_min, cot_max)


@xw.func
def SC_ec2_2004shear_alpha_cw(Ned, Ac, fcd):
    """Calculate factor that affects the maximum shear resistance of the
    concrete based on the prestress.

    EN 1992-1-1 (2005). Eq. (6.11N)

    Args:
    NEd (float): The normal force in the cross-section due to loading or
        prestress (NEd > 0 for compression) in N.
    Ac (float): The cross-sectional area of the concrete in mm2.
    fcd (float): The design strength of the concrete in MPa.

    Returns:
    float: Factor that affects the maximum shear resistance of the concrete
    based on the level of prestress.

    Raises:
    ValueError: The applied prestress exceeds the concrete design strength.
    """
    return structuralcodes.codes.ec2_2004.shear.alpha_cw(Ned, Ac, fcd)


@xw.func
def SC_ec2_2004shear_v(fck):
    """Calculate a strength redcution factor for concrete cracked by shear
    forces.

    EN 1992-1-1 (2005), Eq. (6.6N)

    Args:
    fck (float): The characteristic compressive strength in MPa.

    Returns:
    float: A concrete reduction factor to account for concrete cracked by
    shear forces.
    """
    return structuralcodes.codes.ec2_2004.shear.v(fck)


@xw.func
def SC_ec2_2004shear_v1(fck):
    """Calculate a strength redcution factor for concrete cracked by shear
    forces.

    EN 1992-1-1 (2005), Eq. (6.10N)

    Args:
    fck (float): The characteristic compressive strength in MPa.

    Returns:
    float: A concrete reduction factor to account for concrete cracked by
    shear forces.
    """
    return structuralcodes.codes.ec2_2004.shear.v1(fck)


@xw.func
def SC_ec2_2004shear_vmin(fck, d):
    """Compute the minimum shear resistance of the concrete.

    EN 1992-1-1 (2005), Eq. (6.3)

    Args:
    fck (float): The characteristic compressive strength in MPa.
    d (float): The effective depth of the cross-section in mm.

    Returns:
    float: The minimal shear stress resistance of the concrete in MPa.
    """
    return structuralcodes.codes.ec2_2004.shear.vmin(fck, d)


if __name__ == '__main__':
    xw.Book('pruebaxlwings.xlsm').set_mock_caller()
    main()

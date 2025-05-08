import structuralcodes.codes.ec2_2023
import xlwings as xw


def main():
    wb = xw.Book.caller()
    sheet = wb.sheets[0]
    if sheet['A1'].value == 'Hello xlwings!':
        sheet['A1'].value = 'Bye xlwings!'
    else:
        sheet['A1'].value = 'Hello xlwings!'


@xw.func
def ec2_2023_alpha_c(fcm_28: float):
    """Returns the coefficient to obtain the tangent modulus of elasticity from
    the secant modulus of elasticity.

    EN 1992-1-1, Table 5.1

    Args:
    fcm_28 (float): The characteristic compressive strength at 28 days
        in MPa.

    Returns:
    float: Coefficient alphac so that Ec=alphac*Ecm,28.
    """
    return structuralcodes.codes.ec2_2023.alpha_c(fcm_28)


@xw.func
def ec2_2023_fcm(fck: float, delta_f: float = 8.0):
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
    return structuralcodes.codes.ec2_2023.fcm(fck, delta_f)


@xw.func
def ec2_2023_fctm(fck: float):
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength.

    EN 1992-1-1:2023, Table 5.1.

    Args:
    fck (float): The characteristic compressive strength in MPa.

    Returns:
    float: The mean tensile strength in MPa.
    """
    return structuralcodes.codes.ec2_2023.fctm(fck)


@xw.func
def ec2_2023_fctk_5(fctm: float):
    """Compute the 5% mean concrete tensile strength fractile.

    EN 1992-1-1:2023, Table 5.1.

    Args:
    fctm (float): The mean concrete tensile strength in MPa.

    Returns:
    float: The 5% mean concrete tensile strength fractile in MPa.
    """
    return structuralcodes.codes.ec2_2023.fctk_5(fctm)


@xw.func
def ec2_2023_fctk_95(fctm: float):
    """Compute the 95% mean concrete tensile strength fractile.

    EN 1992-1-1:2023, Table 5.1.

    Args:
    fctm (float): The mean concrete tensile strength in MPa.

    Returns:
    float: The 5% mean concrete tensile strength fractile in MPa.
    """
    return structuralcodes.codes.ec2_2023.fctk_95(fctm)


@xw.func
def ec2_2023_Ecm(fcm: float, kE: float = 9500):
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
    return structuralcodes.codes.ec2_2023.Ecm(fcm, kE)


@xw.func
def ec2_2023_hn(Ac: float, u: float):
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
    return structuralcodes.codes.ec2_2023.hn(Ac, u)


@xw.func
def ec2_2023_phi_correction_factor(fck: float, A_exponent: float):
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
    return structuralcodes.codes.ec2_2023.phi_correction_factor(
        fck, A_exponent
    )


@xw.func
def ec2_2023_eta_cc(fck: float, fck_ref: float = 40):
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
    return structuralcodes.codes.ec2_2023.eta_cc(fck, fck_ref)


@xw.func
def ec2_2023_fcd(fck: float, eta_cc: float, k_tc: float, gamma_c: float):
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
    return structuralcodes.codes.ec2_2023.fcd(fck, eta_cc, k_tc, gamma_c)


@xw.func
def ec2_2023_fctd(fctk_5: float, k_tt: float, gamma_c: float):
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
    return structuralcodes.codes.ec2_2023.fctd(fctk_5, k_tt, gamma_c)


@xw.func
def ec2_2023_eps_c1(fcm: float):
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
    return structuralcodes.codes.ec2_2023.eps_c1(fcm)


@xw.func
def ec2_2023_eps_cu1(fcm: float):
    """Computes the strain at concrete failure of concrete.

    EN 1992-1-1:2023, Eq. (5.10).

    Args:
    fcm (float): The mean strength of concrete in MPa.

    Returns:
    float: The maximum strength at failure of concrete.

    Raises:
    ValueError: If fcm is less than 12+8MPa.
    """
    return structuralcodes.codes.ec2_2023.eps_cu1(fcm)


@xw.func
def ec2_2023_k_sargin(
    Ecm: float,
    fcm: float,
    eps_c1: float,
):
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
    return structuralcodes.codes.ec2_2023.k_sargin(
        Ecm,
        fcm,
        eps_c1,
    )


@xw.func
def ec2_2023_eps_c2():
    """The strain at maximum compressive stress of concrete for the
    parabolic-rectangular law.

    EN 1992-1-1:2023, Eq. 8.4

    Returns:
    float: The strain at maximum compressive stress, absolute value, no
    unit.
    """
    return structuralcodes.codes.ec2_2023.eps_c2()


@xw.func
def ec2_2023_eps_cu2():
    """The ultimate strain of the parabolic-rectangular law.

    EN 1992-1-1:2023, Eq. 8.4

    Returns:
    float: The ultimate strain, absolute value, no unit.
    """
    return structuralcodes.codes.ec2_2023.eps_cu2()


@xw.func
def ec2_2023_n_parabolic_rectangular():
    """The exponent in the parabolic-rectangular law.

    EN 1992-1-1:2023, Eq. 8.4
    Returns:
    float: The exponent n, absolute value, no unit.
    """
    return structuralcodes.codes.ec2_2023.n_parabolic_rectangular()


@xw.func
def ec2_2023_alpha_c_th():
    """Returns the linear coefficient of thermal expansion in 1/C� for
    concrete.

    EN 1992-1-1:2023, 5.1.6-6.

    Returns:
    float: The linear coefficient of thermal expansion in 1/C� for
    concrete.
    """
    return structuralcodes.codes.ec2_2023.alpha_c_th()


@xw.func
def ec2_2023_Es():
    """Returns the value of the modulus of elasticity for weldable reinforcing
    steel.

    EN 1992-1-1:2023, 5.2.4-3.

    Returns:
    float: Modulus of elasticity in MPa.
    """
    return structuralcodes.codes.ec2_2023.Es()


@xw.func
def ec2_2023_alpha_s_th():
    """Returns the linear coefficient of thermal expansion in 1/C� for weldable
    reinforced steel.

    EN 1992-1-1:2023, 5.2.4-5

    Returns:
    float: The linear coefficient of thermal expansion in 1/C� for weldable
    reinforcement steel.
    """
    return structuralcodes.codes.ec2_2023.alpha_s_th()


@xw.func
def ec2_2023_weight_s():
    """Returns the mean unit weight of reinforced steel for the purposes of
    design in kN/m3.

    EN 1992-1-1:2023.2.4-4.

    Returns:
    float: The mean unit weight in kN/m3.
    """
    return structuralcodes.codes.ec2_2023.weight_s()


@xw.func
def ec2_2023_fyd(fyk: float, gamma_s: float):
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
    return structuralcodes.codes.ec2_2023.fyd(fyk, gamma_s)


@xw.func
def ec2_2023_eps_ud(eps_uk: float, gamma_s: float):
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
    return structuralcodes.codes.ec2_2023.eps_ud(eps_uk, gamma_s)


@xw.func
def ec2_2023_sigma_s(
    eps: float, fy: float, k: float, eps_u: float, Es: float = 200000
):
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
    return structuralcodes.codes.ec2_2023.sigma_s(eps, fy, k, eps_u, Es)


@xw.func
def ec2_2023_fpd(fp01k: float, gamma_P: float):
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
    return structuralcodes.codes.ec2_2023.fpd(fp01k, gamma_P)


@xw.func
def ec2_2023_sigma_p(
    eps: float,
    fpy: float,
    fpu: float,
    eps_u: float = 0.035,
    Ep: float = 190000,
):
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
    return structuralcodes.codes.ec2_2023.sigma_p(
        eps,
        fpy,
        fpu,
        eps_u,
        Ep,
    )


@xw.func
def ec2_2023_Ec_eff(fcm: float, phi: float, kE: float = 9500):
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
    return structuralcodes.codes.ec2_2023.Ec_eff(fcm, phi, kE)


@xw.func
def ec2_2023_As_min_y(
    NEd: float, b: float, h: float, fct_eff: float, fyk: float
):
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
    return structuralcodes.codes.ec2_2023.As_min_y(NEd, b, h, fct_eff, fyk)


@xw.func
def ec2_2023_kh(b: float, h: float):
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
    return structuralcodes.codes.ec2_2023.kh(b, h)


@xw.func
def ec2_2023_wk_cal2(
    kw: float, k_1_r: float, srm_cal: float, epssm_epscm: float
):
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
    return structuralcodes.codes.ec2_2023.wk_cal2(
        kw, k_1_r, srm_cal, epssm_epscm
    )


@xw.func
def ec2_2023_k_1_r(h: float, x: float, ay: float):
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
    return structuralcodes.codes.ec2_2023.k_1_r(h, x, ay)


@xw.func
def ec2_2023_epssm_epscm(
    sigma_s: float,
    kt: float,
    fct_eff: float,
    rho_eff: float,
    alphae: float,
    Es: float,
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
    return structuralcodes.codes.ec2_2023.epssm_epscm(
        sigma_s,
        kt,
        fct_eff,
        rho_eff,
        alphae,
        Es,
    )


@xw.func
def ec2_2023_kfl(h: float, xg: float, hceff: float):
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
    return structuralcodes.codes.ec2_2023.kfl(h, xg, hceff)


@xw.func
def ec2_2023_srm_cal(
    c: float,
    kfl_: float,
    kb: float,
    phi: float,
    rho_eff: float,
    kw: float,
    h,
    x: float,
):
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
    return structuralcodes.codes.ec2_2023.srm_cal(
        c,
        kfl_,
        kb,
        phi,
        rho_eff,
        kw,
        h,
        x,
    )


@xw.func
def ec2_2023_wk_cal(
    kw: float,
    h: float,
    xg: float,
    hc_eff: float,
    c: float,
    kb: float,
    phi: float,
    rho_eff: float,
    x: float,
    sigma_s: float,
    kt: float,
    fct_eff: float,
    alphae: float,
    Es: float,
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
    return structuralcodes.codes.ec2_2023.wk_cal(
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
def ec2_2023_delta_simpl(
    delta_loads: float,
    delta_shr: float,
    fck1: float,
    phi1: float,
    b1: float,
    h: float,
    d: float,
    As_1: float,
    Mk: float,
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
    return structuralcodes.codes.ec2_2023.delta_simpl(
        delta_loads,
        delta_shr,
        fck1,
        phi1,
        b1,
        h,
        d,
        As_1,
        Mk,
    )


if __name__ == '__main__':
    xw.Book('pruebaxlwings.xlsm').set_mock_caller()
    main()

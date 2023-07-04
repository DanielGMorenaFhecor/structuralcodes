import typing as t

import structuralcodes.codes.ec2_2004
import structuralcodes.codes.mc2010
import xlwings as xw


def main():
    wb = xw.Book.caller()
    sheet = wb.sheets[0]
    if sheet["A1"].value == "Hello xlwings!":
        sheet["A1"].value = "Bye xlwings!"
    else:
        sheet["A1"].value = "Hello xlwings!"


@xw.func
def ec2_2004_w_max(exposure_class: str, load_combination: str):
    """Computes the recomended value of the maximum crack width.

    EUROCODE 2 1992-1-1:2004, Table (7.1N)

    Args:
        exposure_class (str): The exposure class.
            Possible values: X0, XC1, XC2, XC3, XC4, XD1, XD2, XS1, XS2, XS3
        load_combination (str):
            - f: for frequent load combination
            - qp: for quasi-permanent load combination

    Returns:
        float: The maximum recommended value for the crack width wmax in mm.

    Raises:
        ValueError: if not valid exposure_class or load_combination values."""
    return structuralcodes.codes.ec2_2004.w_max(
        exposure_class, load_combination
    )


@xw.func
def ec2_2004_As_min(
    A_ct: float, sigma_s: float, fct_eff: float, k: float, kc: float
):
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

    Args:
        A_ct (float): is the area of concrete within the tensile zone in mm2.
            The tensile zone is that parg of the section which is calculated
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
        _k (float): is the coefficient which allow for the effect of
            non-uniform self-equilibrating stresses, which lead to a
            reduction of restraint forces. Use 'k_crack_min_steel_area'
            to compute it
            k=1 for webs w<=300mm or flanges widths less than 300mm
            k=0.65 for webs w>=800mm or flanges with widths greater than 800mm
            Intermediate values may be interpolated.
        kc (float): is a coefficient which takes account of the stress
            distribution within the section immediately prior to cracking and
            the change of the lever arm.

    Returns:
        float: the minimm area of reinforcing steel within the tensile
            zone in mm2.

    Raises:
        ValueError: if _k value is not between 0.65 and 1 or kc is not
            larger than 0 and lower than 1."""
    return structuralcodes.codes.ec2_2004.As_min(A_ct, sigma_s, fct_eff, k, kc)


@xw.func
def ec2_2004_k(h: float):
    """Is the coefficient which allow for the effect of
    non-uniform self-equilibrating stresses, which lead to a
    reduction of restraint forces.
    k=1 for webs w<=300mm or flanges widths less than 300mm
    k=0.65 for webs w>=800mm or flanges with widths greater than 800mm

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

    Args:
        h (float): flange length or flange width in mm

    Returns:
        float: k coefficient value

    Raises:
        ValueError: if h is less than 0"""
    return structuralcodes.codes.ec2_2004.k(h)


@xw.func
def ec2_2004_kc_tension():
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm in pure dtension.

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

    Returns:
        float: value of the kc coefficient in pure tension"""
    return structuralcodes.codes.ec2_2004.kc_tension()


@xw.func
def ec2_2004_kc_rect_area(h: float, b: float, fct_eff: float, N_ed: float):
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm for bending+axial combination
    in rectangular sections and webs of box sections and T-sections.

    EUROCODE 2 1992-1-1:2004, Eq. (7.2)

    Args:
        h (float): heigth of the element in mm
        b (float): width of the element in mm
        fct_eff (float): is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct,eff=fct or lower (fct(t)), is cracking
            is expected earlier than 28 days.
        N_ed (str): axial force at the serviceability limit state acting on
            the part of the cross-section under consideration (compressive
            force positive). n_ed should be determined considering the
            characteristic values of prestress and axial forces under the
            relevant combination of actions

    Returns:
        float: value of the kc coefficient

    Raises:
        ValueError: is h or b are less than 0"""
    return structuralcodes.codes.ec2_2004.kc_rect_area(h, b, fct_eff, N_ed)


@xw.func
def ec2_2004_kc_flanges_area(f_cr: float, A_ct: float, fct_eff: float):
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm for bending+axial combination
    in rectangular sections for flanges of box sections and T-sections.

    EUROCODE 2 1992-1-1:2004, Eq. (7.3)

    Args:
        f_cr: is the absolute value in kN of the tensile force within the
            flange immediately prior to cracking due to cracking moment
            calculated with fct,eff
        A_ct (float): is the area of concrete within the tensile zone in mm2.
            The tensile zone is that part of the section which is calculated
            to be in tension just before the formation of the first crack.
        fct_eff (float): is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct,eff=fct or lower (fct(t)), is cracking
            is expected earlier than 28 days.

    Returns:
        float: value of the kc coefficient

    Raises:
        ValueError: is A_ct is less than 0mm2"""
    return structuralcodes.codes.ec2_2004.kc_flanges_area(f_cr, A_ct, fct_eff)


@xw.func
def ec2_2004_xi1(xi: float, phi_p: float, phi_s: float):
    """Computes the adjusted ratio of bond strength taking into account
    the different diameters of prestressing and reinforcing steel.

    EUROCODE 2 1992-1-1:2004, Eq. (7.5)

    Args:
        xi (float): ratio of bond strength of prestressing and reinforcing
            steel, according to Table 6.2 in 6.8.2
        phi_p (float): largest bar diameter in mm of reinforcing steel.
            Equal to 0 if only prestressing is used in control cracking
        phi_s (float): equivalent diameter in mm of tendon acoording
            to 6.8.2

    Returns:
        float: with the value of the ratio

    Raises:
        ValueError: if diameters phi_s or phi_p are lower than 0.
            If ratio of bond strength xi is less than 0.15 or larger than 0.8."""
    return structuralcodes.codes.ec2_2004.xi1(xi, phi_p, phi_s)


@xw.func
def ec2_2004_hc_eff(h: float, d: float, x: float):
    """Returns the effective height of concrete in tension surrounding
    the reinforcement or prestressing tendons.

    EUROCODE 2 1992-1-1:2004, Section (7.3.2-3)

    Args:
        h (float): total depth of the element in mm
        d (float): distance in mm to the level of the steel centroid
        x (float): distance in mm to the zero tensile stress line

    Returns:
        float: the effective height in mm

    Raises:
        ValueError: if any of h, d or x is lower than zero.
        ValueError: if d is greater than h
        ValueError: if x is greater than h"""
    return structuralcodes.codes.ec2_2004.hc_eff(h, d, x)


@xw.func
def ec2_2004_As_min_p(
    A_ct: float,
    sigma_s: float,
    fct_eff: float,
    k: float,
    kc: float,
    Ap: float,
    phi_s: float,
    phi_p: float,
    xi: float,
    delta_s: float,
):
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas in addition with bonded tendons

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

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
        _k (float): is the coefficient which allow for the effect of
            non-uniform self-equilibrating stresses, which lead to a
            reduction of restraint forces. Use 'k_crack_min_steel_area'
            to compute it
            k=1 for webs w<=300mm or flanges widths less than 300mm
            k=0.65 for webs w>=800mm or flanges with widths greater than 800mm
            Intermediate values may be interpolated.
        kc (float): is a coefficient which takes account of the stress
            distribution within the section immediately prior to cracking and
            the change of the lever arm.
        Ap (float): is the area in mm2 of pre or post-tensioned tendons
            within ac_eff
        phi_s (float): largest bar diameter in mm of reinforcing steel.
            Equal to 0 if only prestressing is used in control cracking
        phi_p (float): equivalent diameter in mm of tendon acoording
            to 6.8.2
        chi (float): ratio of bond strength of prestressing and reinforcing
            steel, according to Table 6.2 in 6.8.2
        delta_s (float): stress variation in MPa in prestressing tendons
            from the state of zero strain of the concrete at the same level

    Returns:
        float: the minimm area of reinforcing steel within the tensile
            zone in mm2.

    Raises:
        ValueError: if _k value is not between 0.65 and 1 or kc is not
            larger than 0 and lower than 1. If diameters phi_s or
            phi_p are lower than 0. If ratio of bond xi strength e
            is less than 0.15 or larger than 0.8.
            Is stress variation incr_stress is less than 0."""
    return structuralcodes.codes.ec2_2004.As_min_p(
        A_ct,
        sigma_s,
        fct_eff,
        k,
        kc,
        Ap,
        phi_s,
        phi_p,
        xi,
        delta_s,
    )


@xw.func
def ec2_2004_As_min_2(
    wk: float,
    sigma_s: float,
    fct_eff: float,
    h_cr: float,
    h: float,
    d: float,
    delta_s: float = 0,
    kc: t.Optional[float] = None,
):
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas

    EUROCODE 2 1992-1-1:2004, Table (7.2N), Table (7.3N)

    Args:
        _wk (float): the characteristic crack width value in mm.
        sigma_s (float): the steel stress value in MPa under the relevant
            combination of actions.
        fct_eff (float): is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct,eff=fct or lower (fct(t)), is cracking
            is expected earlier than 28 days.
        h_cr (float): is the depth of the tensile zone immediately prior to
            cracking, considering the characteristic values of prestress and
            axial forces under the quasi-permanent combination of actions.
        h (float): the overall depth of the section in mm.
        d (float): is the effective depth to the centroid of the outer layer
            of the reinforcement.
        delta_s (float, optional): value of prestressed stress in MPa if
            applicable
        kc (float, optional): is a coefficient which takes account of the
            stress distribution within the section immediately prior to
            cracking and the change of the lever arm in a bending section.
            'None' for pure tensile uniform axial section.

    Returns:
        tuple(float, float): with the value of the maximum bar diameters in mm
        in the first position and the maximum bar spacing in mm in the
        second position
    Raises:
        ValueError: if _wk, fct_eff, h_cr, h or d are less than 0
        ValueError: if kc is not between 0 and 1
        ValueError: if combination of wk and stress values are out of scope"""
    return structuralcodes.codes.ec2_2004.As_min_2(
        wk,
        sigma_s,
        fct_eff,
        h_cr,
        h,
        d,
        delta_s,
        kc,
    )


@xw.func
def ec2_2004_alpha_e(Es: float, Ecm: float):
    """Compute the ratio between the steel and mean concrete
    elastic modules.

    EUROCODE 2 1992-1-1:2004, Section 7.3.4-2

    Args:
        Es (float): steel elastic modulus in MPa
        Ecm (float): concrete mean elastic modulus in MPa

    Returns:
        float: ratio between modules
    Raise:
        ValueError: if any of es or ecm is lower than 0."""
    return structuralcodes.codes.ec2_2004.alpha_e(Es, Ecm)


@xw.func
def ec2_2004_rho_p_eff(As_: float, xi1: float, Ap: float, Ac_eff: float):
    """Effective bond ratio between areas

    EUROCODE 2 1992-1-1:2004, Eq. (7.10)

    Args:
        As (float): steel area in mm2
        _xi1 (float): the adjusted ratio of bond according
            to expression (7.5)
        Ap (float): the area in mm2 of post-tensioned tendons in ac_eff
        Ac_eff (float): effective area of concrete in tension surrounding
        the reinforcement or prestressing tendons of depth hc_eff.

    Returns:
        float: with the retio between areas


    Raise:
        ValueError: if any of As, xi1, Ap or Ac_eff is less than 0"""
    return structuralcodes.codes.ec2_2004.rho_p_eff(As_, xi1, Ap, Ac_eff)


@xw.func
def ec2_2004_kt(load_type: str):
    """Returns the kt factor dependent on the load duration for
    the crack width calculation

    Args:
        load_type (str): the load type:
            - 'short' for term loading
            - 'long' for long term loading

    Returns:
        float: with the kt factor

    Raises:
        ValueError: if load_type is not 'short' and not 'long'"""
    return structuralcodes.codes.ec2_2004.kt(load_type)


@xw.func
def ec2_2004_esm_ecm(
    sigma_s: float,
    alpha_e: float,
    rho_p_eff: float,
    kt: float,
    fct_eff: float,
    Es: float,
):
    """Returns the strain difference (esm - ecm) needed to compute the crack
    width. esm is the mean strain in the reinforcement under the relevant
    combination of loads of imposed deformations and taking into account the
    effects of tension stiffening. Only the additional tensile strain beyond
    the state of zero strain of the concrete is considered. ecm is the mean
    strain in the concrete between the cracks.

    EUROCODE 2 1992-1-1:2004, Eq. (7.9)

    Args:
        sigma_s (float): is the stress in MPa in the tension reinforcement
            assuming a cracked section. FOr pretensioned members, s_steel may
            be replaced by increment of s_steel stress variation in
            prestressing tendons from the state of zero strain of the
            concrete at the same level.
        _alpha_e (float): is the ratio Es/Ecm
        _rho_p_eff (float): effective bond ratio between areas given by the
            Eq. (7.10)
        _kt (float): is a factor dependent on the load duration
        fct_eff (float): is the mean value of the tensile strength in MPa
            of the concrete effectvie at the time when the cracks may
            first be expected to occur: fct_eff=fctm or fctm(t) if
            crack is expected earlier than 28 days.
        Es: steel elastic mudulus in MPa

    Returns:
        float: the strain difference between concrete and steel

    Raises:
        ValueError: if any sigma_s, _alpha_e, _rho_p_eff, fct_eff or Es is less
            than 0.
        ValueError: if _kt is not 0.6 and not 0.4"""
    return structuralcodes.codes.ec2_2004.esm_ecm(
        sigma_s,
        alpha_e,
        rho_p_eff,
        kt,
        fct_eff,
        Es,
    )


@xw.func
def ec2_2004_w_spacing(c: float, phi: float):
    """Computes the distance threshold from which the
    maximum crack spacing is constant.

    EUROCODE 2 1992-1-1:2004, Sect. (7.3.4-3)

    Args:
        c (float): cover of  the longitudinal reinforcement in mm
        phi (float): is the bar diameter in mm. Where mixed bar diameters
            used, then it should be replaced for an equivalente bar diameter.

    Returns:
        float: threshold distance in mm

    Raises:
        ValueError: if any of c or phi is less than 0."""
    return structuralcodes.codes.ec2_2004.w_spacing(c, phi)


@xw.func
def ec2_2004_phi_eq(n1: int, n2: int, phi1: float, phi2: float):
    """Computes the equivalent diameter. For a section with n1 bars of
    diameter phi1 and n2 bars of diameter phi2

    EUROCODE 2 1992-1-1:2004, Sect. (7.12)

    Args:
        n1 (int): number of bars with diameter phi1
        n2 (int): number of bars with diameter phi2
        phi1 (float): diameter of n1 bars in mm
        phi2 (float): diamater of n2 bars in mm

    Returns:
        float: the equivalent diameter in mm

    Raises:
        ValueError: if any of n1 or n2 is less than 0
        ValueError: if any of phi1 or phi2 is less than 0
        TypeError: if any of n1 or n2 is not an integer"""
    return structuralcodes.codes.ec2_2004.phi_eq(n1, n2, phi1, phi2)


@xw.func
def ec2_2004_k1(bond_type: str):
    """Get the k1 coefficient which takes account of the bond properties
    of the bounded reinforcement

    EUROCODE 2 1992-1-1:2004, Eq. (7.11-k1)

    Args:
        bond_type (str): the bond property of the reinforcement.
        Possible values:
            - 'bond': for high bond bars
            - 'plane': for bars with an effectively plain surface (e.g.
            prestressing tendons)

    Returns:
        (float): value of the k1 coefficient

    Raises:
        ValueError: if bond_type is neither 'bond' nor 'plane'
        TypeError: if bond_type is not an str"""
    return structuralcodes.codes.ec2_2004.k1(bond_type)


@xw.func
def ec2_2004_k2(epsilon_r: float):
    """Computes a coefficient which takes into account of the
    distribution of strain:

    EUROCODE 2 1992-1-1:2004, Eq. (7.13)

    Args:
        epsilon_r (float): ratio epsilon_2/epsilon_1 where epsilon_1 is
            thre greater and epsilon_2 is the lesser strain at the boundaries
            of the section considererd, assessed on the basis of a cracked
            section. epsilon_r=0 for bending and epsilon_r=1 for pure tension.

    Returns:
        float: the k2 coefficient value.

    Raises:
        ValueError: if epsilon_r is not between 0 and 1."""
    return structuralcodes.codes.ec2_2004.k2(epsilon_r)


@xw.func
def ec2_2004_k3():
    """Returns the k3 coefficient for computing sr_max

    Returns:
        float: value for the coefficient"""
    return structuralcodes.codes.ec2_2004.k3()


@xw.func
def ec2_2004_k4():
    """Returns the k4 coefficient for computing sr_max

    Returns:
        float: value for the coefficient"""
    return structuralcodes.codes.ec2_2004.k4()


@xw.func
def ec2_2004_sr_max_close(
    c: float,
    phi: float,
    rho_p_eff: float,
    k1: float,
    k2: float,
    k3: float,
    k4: float,
):
    """Computes the maximum crack spacing in cases where bonded reinforcement
    is fixed at reasonably close centres within the tension zone
    (w_spacing<=5(c+phi/2)).

    EUROCODE 2 1992-1-1:2004, Eq. (7.11)

    Args:
        c (float): is the cover in mm of the longitudinal reinforcement
        phi (float): is the bar diameter in mm. Where mixed bar diameters
            used, then it should be replaced for an equivalente bar diameter.
        _rho_p_eff (float): effective bond ratio between areas given by the
            Eq. (7.10)
        _k1 (float): coefficient that takes into account the bound properties
            of the bonded reinforcement
        _k2 (float): coefficient that takes into account the distribution of
            of the strain
        _k3 (float): coefficient from the National Annex
        _k4 (float): coefficient from the National Annex

    Returns:
        float: the maximum crack spaing in mm.

    Raises:
        ValueError: if one or more of c, phi, _rho_p_eff, _k3 or _k4
            is lower than zero.
        ValueError: if _k1 is not 0.8 or 1.6
        ValueError: if _k2 is not between 0.5 and 1.0"""
    return structuralcodes.codes.ec2_2004.sr_max_close(
        c,
        phi,
        rho_p_eff,
        k1,
        k2,
        k3,
        k4,
    )


@xw.func
def ec2_2004_sr_max_far(h: float, x: float):
    """Computes the maximum crack spacing in cases where bonded reinforcement
    is fixed at reasonably close centres within the tension zone
    (w_spacing>5(c+phi/2)).

    EUROCODE 2 1992-1-1:2004, Eq. (7.14)

    Args:
        h (float): total depth of the beam in mm
        x (float): distance to non tension area of the element mm

    Returns:
        float: maximum crack spacing in mm

    Raises:
        ValueError: if one of h or x is less than zero.
        ValueError: x is greater than h."""
    return structuralcodes.codes.ec2_2004.sr_max_far(h, x)


@xw.func
def ec2_2004_sr_max_theta(sr_max_y: float, sr_max_z: float, theta: float):
    """Computes the crack spacing sr_max when there is an angle
    between the angle of  principal stress and the direction
    of the reinforcement, for members in two orthogonal directions,
    that is significant (> 15 degrees).

    EUROCODE 2 1992-1-1:2004, Eq. (7.15)

    Args:
        sr_max_y (float): crack spacing in mm in the y-direction.
        sr_max_z (float): crack spacing in mm in the z-direction.
        theta (float): angle in radians between the reinforcement in the
            y-direction and the direction of the principal tensile stress.

    Returns:
        float: the crack spacing in mm.

    Raises:
        ValueError: if sr_max_y or sr_max_z is negative.
        ValueError: if theta is not between 0 and pi/2"""
    return structuralcodes.codes.ec2_2004.sr_max_theta(
        sr_max_y, sr_max_z, theta
    )


@xw.func
def ec2_2004_wk(sr_max: float, esm_ecm: float):
    """Computes the crack width

    EUROCODE 2 1992-1-1:2004, Eq. (7.8)

    Args:
        sr_max (float): the maximum crack length spacing in mm.
        _esm_ecm (float): the difference between the mean strain in the
            reinforcement under relevant combination of loads, including
            the effect of imposed deformations and taking into account
            tension stiffening and the mean strain in the concrete
            between cracks.

    Returns:
        float: crack width in mm.

    Raises:
        ValueError: if any of sr_max or esm_ecm is less than zero."""
    return structuralcodes.codes.ec2_2004.wk(sr_max, esm_ecm)


@xw.func
def mc2010_fcm(fck: float, delta_f: float = 8.0):
    """Compute the mean concrete compressive strength from the characteristic
    strength.

    fib Model Code 2010, Eq. (5.1-1)

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Keyword Args:
        delta_f (float): The difference between the mean and the
        characteristic strength.

    Returns:
        float: The mean compressive strength in MPa."""
    return structuralcodes.codes.mc2010.fcm(fck, delta_f)


@xw.func
def mc2010_fctm(fck: float):
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength.

    fib Model Code 2010, Eqs. (5.1-3a) and (5.1-3b)

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The mean tensile strength in MPa."""
    return structuralcodes.codes.mc2010.fctm(fck)


@xw.func
def mc2010_fctkmin(fctm: float):
    """Compute the lower bound value of the characteristic tensile strength
    from the mean tensile strength.

    fib Model Code 2010, Eq. (5.1-4)

    Args:
        _fctm (float): The mean tensile strength in MPa.

    Returns:
        float: Lower bound of the characteristic tensile strength in MPa."""
    return structuralcodes.codes.mc2010.fctkmin(fctm)


@xw.func
def mc2010_fctkmax(fctm: float):
    """Compute the upper bound value of the characteristic tensile strength
    from the mean tensile strength.

    fib Model Code 2010, Eq. (5.1-5)

    Args:
        _fctm (float): The mean tensile strength in MPa.

    Returns:
        float: Upper bound of the characteristic tensile strength in MPa."""
    return structuralcodes.codes.mc2010.fctkmax(fctm)


@xw.func
def mc2010_Gf(fck: float):
    """Compute tensile fracture energy from characteristic compressive
    strength.

    fib Model Code 2010, Eq. (5.1-9)

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The tensile fracture energy in N/m."""
    return structuralcodes.codes.mc2010.Gf(fck)


if __name__ == "__main__":
    xw.Book("pruebaxlwings.xlsm").set_mock_caller()
    main()

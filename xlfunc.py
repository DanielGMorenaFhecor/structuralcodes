import xlwings as xw

import structuralcodes.codes.ec2_2004
import structuralcodes.codes.mc2010


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
def ec2_2004_crack_min_steel_area(
    a_ct: float, s_steel: float, fct_eff: float, k: float, kc: float
):
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

    Args:
        a_ct (float): is the area of concrete within the tensile zone in mm2.
            The tensile zone is that parg of the section which is calculated
            to be in tension just before the formation of the first crack.
        s_steel (float): is the absolute value of the maximum stress in MPa
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
        ValueError: if k value is not between 0.65 and 1 or kc is not
            larger than 0 and lower than 1."""
    return structuralcodes.codes.ec2_2004.crack_min_steel_area(
        a_ct, s_steel, fct_eff, k, kc
    )


@xw.func
def ec2_2004_k_crack_min_steel_area(h: float):
    """Is the coefficient which allow for the effect of
    non-uniform self-equilibrating stresses, which lead to a
    reduction of restraint forces. Use 'k_crack_min_steel_area'
    to compute it
    k=1 for webs w<=300mm or flanges widths less than 300mm
    k=0.65 for webs w>=800mm or flanges with widths greater than 800mm

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

    Args:
        h (float): flange length or flange width in mm

    Returns:
        float: k coefficient value

    Raises:
        ValueError: if h is less than 0"""
    return structuralcodes.codes.ec2_2004.k_crack_min_steel_area(h)


@xw.func
def ec2_2004_kc_crack_min_steel_area_pure_tension():
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm in pure dtension.

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

    Returns:
        float: value of the kc coefficient in pure tension"""
    return (
        structuralcodes.codes.ec2_2004.kc_crack_min_steel_area_pure_tension()
    )


@xw.func
def ec2_2004_kc_crack_min_steel_area_rectangular(
    h: float, b: float, fct_eff: float, n_ed: float
):
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
        n_ed (str): axial force at the serviceability limit state acting on
            the part of the cross-section under consideration (compressive
            force positive). n_ed should be determined considering the
            characteristic values of prestress and axial forces under the
            relevant combination of actions

    Returns:
        float: value of the kc coefficient

    Raises:
        ValueError: is h or b are less than 0"""
    return structuralcodes.codes.ec2_2004.kc_crack_min_steel_area_rectangular(
        h, b, fct_eff, n_ed
    )


@xw.func
def ec2_2004_kc_crack_min_steel_area_flanges(
    f_cr: float, a_ct: float, fct_eff: float
):
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm for bending+axial combination
    in rectangular sections for flanges of box sections and T-sections.

    EUROCODE 2 1992-1-1:2004, Eq. (7.3)

    Args:
        f_cr: is the absolute value in kN of the tensile force within the
            flange immediately prior to cracking due to cracking moment
            calculated with fct,eff
        a_ct (float): is the area of concrete within the tensile zone in mm2.
            The tensile zone is that part of the section which is calculated
            to be in tension just before the formation of the first crack.
        fct_eff (float): is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct,eff=fct or lower (fct(t)), is cracking
            is expected earlier than 28 days.

    Returns:
        float: value of the kc coefficient

    Raises:
        ValueError: is a_ct is less than 0mm2"""
    return structuralcodes.codes.ec2_2004.kc_crack_min_steel_area_flanges(
        f_cr, a_ct, fct_eff
    )


@xw.func
def ec2_2004_crack_min_steel_area_with_prestresed_tendons(
    a_ct: float,
    s_steel: float,
    fct_eff: float,
    k: float,
    kc: float,
    ap: float,
    d_steel: float,
    d_press: float,
    e: float,
    incr_stress: float,
):
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas in addition with bonded tendons

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

    Args:
        a_ct (float): is the area of concrete within the tensile zone in mm2.
            The tensile zone is that part of the section which is calculated
            to be in tension just before the formation of the first crack.
        s_steel (float): is the absolute value of the maximum stress in MPa
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
            reduction of restraint forces. Use 'k_crack_min_steel_area'
            to compute it
            k=1 for webs w<=300mm or flanges widths less than 300mm
            k=0.65 for webs w>=800mm or flanges with widths greater than 800mm
            Intermediate values may be interpolated.
        kc (float): is a coefficient which takes account of the stress
            distribution within the section immediately prior to cracking and
            the change of the lever arm.
        ac_eff (float): is the effective area in mm2 of concrete in tension
            surrounding or prestressing tendons if depth hc,ef
        ap (float): is the area in mm2 of pre or post-tensioned tendons
            within ac_eff
        d_steel (float): largest bar diameter in mm of reinforcing steel.
            Equal to zero if only prestressing is used in control cracking
        d_press (float): equivalent diameter in mm of tendon acoording
            to 6.8.2
        e (float): ratio of bond strength of prestressing and reinforcing
            steel, according to Table 6.2 in 6.8.2
        incr_stress (float): stress variation in MPa in prestressing tendons
            from the state of zero strain of the concrete at the same level

    Returns:
        float: the minimm area of reinforcing steel within the tensile
            zone in mm2.

    Raises:
        ValueError: if k value is not between 0.65 and 1 or kc is not
            larger than 0 and lower than 1. If diameters d_steel or
            d_press are lower than 0. If ratio of bond strength e
            is less than 0 or larger than 1. If area of tendons ac_eff
            is less than 0. Is stress variation incr_stress is less than 0"""
    return structuralcodes.codes.ec2_2004.crack_min_steel_area_with_prestresed_tendons(
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
    )


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
    xw.Book("borrar.xlsx").set_mock_caller()
    main()

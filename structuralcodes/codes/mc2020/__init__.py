"""The fib Model Code 2020."""

import typing as t

from ._concrete_limit_state_of_cracking import (
    Ac_ef_bar,
    Ac_ef_group,
    As_min_frc,
    Nr_crack,
    Rax,
    beta_TS,
    beta_TS_x,
    beta_TS_y,
    eps_2,
    eps_crack,
    eps_sm_eps_cm,
    eps_sm_eps_cm_restrained,
    eps_sm_eps_cm_theta,
    eps_sm_x_eps_cm_x,
    eps_sm_y_eps_cm_y,
    k1_r,
    k1_r_simpl,
    kfl,
    kfl_simpl,
    max_phi_crack,
    max_st_crack,
    phi_eq,
    phi_p_eq,
    rho_s_ef,
    rho_s_p_ef,
    sigma_s_x,
    sigma_s_y,
    sigma_sr_ef,
    sr_max,
    sr_max_frc,
    sr_max_theta,
    tau_bmp_tau_bms,
    tau_bms,
    theta_reinf,
    up_i,
    wcal,
    wlim_fluid_tightness,
    wlim_frp,
    wlim_prestressed,
    wlim_rfc,
    xi_1,
)

__all__ = [
    'Ac_ef_bar',
    'Ac_ef_group',
    'As_min_frc',
    'Nr_crack',
    'Rax',
    'beta_TS',
    'beta_TS_x',
    'beta_TS_y',
    'eps_2',
    'eps_crack',
    'eps_sm_eps_cm',
    'eps_sm_eps_cm_restrained',
    'eps_sm_eps_cm_theta',
    'eps_sm_x_eps_cm_x',
    'eps_sm_y_eps_cm_y',
    'k1_r',
    'k1_r_simpl',
    'kfl',
    'kfl_simpl',
    'max_phi_crack',
    'max_st_crack',
    'phi_eq',
    'phi_p_eq',
    'rho_s_ef',
    'rho_s_p_ef',
    'sigma_s_x',
    'sigma_s_y',
    'sigma_sr_ef',
    'sr_max',
    'sr_max_frc',
    'sr_max_theta',
    'tau_bmp_tau_bms',
    'tau_bms',
    'theta_reinf',
    'up_i',
    'wcal',
    'wlim_fluid_tightness',
    'wlim_frp',
    'wlim_prestressed',
    'wlim_rfc',
    'xi_1',
]

__title__: str = 'fib Model Code 2020'
__year__: str = '2020'
__materials__: t.Tuple[str] = ('concrete', 'reinforcement')
import numpy as np


def ra_windstr(u: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """

    A Python function to compute wind stress from wind field data with non-linear Cd.
    Adapted from the original work (MATLAB) by Ramkrushn Patel.
    Adapted for public release by: Qin Duan (duanqin23@mails.ucas.ac.cn)

    Usage
    ----------
    Tx, Ty = ra_windstr(u,v)

    Parameters
    ----------
    u : 2D numpy array
        Zonal wind component [m/s].
    v : 2D numpy array
        Meridional wind component [m/s].

    Returns
    -------
    Tx : 2D numpy array
        Zonal wind stress [N/m^2].
    Ty : 2D numpy array
        Meridional wind stress [N/m^2].


    REFERENCE:

    Large, W. G., & Pond, S. (1981). Open ocean momentum flux measurements in moderate to strong winds. J. PHYS. OCEANOGR., 11(3, Mar. 1981), 324-336. https://doi.org/10.1175/1520-0485(1981)011<0324:oomfmi>2.0.co;2

    Trenberth, K. E., W. G. Large, and J. G. Olson, 1990: The Mean Annual Cycle in Global Ocean Wind Stress. J. Phys. Oceanogr., 20, 1742-1760, https://doi.org/10.1175/1520-0485(1990)020<1742:TMACIG>2.0.CO;2.

    Oey, L.-Y., T. Ezer, D.-P. Wang, S.-J. Fan, and X.-Q. Yin (2006), Loop Current warming by Hurricane Wilma, Geophys. Res. Lett., 33, L08613, doi:10.1029/2006GL025873.

    Ramkrushn Patel (2025). Wind Stresses computation (https://ww2.mathworks.cn/matlabcentral/fileexchange/53391-wind-stresses-computation), MATLAB Central File Exchange.

    """

    if not isinstance(u, np.ndarray):
        raise TypeError("u should be np.ndarray")
    if not isinstance(v, np.ndarray):
        raise TypeError("v should be np.ndarray")

    if u.shape != v.shape:
        raise ValueError("u and v must have the same shape")

    rho = 1.2  # kg/m^3, air density
    wspd = np.sqrt(u**2 + v**2)  # wind speed

    # Initialize Cd with NaN
    Cd = np.full_like(wspd, np.nan, dtype=np.float64)

    # Piecewise definition of Cd
    Cd[wspd <= 1] = 0.00218
    mask = (wspd > 1) & (wspd <= 3)
    Cd[mask] = (0.62 + 1.56 / wspd[mask]) * 0.001
    mask = (wspd > 3) & (wspd <= 10)
    Cd[mask] = 0.00114
    mask = (wspd > 10) & (wspd <= 19)
    Cd[mask] = (0.49 + 0.065 * wspd[mask]) * 0.001
    mask = wspd > 19
    Cd[mask] = (1.364 + 0.0234 * wspd[mask] - 0.0002 * wspd[mask] ** 2) * 0.001

    # Compute stresses
    Tx = Cd * rho * wspd * u
    Ty = Cd * rho * wspd * v

    return Tx, Ty

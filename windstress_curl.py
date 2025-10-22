import numpy as np
from windstress import ra_windstr

def ra_windstrcurl(lat, lon, u, v):
    """
    A Python function to compute wind stress curl using finite difference method.
    Adapted from the original work (MATLAB) by Ramkrushn Patel.
    Adapted for public release by: Qin Duan (duanqin23@mails.ucas.ac.cn)
    Requires use with ra_windstr.py

    Usage
    ----------
    WSC_cal = ra_windstrcurl(lat,lon,u,v)
    
    Parameters
    ----------
    lat : 1D numpy array
        Latitude vector [deg].
    lon : 1D numpy array
        Longitude vector [deg].
    u : 2D numpy array, shape (lat,lon)
        Zonal wind component [m/s].
    v : 2D numpy array, shape (lat,lon)
        Meridional wind component [m/s].

    Returns
    -------
    curlZ : 2D numpy array
        Wind stress curl [N/m^3].
    """

    if u.shape != v.shape:
        raise ValueError("u and v must have the same shape")

    rad = np.pi / 180.0
    lt, ln = u.shape

    # wind stresses
    Tx, Ty = ra_windstr(u, v)

    # check latitude spacing consistency
    a = np.diff(lat)
    if not np.allclose(a, a[0]):
        raise ValueError("Latitude difference is not consistent")
    dlat = a[0]
    deltay = dlat * 111176.0  # [m]

    # longitude distance in meters for each grid
    long = np.zeros((lt, ln))
    for ii in range(lt):
        long[ii, :] = lon * 111176.0 * np.cos(lat[ii] * rad)
        # alternatively: lon * 6378137.0 * rad * cos(lat[ii] * rad)

    curlZ = np.full((lt, ln), np.nan)

    # central difference
    for ii in range(1, lt - 1):
        for jj in range(1, ln - 1):
            curlZ[ii, jj] = (Ty[ii, jj + 1] - Ty[ii, jj - 1]) / (2 * (long[ii, jj + 1] - long[ii, jj - 1])) \
                          - (Tx[ii + 1, jj] - Tx[ii - 1, jj]) / (2 * deltay)

    # forward difference at top/left boundaries
    for jj in range(ln - 1):
        curlZ[0, jj] = (Ty[0, jj + 1] - Ty[0, jj]) / (long[0, jj + 1] - long[0, jj]) \
                     - (Tx[1, jj] - Tx[0, jj]) / deltay
    for ii in range(lt - 1):
        curlZ[ii, 0] = (Ty[ii, 1] - Ty[ii, 0]) / (long[ii, 1] - long[ii, 0]) \
                     - (Tx[ii + 1, 0] - Tx[ii, 0]) / deltay
    curlZ[0, -1] = curlZ[0, -2]

    # backward difference at bottom/right boundaries
    for ii in range(1, lt):
        curlZ[ii, -1] = (Ty[ii, -1] - Ty[ii, -2]) / (long[ii, -1] - long[ii, -2]) \
                      - (Tx[ii, -1] - Tx[ii - 1, -1]) / deltay
    for jj in range(1, ln - 1):
        curlZ[-1, jj] = (Ty[-1, jj] - Ty[-1, jj - 1]) / (long[-1, jj] - long[-1, jj - 1]) \
                      - (Tx[-1, jj] - Tx[-2, jj]) / deltay
    curlZ[-1, 0] = curlZ[-1, -2]

    return curlZ
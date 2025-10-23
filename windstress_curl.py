import numpy as np

from windstress import ra_windstr


def ra_windstrcurl(lat: np.ndarray, lon: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """
    A Python function to compute wind stress curl using finite difference method.
    Adapted from the original work (MATLAB) by Ramkrushn Patel.
    Adapted for public release by: Qin Duan (duanqin23@mails.ucas.ac.cn)
    Requires use with func `ra_windstr`

    Usage
    ----------
    curlZ = ra_windstrcurl(lat,lon,u,v)

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

    if isinstance(u, np.ndarray):
        raise TypeError("u should be np.ndarray")
    if isinstance(v, np.ndarray):
        raise TypeError("v should be np.ndarray")
    if isinstance(lon, np.ndarray):
        raise TypeError("lon should be np.ndarray")
    if isinstance(lat, np.ndarray):
        raise TypeError("lat should be np.ndarray")
    if lat.ndim != 1 or lon.ndim != 1:
        raise ValueError(f"the dimension of lat/lon should be 1 but got {lat.ndim}/{lat.ndim}")
    if u.ndim != 2 or v.ndim != 2:
        raise ValueError(f"the dimension of u/v should be 2 but got {u.ndim}/{v.ndim}")
    if u.shape[0] != lat.size or u.shape[1] != lon.size:
        raise ValueError(f"the shape of u shoule be ({lat.size}, {lon.size}) but got {u.shape}")
    if v.shape[0] != lat.size or v.shape[1] != lon.size:
        raise ValueError(f"the shape of v shoule be ({lat.size}, {lon.size}) but got {v.shape}")

    rad = np.pi / 180.0
    latsize, lonsize = u.shape

    # wind stresses
    Tx, Ty = ra_windstr(u, v)

    # check latitude spacing consistency
    a = np.diff(lat)
    if not np.allclose(a, a[0]):
        raise ValueError("Latitude difference is not consistent")
    dlat = a[0]
    dy = dlat * 111_176.0  # [m]

    # longitude distance in meters for each grid
    dx = np.zeros((latsize, lonsize))
    for ii in range(latsize):
        dx[ii, :] = lon * 111_176.0 * np.cos(np.deg2rad(lat[ii]))
        # alternatively: lon * 6378137.0 * rad * cos(lat[ii] * rad)

    curlZ = np.full((latsize, lonsize), np.nan)

    # central difference
    curlZ[1:-1, 1:-1] = (Ty[1:-1, 2:] - Ty[1:-1, :-2]) / (2 * (dx[1:-1, 2:] - dx[1:-1, :-2])) - (
        Tx[2:, 1:-1] - Tx[:-2, 1:-1]
    ) / (2 * dy)

    # forward difference at top/left boundaries
    curlZ[0, :-1] = (Ty[0, 1:] - Ty[0, :-1]) / (dx[0, 1:] - dx[0, :-1]) - (
        Tx[1, :-1] - Tx[0, :-1]
    ) / dy
    curlZ[:-1, 0] = (Ty[:-1, 1] - Ty[:-1, 0]) / (dx[:-1, 1] - dx[:-1, 0]) - (
        Tx[1:, 0] - Tx[:-1, 0]
    ) / dy
    curlZ[0, -1] = curlZ[0, -2]

    # backward difference at bottom/right boundaries
    curlZ[1:, -1] = (Ty[1:, -1] - Ty[1:, -2]) / (dx[1:, -1] - dx[1:, -2]) - (
        Tx[1:, -1] - Tx[:-1, -1]
    )
    curlZ[-1, 1:-1] = (Ty[-1, 1:-1] - Ty[-1, :-2]) / (dx[-1, 1:-1] - dx[-1, :-2]) - (
        Tx[-1, 1:-1] - Tx[-2, 1:-1]
    ) / dy
    curlZ[-1, 0] = curlZ[-1, -2]

    return curlZ

import numpy as np

def calc_sun_angles(lat, lon, stdlon, doy, ftime):
    '''Calculates the Sun Zenith and Azimuth Angles (SZA & SAA).

    Parameters
    ----------
    lat : float
        latitude of the site (degrees).
    long : float
        longitude of the site (degrees).
    stdlng : float
        central longitude of the time zone of the site (degrees).
    doy : float
        day of year of measurement (1-366).
    ftime : float
        time of measurement (decimal hours).

    Returns
    -------
    sza : float
        Sun Zenith Angle (degrees).
    saa : float
        Sun Azimuth Angle (degrees).

    '''

    lat, lon, stdlon, doy, ftime = map(
        np.asarray, (lat, lon, stdlon, doy, ftime))
    # Calculate declination
    declination = 0.409 * np.sin((2.0 * np.pi * doy / 365.0) - 1.39)
    EOT = 0.258 * np.cos(declination) - 7.416 * np.sin(declination) - \
        3.648 * np.cos(2.0 * declination) - 9.228 * np.sin(2.0 * declination)
    LC = (stdlon - lon) / 15.
    time_corr = (-EOT / 60.) + LC
    solar_time = ftime - time_corr
    # Get the hour angle
    w = np.asarray((solar_time - 12.0) * 15.)
    # Get solar elevation angle
    sin_thetha = np.cos(np.radians(w)) * np.cos(declination) * np.cos(np.radians(lat)) + \
                 np.sin(declination) * np.sin(np.radians(lat))
    sun_elev = np.arcsin(sin_thetha)
    # Get solar zenith angle
    sza = np.pi / 2.0 - sun_elev
    sza = np.asarray(np.degrees(sza))
    # Get solar azimuth angle
    cos_phi = np.asarray(
        (np.sin(declination) * np.cos(np.radians(lat)) -
         np.cos(np.radians(w)) * np.cos(declination) * np.sin(np.radians(lat))) /
        np.cos(sun_elev))
    saa = np.zeros(sza.shape)
    saa[w <= 0.0] = np.degrees(np.arccos(cos_phi[w <= 0.0]))
    saa[w > 0.0] = 360. - np.degrees(np.arccos(cos_phi[w > 0.0]))
    return np.asarray(sza), np.asarray(saa)


def calc_omega0_Kustas(LAI, f_C, x_LAD=1, isLAIeff=True):
    ''' Nadir viewing clmping factor

    Estimates the clumping factor forcing equal gap fraction between the real canopy
    and the homogeneous case, after [Kustas1999]_.

    Parameters
    ----------
    LAI : float
        Leaf Area Index, it can be either the effective LAI or the real LAI
        , default input LAI is effective.
    f_C : float
        Apparent fractional cover, estimated from large gaps, means that
        are still gaps within the canopy to be quantified.
    x_LAD : float, optional
        Chi parameter for the ellipsoildal Leaf Angle Distribution function of
        [Campbell1988]_ [default=1, spherical LIDF].
    isLAIeff :  bool, optional
        Defines whether the input LAI is effective or local.

    Returns
    -------
    omega0 : float
        clumping index at nadir.

    References
    ----------
    .. [Kustas1999] William P Kustas, John M Norman, Evaluation of soil and vegetation heat
        flux predictions using a simple two-source model with radiometric temperatures for
        partial canopy cover, Agricultural and Forest Meteorology, Volume 94, Issue 1,
        Pages 13-29, http://dx.doi.org/10.1016/S0168-1923(99)00005-2.
    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998), An introduction to environmental
        biophysics. Springer, New York
        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.
 '''

    # Convert input scalars to numpy array
    LAI, f_C, x_LAD = map(np.asarray, (LAI, f_C, x_LAD))
    theta = np.zeros(LAI.shape)
    # Estimate the beam extinction coefficient based on a ellipsoidal LAD function
    # Eq. 15.4 of Campbell and Norman (1998)
    K_be = np.sqrt(x_LAD**2 + np.tan(theta)**2) / \
        (x_LAD + 1.774 * (x_LAD + 1.182)**-0.733)
    if isLAIeff:
        F = LAI / f_C
    else:  # The input LAI is actually the real LAI
        F = np.array(LAI)
    # Calculate the gap fraction of our canopy
    trans = np.asarray(f_C * np.exp(-K_be * F) + (1.0 - f_C))
    trans[trans <= 0] = 1e-36
    # and then the nadir clumping factor
    omega0 = -np.log(trans) / (F * K_be)
    return omega0


def calc_omega_Kustas(omega0, theta, w_C=1):
    ''' Clumping index at an incidence angle.

    Estimates the clumping index for a given incidence angle assuming randomnly placed canopies.

    Parameters
    ----------
    omega0 : float
        clumping index at nadir, estimated for instance by :func:`calc_omega0_Kustas`.
    theta : float
        incidence angle (degrees).
    w_C :  float, optional
        canopy witdth to height ratio, [default = 1].

    Returns
    -------
    Omega : float
        Clumping index at an incidenc angle.

    References
    ----------
    .. [Kustas1999] William P Kustas, John M Norman, Evaluation of soil and vegetation heat
        flux predictions using a simple two-source model with radiometric temperatures for
        partial canopy cover, Agricultural and Forest Meteorology, Volume 94, Issue 1,
        Pages 13-29, http://dx.doi.org/10.1016/S0168-1923(99)00005-2.
    '''

    w_C = 1.0 / w_C
    omega = omega0 / (omega0 + (1.0 - omega0) *
                      np.exp(-2.2 * (np.radians(theta))**(3.8 - 0.46 * w_C)))
    return omega


def calc_omega_rows(
        LAI,
        f_c0,
        theta=0,
        psi=0,
        D=1,
        x_LAD=1,
        isLAIeff=True):

    ''' Clumping index in row crops.
    Calculates the clumping index for a given incidence angle assuming structured row crops.
    Parameters
    ----------
    LAI : float
        Leaf Area Index, it can be either the effective LAI or the real LAI
        depending on isLAIeff, default input LAI is effective.
    f_c0 : float
        Apparent nadir fractional cover, can be expresses as canopy width/row spacing.
    theta : float, optional
        Incidence angle (degrees), default nadir.
    psi : float, optional
        relative row-sun azimiuth angle
    D :  float, optional
        canopy witdth to height ratio, [default = 1].
    x_LAD : float, optional
        Chi parameter for the ellipsoildal Leaf Angle Distribution function of
        [Campbell1988]_ [default=1, spherical LIDF].
    isLAIeff :  bool, optional
        Defines whether the input LAI is effective or real. [default True]

    Returns
    -------
    omega : float
        clumping index at an incidence angle.

    References
    ----------
    .. [Parry2018] Parry, C. K., H. Nieto, P. Guillevic, N. Agam, W. P. Kustas, J. Alfieri, L. McKee, and A. J. McElrone.
        An intercomparison of radiation partitioning models in vineyard canopies.
        Irrigation Science. Pages 1-14. https://doi.org/10.1007/s00271-019-00621-x.
    '''

    # Convert input scalars in numpy arrays
    # LAI, f_c0, theta, psi, D, x_LAD = map(
    #     np.asarray, (LAI, f_c0, theta, psi, D, x_LAD))
    omega = np.zeros(LAI.shape)
    # Calculate the zenith angle of incidence towards the normal of the row
    # direction
    tan_alpha_x = np.tan(np.radians(theta)) * abs(np.sin(np.radians(psi)))
    # Calculate the fraction that is trasmitted trough vegetation
    f_c = np.asarray(f_c0 * (1.0 + (tan_alpha_x / D)))

    f_c = np.minimum(f_c, 1.0)
    # Estimate the beam extinction coefficient based on a elipsoidal LAD function
    # Eq. 15.4 of Campbell and Norman (1998)
    K_be = np.sqrt(x_LAD**2 + np.tan(np.radians(theta))**2) / \
           (x_LAD + 1.774 * (x_LAD + 1.182)**-0.733)
    if isLAIeff:
        F = LAI / f_c0
    else:
        F = np.asarray(LAI)
    # Calculate the real gap fraction of our canopy
    trans = f_c * np.exp(-K_be * F) + (1.0 - f_c)
    # and then the clumping factor
    omega[trans > 0] = -np.log(trans[trans > 0]) / (F[trans > 0] * K_be[trans > 0])

    return omega



def calc_difuse_ratio(S_dn, sza, press=1013.25, SOLAR_CONSTANT=1320):
    '''Fraction of difuse shortwave radiation.

    Partitions the incoming solar radiation into PAR (visible) and non-PR and
    diffuse and direct beam component of the solar spectrum.

    Parameters
    ----------
    S_dn : float
        Incoming shortwave radiation (W m-2).
    sza : float
        Solar Zenith Angle (degrees).
    Wv : float, optional
        Total column precipitable water vapour (g cm-2), default 1 g cm-2.
    press : float, optional
        atmospheric pressure (mb), default at sea level (1013mb).

    Returns
    -------
    difvis : float
        diffuse fraction in the visible region.
    difnir : float
        diffuse fraction in the NIR region.
    fvis : float
        fration of total visible radiation.
    fnir : float
        fraction of total NIR radiation.

    References
    ----------
    .. [Weiss1985] Weiss and Norman (1985) Partitioning solar radiation into direct and diffuse,
        visible and near-infrared components, Agricultural and Forest Meteorology,
        Volume 34, Issue 2, Pages 205-213,
        http://dx.doi.org/10.1016/0168-1923(85)90020-6.
    '''

    # Convert input scalars to numpy arrays
    S_dn, sza, press = map(np.asarray, (S_dn, sza, press))
    difvis, difnir, fvis, fnir = [np.zeros(S_dn.shape) for i in range(4)]
    fvis = fvis + 0.6
    fnir = fnir + 0.4
    # Calculate potential (clear-sky) visible and NIR solar components
    # Weiss & Norman 1985
    Rdirvis, Rdifvis, Rdirnir, Rdifnir = calc_potential_irradiance_weiss(
        sza, press=press, SOLAR_CONSTANT=SOLAR_CONSTANT)
    # Potential total solar radiation
    potvis = np.asarray(Rdirvis + Rdifvis)
    potvis[potvis <= 0] = 1e-6
    potnir = np.asarray(Rdirnir + Rdifnir)
    potnir[potnir <= 0] = 1e-6
    fclear = S_dn / (potvis + potnir)
    fclear = np.minimum(1.0, fclear)
    # Partition S_dn into VIS and NIR
    fvis = potvis / (potvis + potnir)  # Eq. 7
    fnir = potnir / (potvis + potnir)  # Eq. 8
    fvis = np.maximum(0, fvis)
    fvis = np.minimum(1, fvis)
    fnir = 1.0 - fvis
    # Estimate direct beam and diffuse fractions in VIS and NIR wavebands
    ratiox = np.asarray(fclear)
    ratiox[fclear > 0.9] = 0.9
    dirvis = (Rdirvis / potvis) * (1. - ((.9 - ratiox) / .7)**.6667)  # Eq. 11
    ratiox = np.asarray(fclear)
    ratiox[fclear > 0.88] = 0.88
    dirnir = (Rdirnir / potnir) * \
        (1. - ((.88 - ratiox) / .68)**.6667)  # Eq. 12
    dirvis = np.maximum(0.0, dirvis)
    dirnir = np.maximum(0.0, dirnir)
    dirvis = np.minimum(1, dirvis)
    dirnir = np.minimum(1, dirnir)
    difvis = 1.0 - dirvis
    difnir = 1.0 - dirnir
    return np.asarray(difvis), np.asarray(difnir), np.asarray(fvis), np.asarray(fnir)



def calc_K_be_Campbell(theta, x_LAD=1):
    ''' Beam extinction coefficient

    Calculates the beam extinction coefficient based on [Campbell1998]_ ellipsoidal
    leaf inclination distribution function.

    Parameters
    ----------
    theta : float
        incidence zenith angle (degrees).
    x_LAD : float, optional
        Chi parameter for the ellipsoidal Leaf Angle Distribution function,
        use x_LAD=1 for a spherical LAD.

    Returns
    -------
    K_be : float
        beam extinction coefficient.
    x_LAD: float, optional
        x parameter for the ellipsoidal Leaf Angle Distribution function,
        use x_LAD=1 for a spherical LAD.

    References
    ----------
    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998), An introduction to environmental
        biophysics. Springer, New York
        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.
    '''

    theta = np.radians(theta)
    K_be = np.sqrt(x_LAD**2 + np.tan(theta)**2) / \
        (x_LAD + 1.774 * (x_LAD + 1.182)**-0.733)
    return np.asarray(K_be)


def calc_potential_irradiance_weiss(
        sza,
        press=1013.25,
        SOLAR_CONSTANT=1320,
        fnir_ini=0.5455):
    ''' Estimates the potential visible and NIR irradiance at the surface

    Parameters
    ----------
    sza : float
        Solar Zenith Angle (degrees)
    press : Optional[float]
        atmospheric pressure (mb)

    Returns
    -------
    Rdirvis : float
        Potential direct visible irradiance at the surface (W m-2)
    Rdifvis : float
        Potential diffuse visible irradiance at the surface (W m-2)
    Rdirnir : float
        Potential direct NIR irradiance at the surface (W m-2)
    Rdifnir : float
        Potential diffuse NIR irradiance at the surface (W m-2)

    based on Weiss & Normat 1985, following same strategy in Cupid's RADIN4 subroutine.
    '''

    # Convert input scalars to numpy arrays
    sza, press = map(np.asarray, (sza, press))

    # Set defaout ouput values
    Rdirvis, Rdifvis, Rdirnir, Rdifnir, w = [
        np.zeros(sza.shape) for i in range(5)]

    coszen = np.cos(np.radians(sza))
    # Calculate potential (clear-sky) visible and NIR solar components
    # Weiss & Norman 1985
    # Correct for curvature of atmos in airmas (Kasten and Young,1989)
    i = sza < 90
    airmas = 1.0 / coszen
    # Visible PAR/NIR direct beam radiation
    Sco_vis = SOLAR_CONSTANT * (1.0 - fnir_ini)
    Sco_nir = SOLAR_CONSTANT * fnir_ini
    # Directional trasnmissivity
    # Calculate water vapour absorbance (Wang et al 1976)
    # A=10**(-1.195+.4459*np.log10(1)-.0345*np.log10(1)**2)
    # opticalDepth=np.log(10.)*A
    # T=np.exp(-opticalDepth/coszen)
    # Asssume that most absortion of WV is at the NIR
    Rdirvis[i] = (Sco_vis * np.exp(-.185 * (press[i] / 1313.25) * airmas[i]) -
                  w[i]) * coszen[
                     i]  # Modified Eq1 assuming water vapor absorption
    # Rdirvis=(Sco_vis*exp(-.185*(press/1313.25)*airmas))*coszen
    # #Eq. 1
    Rdirvis = np.maximum(0, Rdirvis)
    # Potential diffuse radiation
    # Eq 3                                      #Eq. 3
    Rdifvis[i] = 0.4 * (Sco_vis * coszen[i] - Rdirvis[i])
    Rdifvis = np.maximum(0, Rdifvis)

    # Same for NIR
    # w=SOLAR_CONSTANT*(1.0-T)
    w = SOLAR_CONSTANT * \
        10 ** (-1.195 + .4459 * np.log10(coszen[i]) - .0345 * np.log10(
        coszen[i]) ** 2)  # Eq. .6
    Rdirnir[i] = (Sco_nir * np.exp(-0.06 * (press[i] / 1313.25)
                                   * airmas[i]) - w) * coszen[i]  # Eq. 4
    Rdirnir = np.maximum(0, Rdirnir)
    # Potential diffuse radiation
    Rdifnir[i] = 0.6 * (Sco_nir * coszen[i] - Rdirvis[i] - w)  # Eq. 5
    Rdifnir = np.maximum(0, Rdifnir)
    Rdirvis, Rdifvis, Rdirnir, Rdifnir = map(
        np.asarray, (Rdirvis, Rdifvis, Rdirnir, Rdifnir))
    return Rdirvis, Rdifvis, Rdirnir, Rdifnir


def calc_spectra_Cambpell(LAI, sza, rho_leaf, tau_leaf, rho_soil, x_LAD=1,
                          LAI_eff=None):
    ''' Canopy spectra

    Estimate canopy spectral using the [Campbell1998]_
    Radiative Transfer Model

    Parameters
    ----------
    LAI : float
        Effective Leaf (Plant) Area Index.
    sza : float
        Sun Zenith Angle (degrees).
    rho_leaf : float, or array_like
        Leaf bihemispherical reflectance
    tau_leaf : float, or array_like
        Leaf bihemispherical transmittance
    rho_soil : float
        Soil bihemispherical reflectance
    x_LAD : float,  optional
        x parameter for the ellipsoildal Leaf Angle Distribution function of
        Campbell 1988 [default=1, spherical LIDF].
    LAI_eff : float or None, optional
        if set, its value is the directional effective LAI
        to be used in the beam radiation, if set to None we assume homogeneous canopies.

    Returns
    -------
    albb : float or array_like
        Beam (black sky) canopy albedo
    albd : float or array_like
        Diffuse (white sky) canopy albedo
    taubt : float or array_like
        Beam (black sky) canopy transmittance
    taudt : float or array_like
        Beam (white sky) canopy transmittance

    References
    ----------
    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998), An introduction to environmental
        biophysics. Springer, New York
        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.
    '''

    # calculate aborprtivity
    amean = 1.0 - rho_leaf - tau_leaf
    del rho_leaf, tau_leaf
    # Calculate canopy beam extinction coefficient
    # Modification to include other LADs
    if isinstance(LAI_eff, type(None)):
        LAI_eff = np.asarray(LAI)
    else:
        LAI_eff = np.asarray(LAI_eff)

    # D I F F U S E   C O M P O N E N T S
    # Integrate to get the diffuse transmitance
    taud = 0
    for angle in range(0, 90, 5):
        akd = calc_K_be_Campbell(angle, x_LAD)  # Eq. 15.4
        taub = np.exp(-akd * LAI)
        taud = taud + taub * np.cos(np.radians(angle)) * \
               np.sin(np.radians(angle)) * np.radians(5)

    taud = 2.0 * taud
    # Diffuse light canopy reflection coefficients  for a deep canopy
    akd = -np.log(taud) / LAI
    del taub, taud
    rcpy = (1.0 - np.sqrt(amean)) / (1.0 + np.sqrt(amean))  # Eq 15.7
    rdcpy = 2.0 * akd * rcpy / (akd + 1.0)  # Eq 15.8
    # Diffuse canopy transmission and albedo coeff for a generic canopy
    expfac = np.sqrt(amean) * akd * LAI
    del akd, LAI
    xnum = (rdcpy * rdcpy - 1.0) * np.exp(-expfac)
    xden = (rdcpy * rho_soil - 1.0) + rdcpy * \
           (rdcpy - rho_soil) * np.exp(-2.0 * expfac)
    taudt = xnum / xden  # Eq 15.11
    fact = ((rdcpy - rho_soil) / (rdcpy * rho_soil - 1.0)) * np.exp(
        -2.0 * expfac)
    del expfac, xnum, xden
    albd = (rdcpy + fact) / (1.0 + rdcpy * fact)  # Eq 15.9
    del rdcpy
    # B E A M   C O M P O N E N T S
    # Direct beam extinction coeff (spher. LAD)
    akb = calc_K_be_Campbell(sza, x_LAD)  # Eq. 15.4
    # Direct beam canopy reflection coefficients for a deep canopy
    rbcpy = 2.0 * akb * rcpy / (akb + 1.0)  # Eq 15.8
    del rcpy, sza, x_LAD
    # Beam canopy transmission and albedo coeff for a generic canopy (visible)
    expfac = np.sqrt(amean) * akb * LAI_eff
    del amean, akb, LAI_eff
    xnum = (rbcpy * rbcpy - 1.0) * np.exp(-expfac)
    xden = (rbcpy * rho_soil - 1.0) + rbcpy * \
           (rbcpy - rho_soil) * np.exp(-2.0 * expfac)
    taubt = xnum / xden  # Eq 15.11
    del xnum, xden
    fact = ((rbcpy - rho_soil) / (rbcpy * rho_soil - 1.0)) * np.exp(
        -2.0 * expfac)
    del expfac
    albb = (rbcpy + fact) / (1.0 + rbcpy * fact)  # Eq 15.9

    return albb, albd, taubt, taudt


def calc_Sn_Campbell(LAI, sza, S_dn_dir, S_dn_dif, fvis, fnir, rho_leaf_vis,
                     tau_leaf_vis, rho_leaf_nir, tau_leaf_nir, rsoilv, rsoiln,
                     x_LAD=1, LAI_eff=None):
    ''' Net shortwave radiation

    Estimate net shorwave radiation for soil and canopy below a canopy using the [Campbell1998]_
    Radiative Transfer Model, and implemented in [Kustas1999]_

    Parameters
    ----------
    LAI : float
        Effective Leaf (Plant) Area Index.
    sza : float
        Sun Zenith Angle (degrees).
    S_dn_dir : float
        Broadband incoming beam shortwave radiation (W m-2).
    S_dn_dif : float
        Broadband incoming diffuse shortwave radiation (W m-2).
    fvis : float
        fration of total visible radiation.
    fnir : float
        fraction of total NIR radiation.
    rho_leaf_vis : float
        Broadband leaf bihemispherical reflectance in the visible region (400-700nm).
    tau_leaf_vis : float
        Broadband leaf bihemispherical transmittance in the visible region (400-700nm).
    rho_leaf_nir : float
        Broadband leaf bihemispherical reflectance in the NIR region (700-2500nm).
    tau_leaf_nir : float
        Broadband leaf bihemispherical transmittance in the NIR region (700-2500nm).
    rsoilv : float
        Broadband soil bihemispherical reflectance in the visible region (400-700nm).
    rsoiln : float
        Broadband soil bihemispherical reflectance in the NIR region (700-2500nm).
    x_LAD : float,  optional
        x parameter for the ellipsoildal Leaf Angle Distribution function of
        Campbell 1988 [default=1, spherical LIDF].
    LAI_eff : float or None, optional
        if set, its value is the directional effective LAI
        to be used in the beam radiation, if set to None we assume homogeneous canopies.

    Returns
    -------
    Sn_C : float
        Canopy net shortwave radiation (W m-2).
    Sn_C_dif : float
        Canopy diffuse shortwave radiation (W m-2).   
        
    Sn_S : float
        Soil net shortwave radiation (W m-2).
        
    Sn_C_par_dir: float
        Canopy PAR direct+diffuse radiation
    
    Sn_C_par_dif: float
        Canopy PAR difuse radiation

    References
    ----------
    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998), An introduction to environmental
        biophysics. Springer, New York
        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.
    .. [Kustas1999] Kustas and Norman (1999) Evaluation of soil and vegetation heat
        flux predictions using a simple two-source model with radiometric temperatures for
        partial canopy cover, Agricultural and Forest Meteorology, Volume 94, Issue 1,
        Pages 13-29, http://dx.doi.org/10.1016/S0168-1923(99)00005-2.
    '''

    rho_leaf = np.array((rho_leaf_vis, rho_leaf_nir))
    tau_leaf = np.array((tau_leaf_vis, tau_leaf_nir))
    rho_soil = np.array((rsoilv, rsoiln))
    albb, albd, taubt, taudt = calc_spectra_Cambpell(LAI,
                                                     sza,
                                                     rho_leaf,
                                                     tau_leaf,
                                                     rho_soil,
                                                     x_LAD=x_LAD,
                                                     LAI_eff=LAI_eff)

    Sn_C = (1.0 - taubt[0]) * (1.0 - albb[0]) * S_dn_dir * fvis + \
           (1.0 - taubt[1]) * (1.0 - albb[1]) * S_dn_dir * fnir + \
           (1.0 - taudt[0]) * (1.0 - albd[0]) * S_dn_dif * fvis + \
           (1.0 - taudt[1]) * (1.0 - albd[1]) * S_dn_dif * fnir
    
    Sn_C_dif=(1.0 - taudt[0]) * (1.0 - albd[0]) * S_dn_dif * fvis + \
             (1.0 - taudt[1]) * (1.0 - albd[1]) * S_dn_dif * fnir

    Sn_S = taubt[0] * (1.0 - rsoilv) * S_dn_dir * fvis + \
           taubt[1] * (1.0 - rsoiln) * S_dn_dir * fnir + \
           taudt[0] * (1.0 - rsoilv) * S_dn_dif * fvis + \
           taudt[1] * (1.0 - rsoiln) * S_dn_dif * fnir
           
    Sn_C_par_dir=(1.0 - taubt[0]) * (1.0 - albb[0]) * S_dn_dir * fvis + \
           (1.0 - taudt[0]) * (1.0 - albd[0]) * S_dn_dif * fvis
    
    Sn_C_par_dif= (1.0 - taudt[0]) * (1.0 - albd[0]) * S_dn_dif * fvis

    return np.asarray(Sn_C), np.asarray(Sn_C_dif),np.asarray(Sn_S), np.array(Sn_C_par_dir),\
            np.array(Sn_C_par_dif)

def LAI_sun_shade(LAI,omega,theta,x_LAD=1):
    
    """
    Calculates sunlit and shade leaf area for clumped canopies [Campbell&Norman,1998]
    
    Parameters
    ----------
    LAI: float
        Leaf area index (m2 leaf m-2 soil)
    omega: float
        clumping factor
    theta : float
        incidence zenith angle (degrees).
    x_LAD : float, optional
        Chi parameter for the ellipsoidal Leaf Angle Distribution function,
        use x_LAD=1 for a spherical LAD.
        
    Returns
    -------
    
    LAI_sun: float
        Sunlit LAI
    LAI_shade: float
        Shaded LAI
    
    References
    ----------
    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998), An introduction to environmental
        biophysics. Springer, New York
        
        
    """
    kbe=calc_K_be_Campbell(theta, x_LAD)
    
    LAI_sun=(1-np.exp(-kbe*omega*LAI))/kbe
    LAI_shade=LAI-LAI_sun
    
    return LAI_sun,LAI_shade


def rint_spheroid(sza,rowdist,plantspace,ratio_rx_rz,htop,hbase,inputparamG,LAI,
                  S_dn_dir,S_dn_dif):
    
    """
    Calculates radiation interception on tree. Assume the canopy is an spheroid
    
    Input:
    -----    
        sza: Solar zentih angle (degrees)
        rowdist: row distance (m)
        plantspace: space between plants (m)
        ratio_rx_rz: relation between horizontal and vertical axis
        htop: top of the canopy
        hbase: base of the canopy (m)
        inputparamG: leaf angle distribution for G calculations
        LAI:Leaf Area Index (m2 leaf m2 soil)
        S_dn_dir : float
            Broadband incoming beam shortwave radiation (W m-2).
        S_dn_dif : float
            Broadband incoming diffuse shortwave radiation (W m-2).
       
        
    Return:
    ------
        iparsun: Intercepted radiation sunlit leaves (micromol m-2 leaf s-1)
        iparshade:Intercepted radiation shade leaves (micromol m-2 leaf s-1)
        LAIsun:LAI sunlit leaves
        LAIshade:LAI shaded leaves
        
    References:
    ----------
        López-Bernal, Á., Morales, A., García-Tejera, O., Testi, L., 
            Orgaz, F., De Melo-Abreu, J. P., and Villalobos, F. J. (2018). 
            OliveCan: A Process-Based Model of Development, Growth and 
            Yield of Olive Orchards. Frontiers in Plant Science 9.
        
        
        
    """ 
    cenit=sza*np.pi/180
    pardr=S_dn_dir
    PARDF=S_dn_dif
    
    rz=(htop-hbase)/2
    rx=rz*ratio_rx_rz
    
    space_per_plant=rowdist*plantspace
    
    v=4/3*np.pi*rx**2*rz#tree volume
    LAD=LAI*space_per_plant/v#Leaf area density m2Leaf m-3canopy volume
    
    R = rx    
    c = ratio_rx_rz
    
    #Qdif (fraction of incoming diffuse radiation intercepted)
    #Integrating in 9 intervals
    QDIF=0
    cenit_dif = 0.078539816
    for c in range(0, 9):
          sza_dif=c*np.pi/180
          g = GFUNCTION(sza_dif,inputparamG)
          u = np.sqrt((np.cos(cenit_dif))**2 + ratio_rx_rz**2 * (np.sin(cenit_dif))** 2)
          b = R * u / np.cos(cenit_dif)
          #Mean interception over PEA(0) (=SA for no overlapping)
          AsX = g * LAD * 1.5 * v / (np.pi * R ** 2 * u)
          t = 2 * (1 - (1 + AsX) * np.exp(-AsX)) / (AsX ** 2)
          #SA for a given reference surface the projected shadow modified by mean interception
          a = np.sqrt(space_per_plant) / 2 / b
    
          if a < 1: 
              f1 = R ** 2 * u / np.cos(cenit_dif) * np.pi
              #f2 = 2 * R ^ 2 * u / Cos(cenit) * (ArcSin(a) - a * Sqr(1 - a ^ 2))
              f2 = R ** 2 * u / np.cos(cenit_dif) * (np.arccos(a) - a * np.sqrt(1 - a ** 2))
        
              SA = (f1 - 2 * f2) * (1 - t) + f2 * (1 - t ** 2)
          else:
              SA = R ** 2 * u / np.cos(cenit_dif) * np.pi * (1 - t)
          
          SA = min(SA, space_per_plant)
          qdifx = SA * np.sin(cenit_dif) * 0.157079632
          QDIF = QDIF + qdifx
          qdifx = 0
          cenit_dif = cenit_dif + 0.157079632
      
    QDIF = QDIF / space_per_plant
        
# Qdir (fraction of incoming direct radiation (over horizontal plane) intercepted)
    g = GFUNCTION(sza,inputparamG)
    u = np.sqrt((np.cos(cenit))**2 + ratio_rx_rz**2 * (np.sin(cenit))**2)
    b = R * u / np.cos(cenit)
# Mean interception over PEA(0) (=SA for no overlapping)
    AsX = g * LAD * 1.5 * v / (np.pi * R **2 * u)
    t = 2 * (1 - (1 + AsX) * np.exp(-AsX)) / (AsX ** 2)

# SA for a given reference surface: the projected shadow modified by mean interception
    a = np.sqrt(space_per_plant) / 2 / b
    if a < 1:
        f1 = R ** 2 * u / np.cos(cenit) * np.pi
        f2 = R ** 2 * u / np.cos(cenit) * (np.arccos(a) - a * np.sqrt(1 - a ** 2))
        SA = (f1 - 2 * f2) * (1 - t) + f2 * (1 - t ** 2)
    else:
        SA = R ** 2 * u / np.cos(cenit) * np.pi * (1 - t)
      

    SA = min(SA, space_per_plant)
    QDIR = SA / space_per_plant
      
    if pardr + PARDF > 0:
          TPAR = pardr + PARDF
          fdif = PARDF / TPAR
    
          RDIF = TPAR * fdif #W/m2 soil (PAR)
          Rdir = TPAR * (1 - fdif) #W/m2 soil (PAR)
    
    
          APARDIF = QDIF * RDIF  #W PAR/m2 soil
          APARDIR = QDIR * Rdir  #W PAR/m2 soil
    
    
    
#          Apar = AparDIR + APARDIF
          
          LAIsun = QDIR * np.cos(cenit) / g
          LAIshade = LAI - LAIsun
    
#          iparsun = APARDIR / LAIsun + APARDIF / LAI #W/m2 leaf
        #1 Watt/m2 PAR = 4.6667 mumol/m2/s PAR
#          iparsun = iparsun * 4.6667 # micromol m-2 leaf s-1
    
#          iparshade = APARDIF / LAI  #W/m2 leaf
#          iparshade = iparshade * 4.6667 #micromol m-2 leaf s-1
          
    return APARDIR,APARDIF,LAIsun,LAIshade

def GFUNCTION(sza,inputparamG):
    
    """
    Calculation of G function based on the procedure proposed by Lemeur,R.(1973).
    
    Inputs: 
    -------    
      sza: Solar zentih angle (degrees)  
      inputparamG:leaf inclination angle distribution function
      
    Return:
    -------        
        G_PROJ: G projection

     Reference: 
     ----------
        Lemeur,   R.,   1973.   A   method   for   simulating   the   direct   
            solar radiation  regime  in  sunflower,  Jerusalem  artichoke,  
            corn  and soybean canopies using actual stand structure data. 
            Agric. For.
            
    """
    sza=sza*np.pi/180
    i = 0
    XG = 0
    aIni=5 * np.pi / 180
    aEnd=75 * np.pi / 180
    aInt=10 * np.pi / 180
#!VA RECORRIENDO LOS DISTINTOS ANGULOS FOLIARES DE 5º A 85º
    for ANG in np.arange(aIni,aEnd,aInt):
        i = i + 1
#CALCULA ANGULO RESPECTO PLANO TRANSVERSAL RAYO SOLAR CON LAS MODIFICACIONES PERTINENTES Y POR CADA ANGULO CENITAL???
        f = np.cos(ANG) * np.cos(sza)
        if ANG + sza > np.pi / 2:
                COG = 1 / (np.tan(ANG)) / (np.tan(sza))
                COG = min(COG, 1)
                SG = np.sqrt(1 - COG ** 2)
                ag = np.arctan(SG / COG)
                a = (1 - 2 * ag / np.pi) * np.cos(ANG) * np.cos(sza)
                f = 2 / np.pi * np.sin(sza) * np.sin(ANG) * SG + a
           
          
        f = np.max(f, 0)
        
           
#VA ACUMULANDO EL VALOR DE G PARA TODOS LOS ANG Y POR CADA CENIT
        XG = XG + f * inputparamG[i]



    G_PROJ = XG
    
    return G_PROJ
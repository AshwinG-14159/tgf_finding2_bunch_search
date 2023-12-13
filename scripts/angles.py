from astropy.io import fits
import astropy.units as u
from astropy.table import Table
import astropy.coordinates as coo
import astropy.constants as con
from astropy.time import Time
import numpy as np
import argparse

def txy(mkffile, trigtime, ra_tran, dec_tran):
    """
    Calculate thetax, thetay using astropy
    Use pitch, roll and yaw information from the MKF file
    """
    # x = -yaw
    # y = +pitch
    # z = +roll

    # Read in the MKF file
    mkfdata = fits.getdata(mkffile, 1)
    # print(f'mkf data: {mkfdata}')
    sel = abs(mkfdata['time'] - trigtime) < 10
    
    try:
        ind = np.where(sel==True)[0][0]
    
    except IndexError:

        print("\nEither the time tagged data you are looking for is not in this file or CZTI might be in SAA. Please check using the script SAA_finder.py\n")
        raise SystemExit

    # Get pitch, roll, yaw
    # yaw is minus x
    try:
        pitch = coo.SkyCoord( np.median(mkfdata['pitch_ra'][sel]) * u.deg, np.median(mkfdata['pitch_dec'][sel]) * u.deg )
        roll = coo.SkyCoord( np.median(mkfdata['roll_ra'][sel]) * u.deg, np.median(mkfdata['roll_dec'][sel]) * u.deg )
        yaw_ra  = (180.0 + np.median(mkfdata['yaw_ra'][sel]) ) % 360
        yaw_dec = -np.median(mkfdata['yaw_dec'][sel])
        minus_yaw = coo.SkyCoord( yaw_ra * u.deg, yaw_dec * u.deg )

        # Transient:
        transient = coo.SkyCoord(ra_tran * u.deg, dec_tran * u.deg)

        #Earth:
        earthx = np.median(mkfdata['posx'][sel]) * u.km
        earthy = np.median(mkfdata['posy'][sel]) * u.km
        earthz = np.median(mkfdata['posz'][sel]) * u.km
        earth = coo.SkyCoord(-earthx, -earthy, -earthz, frame='icrs', representation_type='cartesian')
        earth_czti = np.sqrt(earthz**2 + earthy**2 + earthx**2)
        earth_transient = earth.separation(transient).value

        # Angles from x, y, z axes are:
        ax = minus_yaw.separation(transient)
        ay = pitch.separation(transient)
        az = roll.separation(transient)

        # the components are:
        cx = np.cos(ax.radian) # The .radian is not really needed, but anyway...
        cy = np.cos(ay.radian)
        cz = np.cos(az.radian)

        # Thetax = angle from z axis in ZX plane
        # lets use arctan2(ycoord, xcoord) for this
        thetax = u.rad * np.arctan2(cx, cz)
        thetay = u.rad * np.arctan2(cy, cz)
        phi = np.arctan2(cy,cx)
        phi_new = phi*(u.rad.to(u.deg))
        # Theta and phi of Earths Centre:
        axn = minus_yaw.separation(earth)
        ayn=pitch.separation(earth)
        azn=roll.separation(earth)

        cxn = np.cos(axn.radian)
        cyn = np.cos(ayn.radian)
        czn  = np.cos(azn.radian)
        phin = np.arctan2(cy,cx)
        phi_newn = phin*(u.rad.to(u.deg))
        print("Theta of earth is:",azn.value)
        print("Phi of Earth is:",phi_newn)

        try:
            for i,j in enumerate(phi_new):
                if (j < 0 ):
                    phi_new[i] = 360 + j
                else:
                    phi_new[i] = j
        except TypeError:
            if phi_new<0:
                phi_new = 360.0 + phi_new
            else:
                phi_new = phi_new

        earth_occult_angle = np.arcsin(con.R_earth/earth_czti)

    except RuntimeWarning:

        print("CZTI might be in SAA. Please check using the script SAA_finder.py")
        raise SystemExit

    return az.value, phi_new, thetax.to(u.deg).value, thetay.to(u.deg).value, minus_yaw, pitch, roll, transient, earth, earth_czti, earth_transient, earth_occult_angle.to(u.deg).value

#------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("mkf_file", help=".mkf file to be used for processing", type=str)
    parser.add_argument("ra", help="Source Right Ascension (eg 127.5d, 8h30m, 127d30m)", type=str)
    parser.add_argument("dec", help="Source Declination (eg -27.5d, -27d30m)", type=str)
    parser.add_argument("ut", help="Time in UT (eg 2015-11-30T09:08:49, 2015-11-30 09:08:49)", type=str)
    args = parser.parse_args()
    
    # time in spacecraft seconds
    mjd_t0 = Time(55197.00, format="mjd")
    trigtime = Time(args.ut)
    czti_time = (trigtime - mjd_t0).sec
    print(f"czti time: {czti_time}")

    # RA, dec in degrees
    try:
        ra = coo.Angle(args.ra)
    except u.UnitsError:
        ra = coo.Angle(args.ra, unit=u.deg)

    try:
        dec = coo.Angle(args.dec)
    except u.UnitsError:
        dec = coo.Angle(args.dec, unit=u.deg)

    # Now call the function:
    transient_theta, transient_phi, transient_thetax, transient_thetay, coo_x, coo_y, coo_z, coo_transient, earth, earth_czti, earth_transient, earth_occult_angle = txy(args.mkf_file, czti_time, ra.deg, dec.deg)

    tab = Table(names=['Time in UT','Time in CZTI seconds','Theta_x','Theta_y','Theta','Phi','Earth_transient_angle','Nominal theta', 'Nominal phi'],dtype=('S32','S32','S32','S32','S32','S32','S32', 'S32', 'S32'))
    tab.add_row(['{trigtime}'.format(trigtime=trigtime.iso),'{trigtime:0.1f}'.format(trigtime=czti_time),'{tx:0.6f}'.format(tx=transient_thetax),'{ty:0.6f}'.format(ty=transient_thetay),'{th:0.6f}'.format(th=transient_theta),'{ph:0.6f}'.format(ph=transient_phi),'{occ:0.6f}'.format(occ=earth_transient), '{ra_h:2.0f}:{ra_m:2.0f}:{ra_s:4.1f} ({ra_deg:0.3f})'.format(ra_h=coo_z.ra.hms[0], ra_m=coo_z.ra.hms[1], ra_s=coo_z.ra.hms[2], ra_deg=coo_z.ra.deg), '{dec_h:2.0f}:{dec_m:2.0f}:{dec_s:4.1f} ({dec_deg:0.3f})'.format(dec_h=coo_z.dec.dms[0], dec_m=abs(coo_z.dec.dms[1]), dec_s=abs(coo_z.dec.dms[2]), dec_deg=coo_z.dec.deg)])
    tab.write('table_angles.txt',format='ascii',overwrite=True)


    data = Table(names=['Parameter', 'Value'], dtype=['S32', 'S32'])
    data.add_row(["MKF file", "{mkf_file}".format(mkf_file=args.mkf_file)])
    data.add_row(["Time in UT", "{trigtime}".format(trigtime=trigtime.isot)])
    data.add_row(["Time in CZTI seconds", "{trigtime:0.1f}".format(trigtime=czti_time)])
    data.add_row(["Nominal CZTI RA ", "{ra_h:2.0f}h {ra_m:2.0f}m {ra_s:4.1f}s ({ra_deg:0.3f})".format(ra_h=coo_z.ra.hms[0], ra_m=coo_z.ra.hms[1], ra_s=coo_z.ra.hms[2], ra_deg=coo_z.ra.deg)])
    data.add_row(["Nominal CZTI Dec ", "{dec_h:2.0f}d {dec_m:2.0f}m {dec_s:4.1f}s ({dec_deg:0.3f})".format(dec_h=coo_z.dec.dms[0], dec_m=abs(coo_z.dec.dms[1]), dec_s=abs(coo_z.dec.dms[2]), dec_deg=coo_z.dec.deg)])
    data.add_row(["Transient RA", "{ra_h:2.0f}h {ra_m:2.0f}m {ra_s:4.1f}s ({ra_deg:0.3f})".format(ra_h=ra.hms[0], ra_m=ra.hms[1], ra_s=ra.hms[2], ra_deg=ra.deg)])
    data.add_row(["Transient Dec", "{dec_h:2.0f}d {dec_m:2.0f}m {dec_s:4.1f}s ({dec_deg:0.3f})".format(dec_h=dec.dms[0], dec_m=abs(dec.dms[1]), dec_s=abs(dec.dms[2]), dec_deg=dec.deg)])
    data.add_row(["Theta_x", "{tx:3.3f}".format(tx=transient_thetax)])
    data.add_row(["Theta_y", "{ty:3.3f}".format(ty=transient_thetay)])
    data.add_row(["Theta", "{th:3.3f}".format(th=transient_theta)])
    data.add_row(["Phi", "{ph:3.3f}".format(ph=transient_phi)])
    data.add_row(["Earth-transient angle", "{ang:3.3f}".format(ang=earth_transient)])
    data.add_row(["The earth occultation angle", "{ang:3.4f}".format(ang=earth_occult_angle)])
    data.add_row(["Earth_x", "{dis:3.3f}".format(dis=earth.x)])
    data.add_row(["Earth_y", "{dis:3.3f}".format(dis=earth.y)])
    data.add_row(["Earth_z", "{dis:3.3f}".format(dis=earth.z)])
    data.add_row(["The distance of earth from CZTI", "{dis:5.4f}".format(dis=earth_czti)])
    

    if (earth_transient>earth_occult_angle):
        print("Source is visible\n")
        data.add_row(["Is source visible", "Yes"])

    else:
        print("Source is Earth Occulted\n")
        data.add_row(["Is source visible", "No"])

    data.write("Angles.txt", delimiter='\t', format='ascii', overwrite=True)
    print(data)
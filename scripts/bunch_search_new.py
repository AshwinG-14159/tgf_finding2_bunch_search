import subprocess as subp
import argparse
import glob
import numpy as np
import pandas as pd
from astropy.time import Time
import pexpect
import os
import requests

import csv
import os
from astropy.io import fits
# import matplotlib
# matplotlib.use('QtAgg')
import matplotlib.pyplot as plt
from astropy.table import Table

import astropy.units as u

import time
from scipy.ndimage import label

import argparse
import angles
from scipy.interpolate import interp1d

import astropy.coordinates as coo

import glob

logfile = "../../big_data/extra_logs/log.txt"

def longest_continuous_zeros(arr):
    zero_mask = (arr == 0).astype(int)
    # print(arr, zero_mask)
    zero_groups = np.split(zero_mask, np.where(np.diff(zero_mask) != 0)[0] + 1)
    longest_zeros = max(zero_groups, key=len, default=np.array([]))
    return len(longest_zeros)


def log(text, file = logfile, renew = False):

    if renew:
        f = open(file, 'w')
    else:
        f = open(file, 'a')

    f.write(f"{text}\n")
    f.close()

path_to_orbitinfo = '../data/orbitinfo.csv'


gap = 1
total_attempts = 1000
version = 4
set_of_rows = []


with open(path_to_orbitinfo, 'r') as f:
    r = csv.reader(f)
    count=0
    for row in r:

        if(count<1):
            count+=1
            continue

        if row[0][-6:]!='level2' and row[0][:7] == '2023062':
            if(count%gap!=0):
                count+=1
                continue
            set_of_rows.append(row)
            count+=1
            if(count/gap>total_attempts): # leave 10 between any 2 and take a total of 200. Binning = 0.001
                break


start_text1 = f'{len(set_of_rows)} orbits to search through. These are the top {total_attempts} orbits in {path_to_orbitinfo} leaving a gap of {gap} which are valid'

print(start_text1)
log(start_text1, renew = True)

print('###########################################################')
path_data = '/media/czti/CIFT PC Backup/event_and_bunch'
path_mkfs = '../../big_data/mkfs'
path_plots = '../../big_data/plots'
path_results_csv = '../../big_data/results_csv'

pattern_bunch = "*_bunch.fits"
pattern_veto = "*_quad_clean.evt"
pattern_mkf = "*_level2.mkf"

binning = 0.005

my_type = 1
ind_quad_thres = 300
num_quads = 3
total_thresh = 1600
reg_around_peak = [0.1,10]

plot_stuff = False

count_orbits = 0
count_events_total = 0
count_events_with_mkf = 0

start_text2 = f"""parameters: binning={binning}, individual quadrant threshold = {ind_quad_thres}, num quads = {num_quads}, 
      regions to plot for = {reg_around_peak}, type = {my_type}"""
log(start_text2)

print(start_text2)

file_suffix = f"{my_type}_{ind_quad_thres}_{num_quads}_{binning}"

csv_file = open(f"{path_results_csv}/csv_results_{file_suffix}.csv","w")
writer = csv.writer(csv_file)


for row in set_of_rows:
    mkf_present = True
    row_data = row[0].split('_')
    date = row_data[0]
    orbit = row_data[-1]

    bunch_files = glob.glob(f'{path_data}/*{orbit}{pattern_bunch}')
    veto_files = glob.glob(f'{path_data}/*{orbit}{pattern_veto}')
    mkf_files = glob.glob(f'{path_mkfs}/*{orbit}{pattern_mkf}')
    if(mkf_files==[]):
        log(f"{date}, {orbit} mkf missing")
        mkf_present=False
    else:
        mkf = mkf_files[0]

    if(bunch_files==[] or veto_files == []):
        log(f"{date}, {orbit} bunch or veto missing")
        continue
    
    # if(orbit[:3]!='418'):
    #     continue
    count_orbits+=1
    print(f"trying for {date} {orbit}")
    log(f"trying for {date} {orbit}. {bunch_files[0]} {veto_files[0]} will be used")

    with fits.open(veto_files[0]) as hdul:
        d5 = Table(hdul[5].data)
        time_arr = np.array(d5['Time'])
        veto_arr = np.array(d5['VetoSpec'])
        quad_arr = np.array(d5['QuadID'])
        time_len = time_arr.shape[0]
        mark = time_len//4
        t1_veto = time_arr[:mark]
        t2_veto = time_arr[mark:2*mark]
        t3_veto = time_arr[2*mark:3*mark]
        t4_veto = time_arr[3*mark:4*mark]

        v1 = veto_arr[:mark]
        v2 = veto_arr[mark:2*mark]
        v3 = veto_arr[2*mark:3*mark]
        v4 = veto_arr[3*mark:4*mark]

        v = np.concatenate((v1[np.newaxis,:,:],v2[np.newaxis,:,:],v3[np.newaxis,:,:],v4[np.newaxis,:,:]), axis=0)
        total_v = np.sum(v, axis=2)




    with fits.open(bunch_files[0]) as hdul:
        d1 = Table(hdul[1].data)
        d2 = Table(hdul[2].data)
        d3 = Table(hdul[3].data)
        d4 = Table(hdul[4].data)
        head = hdul[0].header
        tstart = head['TSTARTI'] + head['TSTARTF']
        tstop = head['TSTOPI'] + head['TSTOPF']

    t1 = np.array(d1['Time'])
    t2 = np.array(d2['Time'])
    t3 = np.array(d3['Time'])
    t4 = np.array(d4['Time'])

    c1 = (t1<tstop-10) & (t1>tstart+10)
    c2 = (t2<tstop-10) & (t2>tstart+10)
    c3 = (t3<tstop-10) & (t3>tstart+10)
    c4 = (t4<tstop-10) & (t4>tstart+10)

    c1_veto = (t1_veto<tstop-10) & (t1_veto>tstart+10)
    c2_veto = (t2_veto<tstop-10) & (t2_veto>tstart+10)
    c3_veto = (t3_veto<tstop-10) & (t3_veto>tstart+10)
    c4_veto = (t4_veto<tstop-10) & (t4_veto>tstart+10)

    t1_veto = t1_veto[c1_veto]
    t2_veto = t2_veto[c2_veto]
    t3_veto = t3_veto[c3_veto]
    t4_veto = t4_veto[c4_veto]

    v1 = total_v[0][c1_veto]
    v2 = total_v[1][c2_veto]
    v3 = total_v[2][c3_veto]
    v4 = total_v[3][c4_veto]


    f1 = np.array(d1['NumEvent'], dtype=int)[c1]
    f2 = np.array(d2['NumEvent'], dtype=int)[c2]
    f3 = np.array(d3['NumEvent'], dtype=int)[c3]
    f4 = np.array(d4['NumEvent'], dtype=int)[c4]

    t1 = t1[c1]
    t2 = t2[c2]
    t3 = t3[c3]
    t4 = t4[c4]

    gap = tstop - tstart

    num_bins = int(gap/binning + 1)
    desired_time_bins = np.linspace(tstart, tstop, num_bins)

    a1, d1 = np.histogram(t1, bins=desired_time_bins, weights=f1)
    a2, _ = np.histogram(t2, bins=desired_time_bins, weights=f2)
    a3, _ = np.histogram(t3, bins=desired_time_bins, weights=f3)
    a4, _ = np.histogram(t4, bins=desired_time_bins, weights=f4)


    desired_time_bins = desired_time_bins[:-1]
    a_arr = np.array([a1,a2,a3,a4])


    if(my_type == 1):
        vals = a_arr>ind_quad_thres
        vals = np.sum(vals, axis=0)
        cond = vals>=num_quads
        
        d = desired_time_bins[cond]


        for timestamp in d:
            print(f"peak at {timestamp}")
            log(f"time - {d}")
            written_to_csv = False

            for reg in reg_around_peak:

                # Veto prep:
                c1_vetox = (t1_veto<timestamp+5*reg) & (t1_veto>timestamp-5*reg)
                c2_vetox = (t2_veto<timestamp+5*reg) & (t2_veto>timestamp-5*reg)
                c3_vetox = (t3_veto<timestamp+5*reg) & (t3_veto>timestamp-5*reg)
                c4_vetox = (t4_veto<timestamp+5*reg) & (t4_veto>timestamp-5*reg)

                t1_vetox = t1_veto[c1_vetox]
                t2_vetox = t2_veto[c2_vetox]
                t3_vetox = t3_veto[c3_vetox]
                t4_vetox = t4_veto[c4_vetox]

                v1x = v1[c1_vetox]
                v2x = v2[c2_vetox]
                v3x = v3[c3_vetox]
                v4x = v4[c4_vetox]


                cx = (desired_time_bins>timestamp-reg) & (desired_time_bins<timestamp+reg)
                timesx = desired_time_bins[cx]
                a1x = a_arr[0][cx]
                a2x = a_arr[1][cx]
                a3x = a_arr[2][cx]
                a4x = a_arr[3][cx]

                region_around = np.where(desired_time_bins==timestamp)[0][0]


                pre1 = a_arr[0][(desired_time_bins>timestamp-0.5) & (desired_time_bins<timestamp)]
                pre2 = a_arr[1][(desired_time_bins>timestamp-0.5) & (desired_time_bins<timestamp)]
                pre3 = a_arr[2][(desired_time_bins>timestamp-0.5) & (desired_time_bins<timestamp)]
                pre4 = a_arr[3][(desired_time_bins>timestamp-0.5) & (desired_time_bins<timestamp)]
                l1 = binning * longest_continuous_zeros(pre1)
                l2 = binning * longest_continuous_zeros(pre2)
                l3 = binning * longest_continuous_zeros(pre3)
                l4 = binning * longest_continuous_zeros(pre4)


                if(l1>0.25 and l2>0.25 and l3>0.25 and l4>0.25):
                    continue
                print(f"will be plotted")
                log(f"will be plotted")

                if not written_to_csv:
                    count_events_total+=1
                    if(not mkf_present):
                        row_to_write = [date, orbit, timestamp, -1, -1]
                        writer.writerow(row_to_write)
                        csv_file.flush()
                    else:
                        count_events_with_mkf+=1
                        czti_time = timestamp
                        with fits.open(mkf) as hdul:
                            h = hdul[0].header
                            ra_pnt = f'{h["RA_PNT"]}d'
                            dec_pnt = f'{h["DEC_PNT"]}d'
                    
                        try:
                            ra_pnt = coo.Angle(ra_pnt)
                        except u.UnitsError:
                            ra_pnt = coo.Angle(ra_pnt, unit=u.deg)

                        try:
                            dec_pnt = coo.Angle(dec_pnt)
                        except u.UnitsError:
                            dec_pnt = coo.Angle(dec_pnt, unit=u.deg)

                        transient_theta, transient_phi, transient_thetax, transient_thetay, coo_x, coo_y, coo_z, coo_transient, earth, earth_czti, earth_transient, earth_occult_angle, phi_newn = angles.txy(mkf, czti_time, ra_pnt.deg, dec_pnt.deg)
                        
                        row_to_write = [date, orbit, timestamp, earth_transient, phi_newn]
                        writer.writerow(row_to_write)
                        csv_file.flush()


                if(plot_stuff):
                    fig, axarr = plt.subplots(5, sharex=True, figsize=(7, 9))
                    fig.subplots_adjust(hspace=0.5)

                    axarr[4].plot(t1_vetox, v1x, label='A', drawstyle='steps')
                    axarr[4].plot(t2_vetox, v2x, label='B', drawstyle='steps')
                    axarr[4].plot(t3_vetox, v3x, label='C', drawstyle='steps')
                    axarr[4].plot(t4_vetox, v4x, label='D', drawstyle='steps')
                    axarr[4].set_title('Veto')
                    axarr[4].grid(True)
                    axarr[4].set_ylabel(f'Veto Counts')
                    axarr[4].legend()



                    axarr[0].plot(timesx, a1x)
                    axarr[0].set_title('quad A')
                    axarr[0].set_xlim(timestamp-reg, timestamp+reg)
                    axarr[0].grid(True)
                    axarr[0].set_ylabel(f'Bunch Counts')
                    axarr[0].axhline(y=ind_quad_thres, color='r', linestyle='--', label='cutoff', alpha=0.5)

                    axarr[1].plot(timesx, a2x)
                    axarr[1].set_title('quad B')
                    axarr[1].grid(True)
                    axarr[1].set_ylabel(f'Bunch Counts')
                    axarr[1].axhline(y=ind_quad_thres, color='r', linestyle='--', label='cutoff', alpha=0.5)

                    axarr[2].plot(timesx, a3x)
                    axarr[2].set_title('quad C')
                    axarr[2].grid(True)
                    axarr[2].set_ylabel(f'Bunch Counts')
                    axarr[2].axhline(y=ind_quad_thres, color='r', linestyle='--', label='cutoff', alpha=0.5)

                    axarr[3].plot(timesx, a4x)
                    axarr[3].set_title('quad D')
                    axarr[3].grid(True)
                    axarr[3].set_ylabel(f'Bunch Counts')
                    axarr[3].axhline(y=ind_quad_thres, color='r', linestyle='--', label='cutoff', alpha=0.5)

                    plt.xlabel(f'Time(Binning = {binning})')
                    plt.suptitle(f"Bunch and Veto Quads - Eyeballed Veto - Orbit {orbit}")


                    plt.savefig(f'{path_plots}/{file_suffix}/{row[0][-5:]}_{binning}_{timestamp}_{reg}_{ind_quad_thres}_quads.png')

                    plt.cla()

                    fig, axarr = plt.subplots(2, sharex=True, figsize=(7, 9))
                    fig.subplots_adjust(hspace=0.5)

                    axarr[0].plot(t1_vetox, v1x, label='A', drawstyle='steps')
                    axarr[0].plot(t2_vetox, v2x, label='B', drawstyle='steps')
                    axarr[0].plot(t3_vetox, v3x, label='C', drawstyle='steps')
                    axarr[0].plot(t4_vetox, v4x, label='D', drawstyle='steps')
                    axarr[0].set_title('Veto')
                    axarr[0].grid(True)
                    axarr[0].set_ylabel(f'Veto Counts')
                    axarr[0].legend()

                    a_sum = np.sum(np.array([a1x,a2x,a3x,a4x]), axis=0)


                    axarr[1].plot(timesx, a_sum)
                    axarr[1].set_title('4 Quadrant Summed Bunch Data')
                    axarr[1].set_ylabel(f'Bunch Counts')
                    axarr[1].set_xlim(timestamp-reg, timestamp+reg)
                    axarr[1].grid(True)

                    plt.xlabel(f'Time(Binning = {binning})')
                    plt.suptitle(f"4 quadrant Summed Bunch Data and Veto Quads Eyeballed Veto - Orbit {orbit}")


                    plt.savefig(f'{path_plots}/{file_suffix}/{row[0][-5:]}_{binning}_{timestamp}_{reg}_{ind_quad_thres}_sums.png')
                    plt.cla()

                    plt.close('all')


csv_file.close()











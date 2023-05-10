# -*- coding: utf-8 -*-
"""
Created on Tue May  9 10:13:04 2023

Continous Wave Calibration Read In Functions

Functions for Reading in and integrating the power meter as well as the data
produced by Phillip Black's L-BASS.sh analysis program.

Contact: jordan.norris@student.manchester.ac.uk
        jmnorris98@gmail.com 

@author: Jordan Norris
"""

import numpy as np
from astropy.time import Time
import csv
from datetime import datetime, timedelta, date, time
import os

def zip_rpg_and_pm(rpg_file_path, pm_file_path, starttime, duration, inputhorn, 
                   powermeterID, averaging_time):
    """
    This function takes in the filepath of the various file locations as well
    as the start time, observing duration, input horn and power meter ID as
    well as averaging time and returns minute averaged values with corresponding
    times, CW powers for each spectrometer input for the input horn.

    Parameters
    ----------
    rpg_file_path : string
        String of the filepath for the Phillip Black produced code (may be the
        temp file on moonhut when integrating into his package).
    pm_file_path : string
        String of the file location containing the Power Meter .csv files.
    starttime : Datetime Object
        Datetime Object Corresponding to the start of the observation.
    duration : TimeDelta Object
        Timedelta object corresponding to the lenght of the duration.
    inputhorn : string
        Use 'W' for west horn and 'E' for the east horn.
    powermeterID : integer
        Use the integer value of the power meter ID corresponding to which was
        used for the respective horn.
    averaging_time : float/integer, optional
        Integer or float value of the time desired for the data to be averaged
        over in seconds.

    Returns
    -------
    out_time : DateTime Array
        Datetime array of the corresponding time of the beginning of each
        instance of averaged data.
    pm_cw : Array
        Array of integers corresponding to the averaged PowerMeter
        measured CW signal.
    
    cw_output_1
        Array of CW powers calculated for each averaging time for one
        spectrometer output.
    cw_output_2
        Array of CW powers calculated for each averaging time for the other
        spectrometer output.
    rpg_output_1
        The CW subtracted integrated signal for one spectometer output.
    rpg_output_1
        The CW subtracted integrated signal for the other spectometer output.

    """
    endtime = starttime + duration
    
    compiled_pm_powers, compiled_pm_times = PowerMeterCompiler(starttime, duration, pm_file_path, powermeterID)
    #Temporary BST Fix
    #compiled_pm_times = np.array([i-timedelta(hours=1) for i in compiled_pm_times])
    averaging_time = timedelta(seconds=averaging_time)
    rpg_times, a1p1, a2p2, a1p2, a2p1 = RPGCompiler(rpg_file_path)
    if inputhorn == 'W':
        a1p1_total_times, a1p1_cw_power_array,\
        a1p1_total_minus_cw_array, a1p1_cw_indices = CW_Extraction(rpg_times, a1p1)
        
        a2p2_total_times, a2p2_cw_power_array,\
        a2p2_total_minus_cw_array, a2p2_cw_indices = CW_Extraction(rpg_times, a2p2)
    
    
    
    
        out_time, a1p1_rpg_cw, a2p2_rpg_cw, pm_cw = np.array([]), np.array([]), np.array([]), np.array([])
        out_a1p1_wo_cw, out_a2p2_wo_cw = np.array([]), np.array([])
    
    
        while starttime < endtime:
            indices_pm = np.where(np.logical_and(compiled_pm_times>=starttime, 
                                                 compiled_pm_times<starttime+averaging_time))
            small_pm_powers = compiled_pm_powers[indices_pm]
            
            pm_avgpower = np.average(small_pm_powers)
            pm_cw = np.append(pm_cw, pm_avgpower)
        
            #rmsd = np.std(small_powers)
            indices_a1p1 = np.where(np.logical_and(a1p1_total_times>=starttime, 
                                                   a1p1_total_times<starttime+averaging_time))
            indices_a2p2 = np.where(np.logical_and(a2p2_total_times>=starttime, 
                                                   a2p2_total_times<starttime+averaging_time))
            
            
            small_a1p1_cw = a1p1_cw_power_array[indices_a1p1]
            small_a2p2_cw = a2p2_cw_power_array[indices_a2p2]
            a1p1_cw_avgpower = np.average(small_a1p1_cw)
            a2p2_cw_avgpower = np.average(small_a2p2_cw)
            a1p1_rpg_cw = np.append(a1p1_rpg_cw, a1p1_cw_avgpower)
            a2p2_rpg_cw = np.append(a2p2_rpg_cw, a2p2_cw_avgpower)
            # print('a1p1_cw')
            # print(small_a1p1_cw)
            small_a1p1_wo_cw = a1p1_total_minus_cw_array[indices_a1p1]
            small_a2p2_wo_cw = a2p2_total_minus_cw_array[indices_a2p2]
            
            av_cw_sub_a1p1 = np.average(small_a1p1_wo_cw)
            av_cw_sub_a2p2 = np.average(small_a2p2_wo_cw)
            out_a1p1_wo_cw = np.append(out_a1p1_wo_cw, av_cw_sub_a1p1)
            out_a2p2_wo_cw = np.append(out_a2p2_wo_cw, av_cw_sub_a2p2)
            
            # small_tot = total_power_array[indices_rpg]
            # av_tot = np.average(small_tot)
            # out_total_power = np.append(out_total_power, av_tot)
            
            out_time = np.append(out_time, starttime)
            
            starttime += averaging_time
                    
                    
        return out_time, pm_cw, a1p1_rpg_cw, a2p2_rpg_cw, out_a1p1_wo_cw, out_a2p2_wo_cw
    else:
        a1p2_total_times, a1p2_cw_power_array,\
        a1p2_total_minus_cw_array, a1p2_cw_indices = CW_Extraction(rpg_times, a1p2)
        
        a2p1_total_times, a2p1_cw_power_array,\
        a2p1_total_minus_cw_array, a2p1_cw_indices = CW_Extraction(rpg_times, a2p1)
    
    
    
    
        out_time, a1p2_rpg_cw, a2p1_rpg_cw, pm_cw = np.array([]), np.array([]), np.array([]), np.array([])
        out_a1p2_wo_cw, out_a2p1_wo_cw = np.array([]), np.array([])
    
    
        while starttime < endtime:
            indices_pm = np.where(np.logical_and(compiled_pm_times>=starttime, 
                                                 compiled_pm_times<starttime+averaging_time))
            small_pm_powers = compiled_pm_powers[indices_pm]
            # print(starttime)
            # print(starttime+averaging_time)
            pm_avgpower = np.average(small_pm_powers)
            pm_cw = np.append(pm_cw, pm_avgpower)
        
            #rmsd = np.std(small_powers)
            indices_a1p2 = np.where(np.logical_and(a1p2_total_times>=starttime, 
                                                   a1p2_total_times<starttime+averaging_time))
            
            indices_a2p1 = np.where(np.logical_and(a2p1_total_times>=starttime, 
                                                   a2p1_total_times<starttime+averaging_time))
            
            
            small_a1p2_cw = a1p2_cw_power_array[indices_a1p2]
            
            small_a2p1_cw = a2p1_cw_power_array[indices_a2p1]
            a1p2_cw_avgpower = np.average(small_a1p2_cw)
            a2p1_cw_avgpower = np.average(small_a2p1_cw)
            a1p2_rpg_cw = np.append(a1p2_rpg_cw, a1p2_cw_avgpower)
            a2p1_rpg_cw = np.append(a2p1_rpg_cw, a2p1_cw_avgpower)
            
            small_a1p2_wo_cw = a1p2_total_minus_cw_array[indices_a1p2]
            small_a2p1_wo_cw = a2p1_total_minus_cw_array[indices_a2p1]
            
            av_cw_sub_a1p2 = np.average(small_a1p2_wo_cw)
            av_cw_sub_a2p1 = np.average(small_a2p1_wo_cw)
            out_a1p2_wo_cw = np.append(out_a1p2_wo_cw, av_cw_sub_a1p2)
            out_a2p1_wo_cw = np.append(out_a2p1_wo_cw, av_cw_sub_a2p1)
            
            # small_tot = total_power_array[indices_rpg]
            # av_tot = np.average(small_tot)
            # out_total_power = np.append(out_total_power, av_tot)
            
            out_time = np.append(out_time, starttime)
            
            starttime += averaging_time
                    
                    
        return out_time, pm_cw, a1p2_rpg_cw, a2p1_rpg_cw, out_a1p2_wo_cw, out_a2p1_wo_cw
    
def RPGCompiler(rpgdatapath):
    """
    Reads in the data for the Phillip Black produced frequency-power arrays
    from a1p1_binned.npy, located in the Temp folder in moonhut.
    
    Need: a1p1_binned.npy, a2p2_binned.npy, a1p2_binned.npy, a2p1_binned.npy
    and obshdr.npy in order to read in the data and assign the correct times.
    These are generated when running LBASS.sh on moonhut through a vnc or other
    wise and are found in the temp folder in the pblack directory with the other
    LBASS.sh data.

    Parameters
    ----------
    rpgdatapath : string
        string of the the filepath for the files containing the various .npy
        files genereted by LBASS.sh (may be the temp folder in moonhut).

    Returns
    -------
    datetimes_array : datetime array
        Array of datetime objects corresponding to the time and date of each
        minute of data.
    a1p1_powers : array
        Array of power values where horizontal elements correspond to the
        seperate frequency channels and vertically are the minute by minute data.
    a2p2_powers : array
        Array of power values where horizontal elements correspond to the
        seperate frequency channels and vertically are the minute by minute data..
    a1p2_powers : array
        Array of power values where horizontal elements correspond to the
        seperate frequency channels and vertically are the minute by minute data..
    a2p1_powers : array
        Array of power values where horizontal elements correspond to the
        seperate frequency channels and vertically are the minute by minute data..

    """
    #user_inputs = np.load(rpgdatapath+'inputs.npy')
    a1p1_binned = np.load(rpgdatapath+'a1p1_binned.npy')
    a2p2_binned = np.load(rpgdatapath+'a2p2_binned.npy')
    a1p2_binned = np.load(rpgdatapath+'a1p2_binned.npy')
    a2p1_binned = np.load(rpgdatapath+'a2p1_binned.npy')
    obs_hdr = np.load(rpgdatapath+'obshdr.npy')
    #run_start_date = str(np.load(rpgdatapath+'run_start_date.npy', allow_pickle=True))
    #run_start_date = datetime.strptime(run_start_date, "%Y-%m-%d")
    MJD = Time(obs_hdr[0,8], format='mjd', scale='utc', precision=9)
    
    print('Run Start Date')
    print(MJD)
    run_start_date = MJD.utc.datetime
    print(run_start_date)
    times_since_midnight_ds = [timedelta(seconds=i) for i in a1p1_binned[:,0]]
    
    datetimes_array = np.array([])
    for time_sm in times_since_midnight_ds:
        #print(time_sm)
        date_time_obj = run_start_date + time_sm
        #print(date_time_obj)
        datetimes_array = np.append(datetimes_array, date_time_obj)
    
    a1p1_powers = a1p1_binned[:, 3:]
    a2p2_powers = a2p2_binned[:, 3:]
    a1p2_powers = a1p2_binned[:, 3:]
    a2p1_powers = a2p1_binned[:, 3:]
    
    #print(datetimes_array)
    return datetimes_array, a1p1_powers, a2p2_powers, a1p2_powers, a2p1_powers
    
def CW_Extraction(datetimes_array, powers_array):
    """
    Extracts the CW signal from the frequency power time file by locating the
    largest signal for each minute in time, isolating and then subtracting the
    signal detected in channels with no CW power from it to produce just the CW
    powers.

    Parameters
    ----------
    datetimes_array : datetime array
        Datetime array corresponding to the averaged spectrometer data.
    powers_array : array
        The array of powers detected for each frequency channel for each minute
        in time.

    Returns
    -------
    total_times : datetime array
        Datetime object corresponding to each measurement.
    cw_power_array : array
        array of CW measurements using the side channel subtraction method.
    total_minus_cw_array : array
        integrated power across the entire bandwidth (may need changing
        to 1400 - 1427 MHz) for each minute with the CW signal subtracted.
    cw_indices : array
        The indices corresponding to the frequency channel in which the CW
        signal is dominant.

    """
    cw_indices = np.array([])
    cw_power_array = np.array([])
    total_minus_cw_array = np.array([])
    total_times = np.array([])
    
    for i in range(0, len(powers_array[:,0])): #row by row i
        max_index = np.where(powers_array[i,:] == np.max(powers_array[i,:]))[0]
        print('Max Index:' +str(max_index))
        cw_row_indices = np.array([max_index-2,max_index-1, 
                                   max_index, max_index+1, max_index+2])
        
        cw_channels = powers_array[i, cw_row_indices]
        #print('Original:')
        #print(cw_channels)
        background_indeces = np.array([max_index-3, max_index-4, max_index-5,
                                       max_index+3, max_index+4, max_index+5])
        
        background_channels = powers_array[i, background_indeces]
        #print(background_channels)
        background_average = np.average(background_channels)
        #print(background_average)
        
        
        cw_channels = np.array([i-background_average for i in cw_channels])
        #print('Background Subtracted:')
        #print(cw_channels)
        cw_signal = np.sum(cw_channels)
        #print(cw_signal)
        sig_wo_cw = np.sum(powers_array[i,:]) - cw_signal
        
        cw_indices = np.append(cw_indices, max_index)
        cw_power_array = np.append(cw_power_array, cw_signal)
        total_times = np.append(total_times, datetimes_array[i])
        total_minus_cw_array = np.append(total_minus_cw_array, sig_wo_cw)
    
    return total_times, cw_power_array, total_minus_cw_array, cw_indices
    
def PowerMeterCompiler(starttime, duration, powermeterfile, powermeterID):
    """
    Function Compiles PowerMeter Data From Several .csv files produced by TAP.py
    Or for the regular PowerXpert output .csv files. 
    Parameters
    ----------
    starttime : datetime object
        corresponds to the start of the desired observing period
    duration : timedelta object
        correpsonds to the amount of time of the observation
    powermeterfile : string
        filepath of the folder containing the .csv files
        
    powermeterID : integer
        serial number of the desired power meter (the number leading the date in
        the .csv files)
    

    Returns
    -------
    compiled_pm_powers : numpy array
        array of floats corresponding to the raw powers
        
    compiled_pm_times : numpy array
        array of corresponding datetime objects
    """
    buffer = timedelta(minutes = 60) # change this to suit file size / may be 
    #change buffer to suit size of .csv files
    endtime = starttime + duration
    
    filecontents = np.array(os.listdir(powermeterfile))
    
    #print(starttime)
    timelist = np.array([convert_filename_to_datetime(csvname) for csvname in filecontents])
    
    
    filename_indices_times = np.where(np.logical_and(timelist>=starttime-buffer, 
                                                       timelist<endtime+buffer))
    #print(filename_indices_times)
    
    files_to_read = filecontents[filename_indices_times]
    
    
    serial_number_list = np.array([getPowerMeterID(powermeterfile+'/'+csvname) for csvname in files_to_read])
    
    files_to_read = files_to_read[np.where(serial_number_list == powermeterID)]
    
    compiled_pm_powers, compiled_pm_times = np.array([]), np.array([])
    
    for csvname in files_to_read:
        mw_powers, pm_time_array = PowerMeterFileReader(powermeterfile+'/'+csvname)
        compiled_pm_powers = np.append(compiled_pm_powers, mw_powers)
        compiled_pm_times = np.append(compiled_pm_times, pm_time_array)
        
    
    
    observation_indices = np.where(np.logical_and(compiled_pm_times>=starttime,
                                                  compiled_pm_times<endtime))
    
    compiled_pm_powers = compiled_pm_powers[observation_indices]
    compiled_pm_times = compiled_pm_times[observation_indices]
    
    return compiled_pm_powers, compiled_pm_times

def PowerMeterFileReader(anritsufilepath):
    """
    Reads in the Antritsu power meter data and then converts the powers into
    linear mW as well as adding proper time formats to each measurement so that
    they may be used in conjunction with the RPG data. The date is taken from
    RPG data assuming both measurements were initialised on the same day
    -------
    anritsufilepath: string
        string corresponding to the filepath of the .csv file that needs to be
        read in.
    
    Returns
    -------
    mw_powers : array
        A numpy array of power data in mW
    pm_time_array : array
        A numpy array of datetime data corresponding to the powers.

    """
    dbm_powers = np.genfromtxt(anritsufilepath ,delimiter=',',
                          skip_header=7, usecols=1)
    mw_powers = dBmTOmW(dbm_powers) #power array in milliWatts
    
    del dbm_powers
    
    times = []
    with open(anritsufilepath) as file:
        data = csv.reader(file)
        for row in data:
            times.append(str(row[0]))
    times = np.array(times[7:]) #gathers times as aray of strings
    
    #split anritsufilepath string with /
    filename = anritsufilepath.split('/')[-1]
    filename = filename.split('_')
    
    # calendardate = obsheader[0,7]
    # calendardate = datetime.strptime(calendardate, "isot")
    # print(calendardate)
    obs_date = datetime(int(filename[1]), int(filename[2]), int(filename[3]))
    
    pm_time_array = []
    firstline = True
    
    for i in range(0, np.size(times)):
        if firstline == True:
            timeformat = convertHHMMSStoFullDate(obs_date, times[i])
            pm_time_array.append(timeformat)
            firstline = False
        else:
            if convertHHMMSStoFullDate(obs_date, times[i]).hour < convertHHMMSStoFullDate(obs_date, times[i-1]).hour:
                obs_date = obs_date + timedelta(days=1)
                timeformat = convertHHMMSStoFullDate(obs_date, times[i])
                pm_time_array.append(timeformat)
            else:
                timeformat = convertHHMMSStoFullDate(obs_date, times[i])
                pm_time_array.append(timeformat)
    pm_time_array = np.array(pm_time_array)
    #Temp BST Fix
    #pm_time_array = np.array([i-timedelta(hours=1) for i in pm_time_array])
    
    #pm_time_arrray = pm_time_array
            
    return mw_powers, pm_time_array

def dBmTOmW(x):
    return 10**(x/10)

def convertHHMMSStoFullDate(calendardate, timestring):
    """
    Adds the time given from the .csv file onto the calander date that the
    code is on at that instance (due to calander date not being included
    in the timecode for each measurement)

    Parameters
    ----------
    calendardate : datetime object
        A datetime object corresponding to 00:00 on the calender date.
    timestring : string
        The string corresonding to HH:MM:SS read in from the .csv file.

    Returns
    -------
    timeformat : datetime object
        Datetime object corresponding to the date and time for a given
        power meter measurement.

    """
    t = datetime.strptime(timestring, "%H:%M:%S")
    dt = timedelta(hours=t.hour, minutes=t.minute, seconds=t.second)
    timeformat = calendardate + dt
    return timeformat

def getPowerMeterID(filepath):
    """
    Retreives the powermeter ID from the header of the power meter .csv file.

    Parameters
    ----------
    filepath : String
        Filepath of the .csv file.

    Returns
    -------
    number : Integer
        Power Meter ID as an integer.

    """
    with open(filepath) as file:
        data = csv.reader(file)
        
        for row in data:
            if row[0] == 'Sensor Serial No':
                number = int(row[1])
                break
            else:
                pass
        return number

def convert_filename_to_datetime(csvname):
    """
    Converts the filename of a power meter file to a datetime object.

    Parameters
    ----------
    csvname : String
        Filename produced from PowerXpert in the unaltered format.

    Returns
    -------
    csv_final_datetime : Datetime object
        Datetime object corresponding to the start of data logging for a given
        file.

    """
    csvname = str(csvname)
    csvname = csvname.split('_')
    csv_time = csvname[-1].split('.')[0]
    
    if len(csv_time) == 5:
        csv_time = timedelta(hours=int(csv_time[0]),minutes=int(csv_time[1:3]), 
                             seconds=int(csv_time[3:]))
    else:
        csv_time = timedelta(hours=int(csv_time[:2]),minutes=int(csv_time[2:4]), 
                             seconds=int(csv_time[4:]))
    
    csv_date = datetime(int(csvname[1]), int(csvname[2]), int(csvname[3]))
    
    csv_final_datetime = csv_date + csv_time
    return csv_final_datetime
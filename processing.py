#!/usr/bin/python3
# This is a script with functions needed for the data reduction, use run_process.py to actually analyse
# Updated from B. Quici script By K.Ross 19/5/21

import os
from astropy.io import fits 
import astropy.units as u 
from astropy.table import Table 
from matplotlib import rcParams
from casatasks import (
    flagmanager,
    flagdata,
    importatca,
    mstransform,
    listobs,
    setjy,
    gaincal,
    bandpass,
    fluxscale,
    applycal,
    tclean,
    rmtables,
    impbcor,
    split,
    uvmodelfit,
    exportfits,
)
import numpy as np
from casaplotms import plotms
import matplotlib.pyplot as plt
from casatools import image as IA
import logging 
from argparse import ArgumentParser
import datetime 

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(module)s:%(levelname)s:%(lineno)d %(message)s")
logger.setLevel(logging.INFO)


rcParams['font.family'] = 'serif'
s, dt, axlab = 8, 1.1, 12.
plt.rcParams["xtick.major.size"] = s
plt.rcParams["xtick.minor.size"] = s
plt.rcParams["ytick.major.size"] = s
plt.rcParams["ytick.minor.size"] = s
plt.rcParams["xtick.major.width"] = dt
plt.rcParams["xtick.minor.width"] = dt
plt.rcParams["ytick.major.width"] = dt
plt.rcParams["ytick.minor.width"] = dt
plt.rcParams["xtick.direction"] = 'in'
plt.rcParams["ytick.direction"] = 'in'
plt.rcParams["xtick.major.pad"] = 5.
plt.rcParams["figure.figsize"] = [10., 4.5]

def make_ms(files, visname):
    importatca(
        vis=visname, files=files, options="birdie,noac", edge=4
    )
    return 


def flag_ms(visname):
    logger.debug("Flagging antennas affected by shadowing...")

    flagmanager(vis=visname, mode="save", versionname="before_online_flagging")
    logger.debug("Flagging antennas affected by shadowing...")
    flagdata(vis=visname, mode="shadow", tolerance=0.0, flagbackup=False)
    logger.debug("Flagging visibilities with zero amplitudes...")
    flagdata(vis=visname, mode="clip", clipzeros=True, flagbackup=False)
    logger.debug("Quacking visibilities ...")
    flagdata(
        vis=visname, mode="quack", quackinterval=5.0, quackmode="beg", flagbackup=False
    )
    flagmanager(vis=visname, mode="save", versionname="after_online_flagging")
    flagdata(
        vis=visname,
        mode="tfcrop",
        datacolumn="data",
        action="apply",
        display="report",
        flagbackup=True,
        extendpols=True,
        correlation="",
        flagdimension="freqtime",
        growtime=95.0,
        growfreq=95.0,
        timecutoff=4.0,
        timefit="line",
        freqfit="poly",
        maxnpieces=5,
        combinescans=False,
        ntime="scan",
        extendflags=False,
    )
    logger.debug("Extending flags to all correlations")
    flagdata(
        vis=visname,
        mode="extend",
        action="apply",
        display="report",
        flagbackup=False,
        extendpols=True,
        correlation="",
        growtime=95.0,
        growfreq=95.0,
        growaround=True,
        flagneartime=False,
        flagnearfreq=False,
        combinescans=False,
        ntime="scan",
    )
    return

def split_ms(visname, msname, field="", spw="", n_spw=1, antenna="", scan = "", datacolumn="corrected",listfile=""):
    os.system(f"rm -r {msname}")
    os.system(f"rm -r {msname}.flagversions")
    os.system("rm -r *.last")
    try: 
        mstransform(
            vis=visname,
            outputvis=msname,
            field=field,
            nspw=n_spw,
            regridms=True,
            antenna=antenna,
            scan=scan,
            datacolumn=datacolumn,
            spw=spw
        )
    except: 
        logger.warning("Unable to run mstransform! Check logs ")


    if listfile == "":
        logger.debug("Not listing ms details")
    else: 
        listobs(
            vis=msname,
            listfile=f"{listfile}.dat",
            overwrite=True,
        )
    flagmanager(vis=msname, mode="save", versionname="after_transform")
    return

def calibrate_ms(msname, sec, calfile, pri = "1934_cal_cx", ref = "CA04", solint="60s"):
    if os.path.exists(f"{calfile}.G2"):
        logger.debug("Already found cal file? Not going to overwrite it. Make sure generating correct caltables ")
    else: 
        setjy(
            vis=msname,
            field=pri,
            scalebychan=True,
            standard="Perley-Butler 2010",
            usescratch=True,
        )
        logger.debug(f"Performing gain calibration on {pri}")
        try: 
            gaincal(
                vis=msname,
                caltable=f"{calfile}.G0",
                field=pri,
                refant=ref,
                gaintype="G",
                calmode="p",
                parang=True,
                # minblperant=3,
                solint=solint,
            )
            logger.debug(f"Performing bandpass calibration on {pri}")
            bandpass(
                vis=msname,
                caltable=f"{calfile}.B0",
                field=pri,
                refant=ref,
                solnorm=True,
                solint="inf",
                bandtype="B",
                gaintable=[f"{calfile}.G0"],
                parang=True,
            )
            logger.debug(f"Determining gains on {sec}")
            gaincal(
                vis=msname,
                caltable=f"{calfile}.G1",
                field=f"{pri},{sec}",
                refant=ref,
                gaintype="G",
                calmode="ap",
                parang=True,
                solint="120s",
                gaintable=[f"{calfile}.B0"],
            )
            bandpass(
                vis=msname,
                caltable=f"{calfile}.B1",
                field=pri,
                refant=ref,
                solnorm=True,
                solint="120s",
                bandtype="B",
                gaintable=[f"{calfile}.G1"],
                parang=True,
            )
            logger.debug(f"Deriving gain calibration using {pri}")
            gaincal(
                vis=msname,
                caltable=f"{calfile}.G2",
                field=pri,
                refant=ref,
                gaintype="G",
                calmode="ap",
                parang=True,
                solint=solint,
                # minblperant=3,
                gaintable=[f"{calfile}.B1"],
            )
            logger.debug(f"Deriving gain calibration using {sec}")
            gaincal(
                vis=msname,
                caltable=f"{calfile}.G2",
                field=sec,
                refant=ref,
                gaintype="G",
                calmode="ap",
                parang=True,
                solint=solint,
                # minblperant=3,
                gaintable=[f"{calfile}.B1"],
                append=True,
            )
            logger.debug(
                "Correcting the flux scale using comparison between the primary and secondary calibrator."
            )
            fluxscale(
                vis=msname,
                caltable=f"{calfile}.G2",
                fluxtable=f"{calfile}.F0",
                reference=pri,
            )
        except:
            logger.warning("Couldn't make calfiles!!! ")

    flagmanager(vis=msname, mode="save", versionname="before_applycal")
    return


def applycal_ms(calfile, msname, sec, tar, pri = "1934_cal_cx"):
    applycal(
        vis=msname,
        gaintable=[
            f"{calfile}.B1",
            f"{calfile}.F0",
        ],
        gainfield=[pri, pri, pri],
        field=f"{pri}",
        parang=True,
        flagbackup=False,
    )
    applycal(
        vis=msname,
        gaintable=[
            f"{calfile}.B1",
            f"{calfile}.F0",
        ],
        gainfield=[pri, pri, sec],
        field=f"{sec},{tar}",
        parang=True,
        flagbackup=False,
    )
    applycal(
        vis=msname,
        gaintable=[
            f"{calfile}.B1",
            f"{calfile}.F0",
        ],
        gainfield=[pri, pri, sec],
        field=tar,
        parang=True,
        flagbackup=False,
    )
    return

def flag_postcal(msname, sec, tar, pri="1934_cal_cx"):

    flagmanager(vis=msname, mode="save", versionname="before_rflag")
    flagdata(
        vis=msname,
        mode="rflag",
        field=pri,
        datacolumn="corrected",
        action="apply",
        display="report",
        correlation="ABS_ALL",
        timedevscale=3.0,
        freqdevscale=3.0,
        winsize=3,
        combinescans=True,
        ntime="9999999min",
        extendflags=False,
        flagbackup=False,
    )
    flagdata(
        vis=msname,
        mode="rflag",
        field=sec,
        datacolumn="corrected",
        action="apply",
        display="report",
        correlation="ABS_ALL",
        timedevscale=3.0,
        freqdevscale=3.0,
        winsize=3,
        combinescans=True,
        ntime="9999999min",
        extendflags=False,
        flagbackup=False,
    )
    flagdata(
        vis=msname,
        mode="extend",
        field=pri + "," + sec,
        action="apply",
        display="report",
        flagbackup=False,
        extendpols=True,
        correlation="",
        growtime=95.0,
        growfreq=95.0,
        growaround=True,
        flagneartime=False,
        flagnearfreq=False,
        combinescans=True,
        ntime="9999999min",
    )

    flagdata(
        vis=msname,
        mode="rflag",
        field=tar,
        datacolumn="corrected",
        action="apply",
        display="report",
        correlation="ABS_ALL",
        timedevscale=3.0,
        freqdevscale=3.0,
        winsize=3,
        combinescans=True,
        ntime="9999999min",
        extendflags=False,
        flagbackup=False,
    )

    return 




def slefcal_ms(calfile, srcms, tar, self_round="pcal0",solint="120s",minblperant=4,combine="",spwmap="",applymode="calonly",calmode="p",field="0",gaintable=[""]):
    logger.debug(
        "+ + + + + + + + + + + + + + + + +\n+  Self Cal +\n+ + + + + + + + + + + + + + + + +"
    )
    gaincal(
        vis=srcms,
        caltable=f"{calfile}_{self_round}.cal",
        field=field,
        combine=combine,
        gaintype="G",
        calmode=calmode,
        solint=solint,
        gaintable=gaintable,
        minblperant=minblperant,
    )
    if gaintable[0] == "":
        gaintable = f"{calfile}_{self_round}.cal"
    else:
        gaintable.append(f"{calfile}_{self_round}.cal")
        print(gaintable)
    applycal(
        vis=srcms,
        gaintable=gaintable,
        parang=True,
        field=field,
        spwmap=spwmap,
        applymode=applymode,
        flagbackup=False,
    )
    flagmanager(vis=srcms, mode="save", versionname=f"post {self_round}")
    
    return




if __name__ == "__main__":
    parser = ArgumentParser(
        description="Script to go through basic processing of continuum ATCA data"
    )
    parser.add_argument(
        '--project',
        type=str,
        default=".",
        help="Path to project directory containing data and where most things will be saved (default= ./)"
    )
    parser.add_argument(
        '--target',
        type=str,
        help="String for field name of target to be processed, no default must specify"
    )
    parser.add_argument(
        "--sec",
        type=str,
        help="string field name for secondary calibrator, no default "
    )
    parser.add_argument(
        '--pri',
        type=str,
        default="1934_cal_cx",
        help="String with field name for primary bandpass calibration. (default= 1934_cal_cx)"
    )
    parser.add_argument(
        '--spw',
        type=str,
        default="0",
        help="String for spw number from main ms file to extract and process, (default='0')"
    )
    parser.add_argument(
        "--nspw",
        tpye=int,
        default=1,
        help="Number of spws to split into, default=1"
    )
    parser.add_argument(
        "--ref",
        type=str,
        default="CA04",
        help="Ref antenna for calibration, default=CA04"
    )



    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        default=False,
        help='Enable extra logging'
    )



    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    

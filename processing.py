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
    try:
        importatca(
            vis=visname, files=files, options="birdie,noac", edge=4
        )
    except:
        logger.warning(f"Unable to read in ATCA data??? Something wrong... ")
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
    setjy(
        vis=msname,
        field=pri,
        scalebychan=True,
        standard="Perley-Butler 2010",
        usescratch=True,
    )
    if os.path.exists(f"{calfile}.G0"):
        logger.debug(f"Found {calfile}.G0, skipping")
    else: 
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
        except:
            logger.warning(f"Couldnt do first cal? ")
            return 
    if os.path.exists(f"{calfile}.B0"):
        logger.debug(f"Found {calfile}.B0, skipping")
    else: 
        logger.debug(f"Performing bandpass calibration on {pri}")
        try:
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
        except: 
            logger.warning("Issue with cal? Couldn't do bandpass?")
            return
    if os.path.exists(f"{calfile}.G1"):
        logger.debug(f"Found {calfile}.G1, skipping")
    else: 
        try:
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
        except: 
            logger.warning("Issue with cal? COuldn't make .G1")
            return 
    if os.path.exists(f"{calfile}.B1"):
        logger.debug(f"Found {calfile}.B1, skipping")
    else: 
        logger.debug(f"Performing gain calibration on {pri}")
        try:
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
        except:
            logger.warning("Issue with cal? Couldn't make .B1")
    if os.path.exists(f"{calfile}.G2"):
        logger.debug(f"Found {calfile}.G2, skipping")
    else: 
        try:
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
        except:
            logger.warning("Issue with cal? Coudln't do secondary gain?")
            return 
    try:
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
        logger.warning("Couldn't do final flux scale!  ")
        return 
    
    logger.debug("Completed making all cal files ")
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
    if args.plots is True: 
        plotms(vis=msname, plotfile=f"{calfile}_{pri}_amppostcal.png",xaxis='frequency', yaxis='amp', ydatacolumn='corrected', field=pri,correlation="XX", showgui=False, overwrite=True)
        plotms(vis=msname, plotfile=f"{calfile}_{pri}_phasepostcal.png",xaxis='frequency', yaxis='phase', ydatacolumn='corrected', field=pri,correlation="XX", showgui=False, overwrite=True)
        plotms(vis=msname, plotfile=f"{calfile}_{tar}_phasepostcal.png",xaxis='frequency', yaxis='phase', ydatacolumn='corrected', field=tar,correlation="XX", showgui=False, overwrite=True)
        plotms(vis=msname, plotfile=f"{calfile}_{tar}_amppostcal.png",xaxis='frequency', yaxis='amp', ydatacolumn='corrected', field=tar,correlation="XX", showgui=False, overwrite=True)
    return

def flag_postcal(msname, sec, tar, calfile, pri="1934_cal_cx"):

    flagmanager(vis=msname, mode="save", versionname="before_rflag")
    flagdata(
        vis=msname,
        mode="tfcrop",
        field=f"{pri},{sec},{tar}",
        datacolumn="corrected",
        # growfreq=80,
        # growaround=True,
        # flagnearfreq=True,
    )
    flagdata(
        vis=msname,
        mode="extend",
        field=f"{pri},{sec},{tar}",
        # growfreq=80,
        # growaround=True,
        # flagnearfreq=True,
    )
    if args.plots is True: 
        plotms(vis=msname, plotfile=f"{calfile}_{pri}_amppostcalflag.png",xaxis='frequency', yaxis='amp', ydatacolumn='corrected', field=pri,correlation="XX", showgui=False, overwrite=True)
        plotms(vis=msname, plotfile=f"{calfile}_{pri}_phasepostcalflag.png",xaxis='frequency', yaxis='phase', ydatacolumn='corrected', field=pri,correlation="XX", showgui=False, overwrite=True)
        plotms(vis=msname, plotfile=f"{calfile}_{tar}_phasepostcalflag.png",xaxis='frequency', yaxis='phase', ydatacolumn='corrected', field=tar,correlation="XX", showgui=False, overwrite=True)
        plotms(vis=msname, plotfile=f"{calfile}_{tar}_amppostcalflag.png",xaxis='frequency', yaxis='amp', ydatacolumn='corrected', field=tar,correlation="XX", showgui=False, overwrite=True)
    return 



def imgmfs_ms(imagems, imagename, field="0", spw="", uvrange="", niter=3000, nterms=1, mode="mfs", antenna="",weighting="briggs", robust=-0.5, imsize=2250, cell="0.5arcsec", interactive=True, scan="", datacolumn="corrected"):
    print("Initiating interactive cleaning")

    tclean(
        vis=imagems,
        imagename=f"{imagename}",
        specmode=mode,
        nterms=nterms,
        niter=niter,
        imsize=imsize,
        cell=cell,
        weighting=weighting,
        robust=robust,
        antenna=antenna,
        spw=spw,
        field=field,
        interactive=interactive,
        savemodel="modelcolumn",
        pbcor=False,
        uvrange=uvrange,
        scan=scan,
        datacolumn=datacolumn,
    )
    tclean(
        vis=imagems,
        imagename=f"{imagename}",
        specmode=mode,
        nterms=nterms,
        niter=0,
        spw=spw,
        field=field,
        imsize=imsize,
        cell=cell,
        weighting=weighting,
        robust=robust,
        antenna=antenna,
        interactive=interactive,
        savemodel="modelcolumn",
        pbcor=False,
        calcres=False,
        calcpsf=False,
        uvrange=uvrange,
        datacolumn=datacolumn,
    )
    return



def slefcal_ms(calfile, srcms, tar, selfround="0",solint="60s",minblperant=4,combine="",spwmap="",applymode="",calmode="p",gaintable=[""],solnorm=False):
    logger.debug("Running selfcal")
    gaincal(
        vis=srcms,
        caltable=f"{calfile}_self{selfround}.cal",
        field=tar,
        combine=combine,
        gaintype="G",
        calmode=calmode,
        solint=solint,
        gaintable=gaintable,
        minblperant=minblperant,
        solnorm=solnorm
    )
    if gaintable[0] == "":
        gaintable = [f"{calfile}_self{selfround}.cal"]
    else:
        gaintable.append(f"{calfile}_self{selfround}.cal")
        logger.debug(f"{gaintable}")
    applycal(
        vis=srcms,
        gaintable=gaintable,
        parang=True,
        field=tar,
        spwmap=spwmap,
        applymode=applymode,
        flagbackup=False,
    )
    flagmanager(vis=srcms, mode="save", versionname=f"post {selfround}")
    
    return gaintable 


# TODO: Currnetly only calibrates assuming 1 spw but gives options for it to split to multiple...

if __name__ == "__main__":
    parser = ArgumentParser(
        description="Script to go through basic processing of continuum ATCA data"
    )
    parser.add_argument(
        '--cont',
        default=True,
        help="Option to search for products at each step first and then use them instead of starting from scratch. Default=True"
    )
    parser.add_argument(
        '--applycal',
        default=True,
        help="Apply the cal solutions? This is regardless of continue, is set to false, will override continue and skip the apply step anyway "
    )
    parser.add_argument(
        "--plots",
        default=True,
        help="Plot things like cal phase and amplitude of cals. Default = True "
    )
    parser.add_argument(
        "--flag6",
        default=True, 
        help="Flag antenna 6 altogether for compact configs of ATCA, default=True"
    )
    parser.add_argument(
        '--dir',
        type=str,
        default=".",
        help="Path to project directory containing data and where most things will be saved (default= ./)"
    )
    parser.add_argument(
        '--project',
        type=str,
        default="C3487",
        help="The project code from ATCA to use for naming and for searching for all data files to make ms. default=C3487"
    )
    parser.add_argument(
        "--visname",
        type=str,
        default="c3487_day0_cx.ms",
        help="The visname to do initial flagging and from which the target ms is split default = c3487_day0_cx.ms"
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
        '--band',
        type=str,
        default="cx",
        help="String for spw number from main ms file to extract and process, (default='cx')"
    )
    parser.add_argument(
        "--nspw",
        type=int,
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
        
    

    tar=args.target
    sec=args.sec
    dir=args.dir
    pri=args.pri
    band=args.band

# TODO: Check this actually works? 
    # if args.spw == "0,1":
    #     band = "cx"
    # elif args.spw =="2,3":
    #     band = "l"
    # elif args.spw == "2":
    #     band = "L"

    selfround = "0"
    msname = f"{dir}/{tar}_{band}.ms"
    tarms = f"{dir}/{tar}_{band}_cal.ms"
    visname = f"{dir}/{args.visname}"
    calfile = f"{dir}/{sec}_cal_{band}"
    imagename = f"{dir}/{tar}_{band}_spw"

    files=[]
    for f in os.listdir():
        if f.endswith(f".{args.project}"):
            files.append(f)

    # Read in the ATCA data and make initial visname 
    if args.cont is True:
        logger.warning("Continue is on, checking if files exist before running make_ms")
        if os.path.exists(f"{visname}"):
            logger.debug(f"Found vis, not reading in the ATCA data")
        else: 
            make_ms(files, visname)    
    else: 
        logger.warning(f"Found visname and continue is on so deleting ms")
        os.system(f"rm -r {visname}")
        os.system(f"rm -r {visname}.flagversions") 
        make_ms(files, visname)
    
    if args.flag6 is True:
        flagdata(vis=visname, mode="manual", antenna="5")

    # Flagging data 
    if args.cont is True: 
        logger.debug(f"Continue on: checking for flag versions before flagging ")
        if os.path.exists(f"{visname}.flagversions"):
            logger.debug(f"Found flag version for visname, not running the flagging")
        else: 
            flag_ms(visname)
    else: 
        flag_ms(visname)


    # # Splitting main MS to have just target ms 
    # if args.cont is True:
    #     logger.debug(f"Contunue on, checking if there is the target ms before splitting")
    #     if os.path.exists(f"{msname}"):
    #         logger.debug("Found msname, skipping split")
    #     else: 
    #         split_ms(visname, msname, field=f"{pri},{sec},{tar}", spw=args.spw, n_spw=args.nspw, datacolumn="data",listfile=f"listobs_{msname}.dat")
    # else: 
    #     if os.path.exists(f"{msname}"):
    #         os.system(f"rm -r {msname}")
    #         os.system(f"rm -r {msname}.flagversions")
    #         split_ms(visname, msname, field=f"{pri},{sec},{tar}", spw=args.spw, n_spw=args.nspw, datacolumn="data",listfile=f"listobs_{msname}.dat")
    #         logger.debug(f"remade split version")
    


    # Generating calibration solutions 
    if args.cont is True: 
        logger.debug(f"Continue on: Looking for calibration files before trying calibation")
        if os.path.exists(f"{calfile}.G2"):
            logger.debug("Found last calfile needed so skipping generating cal solutions")
        else: 
            calibrate_ms(visname, sec, calfile, pri=pri,ref = args.ref)
    else: 
        os.system(f"rm -r {calfile}*")
        calibrate_ms(visname, sec, calfile, pri=pri,ref = args.ref)



    # Applying cal solutions 
    if args.applycal is True: 
        logger.debug(f"Apply on: Applying solutions now ")
        applycal_ms(calfile, visname, sec, tar, pri=pri)
    else: 
        logger.warning(f"NOT APPLYING CAL! ")

    flag_postcal(visname, sec, tar, calfile, pri=pri)

    # Running first initial imaging 
    if args.cont is True: 
        logger.debug(f"Continue on: looking for images before starting again")
        if os.path.exists(f"{imagename}0.image"):
            logger.debug(f"Found image, not reimaging")
        else:
            temp_imname = f"{imagename}0" 
            imgmfs_ms(visname, temp_imname, field=tar, spw="0")
            temp_imname = f"{imagename}1" 
            imgmfs_ms(visname, temp_imname, field=tar, spw="1")
    else: 
        os.system(f"rm -r {imagename}*")
        temp_imname = f"{imagename}0"
        imgmfs_ms(visname, temp_imname, field=tar, spw="0")
        temp_imname = f"{imagename}1"
        imgmfs_ms(visname, temp_imname, field=tar, spw="1")


# # TODO: change self so that it just keeps going till it reaches the number specified 
#     # Splitting target MS to have just calibrated target ms 
#     if args.cont is True:
#         logger.debug(f"Contunue on, checking if there is the calibrated target ms before splitting")
#         if os.path.exists(f"{msname}"):
#             logger.debug("Found tar ms, skipping split")
#         else: 
#             split_ms(visname, msname, field=f"{tar}", spw=args.spw, n_spw=args.nspw, datacolumn="corrected")
#     else: 
#         if os.path.exists(f"{msname}"):
#             os.system(f"rm -r {msname}")
#             os.system(f"rm -r {msname}.flagversions")
#             split_ms(visname, msname, field=f"{tar}", spw=args.spw, n_spw=args.nspw, datacolumn="corrected")
#             logger.debug(f"remade split version")

#     # Running first round of selfcal
#     selfround="1"
#     temp_imname = f"{imagename}{selfround}"
#     if args.cont is True: 
#         logger.debug(f"Continue on: checking for selfcal round{selfround}")
#         if os.path.exists(f"{calfile}_self{selfround}.cal"):
#             logger.debug(f"Found first round cal solution, NOT APPLYING!")
#             gaintable = [f"{calfile}_self{selfround}.cal"]
#         else: 
#             gaintable = slefcal_ms(calfile, msname, tar, selfround=selfround)
#     else: 
#         gaintable = slefcal_ms(calfile, msname, tar, selfround=selfround)
    
#     # Running first self imaging 
#     if args.cont is True: 
#         logger.debug(f"Continue on: looking for images before starting again")
#         if os.path.exists(f"{imagename}{selfround}.image"):
#             logger.debug(f"Found image, not reimaging")
#         else:
#             temp_imname = f"{imagename}{selfround}" 
#             imgmfs_ms(msname, temp_imname, field=tar)
#     else: 
#         os.system(f"rm -r {imagename}{selfround}*")
#         temp_imname = f"{imagename}{selfround}"
#         imgmfs_ms(msname, temp_imname, field=tar)    

#     # Running second round of selfcal
#     selfround="2"
#     temp_imname = f"{imagename}{selfround}"
#     if args.cont is True: 
#         logger.debug(f"Continue on: checking for selfcal round{selfround}")
#         if os.path.exists(f"{calfile}_self{selfround}.cal"):
#             logger.debug(f"Found first round cal solution, NOT APPLYING!")
#             gaintable = [f"{calfile}_self1.cal", f"{calfile}_self{selfround}.cal"]
#         else: 
#             gaintable = slefcal_ms(calfile, msname, tar, selfround=selfround, solint="30s", gaintable=gaintable)
#     else: 
#         gaintable = slefcal_ms(calfile, msname, tar, selfround=selfround, solint="30s", gaintable=gaintable)
    
#     # Running second self imaging 
#     if args.cont is True: 
#         logger.debug(f"Continue on: looking for images before starting again")
#         if os.path.exists(f"{imagename}{selfround}.image"):
#             logger.debug(f"Found image, not reimaging")
#         else:
#             temp_imname = f"{imagename}{selfround}" 
#             imgmfs_ms(msname, temp_imname, field=tar)
#     else: 
#         os.system(f"rm -r {imagename}{selfround}*")
#         temp_imname = f"{imagename}{selfround}"
#         imgmfs_ms(msname, temp_imname, field=tar)    


#     # Running third round of selfcal
#     selfround="3"
#     temp_imname = f"{imagename}{selfround}"
#     if args.cont is True: 
#         logger.debug(f"Continue on: checking for selfcal round{selfround}")
#         if os.path.exists(f"{calfile}_self{selfround}.cal"):
#             logger.debug(f"Found first round cal solution, NOT APPLYING!")
#             gaintable = [f"{calfile}_self1.cal", f"{calfile}_self2.cal", f"{calfile}_self{selfround}.cal"]
#         else: 
#             gaintable = slefcal_ms(calfile, msname, tar, selfround=selfround, solint="int", gaintable=gaintable)
#     else: 
#         gaintable = slefcal_ms(calfile, msname, tar, selfround=selfround, solint="int", gaintable=gaintable)
    
#     # Running third self imaging 
#     if args.cont is True: 
#         logger.debug(f"Continue on: looking for images before starting again")
#         if os.path.exists(f"{imagename}{selfround}.image"):
#             logger.debug(f"Found image, not reimaging")
#         else:
#             temp_imname = f"{imagename}{selfround}" 
#             imgmfs_ms(msname, temp_imname, field=tar)
#     else: 
#         os.system(f"rm -r {imagename}{selfround}*")
#         temp_imname = f"{imagename}{selfround}"
#         imgmfs_ms(msname, temp_imname, field=tar)    



#     # Running final amp and phase round of selfcal
#     selfround="4"
#     temp_imname = f"{imagename}{selfround}"
#     if args.cont is True: 
#         logger.debug(f"Continue on: checking for selfcal round{selfround}")
#         if os.path.exists(f"{calfile}_self{selfround}.cal"):
#             logger.debug(f"Found first round cal solution, NOT APPLYING!")
#             gaintable = [f"{calfile}_self1.cal", f"{calfile}_self2.cal", f"{calfile}_self{selfround}.cal"]
#         else: 
#             gaintable = slefcal_ms(calfile, msname, tar, selfround=selfround, solint="inf", gaintable=gaintable,solnorm=True, calmode="ap")
#     else: 
#         gaintable = slefcal_ms(calfile, msname, tar, selfround=selfround, solint="inf", gaintable=gaintable,solnorm=True, calmode="ap")
    
#     # Running third self imaging 
#     if args.cont is True: 
#         logger.debug(f"Continue on: looking for images before starting again")
#         if os.path.exists(f"{imagename}{selfround}.image"):
#             logger.debug(f"Found image, not reimaging")
#         else:
#             temp_imname = f"{imagename}{selfround}" 
#             imgmfs_ms(msname, temp_imname, field=tar)
#     else: 
#         os.system(f"rm -r {imagename}{selfround}*")
#         temp_imname = f"{imagename}{selfround}"
#         imgmfs_ms(msname, temp_imname, field=tar)  
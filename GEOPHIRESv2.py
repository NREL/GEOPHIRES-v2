# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 10:34:04 2017

@author: kbeckers
"""

# GEOPHIRES v2.0
# build date: December 9th 2018
# https://github.com/kfbeckers/GEOPHIRES

# import functions
import math
import datetime
import numpy as np
import time
from mpmath import *
import os
import sys


# user-defined functions
def densitywater(Twater):
    T = Twater + 273.15
    rhowater = (.7983223 + (
                1.50896E-3 - 2.9104E-6 * T) * T) * 1E3  # water density correlation as used in Geophires v1.2 [kg/m3]
    return rhowater;


def viscositywater(Twater):
    muwater = 2.414E-5 * np.power(10, 247.8 / (
                Twater + 273.15 - 140))  # accurate to within 2.5% from 0 to 370 degrees C [Ns/m2]
    # xp = np.linspace(5,150,30)
    # fp = np.array([1519.3, 1307.0, 1138.3, 1002.0, 890.2, 797.3, 719.1, 652.7, 596.1, 547.1, 504.4, 467.0, 433.9, 404.6, 378.5, 355.1, 334.1, 315.0, 297.8, 282.1, 267.8, 254.4, 242.3, 231.3, 221.3, 212.0, 203.4, 195.5, 188.2, 181.4])
    # muwater = np.interp(Twater,xp,fp)
    return muwater;


def heatcapacitywater(Twater):
    Twater = (Twater + 273.15) / 1000
    A = -203.6060
    B = 1523.290
    C = -3196.413
    D = 2474.455
    E = 3.855326
    cpwater = (A + B * Twater + C * Twater ** 2 + D * Twater ** 3 + E / (
                Twater ** 2)) / 18.02 * 1000  # water specific heat capacity in J/kg-K
    return cpwater;


def vaporpressurewater(Twater):
    if Twater < 100:
        A = 8.07131
        B = 1730.63
        C = 233.426
    else:
        A = 8.14019
        B = 1810.94
        C = 244.485
    vaporpressurewater = 133.322 * (
                10 ** (A - B / (C + Twater))) / 1000  # water vapor pressure in kPa using Antione Equation
    return vaporpressurewater;


def run():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # specify path of input file
    fname = os.path.join('Examples', 'example4.txt')

    tic = time.time()

    # read input data (except temperature profile from reservoir)
    try:
        with open(fname, encoding='UTF-8') as f:
            content = f.readlines()
    except:
        print("Error: GEOPHIRES could not read input file (" + fname + ") and will abort simulation.")
        sys.exit()

    # potential other parameters to be read in
    #   external reservoir output filename
    #   tough2file name
    #   distance for surface pipe

    # enduseoption
    # enduseoption = 1: electricity
    # enduseoption = 2: direct-use heat
    # enduseoption = 3: cogen topping cycle
    # enduseoption = 4: cogen bottoming cycle
    # enduseoption = 5: cogen split of mass flow rate
    try:
        enduseoption = int(
            content[[i for i, s in enumerate(content) if 'End-Use Option,' in s][0]].split(',')[1].strip('\n'))
        if not (enduseoption in [1, 2, 31, 32, 41, 42, 51, 52]):
            enduseoption = 1
            print(
                "Warning: Provided end-use option is not 1, 2, 31, 32, 41, 42, 51, or 52. GEOPHIRES will assume default end-use option (1: electricity)")
    except:
        enduseoption = 1
        print(
            "Warning: No valid end-use option provided. GEOPHIRES will assume default end-use option (1: electricity)")

    # pptype: power plant type
    # pptype = 1: Subcritical ORC
    # pptype = 2: Supercritical ORC
    # pptype = 3: Single-Flash
    # pptype = 4: Double-Flash
    if enduseoption in [1, 31, 32, 41, 42, 51, 52]:
        try:
            pptype = int(
                content[[i for i, s in enumerate(content) if 'Power Plant Type,' in s][0]].split(',')[1].strip('\n'))
            if not (pptype in [1, 2, 3, 4]):
                pptype = 1
                print(
                    "Warning: Provided power plant type is not 1, 2, 3 or 4. GEOPHIRES will assume default power plant type (1: subcritical ORC)")
        except:
            pptype = 1
            print(
                "Warning: No valid power plant type provided. GEOPHIRES will assume default power plant type (1: subcritical ORC)")

    # pumpeff: pump efficiency (-)
    try:
        pumpeff = float(
            content[[i for i, s in enumerate(content) if 'Circulation Pump Efficiency,' in s][0]].split(',')[1].strip(
                '\n'))
        if pumpeff < 0.1 or pumpeff > 1:
            pumpeff = 0.75
            print(
                "Warning: Provided circulation pump efficiency outside of range 0.1-1. GEOPHIRES will assume default circulation pump efficiency (0.75)")
    except:
        pumpeff = 0.75
        print(
            "Warning: No valid circulation pump efficiency provided. GEOPHIRES will assume default circulation pump efficiency (0.75)")

    # utilfactor: utilization factor (-)
    try:
        utilfactor = float(
            content[[i for i, s in enumerate(content) if 'Utilization Factor,' in s][0]].split(',')[1].strip('\n'))
        if utilfactor < 0.1 or utilfactor > 1:
            utilfactor = 0.9
            print(
                "Warning: Provided utilization factor outside of range 0.1-1. GEOPHIRES will assume default utilization factor (0.9)")
    except:
        utilfactor = 0.9
        print("Warning: No valid utilization factor provided. GEOPHIRES will assume default utilization factor (0.9)")

    # enduseefficiencyfactor: end-use efficiency for direct-use heat component [-]
    if enduseoption in [2, 31, 32, 41, 42, 51, 52]:
        try:
            enduseefficiencyfactor = float(
                content[[i for i, s in enumerate(content) if 'End-Use Efficiency Factor,' in s][0]].split(',')[1].strip(
                    '\n'))
            if enduseefficiencyfactor < 0.1 or enduseefficiencyfactor > 1:
                enduseefficiencyfactor = 0.9
                print(
                    "Warning: Provided end-use efficiency factor outside of range 0.1-1. GEOPHIRES will assume default end-use efficiency factor (0.9)")
        except:
            enduseefficiencyfactor = 0.9
            print(
                "Warning: No valid end-use efficiency factor provided. GEOPHIRES will assume default end-use efficiency factor (0.9)")

    # chpfraction: fraction of flow rate going to direct-use heat application  (only used in CHP parallel cycle)
    if enduseoption in [51, 52]:
        try:
            chpfraction = float(
                content[[i for i, s in enumerate(content) if 'CHP Fraction,' in s][0]].split(',')[1].strip('\n'))
            if chpfraction < 0.0001 or chpfraction > 0.9999:
                chpfraction = 0.5
                print(
                    "Warning: Provided CHP fraction outside of range 0.0001-0.9999. GEOPHIRES will assume default CHP fraction (0.5)")
        except:
            chpfraction = 0.5
            print("Warning: No valid CHP fraction provided. GEOPHIRES will assume default CHP fraction (0.5)")

    # Tinj: injection temperature (C)
    try:
        Tinj = float(
            content[[i for i, s in enumerate(content) if 'Injection Temperature,' in s][0]].split(',')[1].strip('\n'))
        if Tinj < 0 or Tinj > 200:
            Tinj = 70
            print(
                "Warning: Provided injection temperature outside range of 0-200. GEOPHIRES will assume default injection temperature (70 deg.C)")
    except:
        Tinj = 70
        print(
            "Warning: No valid injection temperature provided. GEOPHIRES will assume default injection temperature (70 deg.C)")

    # Tmax: Maximum allowable Reservoir Temperature (C)
    try:
        Tmax = float(
            content[[i for i, s in enumerate(content) if 'Maximum Temperature,' in s][0]].split(',')[1].strip('\n'))
        if Tmax < 50 or Tmax > 1000:
            Tmax = 400
            print(
                "Warning: Provided maximum temperature outside of range 50-1000. GEOPHIRES will assume default maximum temperature (400 deg.C)")
    except:
        Tmax = 400
        print(
            "Warning: No valid maximum temperature provided. GEOPHIRES will assume default maximum temperature (400 deg.C)")

    # Tchpbottom: power plant entering temperature in the CHP Bottom cycle (in deg.C)
    if enduseoption in [41, 42]:
        try:
            Tchpbottom = float(
                content[[i for i, s in enumerate(content) if 'CHP Bottoming Entering Temperature,' in s][0]].split(',')[
                    1].strip('\n'))
            if Tchpbottom < Tinj or Tchpbottom > Tmax:
                Tchpbottom = 150
                print(
                    "Warning: Provided CHP bottoming entering temperature outside of range Tinj-Tmax. GEOPHIRES will assume default CHP bottom temperature (150 deg.C)")
        except:
            Tchpbottom = 150
            print(
                "Warning: Provided CHP bottoming entering temperature outside of range Tinj-Tmax. GEOPHIRES will assume default CHP bottom temperature (150 deg.C)")

    # Tsurf: surface temperature used for calculating bottomhole temperature (in deg.C)
    try:
        Tsurf = float(
            content[[i for i, s in enumerate(content) if 'Surface Temperature,' in s][0]].split(',')[1].strip('\n'))
        if Tsurf < -50 or Tsurf > 50:
            Tsurf = 15
            print(
                "Warning: Provided surface temperature outside of range -50 to 50. GEOPHIRES will assume default surface temperature (15 deg.C)")
    except:
        Tsurf = 15
        print(
            "Warning: No valid surface temperature provided. GEOPHIRES will assume default surface temperature (15 deg.C)")

    # Tenv: ambient temperature (in deg.C)
    if enduseoption in [1, 31, 32, 41, 42, 51, 52]:
        try:
            Tenv = float(
                content[[i for i, s in enumerate(content) if 'Ambient Temperature,' in s][0]].split(',')[1].strip('\n'))
            if Tenv < -50 or Tenv > 50:
                Tenv = 15
                print(
                    "Warning: Provided ambient temperature outside of range -50 to 50. GEOPHIRES will assume default ambient temperature (15 deg.C)")
        except:
            Tenv = 15
            print(
                "Warning: No valid ambient temperature provided. GEOPHIRES will assume default ambient temperature (15 deg.C)")

    # resoption: Reservoir Option
    #   resoption = 1  Multiple parallel fractures model (LANL)
    #   resoption = 2  Volumetric block model (1D linear heat sweep model (Stanford))
    #   resoption = 3  Drawdown parameter model (Tester)
    #   resoption = 4  Thermal drawdown percentage model (GETEM)
    #   resoption = 5  Generic user-provided temperature profile
    #   resoption = 6  TOUGH2 is called
    try:
        resoption = int(
            content[[i for i, s in enumerate(content) if 'Reservoir Model,' in s][0]].split(',')[1].strip('\n'))
    except:
        resoption = 4
        print(
            "Warning: Parameter 'Reservoir Model' not found. GEOPHIRES will run default reservoir model (Thermal Drawdown Percentage Model)")
    if not (resoption in [1, 2, 3, 4, 5, 6]):
        print(
            "Warning: Selected Reservoir Model not valid. GEOPHIRES will run default reservoir model (Thermal Drawdown Percentage Model)")
        resoption = 4

    # drawdp: Drawdown parameter
    #   used in both resopt 3 and 4
    #   if resoption = 3: drawdp is in units of kg/s/m2
    #   if resoption = 4: drawdp is in units of 1/year
    if resoption == 3 or resoption == 4:
        try:
            drawdp = float(
                content[[i for i, s in enumerate(content) if 'Drawdown Parameter,' in s][0]].split(',')[1].strip('\n'))
            if drawdp < 0 or drawdp > 0.2:
                if resoption == 3:
                    drawdp = 0.0001
                    print(
                        "Warning: Provided drawdown parameter outside of range 0-0.2. GEOPHIRES will assume default drawdown parameter (0.0001 kg/s/m2) for reservoir model 3")
                elif resoption == 4:
                    drawdp = 0.005
                    print(
                        "Warning: Provided drawdown parameter outside of range 0-0.2. GEOPHIRES will assume default drawdown parameter (0.5 %/year) for reservoir model 4")
        except:
            if resoption == 3:
                drawdp = 0.0001
                print(
                    "Warning: No valid drawdown parameter found. GEOPHIRES will assume default drawdown parameter (0.0001 kg/s/m2) for reservoir model 3")
            elif resoption == 4:
                drawdp = 0.005
                print(
                    "Warning: No valid drawdown parameter found. GEOPHIRES will assume default drawdown parameter (0.5 %/year) for reservoir model 4")

    # read file name of reservoir output in case reservoir model 5 is selected
    if resoption == 5:
        try:
            filenamereservoiroutput = \
            content[[i for i, s in enumerate(content) if 'Reservoir Output File Name,' in s][0]].split(',')[1].strip(
                '\n')
        except:
            filenamereservoiroutput = 'ReservoirOutput.txt'
            print(
                "Warning: No valid file name reservoir output found. GEOPHIRES will assume default reservoir output file name (ReservoirOutput.txt)")

    # read TOUGH2 file name if reservoir model 6 is selected. If written 'Doublet', GEOPHIRES will run built-in TOUGH2 doublet model.
    if resoption == 6:
        try:
            tough2modelfilename = \
            content[[i for i, s in enumerate(content) if 'TOUGH2 Model/File Name,' in s][0]].split(',')[1].strip('\n')
        except:
            tough2modelfilename = 'Doublet'
            print(
                "Warning: No valid TOUGH2 model or file name provided. GEOPHIRES will assume default built-in TOUGH2 model (Doublet).")
        if tough2modelfilename == 'Doublet':
            usebuiltintough2model = 1
        else:
            usebuiltintough2model = 0

    # depth: Measured depth of the well (provided in km by user and converted here to m).
    try:
        depth = float(
            content[[i for i, s in enumerate(content) if 'Reservoir Depth,' in s][0]].split(',')[1].strip('\n'))
        if depth < 0.1 or depth > 15:
            depth = 3.
            print(
                "Warning: Provided reservoir depth outside of range 0.1-15. GEOPHIRES will assume default reservoir depth (3 km)")
    except:
        depth = 3.
        print("Warning: No reservoir depth found. GEOPHIRES will assume default reservoir depth (3 km)")
    depth = depth * 1000

    # numseg: number of segments
    try:
        numseg = int(
            content[[i for i, s in enumerate(content) if 'Number of Segments,' in s][0]].split(',')[1].strip('\n'))
        if not (numseg in [1, 2, 3, 4]):
            print(
                "Warning: Provided number of segments outside of range 1-4. GEOPHIRES will assume default number of segments (1)")
    except:
        numseg = 1
        print("Warning: No valid number of segments provided. GEOPHIRES will assume default number of segments (1)")

    # gradient(i): geothermal gradient of layer i (provided in C/km and converted to C/m)
    # layerthickness(i): thickness of layer i (provided in km and converted to m)
    gradient = [0, 0, 0, 0];
    layerthickness = [0, 0, 0, 0];
    try:
        gradient[0] = float(
            content[[i for i, s in enumerate(content) if 'Gradient 1,' in s][0]].split(',')[1].strip('\n')) / 1000
        if gradient[0] < 0 or gradient[0] > 0.5:
            print(
                "Warning: Provided geothermal gradient for layer 1 outside of range 0-500. GEOPHIRES will assume default geothermal gradient (50 deg.C/km)")
            gradient[0] = 50. / 1000
    except:
        gradient[0] = 50. / 1000
        print(
            "Warning: No valid geothermal gradient for layer 1 provided. GEOPHIRES will assume default geothermal gradient (50 deg.C/km)")

    if numseg > 1:
        try:
            gradient[1] = float(
                content[[i for i, s in enumerate(content) if 'Gradient 2,' in s][0]].split(',')[1].strip('\n')) / 1000
            if gradient[1] < 0 or gradient[1] > 0.5:
                print(
                    "Warning: Provided geothermal gradient for layer 2 outside of range 0-500. GEOPHIRES will assume default geothermal gradient (50 deg.C/km)")
                gradient[1] = 50. / 1000
        except:
            gradient[1] = 50. / 1000
            print(
                "Warning: No valid geothermal gradient for layer 2 provided. GEOPHIRES will assume default geothermal gradient (50 deg.C/km)")
        try:
            layerthickness[0] = float(
                content[[i for i, s in enumerate(content) if 'Thickness 1,' in s][0]].split(',')[1].strip('\n')) * 1000
            if layerthickness[0] < 10 or layerthickness[0] > 100000:
                print(
                    "Warning: Provided thickness for layer 1 outside of range 0.01-100. GEOPHIRES will assume default layer thickness (2 km)")
                layerthickness[0] = 2. * 1000
        except:
            layerthickness[0] = 2. * 1000
            print(
                "Warning: No valid thickness for layer 1 provided. GEOPHIRES will assume default layer thickness (2 km)")

    if numseg > 2:
        try:
            gradient[2] = float(
                content[[i for i, s in enumerate(content) if 'Gradient 3,' in s][0]].split(',')[1].strip('\n')) / 1000
            if gradient[2] < 0 or gradient[2] > 0.5:
                print(
                    "Warning: Provided geothermal gradient for layer 3 outside of range 0-500. GEOPHIRES will assume default geothermal gradient (50 deg.C/km)")
                gradient[2] = 50. / 1000
        except:
            gradient[2] = 50. / 1000
            print(
                "Warning: No valid geothermal gradient for layer 3 provided. GEOPHIRES will assume default geothermal gradient (50 deg.C/km)")
        try:
            layerthickness[1] = float(
                content[[i for i, s in enumerate(content) if 'Thickness 2,' in s][0]].split(',')[1].strip('\n')) * 1000
            if layerthickness[1] < 10 or layerthickness[1] > 100000:
                print(
                    "Warning: Provided thickness for layer 2 outside of range 0.01-100. GEOPHIRES will assume default layer thickness (2 km)")
                layerthickness[1] = 2. * 1000
        except:
            layerthickness[1] = 2. * 1000
            print(
                "Warning: No valid thickness for layer 2 provided. GEOPHIRES will assume default layer thickness (2 km)")

    if numseg > 3:
        try:
            gradient[3] = float(
                content[[i for i, s in enumerate(content) if 'Gradient 4,' in s][0]].split(',')[1].strip('\n')) / 1000
            if gradient[3] < 0 or gradient[3] > 0.5:
                print(
                    "Warning: Provided geothermal gradient for layer 4 outside of range 0-500. GEOPHIRES will assume default geothermal gradient (50 deg.C/km)")
                gradient[3] = 50. / 1000
        except:
            gradient[3] = 50. / 1000
            print(
                "Warning: No valid geothermal gradient for layer 4 provided. GEOPHIRES will assume default geothermal gradient (50 deg.C/km)")
        try:
            layerthickness[2] = float(
                content[[i for i, s in enumerate(content) if 'Thickness 3,' in s][0]].split(',')[1].strip('\n')) * 1000
            if layerthickness[2] < 10 or layerthickness[2] > 100000:
                print(
                    "Warning: Provided thickness for layer 3 outside of range 0.01-100. GEOPHIRES will assume default layer thickness (2 km)")
                layerthickness[2] = 2. * 1000
        except:
            layerthickness[2] = 2. * 1000
            print(
                "Warning: No valid thickness for layer 3 provided. GEOPHIRES will assume default layer thickness (2 km)")

    # set thickness of bottom segment to large number to override lower, unused segments
    layerthickness[numseg - 1] = 100000
    # convert 0 C/m gradients to very small number, avoids divide by zero errors later
    gradient = [1e-6 if x == 0 else x for x in gradient]

    # nprod: number of production wells
    # ninj: number of injection wells
    try:
        nprod = float(
            content[[i for i, s in enumerate(content) if 'Number of Production Wells,' in s][0]].split(',')[1].strip(
                '\n'))
        if not (nprod in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]):
            print(
                "Warning: Provided number of production wells is outside range 1-20. GEOPHIRES will assume default number of production wells (2)")
            nprod = 2
    except:
        print(
            "Warning: No valid number of production wells provided. GEOPHIRES will assume default number of production wells (2)")
        nprod = 2
    try:
        ninj = float(
            content[[i for i, s in enumerate(content) if 'Number of Injection Wells,' in s][0]].split(',')[1].strip(
                '\n'))
        if not (ninj in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]):
            print(
                "Warning: Provided number of injection wells is outside range 1-20. GEOPHIRES will assume default number of injection wells (2)")
            ninj = 2
    except:
        print(
            "Warning: No valid number of injection wells provided. GEOPHIRES will assume default number of injection wells (2)")
        ninj = 2

    # prodwelldiam: production well diameter (input as inch and converted to m)
    # injwelldiam: injection well diameter (input as inch and converted to m)
    try:
        prodwelldiam = float(
            content[[i for i, s in enumerate(content) if 'Production Well Diameter,' in s][0]].split(',')[1].strip(
                '\n')) * 0.0254
        if prodwelldiam / 0.0254 < 1 or prodwelldiam / 0.0254 > 30:
            prodwelldiam = 8 * 0.0254
            print(
                "Warning: Provided production well diameter is outside range 1-30. GEOPHIRES will assume default production well diameter (8 inch)")
    except:
        prodwelldiam = 8 * 0.0254
        print(
            "Warning: No valid production well diameter provided. GEOPHIRES will assume default production well diameter (8 inch)")

    try:
        injwelldiam = float(
            content[[i for i, s in enumerate(content) if 'Injection Well Diameter,' in s][0]].split(',')[1].strip(
                '\n')) * 0.0254
        if injwelldiam / 0.0254 < 1 or injwelldiam / 0.0254 > 30:
            injwelldiam = 8 * 0.0254
            print(
                "Warning: Provided injection well diameter is outside range 1-30. GEOPHIRES will assume default injection well diameter (8 inch)")
    except:
        injwelldiam = 8 * 0.0254
        print(
            "Warning: No valid injection well diameter provided. GEOPHIRES will assume default injection well diameter (8 inch)")

    # rameyoptionprod
    # rameyoptionprod = 0: use tempdrop to calculate production well temperature drop
    # rameyoptionprod = 1: use Ramey model to calculate production well temperature drop
    try:
        rameyoptionprod = int(
            content[[i for i, s in enumerate(content) if 'Ramey Production Wellbore Model,' in s][0]].split(',')[
                1].strip('\n'))
        if not (rameyoptionprod in [0, 1]):
            rameyoptionprod = 1
            print(
                "Warning: Selected Ramey Production Wellbore Model parameter not valid. GEOPHIRES will assume default production wellbore model (Ramey model active)")
    except:
        rameyoptionprod = 1
        print(
            "Warning: No valid Ramey Production Wellbore Model parameter provided. GEOPHIRES will assume default productino wellbore model (Ramey model active)")

    # tempdropprod: temperature drop in production well in deg. C (if Ramey model is not used)
    if rameyoptionprod == 0:
        try:
            tempdropprod = float(
                content[[i for i, s in enumerate(content) if 'Production Wellbore Temperature Drop,' in s][0]].split(
                    ',')[1].strip('\n'))
            if tempdropprod < -5 or tempdropprod > 50:
                print(
                    "Warning: Provided production wellbore temperature drop outside of range -5 to 50. GEOPHIRES will assume default production wellbore temperature drop (5deg.C)")
                tempdropprod = 5
        except:
            tempdropprod = 5
            print(
                "Warning: No valid production wellbore temperature drop provided. GEOPHIRES will assume default production wellbore temperature drop (5deg.C)")

    try:
        tempgaininj = float(
            content[[i for i, s in enumerate(content) if 'Injection Wellbore Temperature Gain,' in s][0]].split(',')[
                1].strip('\n'))
        if tempgaininj < -5 or tempgaininj > 50:
            print(
                "Warning: Provided injection wellbore temperature gain outside of range -5 to 50. GEOPHIRES will assume default injection wellbore temperature gain (0deg.C)")
            tempgaininj = 0
    except:
        tempgaininj = 0
        print(
            "Warning: No valid injection wellbore temperature gain provided. GEOPHIRES will assume default injection wellbore temperature gain (0deg.C)")

    # prodwellflowrate: flow rate per production well (kg/s)
    try:
        prodwellflowrate = float(
            content[[i for i, s in enumerate(content) if 'Production Flow Rate per Well,' in s][0]].split(',')[1].strip(
                '\n'))
        if prodwellflowrate < 1 or prodwellflowrate > 500:
            prodwellflowrate = 50
            print(
                "Warning: Provided production wellbore flow rate is outside of range 1-500. GEOPHIRES will assume default flow rate per production well (50 kg/s)")
    except:
        prodwellflowrate = 50
        print(
            "Warning: No valid production wellbore flow rate is provided. GEOPHIRES will assume default flow rate per production well (50 kg/s)")

    # resvoloption: Rock mass volume option
    #   resvoloption = 1  Specify fracnumb, fracsep
    #   resvoloption = 2  specify resvol, fracsep
    #   resvoloption = 3  Specify resvol, fracnumb
    #   resvoloption = 4: Specify resvol only (sufficient for reservoir models 3, 4, 5 and 6)
    try:
        resvoloption = int(
            content[[i for i, s in enumerate(content) if 'Reservoir Volume Option,' in s][0]].split(',')[1].strip('\n'))
        if not resvoloption in [1, 2, 3, 4]:
            if resoption in [1, 2]:
                resvoloption = 3
                print(
                    "Warning: Reservoir volume option should be 1, 2 or 3. GEOPHIRES will assume default reservoir volume option (3)")
            else:
                resvoloption = 4
                print(
                    "Warning: Reservoir volume option should be 1, 2, 3, or 4. GEOPHIRES will assume default reservoir volume option (4)")
    except:
        if resoption in [1, 2]:
            resvoloption = 3
            print(
                "Warning: No valid reservoir volume option provided. GEOPHIRES will assume default reservoir volume option (3)")
        else:
            resvoloption = 4
            print(
                "Warning: No valid reservoir volume option provided. GEOPHIRES will assume default reservoir volume option (4)")

    if resvoloption == 4 and resoption in [1, 2]:
        resvoloption = 3
        print(
            "Warning: If user-selected reservoir model is 1 or 2, then user-selected reservoir volume option cannot be 4 but should be 1, 2, or 3. GEOPHIRES will assume reservoir volume option 3.")

    if resoption in [1, 2] or resvoloption in [1, 2, 3]:  # the first two reservoir models require fracture geometry
        # fracshape: Shape of fractures
        #   fracshape = 1  Circular fracture with known area
        #   fracshape = 2  Circular fracture with known diameter
        #   fracshape = 3  Square fracture
        #   fracshape = 4  Rectangular fracture
        try:
            fracshape = int(
                content[[i for i, s in enumerate(content) if 'Fracture Shape,' in s][0]].split(',')[1].strip('\n'))
            if not (fracshape in [1, 2, 3, 4]):
                fracshape = 1
                print(
                    "Warning: Provided fracture shape should be 1, 2, 3, or 4. GEOPHIRES will assume default fracture shape (1)")
        except:
            fracshape = 1
            print("Warning: No valid fracture shape provided. GEOPHIRES will assume default fracture shape (1)")

        # fracarea: Effective heat transfer area per fracture (m2) (required if fracshape = 1)
        if fracshape == 1:
            try:
                fracarea = float(
                    content[[i for i, s in enumerate(content) if 'Fracture Area,' in s][0]].split(',')[1].strip('\n'))
                if fracarea < 1 or fracarea > 100000000:
                    fracarea = 250000
                    print(
                        "Warning: Provided fracture area outside of range 1-100000000. GEOPHIRES will assume default fracture area (250,000 m2)")
            except:
                fracarea = 250000
                print(
                    "Warning: No valid fracture area provided. GEOPHIRES will assume default fracture area (250,000 m2)")

        # fracheight: Height of fracture = well separation (m)
        if fracshape in [2, 3, 4]:
            try:
                fracheight = float(
                    content[[i for i, s in enumerate(content) if 'Fracture Height,' in s][0]].split(',')[1].strip('\n'))
                if fracheight < 1 or fracheight > 10000:
                    fracheight = 500
                    print(
                        "Warning: Provided fracture height outside of range 1-10000. GEOPHIRES will assume default fracture height (500 m)")
            except:
                fracheight = 500
                print(
                    "Warning: No valid fracture height provided. GEOPHIRES will assume default fracture height (500 m)")

        # fracwidth: Width of fracture (m)
        if fracshape == 4:
            try:
                fracwidth = float(
                    content[[i for i, s in enumerate(content) if 'Fracture Width,' in s][0]].split(',')[1].strip('\n'))
                if fracwidth < 1 or fracwidth > 10000:
                    fracwidth = 500
                    print(
                        "Warning: Provided fracture width outside of range 1-10000. GEOPHIRES will assume default fracture width (500 m)")
            except:
                fracwidth = 500
                print("Warning: No valid fracture width provided. GEOPHIRES will assume default fracture width (500 m)")

        # calculate fracture geometry:
        # fracshape = 1: calculate diameter of circular fracture
        # fracshape = 2: calculate area of circular fracture
        # fracshape = 3: calculate area of square fracture
        # fracshape = 4: calculate area of rectangular fracture
        if fracshape == 1:
            fracheight = math.sqrt(4 / math.pi * fracarea)
            fracwidth = fracheight
        elif fracshape == 2:
            fracwidth = fracheight
            fracarea = math.pi / 4 * fracheight * fracheight
        elif fracshape == 3:
            fracwidth = fracheight
            fracarea = fracheight * fracwidth
        elif fracshape == 4:
            fracarea = fracheight * fracwidth

    # fracnumb: number of fractures
    if resvoloption in [1, 3]:
        try:
            fracnumb = int(
                content[[i for i, s in enumerate(content) if 'Number of Fractures,' in s][0]].split(',')[1].strip('\n'))
            if not (fracnumb in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]):
                fracnumb = 10
                print(
                    "Warning: Provided number of fractures outside of range 1-20. GEOPHIRES will assume default number of fractures (10)")
        except:
            fracnumb = 10
            print(
                "Warning: No valid number of fractures provided. GEOPHIRES will assume default number of fractures (10)")

    # fracsep: fracture separation [m]
    if resvoloption in [1, 2]:
        try:
            fracsep = float(
                content[[i for i, s in enumerate(content) if 'Fracture Separation,' in s][0]].split(',')[1].strip('\n'))
            if fracsep < 1 or fracsep > 10000:
                print(
                    "Warning: Provided fracture separation outside of range 1-10000. GEOPHIRES will assume default fracture separation (50 m)")
                fracsep = 50
        except:
            fracsep = 50
            print(
                "Warning: No valid fracture separation provided. GEOPHIRES will assume default fracture separation (50 m)")

    # resvol: reservoir volume [m^3]
    if resvoloption in [2, 3, 4]:
        try:
            resvol = float(
                content[[i for i, s in enumerate(content) if 'Reservoir Volume,' in s][0]].split(',')[1].strip('\n'))
            if resvol < 10 or resvol > 10000 * 10000 * 10000:
                print(
                    "Warning: Provided reservoir volume outside of range 10-1E12. GEOPHIRES will assume default reservoir volume (1.25E8 m3)")
                resvol = 500. * 500 * 500
        except:
            resvol = 500. * 500 * 500
            print(
                "Warning: No valid reservoir volume provided. GEOPHIRES will assume default reservoir volume (1.25E8 m3)")

    # calculate reservoir geometry:
    # resvoloption = 1: calculate volume of fractured rock mass
    # resvoloption = 2: calculate number of fractures
    # resvoloption = 3: calculate fracture separation
    if resvoloption == 1:
        resvol = (fracnumb - 1) * fracarea * fracsep
    elif resvoloption == 2:
        fracnumb = resvol / fracarea / fracsep + 1
    elif resvoloption == 3:
        fracsep = resvol / fracarea / (fracnumb - 1)

    # waterloss: fraction of water lost = (total geofluid lost)/(total geofluid produced)
    try:
        waterloss = float(
            content[[i for i, s in enumerate(content) if 'Water Loss Fraction,' in s][0]].split(',')[1].strip('\n'))
        if waterloss < 0 or waterloss > 0.99:
            waterloss = 0
            print(
                "Warning: Provided water loss fraction outside of range 0-0.99. GEOPHIRES will assume default water loss fraction (0)")
    except:
        waterloss = 0
        print("Warning: No valid water loss fraction provided. GEOPHIRES will assume default water loss fraction (0)")

    impedancemodelallowed = 1
    productionwellpumping = 1
    setinjectionpressurefixed = 0

    if enduseoption == 1:
        if pptype in [3, 4]:  # simple single- or double-flash power plant assumes no production well pumping
            impedancemodelallowed = 0
            productionwellpumping = 0
            setinjectionpressurefixed = 1
    elif enduseoption in [31, 32]:
        if pptype in [3,
                      4]:  # co-generation topping cycle with single- or double-flash power plant assumes no production well pumping
            impedancemodelallowed = 0
            productionwellpumping = 0
            setinjectionpressurefixed = 1
    elif enduseoption in [41, 42]:
        if pptype in [3,
                      4]:  # co-generation bottoming cycle with single- or double-flash power plant assumes production well pumping
            impedancemodelallowed = 0
            setinjectionpressurefixed = 1
    elif enduseoption in [51, 52]:
        if pptype in [3,
                      4]:  # co-generation parallel cycle with single- or double-flash power plant assumes production well pumping
            impedancemodelallowed = 0
            setinjectionpressurefixed = 1

    impedancemodelused = 0
    if impedancemodelallowed == 1:
        try:
            # impedance: impedance per wellpair (input as GPa*s/m^3 and converted to KPa/kg/s (assuming 1000 for density; density will be corrected for later))
            impedance = float(
                content[[i for i, s in enumerate(content) if 'Reservoir Impedance,' in s][0]].split(',')[1].strip(
                    '\n')) * 1E6 / 1E3
            impedancemodelused = 1
            if impedance < 0.0001 * 1000 or impedance > 10000:
                impedance = 0.1 * 1E6 / 1E3
                print(
                    "Warning: Provided reservoir impedance outside of range 0.0001-1000. GEOPHIRES will assume default reservoir impedance (0.1 GPa*s/m3)")
        except:
            impedancemodelused = 0

    if impedancemodelallowed == 0 or impedancemodelused == 0:
        try:
            # reservoir hydrostatic pressure [kPa]
            Phydrostatic = float(
                content[[i for i, s in enumerate(content) if 'Reservoir Hydrostatic Pressure,' in s][0]].split(',')[
                    1].strip('\n'))
            usebuiltinhydrostaticpressurecorrelation = 0
            if Phydrostatic < 100 or Phydrostatic > 100000:
                usebuiltinhydrostaticpressurecorrelation = 1
                print(
                    "Warning: Provided reservoir hydrostatic pressure outside of range 100-100000 kPa. GEOPHIRES will assume built-in reservoir hydrostatic pressure correlation")
        except:
            usebuiltinhydrostaticpressurecorrelation = 1
            print(
                "Warning: No valid reservoir hydrostatic pressure provided. GEOPHIRES will assume built-in reservoir hydrostatic pressure correlation")

        try:
            # injectivity index [kg/s/bar]
            II = float(
                content[[i for i, s in enumerate(content) if 'Injectivity Index,' in s][0]].split(',')[1].strip('\n'))
            if II < 0.01 or II > 10000:
                II = 10
                print(
                    "Warning: Provided injectivity index outside of range 0.01-10000. GEOPHIRES will assume default injectivity index (10 kg/s/bar)")
        except:
            II = 10
            print(
                "Warning: No valid injectivity index provided. GEOPHIRES will assume default injectivity index (10 kg/s/bar)")

        if productionwellpumping == 1:
            try:
                # productivity index [kg/s/bar]
                PI = float(
                    content[[i for i, s in enumerate(content) if 'Productivity Index,' in s][0]].split(',')[1].strip(
                        '\n'))
                if PI < 0.01 or PI > 10000:
                    PI = 10
                    print(
                        "Warning: Provided productivity index outside of range 0.01-10000. GEOPHIRES will assume default productivity index (10 kg/s/bar)")
            except:
                PI = 10
                print(
                    "Warning: No valid productivity index provided. GEOPHIRES will assume default productivity index (10 kg/s/bar)")

            try:
                # production wellhead pressure [kPa]
                ppwellhead = float(
                    content[[i for i, s in enumerate(content) if 'Production Wellhead Pressure,' in s][0]].split(',')[
                        1].strip('\n'))
                usebuiltinppwellheadcorrelation = 0
                if ppwellhead < 0 or ppwellhead > 10000:
                    usebuiltinppwellheadcorrelation = 1
                    print(
                        "Warning: Provided production wellhead pressure outside of range 0-10000 kPa. GEOPHIRES will calculate production wellhead pressure using built-in correlation")
            except:
                usebuiltinppwellheadcorrelation = 1
                print(
                    "Warning: No valid production wellhead pressure provided. GEOPHIRES will calculate production wellhead pressure using built-in correlation")

        try:
            # plant outlet pressure [kPa]
            Pplantoutlet = float(
                content[[i for i, s in enumerate(content) if 'Plant Outlet Pressure,' in s][0]].split(',')[1].strip(
                    '\n'))
            usebuiltinoutletplantcorrelation = 0
            if Pplantoutlet < 0 or Pplantoutlet > 10000:
                if setinjectionpressurefixed == 1:
                    Pplantoutlet = 100
                    print(
                        "Warning: Provided plant outlet pressure outside of range 0-10000. GEOPHIRES will assume default plant outlet pressure (100 kPa)")
                else:
                    usebuiltinoutletplantcorrelation = 1
                    print(
                        "Warning: Provided plant outlet pressure outside of range 0-10000 kPa. GEOPHIRES will calculate plant outlet pressure based on production wellhead pressure and surface equipment pressure drop of 10 psi")
        except:
            if setinjectionpressurefixed == 1:
                usebuiltinoutletplantcorrelation = 0
                Pplantoutlet = 100
                print(
                    "Warning: No valid plant outlet pressure provided. GEOPHIRES will assume default plant outlet pressure (100 kPa)")
            else:
                usebuiltinoutletplantcorrelation = 1
                print(
                    "Warning: No valid plant outlet pressure provided. GEOPHIRES will calculate plant outlet pressure based on production wellhead pressure and surface equipment pressure drop of 10 psi")

    # impedance: impedance per wellpair (input as GPa*s/m^3 and converted to KPa/kg/s (assuming 1000 for density))
    # try:
    #    impedance = float(content[[i for i, s in enumerate(content) if 'Reservoir Impedance,' in s][0]].split(',')[1].strip('\n'))*1E6/1E3
    #    if impedance < 0.0001*1000 or impedance > 10000:
    #        impedance = 0.1*1E6/1E3
    #        print("Warning: Provided reservoir impedance outside of range 0.0001-1000. GEOPHIRES will assume default reservoir impedance (0.1 GPa*s/m3)")
    # except:
    #    impedance = 0.1*1E6/1E3
    #    print("Warning: No valid reservoir impedance provided. GEOPHIRES will assume default reservoir impedance (0.1 GPa*s/m3)")

    # maxdrawdown: maximum allowable drawdown before redrilling (only works with built in reservoir models)
    if resoption in [1, 2, 3, 4]:
        try:
            maxdrawdown = float(
                content[[i for i, s in enumerate(content) if 'Maximum Drawdown,' in s][0]].split(',')[1].strip('\n'))
            if maxdrawdown < 0 or maxdrawdown > 1:
                maxdrawdown = 1
                print(
                    "Warning: Provided maximum drawdown outside of range 0-1. GEOPHIRES will assume default maximum drawdown (1)")
        except:
            maxdrawdown = 1
            print("Warning: No valid maximum drawdown provided. GEOPHIRES will assume default maximum drawdown (1)")

    # cprock: reservoir heat capacity (in J/kg/K)
    try:
        cprock = float(
            content[[i for i, s in enumerate(content) if 'Reservoir Heat Capacity,' in s][0]].split(',')[1].strip('\n'))
        if cprock < 100 or cprock > 10000:
            cprock = 1000
            print(
                "Warning: Provided reservoir heat capacity outside of range 100-10000. GEOPHIRES will assume default reservoir heat capacity (1000 J/kg/K)")
    except:
        cprock = 1000
        print(
            "Warning: No valid reservoir heat capacity provided. GEOPHIRES will assume default reservoir heat capacity (1000 J/kg/K)")

    # rhorock: reservoir density (in kg/m3)
    try:
        rhorock = float(
            content[[i for i, s in enumerate(content) if 'Reservoir Density,' in s][0]].split(',')[1].strip('\n'))
        if rhorock < 100 or rhorock > 20000:
            rhorock = 2700
            print(
                "Warning: Provided reservoir density outside of range 100-10000. GEOPHIRES will assume default reservoir density (2700 J/kg/K)")
    except:
        rhorock = 2700
        print(
            "Warning: No valid reservoir density provided. GEOPHIRES will assume default reservoir density (2700 J/kg/K)")

    # krock: reservoir thermal conductivity (in W/m/K)
    if rameyoptionprod == 1 or resoption in [1, 2, 3] or (resoption == 6 and usebuiltintough2model == 1):
        try:
            krock = float(
                content[[i for i, s in enumerate(content) if 'Reservoir Thermal Conductivity,' in s][0]].split(',')[
                    1].strip('\n'))
            if krock < 0.01 or krock > 100:
                krock = 3
                print(
                    "Warning: Provided reservoir thermal conductivity outside of range 0.01-100. GEOPHIRES will assume default reservoir thermal conductivity (3 W/m/K)")
        except:
            krock = 3
            print(
                "Warning: No valid reservoir thermal conductivity provided. GEOPHIRES will assume default reservoir thermal conductivity (3 W/m/K)")

    # porrock: reservoir porosity (-)
    if resoption == 2 or (resoption == 6 and usebuiltintough2model == 1):
        try:
            porrock = float(
                content[[i for i, s in enumerate(content) if 'Reservoir Porosity,' in s][0]].split(',')[1].strip('\n'))
            if porrock < 0.001 or porrock > 0.99:
                porrock = 0.04
                print(
                    "Warning: Provided reservoir porosity outside of range 0.001-0.99. GEOPHIRES will assume default reservoir porosity (0.04)")
        except:
            porrock = 0.04
            print(
                "Warning: No valid reservoir porosity provided. GEOPHIRES will assume default reservoir porosity (0.04)")

    # permrock: reservoir permeability (m2)
    if resoption == 6 and usebuiltintough2model == 1:
        try:
            permrock = float(
                content[[i for i, s in enumerate(content) if 'Reservoir Permeability,' in s][0]].split(',')[1].strip(
                    '\n'))
            if permrock < 1E-20 or permrock > 1E-5:
                permrock = 1E-13
                print(
                    "Warning: Provided reservoir permeability outside of range 1E-20 to 1E-5. GEOPHIRES will assume default reservoir permeability (1E-13 m^2)")
        except:
            permrock = 1E-13
            print(
                "Warning: No valid reservoir permeability provided. GEOPHIRES will assume default reservoir permeability (1E-13 m^2)")

    # resthickness: reservoir thickness (m)
    if resoption == 6 and usebuiltintough2model == 1:
        try:
            resthickness = float(
                content[[i for i, s in enumerate(content) if 'Reservoir Thickness,' in s][0]].split(',')[1].strip('\n'))
            if resthickness < 10 or resthickness > 10000:
                resthickness = 250
                print(
                    "Warning: Provided reservoir thickness outside of range 10-10000. GEOPHIRES will assume default reservoir thickness (250 m)")
        except:
            resthickness = 250
            print(
                "Warning: No valid reservoir thickness provided (necessary for using built-in TOUGH2 model). GEOPHIRES will assume default reservoir thickness (250 m)")

    # reswidth: reservoir width (m)
    if resoption == 6 and usebuiltintough2model == 1:
        try:
            reswidth = float(
                content[[i for i, s in enumerate(content) if 'Reservoir Width,' in s][0]].split(',')[1].strip('\n'))
            if reswidth < 10 or reswidth > 10000:
                reswidth = 500
                print(
                    "Warning: Provided reservoir width outside of range 10-10000. GEOPHIRES will assume default reservoir width (500 m)")
        except:
            reswidth = 500
            print(
                "Warning: No valid reservoir width provided (necessary for using built-in TOUGH2 model). GEOPHIRES will assume default reservoir width (500 m)")

    # wellsep: well separation (m)
    if resoption == 6 and usebuiltintough2model == 1:
        try:
            wellsep = float(
                content[[i for i, s in enumerate(content) if 'Well Separation,' in s][0]].split(',')[1].strip('\n'))
            if wellsep < 10 or wellsep > 10000:
                wellsep = 1000
                print(
                    "Warning: Provided well seperation outside of range 10-10000. GEOPHIRES will assume default well seperation (1000 m)")
        except:
            wellsep = 1000
            print(
                "Warning: No valid well seperation provided (necessary for using built-in TOUGH2 model). GEOPHIRES will assume default well seperation (1000 m)")

    # plantlifetime: plant lifetime (years)
    try:
        plantlifetime = int(
            content[[i for i, s in enumerate(content) if 'Plant Lifetime,' in s][0]].split(',')[1].strip('\n'))
        if not (plantlifetime in list(range(1, 101))):
            plantlifetime = 30
            print(
                "Warning: Provided plant lifetime outside of range 1-100. GEOPHIRES will assume default plant lifetime (30 years)")
    except:
        plantlifetime = 30
        print("Warning: No valid plant lifetime provided. GEOPHIRES will assume default plant lifetime (30 years)")

    # econmodel
    # econmodel = 1: use Fixed Charge Rate Model (requires an FCR)
    # econmodel = 2: use standard LCOE/LCOH calculation as found on wikipedia (requries an interest rate).
    # econmodel = 3: use Bicycle LCOE/LCOH model (requires several financial input parameters)
    try:
        econmodel = int(
            content[[i for i, s in enumerate(content) if 'Economic Model,' in s][0]].split(',')[1].strip('\n'))
        if not (econmodel in [1, 2, 3]):
            econmodel = 2
            print(
                "Warning: Provided economic model should be 1, 2, or 3. GEOPHIRES will assume default economic model (2)")
    except:
        econmodel = 2
        print("Warning: No valid economic model provided. GEOPHIRES will assume default economic model (2)")

    # FCR: fixed charge rate required if econmodel = 1
    if econmodel == 1:
        try:
            FCR = float(
                content[[i for i, s in enumerate(content) if 'Fixed Charge Rate,' in s][0]].split(',')[1].strip('\n'))
            if FCR < 0 or FCR > 1:
                FCR = 0.1
                print(
                    "Warning: Provided fixed charge rate is outside of range 0-1. GEOPHIRES will assume default fixed charge rate (0.1)")
        except:
            FCR = 0.1
            print("Warning: No valid fixed charge rate provided. GEOPHIRES will assume default fixed charge rate (0.1)")

    # discountrate: discount rate required if econmodel = 2
    if econmodel == 2:
        try:
            discountrate = float(
                content[[i for i, s in enumerate(content) if 'Discount Rate,' in s][0]].split(',')[1].strip('\n'))
            if discountrate < 0 or discountrate > 1:
                discountrate = 0.07
                print(
                    "Warning: Provided discount rate is outside of range 0-1. GEOPHIRES will assume default discount rate (0.07)")
        except:
            discountrate = 0.07
            print("Warning: No valid discount rate provided. GEOPHIRES will assume default discount rate (0.07)")

    # a whole bunch of BICYCLE parameters provided if econmodel = 3
    if econmodel == 3:
        # bicycle parameters
        # FIB: fraction of investment in bonds (-)
        try:
            FIB = float(
                content[[i for i, s in enumerate(content) if 'Fraction of Investment in Bonds,' in s][0]].split(',')[
                    1].strip('\n'))
            if FIB < 0 or FIB > 1:
                FIB = 0.5
                print(
                    "Warning: Provided fraction of investment in bonds is outside of range 0-1. GEOPHIRES will assume default fraction of investment in bonds (0.5)")
        except:
            FIB = 0.5
            print(
                "Warning: No valid fraction of investment in bonds provided. GEOPHIRES will assume default fraction of investment in bonds (0.5)")

        # BIR: inflated bonds interest rate (-)
        try:
            BIR = float(
                content[[i for i, s in enumerate(content) if 'Inflated Bond Interest Rate,' in s][0]].split(',')[
                    1].strip('\n'))
            if BIR < 0 or BIR > 1:
                BIR = 0.05
                print(
                    "Warning: Provided inflated bond interest rate is outside of range 0-1. GEOPHIRES will assume default inflated bond interest rate (0.05)")
        except:
            BIR = 0.05
            print(
                "Warning: No valid inflated bond interest rate provided. GEOPHIRES will assume default inflated bond interest rate (0.05)")

        # EIR: inflated equity interest rate (-)
        try:
            EIR = float(
                content[[i for i, s in enumerate(content) if 'Inflated Equity Interest Rate,' in s][0]].split(',')[
                    1].strip('\n'))
            if EIR < 0 or EIR > 1:
                EIR = 0.1
                print(
                    "Warning: Provided inflated equity interest rate is outside of range 0-1. GEOPHIRES will assume default inflated equity interest rate (0.1)")
        except:
            EIR = 0.1
            print(
                "Warning: No valid inflated equity interest rate provided. GEOPHIRES will assume default inflated equity interest rate (0.1)")

        # RINFL: inflation rate (-)
        try:
            RINFL = float(
                content[[i for i, s in enumerate(content) if 'Inflation Rate,' in s][0]].split(',')[1].strip('\n'))
            if RINFL < -0.1 or RINFL > 1:
                RINFL = 0.02
                print(
                    "Warning: Provided inflation rate is outside of range -0.1 to 1. GEOPHIRES will assume default inflation rate (0.02)")
        except:
            RINFL = 0.02
            print("Warning: No valid inflation rate provided. GEOPHIRES will assume default inflation rate (0.02)")

        # CTR: combined income tax rate in fraction (-)
        try:
            CTR = float(
                content[[i for i, s in enumerate(content) if 'Combined Income Tax Rate,' in s][0]].split(',')[1].strip(
                    '\n'))
            if CTR < 0 or CTR > 1:
                CTR = 0.3
                print(
                    "Warning: Provided combined income tax rate is outside of range 0 to 1. GEOPHIRES will assume default combined income tax rate (0.3)")
        except:
            CTR = 0.3
            print(
                "Warning: No valid combined income tax rate provided. GEOPHIRES will assume default combined income tax rate (0.3)")

        # GTR: gross revenue tax rate in fraction (-)
        try:
            GTR = float(
                content[[i for i, s in enumerate(content) if 'Gross Revenue Tax Rate,' in s][0]].split(',')[1].strip(
                    '\n'))
            if GTR < 0 or GTR > 1:
                GTR = 0
                print(
                    "Warning: Provided gross revenue tax rate is outside of range 0 to 1. GEOPHIRES will assume default gross revenue tax rate (0)")
        except:
            GTR = 0
            print(
                "Warning: No valid gross revenue tax rate provided. GEOPHIRES will assume default gross revenue tax rate (0)")

        # RITC: investment tax credit rate in fraction (-)
        try:
            RITC = float(
                content[[i for i, s in enumerate(content) if 'Investment Tax Credit Rate,' in s][0]].split(',')[
                    1].strip('\n'))
            if RITC < 0 or RITC > 1:
                RITC = 0
                print(
                    "Warning: Provided investment tax credit rate is outside of range 0 to 1. GEOPHIRES will assume default investment tax credit rate (0)")
        except:
            RITC = 0
            print(
                "Warning: No valid investment tax credit rate provided. GEOPHIRES will assume default investment tax credit rate (0)")

        # PTR: property tax rate in fraction (-)
        try:
            PTR = float(
                content[[i for i, s in enumerate(content) if 'Property Tax Rate,' in s][0]].split(',')[1].strip('\n'))
            if PTR < 0 or PTR > 1:
                PTR = 0
                print(
                    "Warning: Provided property rate is outside of range 0 to 1. GEOPHIRES will assume default property tax rate (0)")
        except:
            PTR = 0
            print("Warning: No valid property tax rate provided. GEOPHIRES will assume default property tax rate (0)")

    # inflrateconstruction: inflation rate during construction (-)
    try:
        inflrateconstruction = float(
            content[[i for i, s in enumerate(content) if 'Inflation Rate During Construction,' in s][0]].split(',')[
                1].strip('\n'))
        if inflrateconstruction < 0 or inflrateconstruction > 1:
            inflrateconstruction = 0
            print(
                "Warning: Provided inflation rate during construction is outside of range 0 to 1. GEOPHIRES will assume default inflation rate during construction (0)")
    except:
        inflrateconstruction = 0
        print(
            "Warning: No valid inflation rate during construction provided. GEOPHIRES will assume default inflation rate during construction (0)")

    # capital cost parameters
    try:  # user can provide total capital cost (M$)
        totalcapcost = float(
            content[[i for i, s in enumerate(content) if 'Total Capital Cost,' in s][0]].split(',')[1].strip('\n'))
        totalcapcostprovided = 1
        if totalcapcost < 0 or totalcapcost > 1000:
            totalcapcostvalid = 0
            print(
                "Warning: Provided total capital cost outside of range 0 to 1000. GEOPHIRES will calculate total capital cost using user-provided costs or built-in correlations for each category.")
        else:
            totalcapcostvalid = 1
    except:
        totalcapcostprovided = 0
        totalcapcostvalid = 0

    # ccwellfixed: well drilling and completion capital cost in M$ (per well)
    try:
        ccwellfixed = float(
            content[[i for i, s in enumerate(content) if 'Well Drilling and Completion Capital Cost,' in s][0]].split(
                ',')[1].strip('\n'))
        ccwellfixedprovided = 1
        if ccwellfixed < 0 or ccwellfixed > 200:
            ccwellfixedvalid = 0
        else:
            ccwellfixedvalid = 1
    except:
        ccwellfixedprovided = 0
        ccwellfixedvalid = 0

    # ccwelladjfactor: adj factor for built-in correlation well drilling and completion cost
    try:
        ccwelladjfactor = float(content[[i for i, s in enumerate(content) if
                                         'Well Drilling and Completion Capital Cost Adjustment Factor,' in s][0]].split(
            ',')[1].strip('\n'))
        ccwelladjfactorprovided = 1
        if ccwelladjfactor < 0 or ccwelladjfactor > 10:
            ccwelladjfactorvalid = 0
        else:
            ccwelladjfactorvalid = 1
    except:
        ccwelladjfactorprovided = 0
        ccwelladjfactorvalid = 0

    if ccwellfixedvalid == 1 and ccwelladjfactorvalid == 1:
        print(
            "Warning: Provided well drilling and completion cost adjustment factor not considered because valid total well drilling and completion cost provided.")
    elif ccwellfixedprovided == 0 and ccwelladjfactorprovided == 0:
        ccwelladjfactor = 1
        print(
            "Warning: No valid well drilling and completion total cost or adjustment factor provided. GEOPHIRES will assume default built-in well drilling and completion cost correlation with adjustment factor = 1.")
    elif ccwellfixedprovided == 1 and ccwellfixedvalid == 0:
        print(
            "Provided well drilling and completion cost outside of range 0-1000. GEOPHIRES will assume default built-in well drilling and completion cost correlation with adjustment factor = 1.")
        ccwelladjfactor = 1
    elif ccwellfixedprovided == 0 and ccwelladjfactorprovided == 1 and ccwelladjfactorvalid == 0:
        print(
            "Provided well drilling and completion cost adjustment factor outside of range 0-10. GEOPHIRES will assume default built-in well drilling and completion cost correlation with adjustment factor = 1.")
        ccwelladjfactor = 1

    # Drilling cost correlation (should be 1, 2, 3, or 4) if no valid fixed well drilling cost is provided
    if ccwellfixedvalid == 0:
        try:
            wellcorrelation = int(
                content[[i for i, s in enumerate(content) if 'Well Drilling Cost Correlation,' in s][0]].split(',')[
                    1].strip('\n'))
            if not (wellcorrelation in [1, 2, 3, 4]):
                wellcorrelation = 1
                print(
                    "Warning: Selected well drilling cost correlation number should be 1, 2, 3 or 4. GEOPHIRES will assume default well drilling cost correlation (1)")
        except:
            wellcorrelation = 1
            print(
                "Warning: No valid well drilling cost correlation number provided. GEOPHIRES will assume default well drilling cost correlation (1)")
    # ccstimfixed: reservoir stimulation cost in M$
    try:
        ccstimfixed = float(
            content[[i for i, s in enumerate(content) if 'Reservoir Stimulation Capital Cost,' in s][0]].split(',')[
                1].strip('\n'))
        ccstimfixedprovided = 1
        if ccstimfixed < 0 or ccstimfixed > 100:
            ccstimfixedvalid = 0
        else:
            ccstimfixedvalid = 1
    except:
        ccstimfixedprovided = 0
        ccstimfixedvalid = 0

    # ccstimadjfactor: adj factor for built-in correlation for reservoir stimulation cost
    try:
        ccstimadjfactor = float(content[[i for i, s in enumerate(content) if
                                         'Reservoir Stimulation Capital Cost Adjustment Factor,' in s][0]].split(',')[
                                    1].strip('\n'))
        ccstimadjfactorprovided = 1
        if ccstimadjfactor < 0 or ccstimadjfactor > 10:
            ccstimadjfactorvalid = 0
        else:
            ccstimadjfactorvalid = 1
    except:
        ccstimadjfactorprovided = 0
        ccstimadjfactorvalid = 0

    if ccstimfixedvalid == 1 and ccstimadjfactorvalid == 1:
        print(
            "Warning: Provided reservoir stimulation cost adjustment factor not considered because valid total reservoir stimulation cost provided.")
    elif ccstimfixedprovided == 0 and ccstimadjfactorprovided == 0:
        ccstimadjfactor = 1
        print(
            "Warning: No valid reservoir stimulation total cost or adjustment factor provided. GEOPHIRES will assume default built-in reservoir stimulation cost correlation with adjustment factor = 1.")
    elif ccstimfixedprovided == 1 and ccstimfixedvalid == 0:
        print(
            "Provided reservoir stimulation cost outside of range 0-100. GEOPHIRES will assume default built-in reservoir stimulation cost correlation with adjustment factor = 1.")
        ccstimadjfactor = 1
    elif ccstimfixedprovided == 0 and ccstimadjfactorprovided == 1 and ccstimadjfactorvalid == 0:
        print(
            "Provided reservoir stimulation cost adjustment factor outside of range 0-10. GEOPHIRES will assume default reservoir stimulation cost correlation with adjustment factor = 1.")
        ccstimadjfactor = 1

    # ccplantfixed: surface plant cost in M$
    try:
        ccplantfixed = float(
            content[[i for i, s in enumerate(content) if 'Surface Plant Capital Cost,' in s][0]].split(',')[1].strip(
                '\n'))
        ccplantfixedprovided = 1
        if ccplantfixed < 0 or ccplantfixed > 1000:
            ccplantfixedvalid = 0
        else:
            ccplantfixedvalid = 1
    except:
        ccplantfixedprovided = 0
        ccplantfixedvalid = 0

    # ccplantadjfactor: adj factor for built-in surface plant cost correlation
    try:
        ccplantadjfactor = float(content[[i for i, s in enumerate(content) if
                                          'Surface Plant Capital Cost Adjustment Factor,' in s][0]].split(',')[1].strip(
            '\n'))
        ccplantadjfactorprovided = 1
        if ccplantadjfactor < 0 or ccplantadjfactor > 10:
            ccplantadjfactorvalid = 0
        else:
            ccplantadjfactorvalid = 1
    except:
        ccplantadjfactorprovided = 0
        ccplantadjfactorvalid = 0
        ccplantadjfactor = 1

    if totalcapcostvalid == 1:
        if ccplantfixedprovided == 1:
            print("Warning: Provided surface plant cost not considered because valid total capital cost provided.")
        if ccplantadjfactorprovided == 1:
            print(
                "Warning: Provided surface plant cost adjustment factor not considered because valid total capital cost provided.")
    else:
        if ccplantfixedvalid == 1 and ccplantadjfactorvalid == 1:
            print(
                "Warning: Provided surface plant cost adjustment factor not considered because valid total surface plant cost provided.")
        elif ccplantfixedprovided == 0 and ccplantadjfactorprovided == 0:
            ccplantadjfactor = 1
            print(
                "Warning: No valid surface plant total cost or adjustment factor provided. GEOPHIRES will assume default built-in surface plant cost correlation with adjustment factor = 1.")
        elif ccplantfixedprovided == 1 and ccplantfixedvalid == 0:
            print(
                "Provided surface plant cost outside of range 0-1000. GEOPHIRES will assume default built-in surface plant cost correlation with adjustment factor = 1.")
            ccplantadjfactor = 1
        elif ccplantfixedprovided == 0 and ccplantadjfactorprovided == 1 and ccplantadjfactorvalid == 0:
            print(
                "Provided surface plant cost adjustment factor outside of range 0-10. GEOPHIRES will assume default surface plant cost correlation with adjustment factor = 1.")
            ccplantadjfactor = 1

    # ccgathfixed: field gathering system network cost in M$
    try:
        ccgathfixed = float(
            content[[i for i, s in enumerate(content) if 'Field Gathering System Capital Cost,' in s][0]].split(',')[
                1].strip('\n'))
        ccgathfixedprovided = 1
        if ccgathfixed < 0 or ccgathfixed > 100:
            ccgathfixedvalid = 0
        else:
            ccgathfixedvalid = 1
    except:
        ccgathfixedprovided = 0
        ccgathfixedvalid = 0

    # ccgathadjfactor: adj factor for built-in field gathering system cost correlation
    try:
        ccgathadjfactor = float(content[[i for i, s in enumerate(content) if
                                         'Field Gathering System Capital Cost Adjustment Factor,' in s][0]].split(',')[
                                    1].strip('\n'))
        ccgathadjfactorprovided = 1
        if ccgathadjfactor < 0 or ccgathadjfactor > 10:
            ccgathadjfactorvalid = 0
        else:
            ccgathadjfactorvalid = 1
    except:
        ccgathadjfactorprovided = 0
        ccgathadjfactorvalid = 0
        ccgathadjfactor = 1

    if totalcapcostvalid == 1:
        if ccgathfixedprovided == 1:
            print(
                "Warning: Provided field gathering system cost not considered because valid total capital cost provided.")
        if ccgathadjfactorprovided == 1:
            print(
                "Warning: Provided field gathering system cost adjustment factor not considered because valid total capital cost provided.")
    else:
        if ccgathfixedvalid == 1 and ccgathadjfactorvalid == 1:
            print(
                "Warning: Provided field gathering system cost adjustment factor not considered because valid total field gathering system cost provided.")
        elif ccgathfixedprovided == 0 and ccgathadjfactorprovided == 0:
            ccgathadjfactor = 1
            print(
                "Warning: No valid field gathering system total cost or adjustment factor provided. GEOPHIRES will assume default built-in field gathering system cost correlation with adjustment factor = 1.")
        elif ccgathfixedprovided == 1 and ccgathfixedvalid == 0:
            print(
                "Provided field gathering system cost outside of range 0-100. GEOPHIRES will assume default built-in field gathering system cost correlation with adjustment factor = 1.")
            ccgathadjfactor = 1
        elif ccgathfixedprovided == 0 and ccgathadjfactorprovided == 1 and ccgathadjfactorvalid == 0:
            print(
                "Provided field gathering system cost adjustment factor outside of range 0-10. GEOPHIRES will assume default field gathering system cost correlation with adjustment factor = 1.")
            ccgathadjfactor = 1

    # ccexplfixed: exploration cost in M$
    try:
        ccexplfixed = float(
            content[[i for i, s in enumerate(content) if 'Exploration Capital Cost,' in s][0]].split(',')[1].strip(
                '\n'))
        ccexplfixedprovided = 1
        if ccexplfixed < 0 or ccexplfixed > 100:
            ccexplfixedvalid = 0
        else:
            ccexplfixedvalid = 1
    except:
        ccexplfixedprovided = 0
        ccexplfixedvalid = 0

    # ccexpladjfactor: adj factor for built-in exploration cost correlation
    try:
        ccexpladjfactor = float(
            content[[i for i, s in enumerate(content) if 'Exploration Capital Cost Adjustment Factor,' in s][0]].split(
                ',')[1].strip('\n'))
        ccexpladjfactorprovided = 1
        if ccexpladjfactor < 0 or ccexpladjfactor > 10:
            ccexpladjfactorvalid = 0
        else:
            ccexpladjfactorvalid = 1
    except:
        ccexpladjfactorprovided = 0
        ccexpladjfactorvalid = 0

    if totalcapcostvalid == 1:
        if ccexplfixedprovided == 1:
            print("Warning: Provided exploration cost not considered because valid total capital cost provided.")
        if ccexpladjfactorprovided == 1:
            print(
                "Warning: Provided exploration cost adjustment factor not considered because valid total capital cost provided.")
    else:
        if ccexplfixedvalid == 1 and ccexpladjfactorvalid == 1:
            print(
                "Warning: Provided exploration cost adjustment factor not considered because valid total exploration cost provided.")
        elif ccexplfixedprovided == 0 and ccexpladjfactorprovided == 0:
            ccexpladjfactor = 1
            print(
                "Warning: No valid exploration total cost or adjustment factor provided. GEOPHIRES will assume default built-in exploration cost correlation with adjustment factor = 1.")
        elif ccexplfixedprovided == 1 and ccexplfixedvalid == 0:
            print(
                "Provided exploration cost outside of range 0-100. GEOPHIRES will assume default built-in exploration cost correlation with adjustment factor = 1.")
            ccexpladjfactor = 1
        elif ccexplfixedprovided == 0 and ccexpladjfactorprovided == 1 and ccexpladjfactorvalid == 0:
            print(
                "Provided exploration cost adjustment factor outside of range 0-10. GEOPHIRES will assume default exploration cost correlation with adjustment factor = 1.")
            ccexpladjfactor = 1

    # pipinglength: surface piping length (-)
    try:
        pipinglength = float(
            content[[i for i, s in enumerate(content) if 'Surface Piping Length,' in s][0]].split(',')[1].strip('\n'))
        if pipinglength < 0 or pipinglength > 100:
            pipinglength = 5
            print(
                "Warning: Provided surface transmission piping length outside of range 0-100. GEOPHIRES will assume default piping length (5km)")
    except:
        pipinglength = 0

    # O&M cost parameters
    # oamtotalfixed: total O&M cost in M$/year
    try:  # user can provide total O&M cost (M$)
        oamtotalfixed = float(
            content[[i for i, s in enumerate(content) if 'Total O&M Cost,' in s][0]].split(',')[1].strip('\n'))
        oamtotalfixedprovided = 1
        if oamtotalfixed < 0 or oamtotalfixed > 100:
            oamtotalfixedvalid = 0
            print(
                "Warning: Provided total annual O&M cost outside of range 0 to 100. GEOPHIRES will calculate total O&M cost using user-provided costs or built-in correlations for each category.")
        else:
            oamtotalfixedvalid = 1
    except:
        oamtotalfixedprovided = 0
        oamtotalfixedvalid = 0

    # oamwellfixed: total wellfield O&M cost in M$/year
    try:
        oamwellfixed = float(
            content[[i for i, s in enumerate(content) if 'Wellfield O&M Cost,' in s][0]].split(',')[1].strip('\n'))
        oamwellfixedprovided = 1
        if oamwellfixed < 0 or oamwellfixed > 100:
            oamwellfixedvalid = 0
        else:
            oamwellfixedvalid = 1
    except:
        oamwellfixedprovided = 0
        oamwellfixedvalid = 0

    # oamwelladjfactor: adj factor to built-in correlation for wellfield O&M cost
    try:
        oamwelladjfactor = float(
            content[[i for i, s in enumerate(content) if 'Wellfield O&M Cost Adjustment Factor,' in s][0]].split(',')[
                1].strip('\n'))
        oamwelladjfactorprovided = 1
        if oamwelladjfactor < 0 or oamwelladjfactor > 10:
            oamwelladjfactorvalid = 0
        else:
            oamwelladjfactorvalid = 1
    except:
        oamwelladjfactorprovided = 0
        oamwelladjfactorvalid = 0

    if oamtotalfixedvalid == 1:
        if oamwellfixedprovided == 1:
            print(
                "Warning: Provided total wellfield O&M cost not considered because valid total annual O&M cost provided.")
        if oamwelladjfactorprovided == 1:
            print(
                "Warning: Provided wellfield O&M cost adjustment factor not considered because valid total annual O&M cost provided.")
    else:
        if oamwellfixedvalid == 1 and oamwelladjfactorvalid == 1:
            print(
                "Warning: Provided wellfield O&M cost adjustment factor not considered because valid total wellfield O&M cost provided.")
        elif oamwellfixedprovided == 0 and oamwelladjfactorprovided == 0:
            oamwelladjfactor = 1
            print(
                "Warning: No valid total wellfield O&M cost or adjustment factor provided. GEOPHIRES will assume default built-in wellfield O&M cost correlation with adjustment factor = 1.")
        elif oamwellfixedprovided == 1 and oamwellfixedvalid == 0:
            print(
                "Provided total wellfield O&M cost outside of range 0-100. GEOPHIRES will assume default built-in wellfield O&M cost correlation with adjustment factor = 1.")
            oamwelladjfactor = 1
        elif oamwellfixedprovided == 0 and oamwelladjfactorprovided == 1 and oamwelladjfactorvalid == 0:
            print(
                "Provided wellfield O&M cost adjustment factor outside of range 0-10. GEOPHIRES will assume default wellfield O&M cost correlation with adjustment factor = 1.")
            oamwelladjfactor = 1

    # oamplantfixed: plant O&M cost in M$/year
    try:
        oamplantfixed = float(
            content[[i for i, s in enumerate(content) if 'Surface Plant O&M Cost,' in s][0]].split(',')[1].strip('\n'))
        oamplantfixedprovided = 1
        if oamplantfixed < 0 or oamplantfixed > 100:
            oamplantfixedvalid = 0
        else:
            oamplantfixedvalid = 1
    except:
        oamplantfixedprovided = 0
        oamplantfixedvalid = 0

    # oamplantadjfactor: adj factor for built-in correlation for plant O&M cost
    try:
        oamplantadjfactor = float(
            content[[i for i, s in enumerate(content) if 'Surface Plant O&M Cost Adjustment Factor,' in s][0]].split(
                ',')[1].strip('\n'))
        oamplantadjfactorprovided = 1
        if oamplantadjfactor < 0 or oamplantadjfactor > 10:
            oamplantadjfactorvalid = 0
        else:
            oamplantadjfactorvalid = 1
    except:
        oamplantadjfactorprovided = 0
        oamplantadjfactorvalid = 0

    if oamtotalfixedvalid == 1:
        if oamplantfixedprovided == 1:
            print(
                "Warning: Provided total surface plant O&M cost not considered because valid total annual O&M cost provided.")
        if oamplantadjfactorprovided == 1:
            print(
                "Warning: Provided surface plant O&M cost adjustment factor not considered because valid total annual O&M cost provided.")
    else:
        if oamplantfixedvalid == 1 and oamplantadjfactorvalid == 1:
            print(
                "Warning: Provided surface plant O&M cost adjustment factor not considered because valid total surface plant O&M cost provided.")
        elif oamplantfixedprovided == 0 and oamplantadjfactorprovided == 0:
            oamplantadjfactor = 1
            print(
                "Warning: No valid surface plant O&M cost or adjustment factor provided. GEOPHIRES will assume default built-in surface plant O&M cost correlation with adjustment factor = 1.")
        elif oamplantfixedprovided == 1 and oamplantfixedvalid == 0:
            print(
                "Provided surface plant O&M cost outside of range 0-100. GEOPHIRES will assume default built-in surface plant O&M cost correlation with adjustment factor = 1.")
            oamplantadjfactor = 1
        elif oamplantfixedprovided == 0 and oamplantadjfactorprovided == 1 and oamplantadjfactorvalid == 0:
            print(
                "Provided surface plant O&M cost adjustment factor outside of range 0-10. GEOPHIRES will assume default surface plant O&M cost correlation with adjustment factor = 1.")
            oamplantadjfactor = 1

    # oamwaterfixed: total water cost in M$/year
    try:
        oamwaterfixed = float(
            content[[i for i, s in enumerate(content) if 'Water Cost,' in s][0]].split(',')[1].strip('\n'))
        oamwaterfixedprovided = 1
        if oamwaterfixed < 0 or oamwaterfixed > 100:
            oamwaterfixedvalid = 0
        else:
            oamwaterfixedvalid = 1
    except:
        oamwaterfixedprovided = 0
        oamwaterfixedvalid = 0

    # oamwateradjfactor: adj factor for built-in correlation for water cost
    try:
        oamwateradjfactor = float(
            content[[i for i, s in enumerate(content) if 'Water Cost Adjustment Factor,' in s][0]].split(',')[1].strip(
                '\n'))
        oamwateradjfactorprovided = 1
        if oamwateradjfactor < 0 or oamwateradjfactor > 10:
            oamwateradjfactorvalid = 0
        else:
            oamwateradjfactorvalid = 1
    except:
        oamwateradjfactorprovided = 0
        oamwateradjfactorvalid = 0

    if oamtotalfixedvalid == 1:
        if oamwaterfixedprovided == 1:
            print("Warning: Provided total water cost not considered because valid total annual O&M cost provided.")
        if oamwateradjfactorprovided == 1:
            print(
                "Warning: Provided water cost adjustment factor not considered because valid total annual O&M cost provided.")
    else:
        if oamwaterfixedvalid == 1 and oamwateradjfactorvalid == 1:
            print(
                "Warning: Provided water cost adjustment factor not considered because valid total water cost provided.")
        elif oamwaterfixedprovided == 0 and oamwateradjfactorprovided == 0:
            oamwateradjfactor = 1
            print(
                "Warning: No valid total water cost or adjustment factor provided. GEOPHIRES will assume default built-in water cost correlation with adjustment factor = 1.")
        elif oamwaterfixedprovided == 1 and oamwaterfixedvalid == 0:
            print(
                "Provided total water cost outside of range 0-100. GEOPHIRES will assume default built-in water cost correlation with adjustment factor = 1.")
            oamwateradjfactor = 1
        elif oamwaterfixedprovided == 0 and oamwateradjfactorprovided == 1 and oamwateradjfactorvalid == 0:
            print(
                "Provided water cost adjustment factor outside of range 0-10. GEOPHIRES will assume default water cost correlation with adjustment factor = 1.")
            oamwateradjfactor = 1

    # elecprice: electricity price (in $/kWh) to calculate pumping cost in case of direct-use or additional revenue stream from electricity sales in co-gen option
    if enduseoption in [2, 32, 42, 52]:
        try:
            elecprice = float(
                content[[i for i, s in enumerate(content) if 'Electricity Rate' in s][0]].split(',')[1].strip('\n'))
            if elecprice < 0 or elecprice > 1:
                elecprice = 0.07
                print(
                    "Warning: Provided electricty rate is outside of range 0-1. GEOPHIRES will assume default electricity rate ($0.07/kWh)")
        except:
            elecprice = 0.07
            print(
                "Warning: No valid electricity rate provided. GEOPHIRES will assume default electricity rate ($0.07/kWh)")

    # heatprice: heat price (in $/kWh) to calculate additional revenue stream from heat sales in co-gen option
    if enduseoption in [31, 41, 51]:
        try:
            heatprice = float(
                content[[i for i, s in enumerate(content) if 'Heat Rate' in s][0]].split(',')[1].strip('\n'))
            if heatprice < 0 or heatprice > 1:
                heatprice = 0.02
                print(
                    "Warning: Provided heat rate is outside of range 0-1. GEOPHIRES will assume default heat rate ($0.02/kWh)")
        except:
            heatprice = 0.02
            print("Warning: No valid heat rate provided. GEOPHIRES will assume default heat rate ($0.02/kWh)")

    # printoutput
    # printoutput = 0: do not print output to console
    # printoutput = 1: print output to console (default)
    try:
        printoutput = int(
            content[[i for i, s in enumerate(content) if 'Print Output to Console' in s][0]].split(',')[1].strip('\n'))
        if not (printoutput in [0, 1]):
            printoutput = 1
            print(
                "Warning: Provided print output option should be 0 or 1. GEOPHIRES will assume default print output option (1)")
    except:
        printoutput = 1
        print("Warning: No valid print output option provided. GEOPHIRES will assume default print output option (1)")

    # number of timesteps per year [1/year]
    try:
        timestepsperyear = int(
            content[[i for i, s in enumerate(content) if 'Time steps per year' in s][0]].split(',')[1].strip('\n'))
        if not (timestepsperyear in list(range(1, 101))):
            timestepsperyear = 4
            print(
                "Warning: Provided number of time steps per year outside of range 1-100. GEOPHIRES will assume default number of time steps per year (4)")
    except:
        timestepsperyear = 4
        print(
            "Warning: No valid number of time steps per year provided. GEOPHIRES will assume default number of time steps per year (4)")

    # some additional preprocessing calculations
    # calculate maximum well depth (m)
    intersecttemperature = [1000., 1000., 1000., 1000.]
    if numseg == 1:
        maxdepth = (Tmax - Tsurf) / gradient[0]
    else:
        maxdepth = 0
        intersecttemperature[0] = Tsurf + gradient[0] * layerthickness[0]
        for i in range(1, numseg - 1):
            intersecttemperature[i] = intersecttemperature[i - 1] + gradient[i] * layerthickness[i]
        layerindex = next(loc for loc, val in enumerate(intersecttemperature) if val > Tmax)
        if layerindex > 0:
            for i in range(0, layerindex):
                maxdepth = maxdepth + layerthickness[i]
            maxdepth = maxdepth + (Tmax - intersecttemperature[layerindex - 1]) / gradient[layerindex]
        else:
            maxdepth = (Tmax - Tsurf) / gradient[0]

    if depth > maxdepth:
        depth = maxdepth

    # calculate initial reservoir temperature
    intersecttemperature = [Tsurf] + intersecttemperature
    totaldepth = np.append(np.array([0.]), np.cumsum(layerthickness))
    temperatureindex = max(loc for loc, val in enumerate(depth > totaldepth) if val == True)
    Trock = intersecttemperature[temperatureindex] + gradient[temperatureindex] * (depth - totaldepth[temperatureindex])

    # calculate average geothermal gradient
    if numseg == 1:
        averagegradient = gradient[0]
    else:
        averagegradient = (Trock - Tsurf) / depth

    # specify time-stepping vectors
    timevector = np.linspace(0, plantlifetime, timestepsperyear * plantlifetime + 1)
    Tresoutput = np.zeros(len(timevector))

    # calculate reservoir water properties
    cpwater = heatcapacitywater(
        Tinj * 0.5 + (Trock * 0.9 + Tinj * 0.1) * 0.5)  # J/kg/K (based on TARB in Geophires v1.2)
    rhowater = densitywater(Tinj * 0.5 + (Trock * 0.9 + Tinj * 0.1) * 0.5)

    # temperature gain in injection wells
    Tinj = Tinj + tempgaininj

    # calculate reservoir temperature output (internal or external)
    #   resoption = 1  Multiple parallel fractures model (LANL)
    #   resoption = 2  Volumetric block model (1D linear heat sweep model (Stanford))
    #   resoption = 3  Drawdown parameter model (Tester)
    #   resoption = 4  Thermal drawdown percentage model (GETEM)
    #   resoption = 5  Generic user-provided temperature profile
    #   resoption = 6  Tough2 is called

    if resoption == 1:

        # convert flowrate to volumetric rate
        q = nprod * prodwellflowrate / rhowater  # m^3/s

        # specify Laplace-space function
        fp = lambda s: (1. / s) * exp(-sqrt(s) * tanh(
            (rhowater * cpwater * (q / fracnumb / fracwidth) * (fracsep / 2.) / (2. * krock * fracheight)) * sqrt(s)))

        # calculate non-dimensional time
        td = (rhowater * cpwater) ** 2 / (4 * krock * rhorock * cprock) * (
                    q / fracnumb / fracwidth / fracheight) ** 2 * timevector * 365. * 24. * 3600

        # calculate non-dimensional temperature array
        Twnd = []
        try:
            for t in range(1, len(timevector)):
                Twnd = Twnd + [float(invertlaplace(fp, td[t], method='talbot'))]
        except:
            print(
                "Error: GEOPHIRES could not execute numerical inverse laplace calculation for reservoir model 1. Simulation will abort.")
            sys.exit()

        Twnd = np.asarray(Twnd)

        # calculate dimensional temperature, add initial rock temperature to beginning of array
        Tresoutput = Trock - (Twnd * (Trock - Tinj))
        Tresoutput = np.append([Trock], Tresoutput)

    elif resoption == 2:

        # specify rock properties
        phi = porrock  # porosity [%]
        h = 500.  # heat transfer coefficient [W/m^2 K]
        shape = 0.2  # ratio of conduction path length
        alpha = krock / (rhorock * cprock)

        # storage ratio
        gamma = (rhowater * cpwater * phi) / (rhorock * cprock * (1 - phi))
        # effective rock radius
        r_efr = 0.83 * (0.75 * (fracsep * fracheight * fracwidth) / math.pi) ** (1. / 3.)
        # Biot number
        Bi = h * r_efr / krock
        # effective rock time constant
        tau_efr = r_efr ** 2. * (shape + 1. / Bi) / (3. * alpha)

        # reservoir dimensions and flow properties
        hl = (fracnumb - 1) * fracsep
        wl = fracwidth
        aave = hl * wl
        u0 = nprod * prodwellflowrate / (rhowater * aave)
        tres = (fracheight * phi) / u0

        # number of heat transfer units
        ntu = tres / tau_efr

        # specify Laplace-space function
        fp = lambda s: (1 / s) * (1 - exp(-(1 + ntu / (gamma * (s + ntu))) * s))

        # calculate non-dimensional temperature array
        Twnd = []
        try:
            for t in range(1, len(timevector)):
                Twnd = Twnd + [float(invertlaplace(fp, timevector[t] * 365. * 24. * 3600. / tres, method='talbot'))]
        except:
            print(
                "Error: GEOPHIRES could not execute numerical inverse laplace calculation for reservoir model 2. Simulation will abort.")
            sys.exit()
        Twnd = np.asarray(Twnd)

        # calculate dimensional temperature, add error-handling for non-sensical temperatures
        Tresoutput = Twnd * (Trock - Tinj) + Tinj
        Tresoutput = np.append([Trock], Tresoutput)
        Tresoutput = np.asarray([Trock if x > Trock or x < Tinj else x for x in Tresoutput])

    elif resoption == 3:

        Tresoutput[0] = Trock
        for i in range(1, len(timevector)): Tresoutput[i] = math.erf(
            1. / drawdp / cpwater * math.sqrt(krock * rhorock * cprock / timevector[i] / (365. * 24. * 3600.))) * (
                                                                        Trock - Tinj) + Tinj

    elif resoption == 4:

        Tresoutput = (1 - drawdp * timevector) * (Trock - Tinj) + Tinj  # this is no longer as in thesis (equation 4.16)

    elif resoption == 5:

        Tresoutput[0] = Trock
        restempfname = filenamereservoiroutput
        try:
            with open(restempfname) as f:
                contentprodtemp = f.readlines()
        except:
            print(
                'Error: GEOPHIRES could not read reservoir output file (' + filenamereservoiroutput + ') and will abort simulation.')
            sys.exit()
        numlines = len(contentprodtemp)
        if numlines != plantlifetime * timestepsperyear + 1:
            print('Error: Reservoir output file (' + filenamereservoiroutput + ') does not have required ' + str(
                plantlifetime * timestepsperyear + 1) + ' lines. GEOPHIRES will abort simulation.')
            sys.exit()
        for i in range(0, numlines):
            Tresoutput[i] = float(contentprodtemp[i].split(',')[1].strip('\n'))

    elif resoption == 6:

        # GEOPHIRES assumes TOUGH2 executable and input file are in same directory as GEOPHIRESv2.py

        # create tough2 input file
        path_to_exe = str('xt2_eos1.exe')
        if not os.path.exists(os.path.join(os.getcwd(), path_to_exe)):
            raise OSError(
                'TOUGH2 executable file does not exist in current working directory. GEOPHIRES will abort simulation.')
        if tough2modelfilename == 'Doublet':
            infile = str('Doublet.dat')
            outfile = str('Doublet.out')
            initialtemp = Trock
            rockthermalcond = krock
            rockheatcap = cprock
            rockdensity = rhorock
            rockpor = porrock
            rockperm = permrock
            reservoirthickness = resthickness
            reservoirwidth = reswidth
            wellseperation = wellsep
            DeltaXgrid = wellseperation / 15
            DeltaYgrid = reservoirwidth / 11
            DeltaZgrid = reservoirthickness / 5
            flowrate = prodwellflowrate

            # convert injection temperature to injection enthalpy
            arraytinj = np.array([1.8, 11.4, 23.4, 35.4, 47.4, 59.4, 71.3, 83.3, 95.2, 107.1, 118.9])
            arrayhinj = np.array([1.0E4, 5.0E4, 1.0E5, 1.5E5, 2.0E5, 2.5E5, 3.0E5, 3.5E5, 4.0E5, 4.5E5, 5.0E5])
            injenthalpy = np.interp(Tinj, arraytinj, arrayhinj)
            # write doublet input file
            f = open(infile, 'w')
            f.write('Doublet\n')
            f.write('MESHMAKER1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n')
            f.write('XYZ\n')
            f.write('	0.\n')
            f.write('NX      17 %9.3f\n' % (DeltaXgrid))
            f.write('NY      11 %9.3f\n' % (DeltaYgrid))
            f.write('NZ       5 %9.3f\n' % (DeltaZgrid))
            f.write('\n')
            f.write('\n')
            f.write('ROCKS----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n')
            f.write('POMED    3%10.1f %9.4f %9.2E %9.2E %9.2E %9.4f %9.2f          \n' % (
            rockdensity, rockpor, rockperm, rockperm, rockperm, rockthermalcond, rockheatcap))
            f.write('       0.0       0.0       2.0       0.0       0.0\n')
            f.write('    3            0.3      0.05\n')
            f.write('    8\n')
            f.write('\n')
            f.write('MULTI----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n')
            f.write('    1    2    2    6\n')
            f.write('START----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n')
            f.write('PARAM----1-MOP* 123456789012345678901234----*----5----*----6----*----7----*----8\n')
            f.write(' 8 19999       5000000000001  03 000   0                                        \n')
            f.write('       0.0 %9.3E 5259490.0       0.0                9.81       4.0       1.0\n' % (
                        plantlifetime * 365 * 24 * 3600))
            f.write('    1.0E-5       1.0                 1.0       1.0          \n')
            f.write('           1000000.0          %10.1f\n' % (initialtemp))
            f.write('                                                                                \n')
            f.write('SOLVR----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n')
            f.write('3  Z1   O0       0.1    1.0E-6\n')
            f.write('\n')
            f.write('\n')
            f.write('GENER----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n')
            f.write('A36 2  012                   0     COM1  %9.3f %9.1f          \n' % (flowrate, injenthalpy))
            f.write('A3616  021                   0     MASS  %9.3f             \n' % (-flowrate))
            f.write('\n')
            f.write('INCON----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n')
            f.write('\n')
            f.write('FOFT ----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n')
            f.write('A36 2\n')
            f.write('A3616\n')
            f.write('\n')
            f.write('GOFT ----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n')
            f.write('A36 2  012\n')
            f.write('A3616  021\n')
            f.write('\n')
            f.write('ENDCY\n')
            f.close()
            print("GEOPHIRES will run TOUGH2 simulation with built-in Doublet model ...")

        else:
            infile = tough2modelfilename
            outfile = str('tough2output.out')
            print(
                "GEOPHIRES will run TOUGH2 simulation with user-provided input file = " + tough2modelfilename + " ...")

        # run TOUGH2 executable
        try:
            os.system('%s < %s > %s' % (path_to_exe, infile, outfile))
        except:
            print("Error: GEOPHIRES could not run TOUGH2 and will abort simulation.")
            sys.exit()

        # read output temperature and pressure
        try:
            fname = 'FOFT'
            with open(fname) as f:
                content = f.readlines()

            NumerOfResults = len(content)
            SimTimes = np.zeros(NumerOfResults)
            ProdPressure = np.zeros(NumerOfResults)
            ProdTemperature = np.zeros(NumerOfResults)
            for i in range(0, NumerOfResults):
                SimTimes[i] = float(content[i].split(',')[1].strip('\n'))
                ProdPressure[i] = float(content[i].split(',')[8].strip('\n'))
                ProdTemperature[i] = float(content[i].split(',')[9].strip('\n'))

            # print(ProdTemperature)
            Tresoutput = np.interp(timevector * 365 * 24 * 3600, SimTimes, ProdTemperature)
        except:
            print(
                "Error: GEOPHIRES could not import production temperature and pressure from TOUGH2 output file (" + fname + ") and will abort simulation.")

        # define function to extract temperature profile from outfile (move up to top of script?)

    # calculate wellbore temperature drop
    ProdTempDrop = 0
    if rameyoptionprod == 0:
        ProdTempDrop = tempdropprod
    elif rameyoptionprod == 1:
        alpharock = krock / (rhorock * cprock)
        framey = np.zeros(len(timevector))
        framey[1:] = -np.log(1.1 * (prodwelldiam / 2.) / np.sqrt(
            4. * alpharock * timevector[1:] * 365. * 24. * 3600. * utilfactor)) - 0.29
        framey[0] = -np.log(1.1 * (prodwelldiam / 2.) / np.sqrt(4. * alpharock * timevector[
            1] * 365. * 24. * 3600. * utilfactor)) - 0.29  # assume outside diameter of casing is 10% larger than inside diameter of production pipe (=prodwelldiam)
        # assume borehole thermal resistance negligible to rock thermal resistance
        rameyA = prodwellflowrate * cpwater * framey / 2 / math.pi / krock
        # this code is only valid so far for 1 gradient and deviation = 0 !!!!!!!!   For multiple gradients, use Ramey's model for every layer

        ProdTempDrop = -((Trock - Tresoutput) - averagegradient * (depth - rameyA) + (
                    Tresoutput - averagegradient * rameyA - Trock) * np.exp(-depth / rameyA))

    ProducedTemperature = Tresoutput - ProdTempDrop

    # redrilling
    redrill = 0
    if resoption < 5:  # only applies to the built-in analytical reservoir models
        indexfirstmaxdrawdown = np.argmax(ProducedTemperature < (1 - maxdrawdown) * ProducedTemperature[0])
        if indexfirstmaxdrawdown > 0:  # redrilling necessary
            redrill = int(np.floor(len(ProducedTemperature) / indexfirstmaxdrawdown))
            ProducedTemperatureRepeatead = np.tile(ProducedTemperature[0:indexfirstmaxdrawdown], redrill + 1)
            ProducedTemperature = ProducedTemperatureRepeatead[0:len(ProducedTemperature)]

    # ------------------------------------------
    # calculate pressure drops and pumping power
    # ------------------------------------------
    # production wellbore fluid conditions [kPa]
    Tprodaverage = Tresoutput - ProdTempDrop / 4.  # most of temperature drop happens in upper section (because surrounding rock temperature is lowest in upper section)
    rhowaterprod = densitywater(Tprodaverage)  # replace with correlation based on Tprodaverage
    muwaterprod = viscositywater(Tprodaverage)  # replace with correlation based on Tprodaverage
    vprod = prodwellflowrate / rhowaterprod / (math.pi / 4. * prodwelldiam ** 2)
    Rewaterprod = 4. * prodwellflowrate / (muwaterprod * math.pi * prodwelldiam)  # laminar or turbulent flow?
    Rewaterprodaverage = np.average(Rewaterprod)
    if Rewaterprodaverage < 2300.:
        f3 = 64. / Rewaterprod
    else:
        relroughness = 1E-4 / prodwelldiam
        f3 = 1. / np.power(-2 * np.log10(relroughness / 3.7 + 5.74 / np.power(Rewaterprod, 0.9)), 2.)
        f3 = 1. / np.power((-2 * np.log10(relroughness / 3.7 + 2.51 / Rewaterprod / np.sqrt(f3))), 2.)
        f3 = 1. / np.power((-2 * np.log10(relroughness / 3.7 + 2.51 / Rewaterprod / np.sqrt(f3))), 2.)
        f3 = 1. / np.power((-2 * np.log10(relroughness / 3.7 + 2.51 / Rewaterprod / np.sqrt(f3))), 2.)
        f3 = 1. / np.power((-2 * np.log10(relroughness / 3.7 + 2.51 / Rewaterprod / np.sqrt(f3))), 2.)
        f3 = 1. / np.power((-2 * np.log10(relroughness / 3.7 + 2.51 / Rewaterprod / np.sqrt(f3))),
                           2.)  # 6 iterations to converge

    # injection well conditions
    Tinjaverage = Tinj
    rhowaterinj = densitywater(Tinjaverage) * np.linspace(1, 1, len(ProducedTemperature))
    muwaterinj = viscositywater(Tinjaverage) * np.linspace(1, 1,
                                                           len(ProducedTemperature))  # replace with correlation based on Tinjaverage
    vinj = nprod / ninj * prodwellflowrate * (1. + waterloss) / rhowaterinj / (math.pi / 4. * injwelldiam ** 2)
    Rewaterinj = 4. * nprod / ninj * prodwellflowrate * (1. + waterloss) / (
                muwaterinj * math.pi * injwelldiam)  # laminar or turbulent flow?
    Rewaterinjaverage = np.average(Rewaterinj)
    if Rewaterinjaverage < 2300.:  # laminar flow
        f1 = 64. / Rewaterinj
    else:  # turbulent flow
        relroughness = 1E-4 / injwelldiam
        f1 = 1. / np.power(-2 * np.log10(relroughness / 3.7 + 5.74 / np.power(Rewaterinj, 0.9)), 2.)
        f1 = 1. / np.power((-2 * np.log10(relroughness / 3.7 + 2.51 / Rewaterinj / np.sqrt(f1))), 2.)
        f1 = 1. / np.power((-2 * np.log10(relroughness / 3.7 + 2.51 / Rewaterinj / np.sqrt(f1))), 2.)
        f1 = 1. / np.power((-2 * np.log10(relroughness / 3.7 + 2.51 / Rewaterinj / np.sqrt(f1))), 2.)
        f1 = 1. / np.power((-2 * np.log10(relroughness / 3.7 + 2.51 / Rewaterinj / np.sqrt(f1))), 2.)
        f1 = 1. / np.power((-2 * np.log10(relroughness / 3.7 + 2.51 / Rewaterinj / np.sqrt(f1))),
                           2.)  # 6 iterations to converge

    if impedancemodelused == 1:  # assumed everything stays liquid throughout
        # injecion well pressure drop [kPa]
        DP1 = f1 * (rhowaterinj * vinj ** 2 / 2) * (depth / injwelldiam) / 1E3  # /1E3 to convert from Pa to kPa

        # reservoir pressure drop [kPa]
        rhowaterreservoir = densitywater(0.1 * Tinj + 0.9 * Tresoutput)  # based on TARB in Geophires v1.2
        DP2 = impedance * nprod * prodwellflowrate * 1000. / rhowaterreservoir

        # production well pressure drop [kPa]
        DP3 = f3 * (rhowaterprod * vprod ** 2 / 2.) * (depth / prodwelldiam) / 1E3  # /1E3 to convert from Pa to kPa

        # buoyancy pressure drop [kPa]
        DP4 = (rhowaterprod - rhowaterinj) * depth * 9.81 / 1E3  # /1E3 to convert from Pa to kPa

        # overall pressure drop
        DP = DP1 + DP2 + DP3 + DP4

        # calculate pumping power [MWe] (approximate)
        PumpingPower = DP * nprod * prodwellflowrate * (1 + waterloss) / rhowaterinj / pumpeff / 1E3

        # in GEOPHIRES v1.2, negative pumping power values become zero (b/c we are not generating electricity)
        PumpingPower = [0. if x < 0. else x for x in PumpingPower]

    else:  # PI and II are used
        # reservoir hydrostatic pressure [kPa]
        if usebuiltinhydrostaticpressurecorrelation == 1:
            CP = 4.64E-7
            CT = 9E-4 / (30.796 * Trock ** (-0.552))
            Phydrostatic = 0 + 1. / CP * (math.exp(
                densitywater(Tsurf) * 9.81 * CP / 1000 * (depth - CT / 2 * averagegradient * depth ** 2)) - 1)

        if productionwellpumping == 1:
            Pexcess = 344.7  # [kPa] = 50 psi. Excess pressure covers non-condensable gas pressure and net positive suction head for the pump
            Pminimum = vaporpressurewater(
                Trock) + Pexcess  # [kPa] is minimum production pump inlet pressure and minimum wellhead pressure
            if usebuiltinppwellheadcorrelation == 1:
                # production wellhead pressure [kPa]
                Pprodwellhead = Pminimum
            else:
                Pprodwellhead = ppwellhead
                if Pprodwellhead < Pminimum:
                    Pprodwellhead = Pminimum
                    print(
                        "Warning: provided production wellhead pressure under minimum pressure. GEOPHIRES will assume minimum wellhead pressure")

            PIkPa = PI / 100  # convert PI from kg/s/bar to kg/s/kPa

            # calculate pumping depth
            pumpdepth = depth + (Pminimum - Phydrostatic + prodwellflowrate / PIkPa) / (
                        f3 * (rhowaterprod * vprod ** 2 / 2.) * (1 / prodwelldiam) / 1E3 + rhowaterprod * 9.81 / 1E3)
            pumpdepthfinal = np.max(pumpdepth)
            if pumpdepthfinal < 0:
                pumpdepthfinal = 0
                print(
                    "Warning: GEOPHIRES calculates negative production well pumping depth. No production well pumps will be assumed")
            elif pumpdepthfinal > 600:
                print(
                    "Warning: GEOPHIRES calculates pump depth to be deeper than 600 m. Verify reservoir pressure, production well flow rate and production well dimensions")

            # calculate production well pumping pressure [kPa]
            DP3 = Pprodwellhead - (Phydrostatic - prodwellflowrate / PIkPa - rhowaterprod * 9.81 * depth / 1E3 - f3 * (
                        rhowaterprod * vprod ** 2 / 2.) * (depth / prodwelldiam) / 1E3)
            # DP3 = [0 if x<0 else x for x in DP3] #set negative values to 0
            PumpingPowerProd = DP3 * nprod * prodwellflowrate / rhowaterprod / pumpeff / 1E3  # [MWe] total pumping power for production wells
            PumpingPowerProd = np.array([0. if x < 0. else x for x in PumpingPowerProd])

        IIkPa = II / 100  # convert II from kg/s/bar to kg/s/kPa

        # necessary injection wellhead pressure [kPa]
        Pinjwellhead = Phydrostatic + prodwellflowrate * (
                    1 + waterloss) * nprod / ninj / IIkPa - rhowaterinj * 9.81 * depth / 1E3 + f1 * (
                                   rhowaterinj * vinj ** 2 / 2) * (depth / injwelldiam) / 1E3

        # plant outlet pressure [kPa]
        if usebuiltinoutletplantcorrelation == 1:
            DPSurfaceplant = 68.95  # [kPa] assumes 10 psi pressure drop in surface equipment
            Pplantoutlet = Pprodwellhead - DPSurfaceplant

        # injection pump pressure [kPa]
        DP1 = Pinjwellhead - Pplantoutlet
        # DP1 = [0 if x<0 else x for x in DP1] #set negative values to 0
        PumpingPowerInj = DP1 * nprod * prodwellflowrate * (
                    1 + waterloss) / rhowaterinj / pumpeff / 1E3  # [MWe] total pumping power for injection wells
        PumpingPowerInj = np.array([0. if x < 0. else x for x in PumpingPowerInj])

        # total pumping power
        if productionwellpumping == 1:
            PumpingPower = PumpingPowerInj + PumpingPowerProd
        else:
            PumpingPower = PumpingPowerInj

        # negative pumping power values become zero (b/c we are not generating electricity)
        PumpingPower = [0. if x < 0. else x for x in PumpingPower]

    # ----------------------------------------------
    # calculate produced electricity/direct-use heat
    # ----------------------------------------------
    if enduseoption == 2:  # direct-use
        HeatExtracted = nprod * prodwellflowrate * cpwater * (
                    ProducedTemperature - Tinj) / 1E6  # heat extracted from geofluid [MWth]
        HeatProduced = HeatExtracted * enduseefficiencyfactor  # useful direct-use heat provided to application [MWth]
    else:
        if (math.floor(enduseoption / 10) == 4):
            TenteringPP = Tchpbottom
        else:
            TenteringPP = ProducedTemperature
        # Availability water (copied from GEOPHIRES v1.0 Fortran Code)
        A = 4.041650
        B = -1.204E-2
        C = 1.60500E-5
        T0 = Tenv + 273.15
        T1 = TenteringPP + 273.15
        T2 = Tenv + 273.15
        Availability = ((A - B * T0) * (T1 - T2) + (B - C * T0) / 2.0 * (T1 ** 2 - T2 ** 2) + C / 3.0 * (
                    T1 ** 3 - T2 ** 3) - A * T0 * np.log(T1 / T2)) * 2.2046 / 947.83  # MJ/kg

        if pptype == 1:  # Subcritical ORC
            if (Tenv < 15.):
                C1 = 2.746E-3
                C0 = -8.3806E-2
                D1 = 2.713E-3
                D0 = -9.1841E-2
                Tfraction = (Tenv - 5.) / 10.
            else:
                C1 = 2.713E-3
                C0 = -9.1841E-2
                D1 = 2.676E-3
                D0 = -1.012E-1
                Tfraction = (Tenv - 15.) / 10.
            etaull = C1 * TenteringPP + C0
            etauul = D1 * TenteringPP + D0
            etau = (1 - Tfraction) * etaull + Tfraction * etauul
            if (Tenv < 15.):
                C1 = 0.0894
                C0 = 55.6
                D1 = 0.0894
                D0 = 62.6
                Tfraction = (Tenv - 5.) / 10.
            else:
                C1 = 0.0894
                C0 = 62.6
                D1 = 0.0894
                D0 = 69.6
                Tfraction = (Tenv - 15.) / 10.
            reinjtll = C1 * TenteringPP + C0
            reinjtul = D1 * TenteringPP + D0
            ReinjTemp = (1. - Tfraction) * reinjtll + Tfraction * reinjtul
        elif pptype == 2:  # Supercritical ORC
            if (Tenv < 15.):
                C2 = -1.55E-5
                C1 = 7.604E-3
                C0 = -3.78E-1
                D2 = -1.499E-5
                D1 = 7.4268E-3
                D0 = -3.7915E-1
                Tfraction = (Tenv - 5.) / 10.
            else:
                C2 = -1.499E-5
                C1 = 7.4268E-3
                C0 = -3.7915E-1
                D2 = -1.55E-5
                D1 = 7.55136E-3
                D0 = -4.041E-1
                Tfraction = (Tenv - 15.) / 10.
            etaull = C2 * TenteringPP ** 2 + C1 * TenteringPP + C0
            etauul = D2 * TenteringPP ** 2 + D1 * TenteringPP + D0
            etau = (1 - Tfraction) * etaull + Tfraction * etauul
            if (Tenv < 15.):
                C1 = 0.02
                C0 = 49.26
                D1 = 0.02
                D0 = 56.26
                Tfraction = (Tenv - 5.) / 10.
            else:
                C1 = 0.02
                C0 = 56.26
                D1 = 0.02
                D0 = 63.26
                Tfraction = (Tenv - 15.) / 10.
            reinjtll = C1 * TenteringPP + C0
            reinjtul = D1 * TenteringPP + D0
            ReinjTemp = (1. - Tfraction) * reinjtll + Tfraction * reinjtul
        elif pptype == 3:  # single-flash
            if (Tenv < 15.):
                C2 = -4.27318E-7
                C1 = 8.65629E-4
                C0 = 1.78931E-1
                D2 = -5.85412E-7
                D1 = 9.68352E-4
                D0 = 1.58056E-1
                Tfraction = (Tenv - 5.) / 10.
            else:
                C2 = -5.85412E-7
                C1 = 9.68352E-4
                C0 = 1.58056E-1
                D2 = -7.78996E-7
                D1 = 1.09230E-3
                D0 = 1.33708E-1
                Tfraction = (Tenv - 15.) / 10.
            etaull = C2 * TenteringPP ** 2 + C1 * TenteringPP + C0
            etauul = D2 * TenteringPP ** 2 + D1 * TenteringPP + D0
            etau = (1. - Tfraction) * etaull + Tfraction * etauul
            if (Tenv < 15.):
                C2 = -1.11519E-3
                C1 = 7.79126E-1
                C0 = -10.2242
                D2 = -1.10232E-3
                D1 = 7.83893E-1
                D0 = -5.17039
                Tfraction = (Tenv - 5.) / 10.
            else:
                C2 = -1.10232E-3
                C1 = 7.83893E-1
                C0 = -5.17039
                D2 = -1.08914E-3
                D1 = 7.88562E-1
                D0 = -1.89707E-1
                Tfraction = (Tenv - 15.) / 10.
            reinjtll = C2 * TenteringPP ** 2 + C1 * TenteringPP + C0
            reinjtul = D2 * TenteringPP ** 2 + D1 * TenteringPP + D0
            ReinjTemp = (1. - Tfraction) * reinjtll + Tfraction * reinjtul
        elif pptype == 4:  # double-flash
            if (Tenv < 15.):
                C2 = -1.200E-6
                C1 = 1.22731E-3
                C0 = 2.26956E-1
                D2 = -1.42165E-6
                D1 = 1.37050E-3
                D0 = 1.99847E-1
                Tfraction = (Tenv - 5.) / 10.
            else:
                C2 = -1.42165E-6
                C1 = 1.37050E-3
                C0 = 1.99847E-1
                D2 = -1.66771E-6
                D1 = 1.53079E-3
                D0 = 1.69439E-1
                Tfraction = (Tenv - 15.) / 10.
            etaull = C2 * TenteringPP ** 2 + C1 * TenteringPP + C0
            etauul = D2 * TenteringPP ** 2 + D1 * TenteringPP + D0
            etau = (1. - Tfraction) * etaull + Tfraction * etauul
            if (Tenv < 15.):
                C2 = -7.70928E-4
                C1 = 5.02466E-1
                C0 = 5.22091
                D2 = -7.69455E-4
                D1 = 5.09406E-1
                D0 = 11.6859
                Tfraction = (Tenv - 5.) / 10.
            else:
                C2 = -7.69455E-4
                C1 = 5.09406E-1
                C0 = 11.6859
                D2 = -7.67751E-4
                D1 = 5.16356E-1
                D0 = 18.0798
                Tfraction = (Tenv - 15.) / 10.
            reinjtll = C2 * TenteringPP ** 2 + C1 * TenteringPP + C0
            reinjtul = D2 * TenteringPP ** 2 + D1 * TenteringPP + D0
            ReinjTemp = (1. - Tfraction) * reinjtll + Tfraction * reinjtul

        # check if reinjectemp (model calculated) >= Tinj (user provided)
        if enduseoption == 1:  # pure electricity
            if np.min(ReinjTemp) < Tinj:
                Tinj = np.min(ReinjTemp)
                print("Warning: injection temperature lowered")
        elif (math.floor(enduseoption / 10) == 3):  # enduseoption = 3: cogen topping cycle
            if np.min(ReinjTemp) < Tinj:
                Tinj = np.min(ReinjTemp)
                print("Warning: injection temperature lowered")
        elif (math.floor(enduseoption / 10) == 4):  # enduseoption = 4: cogen bottoming cycle
            if np.min(ReinjTemp) < Tinj:
                Tinj = np.min(ReinjTemp)
                print("Warning: injection temperature lowered")
        elif (math.floor(enduseoption / 10) == 5):  # enduseoption = 5: cogen split of mass flow rate
            if np.min(ReinjTemp) < Tinj:
                # Tinj = np.min(ReinjTemp)
                print("Warning: injection temperature incorrect but cannot be lowered")
                # chpfraction*Tinj+(1-chpfraction)

        # calculate electricity/heat
        if enduseoption == 1:  # pure electricity
            ElectricityProduced = Availability * etau * nprod * prodwellflowrate
            HeatExtracted = nprod * prodwellflowrate * cpwater * (
                        ProducedTemperature - Tinj) / 1E6  # Heat extracted from geofluid [MWth]
            HeatExtractedTowardsElectricity = HeatExtracted
        elif (math.floor(enduseoption / 10) == 3):  # enduseoption = 3: cogen topping cycle
            ElectricityProduced = Availability * etau * nprod * prodwellflowrate
            HeatExtracted = nprod * prodwellflowrate * cpwater * (
                        ProducedTemperature - Tinj) / 1E6  # Heat extracted from geofluid [MWth]
            HeatProduced = enduseefficiencyfactor * nprod * prodwellflowrate * cpwater * (
                        ReinjTemp - Tinj) / 1E6  # Useful heat for direct-use application [MWth]
            HeatExtractedTowardsElectricity = nprod * prodwellflowrate * cpwater * (
                        ProducedTemperature - ReinjTemp) / 1E6
        elif (math.floor(enduseoption / 10) == 4):  # enduseoption = 4: cogen bottoming cycle
            ElectricityProduced = Availability * etau * nprod * prodwellflowrate
            HeatExtracted = nprod * prodwellflowrate * cpwater * (
                        ProducedTemperature - Tinj) / 1E6  # Heat extracted from geofluid [MWth]
            HeatProduced = enduseefficiencyfactor * nprod * prodwellflowrate * cpwater * (
                        ProducedTemperature - Tchpbottom) / 1E6  # Useful heat for direct-use application [MWth]
            HeatExtractedTowardsElectricity = nprod * prodwellflowrate * cpwater * (Tchpbottom - Tinj) / 1E6
        elif (math.floor(enduseoption / 10) == 5):  # enduseoption = 5: cogen split of mass flow rate
            ElectricityProduced = Availability * etau * nprod * prodwellflowrate * (
                        1. - chpfraction)  # electricity part [MWe]
            HeatExtracted = nprod * prodwellflowrate * cpwater * (
                        ProducedTemperature - Tinj) / 1E6  # Total amount of heat extracted from geofluid [MWth]
            HeatProduced = enduseefficiencyfactor * chpfraction * nprod * prodwellflowrate * cpwater * (
                        ProducedTemperature - Tinj) / 1E6  # useful heat part for direct-use application [MWth]
            HeatExtractedTowardsElectricity = (1. - chpfraction) * nprod * prodwellflowrate * cpwater * (
                        ProducedTemperature - Tinj) / 1E6

        # subtract pumping power for net electricity and  calculate first law efficiency
        if enduseoption == 1 or enduseoption > 2:
            NetElectricityProduced = ElectricityProduced - PumpingPower
            FirstLawEfficiency = NetElectricityProduced / HeatExtractedTowardsElectricity

    # -------------
    # capital costs
    # -------------
    # well costs (using GeoVision drilling correlations). These are calculated whether or not totalcapcostvalid = 1
    if ccwellfixedvalid == 1:
        C1well = ccwellfixed
        Cwell = C1well * (nprod + ninj)
    else:
        if wellcorrelation == 1:  # vertical open-hole, small diameter
            C1well = (
                                 0.3021 * depth ** 2 + 584.9112 * depth + 751368.) * 1E-6  # well drilling and completion cost in M$/well
        elif wellcorrelation == 2:  # deviated liner, small diameter
            C1well = (0.2898 * depth ** 2 + 822.1507 * depth + 680563.) * 1E-6
        elif wellcorrelation == 3:  # vertical open-hole, large diameter
            C1well = (0.2818 * depth ** 2 + 1275.5213 * depth + 632315.) * 1E-6
        elif wellcorrelation == 4:  # deviated liner, large diameter
            C1well = (0.2553 * depth ** 2 + 1716.7157 * depth + 500867.) * 1E-6
        if depth < 500.:
            print("Warning: drilling cost correlation extrapolated for drilling depth < 500 m")
        if depth > 7000.:
            print("Warning: drilling cost correlation extrapolated for drilling depth > 7000 m")
        C1well = ccwelladjfactor * C1well
        Cwell = 1.05 * C1well * (nprod + ninj)  # 1.05 for 5% indirect costs

    # reservoir stimulation costs (M$/injection well). These are calculated whether or not totalcapcostvalid = 1
    if ccstimfixedvalid == 1:
        Cstim = ccstimfixed
    else:
        Cstim = 1.05 * 1.15 * ccstimadjfactor * ninj * 1.25  # 1.15 for 15% contingency and 1.05 for 5% indirect costs

    # field gathering system costs (M$)
    if ccgathfixedvalid == 1:
        Cgath = ccgathfixed
    else:
        # Cgath = ccgathadjfactor*50-6*np.max(HeatExtracted)*1000. (GEOPHIRES v1 correlation)

        if impedancemodelused == 1:
            pumphp = np.max(PumpingPower) * 1341
            numberofpumps = np.ceil(pumphp / 2000)  # pump can be maximum 2,000 hp
            if numberofpumps == 0:
                Cpumps = 0
            else:
                pumphpcorrected = pumphp / numberofpumps
                Cpumps = numberofpumps * 1.5 * ((1750 * (pumphpcorrected) ** 0.7) * 3 * (pumphpcorrected) ** (-0.11))
        else:
            if productionwellpumping == 1:
                prodpumphp = np.max(PumpingPowerProd) / nprod * 1341
                Cpumpsprod = nprod * 1.5 * (1750 * (prodpumphp) ** 0.7 + 5750 * (prodpumphp) ** 0.2 + 10000 + np.max(
                    pumpdepth) * 50 * 3.281)  # see page 46 in user's manual asusming rental of rig for 1 day.
            else:
                Cpumpsprod = 0

            injpumphp = np.max(PumpingPowerInj) * 1341
            numberofinjpumps = np.ceil(injpumphp / 2000)  # pump can be maximum 2,000 hp
            if numberofinjpumps == 0:
                Cpumpsinj = 0
            else:
                injpumphpcorrected = injpumphp / numberofinjpumps
                Cpumpsinj = numberofinjpumps * 1.5 * (1750 * (injpumphpcorrected) ** 0.7) * 3 * (
                    injpumphpcorrected) ** (-0.11)

            Cpumps = Cpumpsinj + Cpumpsprod

        Cgath = 1.15 * ccgathadjfactor * 1.12 * ((
                                                             nprod + ninj) * 750 * 500. + Cpumps) / 1E6  # Based on GETEM 2016 #1.15 for 15% contingency and 1.12 for 12% indirect costs

        # plant costs
    if enduseoption == 2:  # direct-use
        if ccplantfixedvalid == 1:
            Cplant = ccplantfixed
        else:
            Cplant = 1.12 * 1.15 * ccplantadjfactor * 250E-6 * np.max(
                HeatExtracted) * 1000.  # 1.15 for 15% contingency and 1.12 for 12% indirect costs
    else:  # all other options have power plant
        if pptype == 1:  # sub-critical ORC
            MaxProducedTemperature = np.max(TenteringPP)
            if (MaxProducedTemperature < 150.):
                C3 = -1.458333E-3
                C2 = 7.6875E-1
                C1 = -1.347917E2
                C0 = 1.0075E4
                CCAPP1 = C3 * MaxProducedTemperature ** 3 + C2 * MaxProducedTemperature ** 2 + C1 * MaxProducedTemperature + C0
            else:
                CCAPP1 = 2231 - 2 * (MaxProducedTemperature - 150.)
            Cplantcorrelation = CCAPP1 * math.pow(np.max(ElectricityProduced) / 15., -0.06) * np.max(
                ElectricityProduced) * 1000. / 1E6

        elif pptype == 2:  # supercritical ORC
            MaxProducedTemperature = np.max(TenteringPP)
            if (MaxProducedTemperature < 150.):
                C3 = -1.458333E-3
                C2 = 7.6875E-1
                C1 = -1.347917E2
                C0 = 1.0075E4
                CCAPP1 = C3 * MaxProducedTemperature ** 3 + C2 * MaxProducedTemperature ** 2 + C1 * MaxProducedTemperature + C0
            else:
                CCAPP1 = 2231 - 2 * (MaxProducedTemperature - 150.)
            Cplantcorrelation = 1.1 * CCAPP1 * math.pow(np.max(ElectricityProduced) / 15., -0.06) * np.max(
                ElectricityProduced) * 1000. / 1E6  # factor 1.1 to make supercritical 10% more expansive than subcritical

        elif pptype == 3:  # single-flash
            if (np.max(ElectricityProduced) < 10.):
                C2 = 4.8472E-2
                C1 = -35.2186
                C0 = 8.4474E3
                D2 = 4.0604E-2
                D1 = -29.3817
                D0 = 6.9911E3
                PLL = 5.
                PRL = 10.
            elif (np.max(ElectricityProduced) < 25.):
                C2 = 4.0604E-2
                C1 = -29.3817
                C0 = 6.9911E3
                D2 = 3.2773E-2
                D1 = -23.5519
                D0 = 5.5263E3
                PLL = 10.
                PRL = 25.
            elif (np.max(ElectricityProduced) < 50.):
                C2 = 3.2773E-2
                C1 = -23.5519
                C0 = 5.5263E3
                D2 = 3.4716E-2
                D1 = -23.8139
                D0 = 5.1787E3
                PLL = 25.
                PRL = 50.
            elif (np.max(ElectricityProduced) < 75.):
                C2 = 3.4716E-2
                C1 = -23.8139
                C0 = 5.1787E3
                D2 = 3.5271E-2
                D1 = -24.3962
                D0 = 5.1972E3
                PLL = 50.
                PRL = 75.
            else:
                C2 = 3.5271E-2
                C1 = -24.3962
                C0 = 5.1972E3
                D2 = 3.3908E-2
                D1 = -23.4890
                D0 = 5.0238E3
                PLL = 75.
                PRL = 100.
            maxProdTemp = np.max(TenteringPP)
            CCAPPLL = C2 * maxProdTemp ** 2 + C1 * maxProdTemp + C0
            CCAPPRL = D2 * maxProdTemp ** 2 + D1 * maxProdTemp + D0
            b = math.log(CCAPPRL / CCAPPLL) / math.log(PRL / PLL)
            a = CCAPPRL / PRL ** b
            Cplantcorrelation = 0.8 * a * math.pow(np.max(ElectricityProduced), b) * np.max(
                ElectricityProduced) * 1000. / 1E6  # factor 0.75 to make double flash 25% more expansive than single flash

        elif pptype == 4:  # double-flash
            if (np.max(ElectricityProduced) < 10.):
                C2 = 4.8472E-2
                C1 = -35.2186
                C0 = 8.4474E3
                D2 = 4.0604E-2
                D1 = -29.3817
                D0 = 6.9911E3
                PLL = 5.
                PRL = 10.
            elif (np.max(ElectricityProduced) < 25.):
                C2 = 4.0604E-2
                C1 = -29.3817
                C0 = 6.9911E3
                D2 = 3.2773E-2
                D1 = -23.5519
                D0 = 5.5263E3
                PLL = 10.
                PRL = 25.
            elif (np.max(ElectricityProduced) < 50.):
                C2 = 3.2773E-2
                C1 = -23.5519
                C0 = 5.5263E3
                D2 = 3.4716E-2
                D1 = -23.8139
                D0 = 5.1787E3
                PLL = 25.
                PRL = 50.
            elif (np.max(ElectricityProduced) < 75.):
                C2 = 3.4716E-2
                C1 = -23.8139
                C0 = 5.1787E3
                D2 = 3.5271E-2
                D1 = -24.3962
                D0 = 5.1972E3
                PLL = 50.
                PRL = 75.
            else:
                C2 = 3.5271E-2
                C1 = -24.3962
                C0 = 5.1972E3
                D2 = 3.3908E-2
                D1 = -23.4890
                D0 = 5.0238E3
                PLL = 75.
                PRL = 100.
            maxProdTemp = np.max(TenteringPP)
            CCAPPLL = C2 * maxProdTemp ** 2 + C1 * maxProdTemp + C0
            CCAPPRL = D2 * maxProdTemp ** 2 + D1 * maxProdTemp + D0
            b = math.log(CCAPPRL / CCAPPLL) / math.log(PRL / PLL)
            a = CCAPPRL / PRL ** b
            Cplantcorrelation = a * math.pow(np.max(ElectricityProduced), b) * np.max(ElectricityProduced) * 1000. / 1E6

        if ccplantfixedvalid == 1:
            Cplant = ccplantfixed
        else:
            Cplant = 1.12 * 1.15 * ccplantadjfactor * Cplantcorrelation * 1.02  # 1.02 to convert cost from 2012 to 2016 #factor 1.15 for 15% contingency and 1.12 for 12% indirect costs.

    # add direct-use plant cost of co-gen system to Cplant (only of no total ccplant was provided)
    if ccplantfixedvalid == 0:  # 1.15 below for contingency and 1.12 for indirect costs
        if (math.floor(enduseoption / 10) == 3):  # enduseoption = 3: cogen topping cycle
            Cplant = Cplant + 1.12 * 1.15 * ccplantadjfactor * 250E-6 * np.max(
                HeatProduced / enduseefficiencyfactor) * 1000.
        elif (math.floor(enduseoption / 10) == 4):  # enduseoption = 4: cogen bottoming cycle
            Cplant = Cplant + 1.12 * 1.15 * ccplantadjfactor * 250E-6 * np.max(
                HeatProduced / enduseefficiencyfactor) * 1000.
        elif (math.floor(enduseoption / 10) == 5):  # enduseoption = 5: cogen parallel cycle
            Cplant = Cplant + 1.12 * 1.15 * ccplantadjfactor * 250E-6 * np.max(
                HeatProduced / enduseefficiencyfactor) * 1000.

    if totalcapcostvalid == 0:
        # exploration costs (same as in Geophires v1.2) (M$)
        if ccexplfixedvalid == 1:
            Cexpl = ccexplfixed
        else:
            Cexpl = 1.15 * ccexpladjfactor * 1.12 * (
                        1. + C1well * 0.6)  # 1.15 for 15% contingency and 1.12 for 12% indirect costs

        # Surface Piping Length Costs (M$) #assumed $750k/km
        Cpiping = 750 / 1000 * pipinglength

        Ccap = Cexpl + Cwell + Cstim + Cgath + Cplant + Cpiping
    else:
        Ccap = totalcapcost

    # ---------
    # O&M costs
    # ---------
    if oamtotalfixedvalid == 0:
        # labor cost
        if enduseoption == 1:  # electricity
            if np.max(ElectricityProduced) < 2.5:
                Claborcorrelation = 236. / 1E3  # M$/year
            else:
                Claborcorrelation = (589. * math.log(np.max(ElectricityProduced)) - 304.) / 1E3  # M$/year
        else:
            if np.max(HeatExtracted) < 2.5 * 5.:
                Claborcorrelation = 236. / 1E3  # M$/year
            else:
                Claborcorrelation = (589. * math.log(np.max(HeatExtracted) / 5.) - 304.) / 1E3  # M$/year
        Claborcorrelation = Claborcorrelation * 1.1  # 1.1 to convert from 2012 to 2016$ with BLS employment cost index (for utilities in March)

        # plant O&M cost
        if oamplantfixedvalid == 1:
            Coamplant = oamplantfixed
        else:
            Coamplant = oamplantadjfactor * (1.5 / 100. * Cplant + 0.75 * Claborcorrelation)

        # wellfield O&M cost
        if oamwellfixedvalid == 1:
            Coamwell = oamwellfixed
        else:
            Coamwell = oamwelladjfactor * (1. / 100. * (Cwell + Cgath) + 0.25 * Claborcorrelation)

        # water O&M cost
        if oamwaterfixedvalid == 1:
            Coamwater = oamwaterfixed
        else:
            Coamwater = oamwateradjfactor * (
                        nprod * prodwellflowrate * waterloss * utilfactor * 365. * 24. * 3600. / 1E6 * 925. / 1E6)  # here is assumed 1 l per kg maybe correct with real temp. (M$/year) 925$/ML = 3.5$/1,000 gallon

        Coam = Coamwell + Coamplant + Coamwater  # total O&M cost (M$/year)
    else:
        Coam = oamtotalfixed  # total O&M cost (M$/year)

    if redrill > 0:  # account for well redrilling
        Coam = Coam + (Cwell + Cstim) * redrill / plantlifetime

    # ---------------------------------------------
    # Calculate annual electricity/heat production
    # ---------------------------------------------
    HeatkWhExtracted = np.zeros(
        plantlifetime)  # all end-use options have "heat extracted from reservoir" and pumping kWs
    PumpingkWh = np.zeros(plantlifetime)
    for i in range(0, plantlifetime):
        HeatkWhExtracted[i] = np.trapz(HeatExtracted[(0 + i * timestepsperyear):((i + 1) * timestepsperyear) + 1],
                                       dx=1. / timestepsperyear * 365. * 24.) * 1000. * utilfactor
        PumpingkWh[i] = np.trapz(PumpingPower[(0 + i * timestepsperyear):((i + 1) * timestepsperyear) + 1],
                                 dx=1. / timestepsperyear * 365. * 24.) * 1000. * utilfactor
    if enduseoption == 1 or enduseoption > 2:  # all these end-use options have an electricity generation component
        TotalkWhProduced = np.zeros(plantlifetime)
        NetkWhProduced = np.zeros(plantlifetime)
        for i in range(0, plantlifetime):
            TotalkWhProduced[i] = np.trapz(
                ElectricityProduced[(0 + i * timestepsperyear):((i + 1) * timestepsperyear) + 1],
                dx=1. / timestepsperyear * 365. * 24.) * 1000. * utilfactor
            NetkWhProduced[i] = np.trapz(
                NetElectricityProduced[(0 + i * timestepsperyear):((i + 1) * timestepsperyear) + 1],
                dx=1. / timestepsperyear * 365. * 24.) * 1000. * utilfactor
    if enduseoption > 1:  # all those end-use options have a direct-use component
        HeatkWhProduced = np.zeros(plantlifetime)
        for i in range(0, plantlifetime):
            HeatkWhProduced[i] = np.trapz(HeatProduced[(0 + i * timestepsperyear):((i + 1) * timestepsperyear) + 1],
                                          dx=1. / timestepsperyear * 365. * 24.) * 1000. * utilfactor

    # --------------------------------
    # calculate reservoir heat content
    # --------------------------------
    InitialReservoirHeatContent = resvol * rhorock * cprock * (Trock - Tinj) / 1E15  # 10^15 J
    RemainingReservoirHeatContent = InitialReservoirHeatContent - np.cumsum(HeatkWhExtracted) * 3600 * 1E3 / 1E15

    # ---------------------------
    # Calculate LCOE/LCOH
    # ---------------------------
    if econmodel == 1:  # simple FCR model
        if enduseoption == 1:
            Price = (FCR * (1 + inflrateconstruction) * Ccap + Coam) / np.average(NetkWhProduced) * 1E8  # cents/kWh
        elif enduseoption == 2:
            averageannualpumpingcosts = np.average(PumpingkWh) * elecprice / 1E6  # M$/year
            Price = (FCR * (1 + inflrateconstruction) * Ccap + Coam + averageannualpumpingcosts) / np.average(
                HeatkWhProduced) * 1E8  # cents/kWh
            Price = Price * 2.931  # $/Million Btu
        elif enduseoption > 2:  # cogeneration
            if enduseoption % 10 == 1:  # heat sales is additional income revenue stream
                averageannualheatincome = np.average(
                    HeatkWhProduced) * heatprice / 1E6  # M$/year ASSUMING heatprice IS IN $/KWH FOR HEAT SALES
                Price = (FCR * (1 + inflrateconstruction) * Ccap + Coam - averageannualheatincome) / np.average(
                    NetkWhProduced) * 1E8  # cents/kWh
            elif enduseoption % 10 == 2:  # electricity sales is additional income revenue stream
                averageannualelectricityincome = np.average(NetkWhProduced) * elecprice / 1E6  # M$/year
                Price = (FCR * (1 + inflrateconstruction) * Ccap + Coam - averageannualelectricityincome) / np.average(
                    HeatkWhProduced) * 1E8  # cents/kWh
                Price = Price * 2.931  # $/MMBTU
    elif econmodel == 2:  # standard levelized cost model
        discountvector = 1. / np.power(1 + discountrate, np.linspace(0, plantlifetime - 1, plantlifetime))
        if enduseoption == 1:
            Price = ((1 + inflrateconstruction) * Ccap + np.sum(Coam * discountvector)) / np.sum(
                NetkWhProduced * discountvector) * 1E8  # cents/kWh
        elif enduseoption == 2:
            averageannualpumpingcosts = np.average(PumpingkWh) * elecprice / 1E6  # M$/year
            Price = ((1 + inflrateconstruction) * Ccap + np.sum(
                (Coam + PumpingkWh * elecprice / 1E6) * discountvector)) / np.sum(
                HeatkWhProduced * discountvector) * 1E8  # cents/kWh
            Price = Price * 2.931  # $/MMBTU
        elif enduseoption > 2:
            if enduseoption % 10 == 1:  # heat sales is additional income revenue stream
                annualheatincome = HeatkWhProduced * heatprice / 1E6  # M$/year ASSUMING heatprice IS IN $/KWH FOR HEAT SALES
                Price = ((1 + inflrateconstruction) * Ccap + np.sum(
                    (Coam - annualheatincome) * discountvector)) / np.sum(
                    NetkWhProduced * discountvector) * 1E8  # cents/kWh
            elif enduseoption % 10 == 2:  # electricity sales is additional income revenue stream
                annualelectricityincome = NetkWhProduced * elecprice / 1E6  # M$/year
                Price = ((1 + inflrateconstruction) * Ccap + np.sum(
                    (Coam - annualelectricityincome) * discountvector)) / np.sum(
                    HeatkWhProduced * discountvector) * 1E8  # cents/kWh
                Price = Price * 2.931  # $/MMBTU
    elif econmodel == 3:  # bicycle model
        iave = FIB * BIR * (1 - CTR) + (1 - FIB) * EIR  # average return on investment (tax and inflation adjusted)
        CRF = iave / (1 - np.power(1 + iave, -plantlifetime))  # capital recovery factor
        inflationvector = np.power(1 + RINFL, np.linspace(1, plantlifetime, plantlifetime))
        discountvector = 1. / np.power(1 + iave, np.linspace(1, plantlifetime, plantlifetime))
        NPVcap = np.sum((1 + inflrateconstruction) * Ccap * CRF * discountvector)
        NPVfc = np.sum((1 + inflrateconstruction) * Ccap * PTR * inflationvector * discountvector)
        NPVit = np.sum(
            CTR / (1 - CTR) * ((1 + inflrateconstruction) * Ccap * CRF - Ccap / plantlifetime) * discountvector)
        NPVitc = (1 + inflrateconstruction) * Ccap * RITC / (1 - CTR)
        if enduseoption == 1:
            NPVoandm = np.sum(Coam * inflationvector * discountvector)
            NPVgrt = GTR / (1 - GTR) * (NPVcap + NPVoandm + NPVfc + NPVit - NPVitc)
            Price = (NPVcap + NPVoandm + NPVfc + NPVit + NPVgrt - NPVitc) / np.sum(
                NetkWhProduced * inflationvector * discountvector) * 1E8
        elif enduseoption == 2:
            PumpingCosts = PumpingkWh * elecprice / 1E6
            averageannualpumpingcosts = np.average(PumpingkWh) * elecprice / 1E6  # M$/year
            NPVoandm = np.sum((Coam + PumpingCosts) * inflationvector * discountvector)
            NPVgrt = GTR / (1 - GTR) * (NPVcap + NPVoandm + NPVfc + NPVit - NPVitc)
            Price = (NPVcap + NPVoandm + NPVfc + NPVit + NPVgrt - NPVitc) / np.sum(
                HeatkWhProduced * inflationvector * discountvector) * 1E8
            Price = Price * 2.931  # $/MMBTU
        elif enduseoption > 2:
            if enduseoption % 10 == 1:  # heat sales is additional income revenue stream
                annualheatincome = HeatkWhProduced * heatprice / 1E6  # M$/year ASSUMING ELECPRICE IS IN $/KWH FOR HEAT SALES
                NPVoandm = np.sum(Coam * inflationvector * discountvector)
                NPVgrt = GTR / (1 - GTR) * (NPVcap + NPVoandm + NPVfc + NPVit - NPVitc)
                Price = (NPVcap + NPVoandm + NPVfc + NPVit + NPVgrt - NPVitc - np.sum(
                    annualheatincome * inflationvector * discountvector)) / np.sum(
                    NetkWhProduced * inflationvector * discountvector) * 1E8
            elif enduseoption % 10 == 2:  # electricity sales is additional income revenue stream
                annualelectricityincome = NetkWhProduced * elecprice / 1E6  # M$/year
                NPVoandm = np.sum(Coam * inflationvector * discountvector)
                NPVgrt = GTR / (1 - GTR) * (NPVcap + NPVoandm + NPVfc + NPVit - NPVitc)
                Price = (NPVcap + NPVoandm + NPVfc + NPVit + NPVgrt - NPVitc - np.sum(
                    annualelectricityincome * inflationvector * discountvector)) / np.sum(
                    HeatkWhProduced * inflationvector * discountvector) * 1E8
                Price = Price * 2.931  # $/MMBTU
    # ---------------------------------------
    # write results to output file and screen
    # ---------------------------------------
    f = open('HDR.out', 'w')
    f.write('                               *****************\n')
    f.write('                               ***CASE REPORT***\n')
    f.write('                               *****************\n')
    f.write('\n')
    f.write('                           ***SUMMARY OF RESULTS***\n')
    f.write('\n')
    if enduseoption == 1:
        f.write("      End-Use Option = Electricity\n")
    elif enduseoption == 2:
        f.write("      End-Use Option = Direct-Use Heat\n")
    elif enduseoption == 31:  # topping cycle
        f.write("      End-Use Option = Cogeneration Topping Cycle\n")
        f.write("      Heat sales considered as extra income\n")
    elif enduseoption == 32:  # topping cycle
        f.write("      End-Use Option = Cogeneration Topping Cycle\n")
        f.write("      Electricity sales considered as extra income\n")
    elif enduseoption == 41:  # bottoming cycle
        f.write("      End-Use Option = Cogeneration Bottoming Cycle\n")
        f.write("      Heat Sales considered as extra income\n")
    elif enduseoption == 42:  # bottoming cycle
        f.write("      End-Use Option = Cogeneration Bottoming Cycle\n")
        f.write("      Electricity sales considered as extra income\n")
    elif enduseoption == 51:  # cogen split of mass flow rate
        f.write("      End-Use Option = Cogeneration Parallel Cycle\n")
        f.write("      Heat sales considered as extra income\n")
    elif enduseoption == 52:  # cogen split of mass flow rate
        f.write("      End-Use Option = Cogeneration Parallel Cycle\n")
        f.write("      Electricity sales considered as extra income\n")

    if enduseoption == 1 or enduseoption > 2:  # there is an electricity component
        f.write("      Average Net Electricity Production (MWe)          " + "{0:10.2f}".format(
            np.average(NetElectricityProduced)) + "\n")
    if enduseoption > 1:  # there is a direct-use component
        f.write("      Average Direct-Use Heat Production (MWth)         " + "{0:10.2f}".format(
            np.average(HeatProduced)) + "\n")

    if (enduseoption % 10) == 1:  # levelized cost expressed as LCOE
        f.write("      Electricity breakeven price (cents/kWh)           " + "{0:10.2f}".format(Price) + "\n")
    elif (enduseoption % 10) == 2:  # levelized cost expressed as LCOH
        f.write("      Direct-Use heat breakeven price ($/MMBTU)         " + "{0:10.2f}".format(Price) + "\n")

    f.write("      Number of production wells                     " + "{0:10.0f}".format(nprod) + "\n")
    f.write("      Number of injection wells                      " + "{0:10.0f}".format(ninj) + "\n")
    f.write("      Flowrate per production well (kg/s)              " + "{0:10.1f}".format(prodwellflowrate) + "\n")
    f.write("      Well depth (m)                                   " + "{0:10.1f}".format(depth) + "\n")
    if numseg == 1:
        f.write(
            '      Geothermal gradient (deg.C/km)                   ' + "{0:10.1f}".format((gradient[0] * 1E3)) + '\n')
    elif numseg == 2:
        f.write(
            '      Segment 1 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[0] * 1E3)) + '\n')
        f.write('      Segment 1 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[0] / 1E3)) + '\n')
        f.write(
            '      Segment 2 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[1] * 1E3)) + '\n')
    elif numseg == 3:
        f.write(
            '      Segment 1 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[0] * 1E3)) + '\n')
        f.write('      Segment 1 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[0] / 1E3)) + '\n')
        f.write(
            '      Segment 2 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[1] * 1E3)) + '\n')
        f.write('      Segment 2 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[1] / 1E3)) + '\n')
        f.write(
            '      Segment 3 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[2] * 1E3)) + '\n')
    elif numseg == 4:
        f.write(
            '      Segment 1 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[0] * 1E3)) + '\n')
        f.write('      Segment 1 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[0] / 1E3)) + '\n')
        f.write(
            '      Segment 2 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[1] * 1E3)) + '\n')
        f.write('      Segment 2 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[1] / 1E3)) + '\n')
        f.write(
            '      Segment 3 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[2] * 1E3)) + '\n')
        f.write('      Segment 3 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[2] / 1E3)) + '\n')
        f.write(
            '      Segment 4 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[3] * 1E3)) + '\n')

    f.write('\n')
    f.write('\n')
    f.write('                           ***ECONOMIC PARAMETERS***\n')
    f.write('\n')
    if econmodel == 1:
        f.write("      Economic Model = Fixed Charge Rate (FCR) Model\n")
        f.write("      Fixed Charge Rate (FCR) (%)                       " + "{0:10.2f}".format((FCR * 100)) + '\n')
    elif econmodel == 2:
        f.write("      Economic Model = Standard Levelized Cost Model\n")
        f.write("      Interest Rate (%)                                 " + "{0:10.2f}".format(
            (discountrate * 100)) + '\n')
    elif econmodel == 3:
        f.write("      Economic Model  = BICYCLE Model\n")
    f.write('      Accrued financing during construction (%)         ' + "{0:10.2f}".format(
        (inflrateconstruction * 100)) + '\n')
    f.write('      Project lifetime (years)                       ' + "{0:10.0f}".format((plantlifetime)) + '\n')
    f.write('      Capacity factor (%)                              ' + "{0:10.1f}".format((utilfactor) * 100) + '\n')

    f.write('\n')
    f.write('                          ***ENGINEERING PARAMETERS***\n')
    f.write('\n')
    f.write("      Well depth (m)                                   " + "{0:10.1f}".format(depth) + "\n")
    f.write("      Water loss rate (%)                              " + "{0:10.1f}".format(waterloss * 100) + "\n")
    f.write("      Pump efficiency (%)                              " + "{0:10.1f}".format(pumpeff * 100) + "\n")
    f.write("      Injection temperature (deg.C)                    " + "{0:10.1f}".format(Tinj) + "\n")
    if rameyoptionprod == 1:
        f.write("      Production Wellbore heat transmission calculated with Ramey's model\n")
        f.write("      Average production well temperature drop (deg.C) " + "{0:10.1f}".format(
            np.average(ProdTempDrop)) + "\n")
    elif rameyoptionprod == 0:
        f.write("      User-provided production well temperature drop\n")
        f.write("      Constant production well temperature drop (deg.C)" + "{0:10.1f}".format(tempdropprod) + "\n")
    f.write("      Flowrate per production well (kg/s)              " + "{0:10.1f}".format(prodwellflowrate) + "\n")
    f.write(
        "      Injection well casing ID (inches)                  " + "{0:10.3f}".format(injwelldiam / 0.0254) + "\n")
    f.write(
        "      Produciton well casing ID (inches)                 " + "{0:10.3f}".format(prodwelldiam / 0.0254) + "\n")
    f.write("      Number of times redrilling                     " + "{0:10.0f}".format(redrill) + "\n")
    if enduseoption == 1 or enduseoption > 2:
        if pptype == 1:
            f.write("      Power plant type                                        Subcritical ORC\n")
        elif pptype == 2:
            f.write("      Power plant type                                        Supercritical ORC\n")
        elif pptype == 3:
            f.write("      Power plant type                                        Single-Flash\n")
        elif pptype == 4:
            f.write("      Power plant type                                        Double-Flash\n")
    f.write('\n')
    f.write('\n')
    f.write('                         ***RESOURCE CHARACTERISTICS***\n')
    f.write('\n')
    f.write('      Maximum reservoir temperature (deg.C)            ' + "{0:10.1f}".format((Tmax)) + '\n')
    f.write('      Number of segments                             ' + "{0:10.0f}".format((numseg)) + '\n')
    if numseg == 1:
        f.write(
            '      Geothermal gradient (deg.C/km)                   ' + "{0:10.1f}".format((gradient[0] * 1E3)) + '\n')
    elif numseg == 2:
        f.write(
            '      Segment 1 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[0] * 1E3)) + '\n')
        f.write('      Segment 1 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[0] / 1E3)) + '\n')
        f.write(
            '      Segment 2 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[1] * 1E3)) + '\n')
    elif numseg == 3:
        f.write(
            '      Segment 1 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[0] * 1E3)) + '\n')
        f.write('      Segment 1 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[0] / 1E3)) + '\n')
        f.write(
            '      Segment 2 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[1] * 1E3)) + '\n')
        f.write('      Segment 2 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[1] / 1E3)) + '\n')
        f.write(
            '      Segment 3 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[2] * 1E3)) + '\n')
    elif numseg == 4:
        f.write(
            '      Segment 1 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[0] * 1E3)) + '\n')
        f.write('      Segment 1 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[0] / 1E3)) + '\n')
        f.write(
            '      Segment 2 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[1] * 1E3)) + '\n')
        f.write('      Segment 2 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[1] / 1E3)) + '\n')
        f.write(
            '      Segment 3 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[2] * 1E3)) + '\n')
        f.write('      Segment 3 thickness (km)                       ' + "{0:10.0f}".format(
            (layerthickness[2] / 1E3)) + '\n')
        f.write(
            '      Segment 4 geothermal gradient (deg.C/km)         ' + "{0:10.1f}".format((gradient[3] * 1E3)) + '\n')

    f.write('\n')
    f.write('\n')
    f.write('                           ***RESERVOIR PARAMETERS***\n')
    f.write('\n')
    if resoption == 1:
        f.write("      Reservoir Model = Multiple Parallel Fractures Model\n")
    elif resoption == 2:
        f.write("      Reservoir Model = 1-D Linear Heat Sweep Model\n")
    elif resoption == 3:
        f.write("      Reservoir Model = Single Fracture m/A Thermal Drawdown Model\n")
        f.write("      m/A Drawdown Parameter (kg/s/m^2)                       " + "{0:.5f}".format(drawdp) + "\n")
    elif resoption == 4:
        f.write("      Reservoir Model = Annual Percentage Thermal Drawdown Model\n")
        f.write(
            "      Annual Thermal Drawdown (%/year)                        " + "{0:.3f}".format(drawdp * 100) + "\n")
    elif resoption == 5:
        f.write("      Reservoir Model = User-Provided Temperature Profile\n")
    elif resoption == 6:
        f.write("      Reservoir Model = TOUGH2 Simulator\n")

    f.write("      Bottom-hole temperature (deg.C)                   " + "{0:10.2f}".format((Trock)) + '\n')
    if resoption > 3:
        f.write('      Warning: the reservoir dimensions and thermo-physical properties \n')
        f.write('               listed below are default values if not provided by the user.   \n')
        f.write('               They are only used for calculating remaining heat content.  \n')

    if resoption == 1 or resoption == 2:
        if fracshape == 1:
            f.write('      Fracture model = circular fracture with known area \n')
            f.write("      Well seperation = fracture diameter (m)           " + "{0:10.2f}".format(fracheight) + '\n')
        elif fracshape == 2:
            f.write('      Fracture model = circular fracture with known diameter\n')
            f.write("      Well seperation = fracture diameter (m)           " + "{0:10.2f}".format(fracheight) + '\n')
        elif fracshape == 3:
            f.write('      Fracture model = square fracture with known fracture height\n')
            f.write("      Well seperation = fracture height (m)             " + "{0:10.2f}".format(fracheight) + '\n')
        elif fracshape == 4:
            f.write('      Fracture model = rectangular fracture with known fracture height and width\n')
            f.write("      Well seperation = fracture height (m)             " + "{0:10.2f}".format(fracheight) + '\n')
            f.write("      Fracture width (m)                                " + "{0:10.2f}".format(fracwidth) + '\n')
        f.write("      Fracture area (m^2)                            " + "{0:10.0f}".format(fracarea) + '\n')
    if resvoloption == 1:
        f.write('      Reservoir volume calculated with fracture separation and number of fractures as input\n')
    elif resvoloption == 2:
        f.write('      Number of fractures calculated with reservoir volume and fracture separation as input\n')
    elif resvoloption == 3:
        f.write('      Fracture separation calculated with reservoir volume and number of fractures as input\n')
    elif resvoloption == 4:
        f.write('      Reservoir volume provided as input\n')
    if resvoloption in [1, 2, 3]:
        f.write("      Number of fractures                               " + "{0:10.2f}".format(fracnumb) + '\n')
        f.write("      Fracture separation (m)                           " + "{0:10.2f}".format(fracsep) + '\n')
    f.write("      Reservoir volume (m^3)                         " + "{0:10.0f}".format(resvol) + '\n')
    if impedancemodelused == 1:
        f.write(
            "      Reservoir impedance (GPa/m^3/s)                   " + "{0:10.2f}".format((impedance / 1000)) + '\n')
    else:
        f.write("      Reservoir hydrostatic pressure (kPa)              " + "{0:10.2f}".format(Phydrostatic) + '\n')
        f.write("      Plant outlet pressure (kPa)                       " + "{0:10.2f}".format(Pplantoutlet) + '\n')
        if productionwellpumping == 1:
            f.write(
                "      Production wellhead pressure (kPa)                " + "{0:10.2f}".format(Pprodwellhead) + '\n')
            f.write("      Productivity Index (kg/s/bar)                     " + "{0:10.2f}".format(PI) + '\n')
        f.write("      Injectivity Index (kg/s/bar)                      " + "{0:10.2f}".format(II) + '\n')

    f.write("      Reservoir density (kg/m^3)                        " + "{0:10.2f}".format((rhorock)) + '\n')
    if rameyoptionprod == 1 or resoption in [1, 2, 3, 6]:
        f.write("      Reservoir thermal conductivity (W/m/K)            " + "{0:10.2f}".format((krock)) + '\n')
    f.write("      Reservoir heat capacity (J/kg/K)                  " + "{0:10.2f}".format((cprock)) + '\n')
    if resoption == 2 or (resoption == 6 and usebuiltintough2model == 1):
        f.write("      Reservoir porosity (%)                            " + "{0:10.2f}".format((porrock * 100)) + '\n')
    if resoption == 6 and usebuiltintough2model == 1:
        f.write("      Reservoir permeability (m^2)                      " + "{0:10.2E}".format(permrock) + '\n')
        f.write("      Reservoir thickness (m)                           " + "{0:10.2f}".format(resthickness) + '\n')
        f.write("      Reservoir width (m)                               " + "{0:10.2f}".format(reswidth) + '\n')
        f.write("      Well separation (m)                               " + "{0:10.2f}".format(wellsep) + '\n')

    f.write('\n')
    f.write('\n')
    f.write('                          ***CAPITAL COSTS (M$)***\n')
    f.write('\n')
    if totalcapcostvalid == 0:
        f.write('      Drilling and completion costs                     ' + "{0:10.2f}".format((Cwell)) + '\n')
        f.write('      Drilling and completion costs per well            ' + "{0:10.2f}".format(
            (Cwell / (nprod + ninj))) + '\n')
        f.write('      Stimulation costs                                 ' + "{0:10.2f}".format((Cstim)) + '\n')
        f.write('      Surface power plant costs                         ' + "{0:10.2f}".format((Cplant)) + '\n')
        f.write('      Field gathering system costs                      ' + "{0:10.2f}".format((Cgath)) + '\n')
        if pipinglength > 0:
            f.write('      Transmission pipeline cost                        ' + "{0:10.2f}".format((Cpiping)) + '\n')
        f.write(
            '      Total surface equipment costs                     ' + "{0:10.2f}".format((Cplant + Cgath)) + '\n')
        f.write('      Exploration costs                                 ' + "{0:10.2f}".format((Cexpl)) + '\n')
    if totalcapcostvalid == 1 and redrill > 0:
        f.write('      Drilling and completion costs (for redrilling)    ' + "{0:10.2f}".format((Cwell)) + '\n')
        f.write('      Drilling and completion costs per redrilled well  ' + "{0:10.2f}".format(
            (Cwell / (nprod + ninj))) + '\n')
        f.write('      Stimulation costs (for redrilling)                ' + "{0:10.2f}".format((Cstim)) + '\n')
    f.write('      Total capital costs                               ' + "{0:10.2f}".format((Ccap)) + '\n')
    if econmodel == 1:
        f.write('      Annualized capital costs                          ' + "{0:10.2f}".format(
            (Ccap * (1 + inflrateconstruction) * FCR)) + '\n')

    f.write('\n')
    f.write('\n')
    f.write('                ***OPERATING AND MAINTENANCE COSTS (M$/yr)***\n')
    f.write('\n')
    if oamtotalfixedvalid == 0:
        f.write('      Wellfield maintenance costs                       ' + "{0:10.2f}".format((Coamwell)) + '\n')
        f.write('      Power plant maintenance costs                     ' + "{0:10.2f}".format((Coamplant)) + '\n')
        f.write('      Water costs                                       ' + "{0:10.2f}".format((Coamwater)) + '\n')
        if enduseoption == 2:
            f.write('      Average annual pumping costs                      ' + "{0:10.2f}".format(
                (averageannualpumpingcosts)) + '\n')
            f.write('      Total operating and maintenance costs             ' + "{0:10.2f}".format(
                (Coam + averageannualpumpingcosts)) + '\n')
        else:
            f.write('      Total operating and maintenance costs             ' + "{0:10.2f}".format((Coam)) + '\n')

    f.write('\n')
    f.write('\n')
    f.write('                           ***POWER GENERATION RESULTS***\n')
    f.write('\n')
    if enduseoption == 1 or enduseoption > 2:  # there is electricity component
        f.write(
            '      Initial geofluid availability (MWe/(kg/s)         ' + "{0:10.2f}".format((Availability[0])) + '\n')
        f.write('      Initial net power generation (MWe)                ' + "{0:10.2f}".format(
            (NetElectricityProduced[0])) + '\n')
        f.write('      Average net power generation (MWe)                ' + "{0:10.2f}".format(
            np.average(NetElectricityProduced)) + '\n')
        f.write('      Initial pumping power/net installed power (%)     ' + "{0:10.2f}".format(
            (PumpingPower[0] / NetElectricityProduced[0] * 100)) + '\n')
        f.write('      Average Annual Net Electricity Generation (GWh/yr)' + "{0:10.2f}".format(
            np.average(NetkWhProduced / 1E6)) + '\n')
    if enduseoption > 1:  # there is direct-use component
        f.write(
            '      Initial direct-use heat production (MWth)         ' + "{0:10.2f}".format((HeatProduced[0])) + '\n')
        f.write('      Average direct-use heat production (MWth)         ' + "{0:10.2f}".format(
            np.average(HeatProduced)) + '\n')
        f.write('      Average annual heat production (GWh/yr)           ' + "{0:10.2f}".format(
            np.average(HeatkWhProduced / 1E6)) + '\n')

    if impedancemodelused == 1:
        f.write('      Average total geofluid pressure drop (kPa)       ' + "{0:10.1f}".format(np.average(DP)) + '\n')
        f.write('      Average injection well pressure drop (kPa)       ' + "{0:10.1f}".format(np.average(DP1)) + '\n')
        f.write('      Average reservoir pressure drop (kPa)            ' + "{0:10.1f}".format(np.average(DP2)) + '\n')
        f.write('      Average production well pressure drop (kPa)      ' + "{0:10.1f}".format(np.average(DP3)) + '\n')
        f.write('      Average buoyancy correction (kPa)                ' + "{0:10.1f}".format(np.average(DP4)) + '\n')
    else:
        f.write(
            '      Average injection well pump pressure drop (kPa)       ' + "{0:10.1f}".format(np.average(DP1)) + '\n')
        if productionwellpumping == 1:
            f.write('      Average production well pump pressure drop (kPa)      ' + "{0:10.1f}".format(
                np.average(DP3)) + '\n')

    f.write('\n')
    f.write('                                        ******************************\n')
    f.write('                                        *  POWER GENERATION PROFILE  *\n')
    f.write('                                        ******************************\n')
    if enduseoption == 1:  # only electricity
        f.write(
            '  YEAR            THERMAL                    GEOFLUID                    PUMP                    NET                    FIRST LAW\n')
        f.write(
            '                  DRAWDOWN                  TEMPERATURE                  POWER                  POWER                   EFFICIENCY\n')
        f.write(
            '                                              (deg C)                    (MWe)                  (MWe)                       (%)\n')
        for i in range(0, plantlifetime + 1):
            f.write(
                '  {0:2.0f}              {1:8.4f}                   {2:8.2f}                  {3:8.4f}               {4:8.4f}                   {5:8.4f}'.format(
                    i, ProducedTemperature[i * timestepsperyear] / ProducedTemperature[0],
                    ProducedTemperature[i * timestepsperyear], PumpingPower[i * timestepsperyear],
                    NetElectricityProduced[i * timestepsperyear],
                    FirstLawEfficiency[i * timestepsperyear] * 100) + '\n')
    elif enduseoption == 2:  # only direct-use
        f.write(
            '  YEAR            THERMAL                    GEOFLUID                    PUMP                    NET\n')
        f.write(
            '                  DRAWDOWN                  TEMPERATURE                  POWER                   HEAT\n')
        f.write(
            '                                              (deg C)                    (MWe)                  (MWth)\n')
        for i in range(0, plantlifetime + 1):
            f.write(
                '  {0:2.0f}              {1:8.4f}                   {2:8.2f}                  {3:8.4f}               {4:8.4f}'.format(
                    i, ProducedTemperature[i * timestepsperyear] / ProducedTemperature[0],
                    ProducedTemperature[i * timestepsperyear], PumpingPower[i * timestepsperyear],
                    HeatProduced[i * timestepsperyear]) + '\n')
    elif enduseoption > 2:  # both electricity and direct-use
        f.write(
            '  YEAR            THERMAL                    GEOFLUID                    PUMP                    NET                     NET                    FIRST LAW\n')
        f.write(
            '                  DRAWDOWN                  TEMPERATURE                  POWER                  POWER                    HEAT                   EFFICIENCY\n')
        f.write(
            '                                              (deg C)                    (MWe)                  (MWe)                   (MWth)                      (%)\n')
        for i in range(0, plantlifetime + 1):
            f.write(
                '  {0:2.0f}              {1:8.4f}                   {2:8.2f}                  {3:8.4f}               {4:8.4f}                   {5:8.4f}                    {6:8.4f}'.format(
                    i, ProducedTemperature[i * timestepsperyear] / ProducedTemperature[0],
                    ProducedTemperature[i * timestepsperyear], PumpingPower[i * timestepsperyear],
                    NetElectricityProduced[i * timestepsperyear], HeatProduced[i * timestepsperyear],
                    FirstLawEfficiency[i * timestepsperyear] * 100) + '\n')
    f.write('\n')

    f.write('\n')
    f.write('                              ***************************************************************\n')
    f.write('                              *  HEAT AND/OR ELECTRICITY EXTRACTION AND GENERATION PROFILE  *\n')
    f.write('                              ***************************************************************\n')
    if enduseoption == 1:  # only electricity
        f.write(
            '  YEAR             ELECTRICITY                   HEAT                RESERVOIR            PERCENTAGE OF\n')
        f.write(
            '                    PROVIDED                   EXTRACTED            HEAT CONTENT        TOTAL HEAT MINED\n')
        f.write('                   (GWh/year)                  (GWh/year)            (10^15 J)                 (%)\n')
        for i in range(0, plantlifetime):
            f.write(
                '  {0:2.0f}              {1:8.1f}                    {2:8.1f}              {3:8.2f}               {4:8.2f}'.format(
                    i + 1, NetkWhProduced[i] / 1E6, HeatkWhExtracted[i] / 1E6, RemainingReservoirHeatContent[i], (
                                InitialReservoirHeatContent - RemainingReservoirHeatContent[
                            i]) * 100 / InitialReservoirHeatContent) + '\n')
    elif enduseoption == 2:  # only direct-use
        f.write(
            '  YEAR               HEAT                       HEAT                RESERVOIR            PERCENTAGE OF\n')
        f.write(
            '                    PROVIDED                   EXTRACTED            HEAT CONTENT        TOTAL HEAT MINED\n')
        f.write('                   (GWh/year)                  (GWh/year)            (10^15 J)                 (%)\n')
        for i in range(0, plantlifetime):
            f.write(
                '  {0:2.0f}              {1:8.1f}                    {2:8.1f}              {3:8.2f}               {4:8.2f}'.format(
                    i + 1, HeatkWhProduced[i] / 1E6, HeatkWhExtracted[i] / 1E6, RemainingReservoirHeatContent[i], (
                                InitialReservoirHeatContent - RemainingReservoirHeatContent[
                            i]) * 100 / InitialReservoirHeatContent) + '\n')
    elif enduseoption > 2:  # both electricity and direct-use
        f.write(
            '  YEAR               HEAT                   ELECTRICITY                  HEAT                RESERVOIR            PERCENTAGE OF\n')
        f.write(
            '                    PROVIDED                 PROVIDED                  EXTRACTED            HEAT CONTENT        TOTAL HEAT MINED\n')
        f.write(
            '                   (GWh/year)               (GWh/year)                 (GWh/year)            (10^15 J)                 (%)\n')
        for i in range(0, plantlifetime):
            f.write(
                '  {0:2.0f}              {1:8.1f}                 {2:8.1f}                    {3:8.2f}              {4:8.2f}               {5:8.2f}'.format(
                    i + 1, HeatkWhProduced[i] / 1E6, NetkWhProduced[i] / 1E6, HeatkWhExtracted[i] / 1E6,
                    RemainingReservoirHeatContent[i], (InitialReservoirHeatContent - RemainingReservoirHeatContent[
                        i]) * 100 / InitialReservoirHeatContent) + '\n')

    f.write('\n')

    f.close()

    # print results to console screen
    if printoutput == 1:
        print("")
        print("----------------------------")
        print("GEOPHIRES Simulation Results")
        print("----------------------------")
        print("")
        print("1. Simulation Metadata")
        print("----------------------")
        print(" GEOPHIRES Version = 2.0")
        print(" GEOPHIRES Build Date = 2018-01-02")
        currentdate = datetime.datetime.now().strftime("%Y-%m-%d")
        currenttime = datetime.datetime.now().strftime("%H:%M")
        print(" Simulation Date = " + currentdate)
        print(" Simulation Time = " + currenttime)
        print(" Calculation Time = " + "{0:.3f}".format((time.time() - tic)) + " s")

        print("")
        print("2. Summary of Simulation Results")
        print("--------------------------------")
        if printoutput == 1:
            # say what type of end-use option
            if enduseoption == 1:
                print(" End-Use Option = Electricity")
            elif enduseoption == 2:
                print(" End-Use Option = Direct-Use Heat")
            elif enduseoption == 31:  # topping cycle
                print(" End-Use Option = Cogeneration Topping Cycle")
                print(" Heat sales considered as extra income")
            elif enduseoption == 32:  # topping cycle
                print(" End-Use Option = Cogeneration Topping Cycle")
                print(" Electricity sales considered as extra income")
            elif enduseoption == 41:  # bottoming cycle
                print(" End-Use Option = Cogeneration Bottoming Cycle")
                print(" Heat Sales considered as extra income")
            elif enduseoption == 42:  # bottoming cycle
                print(" End-Use Option = Cogeneration Bottoming Cycle")
                print(" Electricity sales considered as extra income")
            elif enduseoption == 51:  # cogen split of mass flow rate
                print(" End-Use Option = Cogeneration Parallel Cycle")
                print(" Heat sales considered as extra income")
            elif enduseoption == 52:  # cogen split of mass flow rate
                print(" End-Use Option = Cogeneration Parallel Cycle")
                print(" Electricity sales considered as extra income")
            # say what type of power plant
            if enduseoption == 1 or enduseoption > 2:
                if pptype == 1:
                    print(" Power Plant Type = Subcritical ORC")
                elif pptype == 2:
                    print(" Power Plant Type = Supercritical ORC")
                elif pptype == 3:
                    print(" Power Plant Type = Single-Flash")
                elif pptype == 4:
                    print(" Power Plant Type = Double-Flash")

            # print(NetElectricityProduced)
            if enduseoption == 1 or enduseoption > 2:
                print(" Average Net Electricity Generation = " + "{0:.2f}".format(
                    np.average(NetElectricityProduced)) + " MWe")
            if enduseoption > 1:
                print(" Average Net Heat Production = " + "{0:.2f}".format(np.average(HeatProduced)) + " MWth")

            # print LCOE/LCOH
            if enduseoption == 1:
                print(" LCOE = " + "{0:.1f}".format((Price)) + " cents/kWh")
            elif enduseoption == 2:
                print(" LCOH = " + "{0:.1f}".format((Price)) + " $/MMBTU")
            elif enduseoption % 10 == 1:  # heat sales is additional income revenuw stream
                print(" LCOE = " + "{0:.1f}".format((Price)) + " cents/kWh")
                print(" Additional average annual revenue from heat sales = " + "{0:.1f}".format(
                    np.average(annualheatincome)) + " M$/year")
            elif enduseoption % 10 == 1:  # electricity sales is additional income revenuw stream
                print(" LCOH = " + "{0:.1f}".format((Price)) + " $/MMBTU")
                print(" Additional average annual revenue from electricity sales = " + "{0:.1f}".format(
                    np.average(annualelectricityincome)) + " M$/year")
            # say what type of economic model is used
            if econmodel == 1:
                print(" Economic Model Used = Fixed Charge Rate (FCR) Model")
                print(" Fixed Charge Rate (FCR) = " + "{0:.2f}".format((FCR * 100)) + "%")
            elif econmodel == 2:
                print(" Economic Model Used = Standard Levelized Cost Model")
                print(" Discount Rate = " + "{0:.2f}".format((discountrate * 100)) + "%")
            elif econmodel == 3:
                print(" Economic Model Used = BICYCLE Model")

            print("")
            print("3. Reservoir Simulation Results")
            print("-------------------------------")
            if resoption == 1:
                print(" Reservoir Model = Multiple Parallel Fractures Model")
            elif resoption == 2:
                print(" Reservoir Model = 1-D Linear Heat Sweep Model")
            elif resoption == 3:
                print(" Reservoir Model = Single Fracture m/A Thermal Drawdown Model")
                print(" m/A Drawdown Parameter = " + "{0:.5f}".format(drawdp) + " kg/s/m^2")
            elif resoption == 4:
                print(" Reservoir Model = Annual Percentage Thermal Drawdown Model")
                print(" Annual Thermal Drawdown = " + "{0:.3f}".format(drawdp * 100) + " %/year")
            elif resoption == 5:
                print(" Reservoir Model = User-Provided Temperature Profile")
            elif resoption == 6:
                print(" Reservoir Model = TOUGH2 EOS1 Simulator")

            print(" Number of Production Wells = " + "{0:.0f}".format((nprod)))
            print(" Number of Injection Wells = " + "{0:.0f}".format((ninj)))
            print(" Number of Times Redrilling = " + "{0:.0f}".format((redrill)))
            print(" Well Depth = " + "{0:.1f}".format((depth)) + " m")
            print(" Flow Rate per Production Well = " + "{0:.0f}".format((prodwellflowrate)) + " kg/s")
            print(" Initial Reservoir Temperature = " + "{0:.1f}".format(Trock) + "°C")
            print(" Maximum Production Temperature = " + "{0:.1f}".format(np.max(ProducedTemperature)) + "°C")
            print(" Average Production Temperature = " + "{0:.1f}".format(np.average(ProducedTemperature)) + "°C")
            print(" Minimum Production Temperature = " + "{0:.1f}".format(np.min(ProducedTemperature)) + "°C")
            print(" Initial Production Temperature = " + "{0:.1f}".format((ProducedTemperature[0])) + "°C")
            print(" Average Reservoir Heat Extraction = " + "{0:.2f}".format(np.average(HeatExtracted)) + " MWth")
            if rameyoptionprod == 1:
                print(" Production Wellbore Heat Transmission Model = Ramey Model")
                print(
                    " Average Production Well Temperature Drop = " + "{0:.1f}".format(np.average(ProdTempDrop)) + "°C")
            elif rameyoptionprod == 0:
                print(" Wellbore Heat Transmission Model = Constant Temperature Drop of " + "{0:.1f}".format(
                    tempdropprod) + "°C")
            if impedancemodelused == 1:
                print(" Total Average Pressure Drop = " + "{0:.1f}".format(np.average(DP)) + " kPa")
                print("    Average Injection Well Pressure Drop = " + "{0:.1f}".format(np.average(DP1)) + " kPa")
                print("    Average Reservoir Pressure Drop = " + "{0:.1f}".format(np.average(DP2)) + " kPa")
                print("    Average Production Well Pressure Drop = " + "{0:.1f}".format(np.average(DP3)) + " kPa")
                print("    Average Buoyancy Pressure Drop = " + "{0:.1f}".format(np.average(DP4)) + " kPa")
            else:
                print(" Average Injection Well Pump Pressure Drop = " + "{0:.1f}".format(np.average(DP1)) + " kPa")
                if productionwellpumping == 1:
                    print(" Average Production Well Pump Pressure Drop = " + "{0:.1f}".format(np.average(DP3)) + " kPa")

            print("")
            print("4. Surface Equipment Simulation Results")
            print("---------------------------------------")
            if enduseoption == 1 or enduseoption > 2:
                print(
                    " Maximum Total Electricity Generation = " + "{0:.2f}".format(np.max(ElectricityProduced)) + " MWe")
                print(" Average Total Electricity Generation = " + "{0:.2f}".format(
                    np.average(ElectricityProduced)) + " MWe")
                print(
                    " Minimum Total Electricity Generation = " + "{0:.2f}".format(np.min(ElectricityProduced)) + " MWe")
                print(" Initial Total Electricity Generation = " + "{0:.2f}".format((ElectricityProduced[0])) + " MWe")
                print(" Maximum Net Electricity Generation = " + "{0:.2f}".format(
                    np.max(NetElectricityProduced)) + " MWe")
                print(" Average Net Electricity Generation = " + "{0:.2f}".format(
                    np.average(NetElectricityProduced)) + " MWe")
                print(" Minimum Net Electricity Generation = " + "{0:.2f}".format(
                    np.min(NetElectricityProduced)) + " MWe")
                print(" Initial Net Electricity Generation = " + "{0:.2f}".format((NetElectricityProduced[0])) + " MWe")
                print(" Average Annual Total Electricity Generation = " + "{0:.2f}".format(
                    np.average(TotalkWhProduced / 1E6)) + " GWh")
                print(" Average Annual Net Electricity Generation = " + "{0:.2f}".format(
                    np.average(NetkWhProduced / 1E6)) + " GWh")
            if enduseoption > 1:
                print(" Maximum Net Heat Production = " + "{0:.2f}".format(np.max(HeatProduced)) + " MWth")
                print(" Average Net Heat Production = " + "{0:.2f}".format(np.average(HeatProduced)) + " MWth")
                print(" Minimum Net Heat Production = " + "{0:.2f}".format(np.min(HeatProduced)) + " MWth")
                print(" Initial Net Heat Production = " + "{0:.2f}".format((HeatProduced[0])) + " MWth")
                print(
                    " Average Annual Heat Production = " + "{0:.2f}".format(np.average(HeatkWhProduced / 1E6)) + " GWh")
            print(" Average Pumping Power = " + "{0:.2f}".format(np.average(PumpingPower)) + " MWe")

            print("")
            print("5. Capital and O&M Costs")
            print("------------------------")
            print(" Total Capital Cost = " + "{0:.2f}".format(Ccap) + " M$")
            if totalcapcostvalid == 0:
                print("   Wellfield Cost = " + "{0:.2f}".format(Cwell) + " M$")
                print("   Surface Plant Cost = " + "{0:.2f}".format(Cplant) + " M$")
                print("   Exploration Cost = " + "{0:.2f}".format(Cexpl) + " M$")
                print("   Field Gathering System Cost = " + "{0:.2f}".format(Cgath) + " M$")
                if pipinglength > 0:
                    print("   Transmission Pipeline Cost = " + "{0:.2f}".format(Cpiping) + " M$")
                print("   Stimulation Cost = " + "{0:.2f}".format(Cstim) + " M$")
            if enduseoption == 2:
                print(" Total O&M Cost = " + "{0:.2f}".format(Coam + averageannualpumpingcosts) + " M$/year")
            else:
                print(" Total O&M Cost = " + "{0:.2f}".format(Coam) + " M$/year")
            if oamtotalfixedvalid == 0:
                print("   Wellfield O&M Cost = " + "{0:.2f}".format(Coamwell) + " M$/year")
                print("   Surface Plant O&M Cost = " + "{0:.2f}".format(Coamplant) + " M$/year")
                print("   Make-Up Water O&M Cost = " + "{0:.2f}".format(Coamwater) + " M$/year")
                if enduseoption == 2:
                    print(
                        "   Average annual pumping costs = " + "{0:.2f}".format(averageannualpumpingcosts) + " M$/year")

            print("")
            print("6. Power Generation Profile")
            print("---------------------------")

    if enduseoption == 1:  # only electricity
        print('  YEAR   THERMAL     GEOFLUID       PUMP      NET      FIRST LAW')
        print('         DRAWDOWN    TEMPERATURE    POWER     POWER    EFFICIENCY')
        print('         (-)         (deg C)        (MWe)     (MWe)    (%)')
        for i in range(0, plantlifetime + 1):
            print('  {0:2.0f}   {1:8.4f}     {2:8.2f}      {3:8.4f}  {4:8.4f}  {5:8.4f}'.format(i, ProducedTemperature[
                i * timestepsperyear] / ProducedTemperature[0], ProducedTemperature[i * timestepsperyear], PumpingPower[
                                                                                                    i * timestepsperyear],
                                                                                                NetElectricityProduced[
                                                                                                    i * timestepsperyear],
                                                                                                FirstLawEfficiency[
                                                                                                    i * timestepsperyear] * 100))
    elif enduseoption == 2:  # only direct-use
        print('  YEAR   THERMAL      GEOFLUID      PUMP      NET')
        print('         DRAWDOWN     TEMPERATURE   POWER     HEAT')
        print('         (-)          (deg C)       (MWe)     (MWth)')
        for i in range(0, plantlifetime + 1):
            print('  {0:2.0f}   {1:8.4f}     {2:8.2f}      {3:8.4f}   {4:8.4f}'.format(i, ProducedTemperature[
                i * timestepsperyear] / ProducedTemperature[0], ProducedTemperature[i * timestepsperyear], PumpingPower[
                                                                                           i * timestepsperyear],
                                                                                       HeatProduced[
                                                                                           i * timestepsperyear]))
    elif enduseoption > 2:  # both electricity and direct-use
        print('  YEAR   THERMAL      GEOFLUID      PUMP      NET       NET       FIRST LAW')
        print('         DRAWDOWN     TEMPERATURE   POWER     POWER     HEAT      EFFICIENCY')
        print('         (-)          (deg C)       (MWe)     (MWe)     (MWth)    (%)')
        for i in range(0, plantlifetime + 1):
            print('  {0:2.0f}    {1:8.4f}    {2:8.2f}      {3:8.4f}   {4:8.4f}  {5:8.4f}  {6:8.4f}'.format(i,
                                                                                                           ProducedTemperature[
                                                                                                               i * timestepsperyear] /
                                                                                                           ProducedTemperature[
                                                                                                               0],
                                                                                                           ProducedTemperature[
                                                                                                               i * timestepsperyear],
                                                                                                           PumpingPower[
                                                                                                               i * timestepsperyear],
                                                                                                           NetElectricityProduced[
                                                                                                               i * timestepsperyear],
                                                                                                           HeatProduced[
                                                                                                               i * timestepsperyear],
                                                                                                           FirstLawEfficiency[
                                                                                                               i * timestepsperyear] * 100))


if __name__ == '__main__':
    run()

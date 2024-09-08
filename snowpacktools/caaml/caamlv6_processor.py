################################################################################
# Copyright 2022 Avalanche Warning Service Tyrol                               #
################################################################################
# This is free software you can redistribute/modify under the terms of the     #
# GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
################################################################################

import numpy as np
from lxml import etree as et
import os

caaml_ns = "{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}"
gml_ns = "{http://www.opengis.net/gml}"

def get_prof_metadata(file):
    '''Returns dictionary of relevant metadata for provided snowprofile.'''

    tree = et.parse(file)
    xroot = tree.getroot()
    
    """Altitude, location, aspect"""
    child_locRef            = xroot.find(caaml_ns + 'locRef')
    child_validElevation    = child_locRef.find(caaml_ns + 'validElevation')
    child_ElevationPosition = child_validElevation.find(caaml_ns + 'ElevationPosition')
    child_position          = child_ElevationPosition.find(caaml_ns + 'position')
    alt                     = float(child_position.text)

    child_pointLocation = child_locRef.find(caaml_ns + 'pointLocation')
    child_Point         = child_pointLocation.find(gml_ns + 'Point') # be careful with namespace (gml)
    child_pos           = child_Point.find(gml_ns + 'pos')
    lon_lat             = child_pos.text

    lon_lat = lon_lat.split(' ')
    lon     = float(lon_lat[0])
    lat     = float(lon_lat[1])

    child_validAspect    = child_locRef.find(caaml_ns + 'validAspect')
    child_AspectPosition = child_validAspect.find(caaml_ns + 'AspectPosition')
    child_positionA      = child_AspectPosition.find(caaml_ns + 'position')
    aspect               = child_positionA.text

    child_validSlopeAngle    = child_locRef.find(caaml_ns + 'validSlopeAngle')
    if child_validSlopeAngle is None:
        slope_angle = 0
    else:
        child_SlopeAnglePosition = child_validSlopeAngle.find(caaml_ns + 'SlopeAnglePosition')
        child_position_angle     = child_SlopeAnglePosition.find(caaml_ns + 'position')
        slope_angle              = int(child_position_angle.text)

    """Date"""
    child_timeRef      = xroot.find(caaml_ns + 'timeRef')
    child_recordTime   =  child_timeRef.find(caaml_ns + 'recordTime')
    child_TimeInstant  = child_recordTime.find(caaml_ns + 'TimeInstant')
    child_timePosition = child_TimeInstant.find(caaml_ns + 'timePosition')
    date               = child_timePosition.text # string

    """Name"""
    child_name = child_locRef.find(caaml_ns + 'name')
    name       = child_name.text

    prof_meta = {'lon'        : lon,
                 'lat'        : lat,
                 'alt'        : alt,
                 'datetime'   : date,
                 'name'       : name,
                 'aspect'     : aspect,
                 'slope_angle': slope_angle}

    return prof_meta


def add_snp_metadata(tree: object):
    """Add relevant metadata to caaml file (section customData:snp)"""
    
    xroot = tree.getroot()

    child_metaData                  = xroot.find(caaml_ns + "metaData")
    child_customData                = child_metaData.find(caaml_ns + "customData")
    if child_customData is None:
        child_customData            = et.SubElement(child_metaData, caaml_ns + "customData")

    child_snp                       = child_customData.find(caaml_ns + "snp")
    if child_snp is None:
        """Add SNP relevant metadata to CAAML file, if not existing"""
        child_snp                       = et.SubElement(child_customData,caaml_ns + 'snp')
        child_CanopyHeight              = et.SubElement(child_snp,caaml_ns + 'CanopyHeight')
        child_CanopyBasalArea           = et.SubElement(child_snp,caaml_ns + 'CanopyBasalArea')
        child_CanopyLAI                 = et.SubElement(child_snp,caaml_ns + 'CanopyLAI')
        child_CanopyDirectThroughfall   = et.SubElement(child_snp,caaml_ns + 'CanopyDirectThroughfall')
        child_SoilAlb                   = et.SubElement(child_snp,caaml_ns + 'SoilAlb')
        child_ErosionLevel              = et.SubElement(child_snp,caaml_ns + 'ErosionLevel')
        child_CanopyHeight.set("uom","m")
        child_CanopyBasalArea.set("uom","m^2")
        
        child_CanopyHeight.text             = "0.0"
        child_CanopyBasalArea.text          = "0.0"
        child_CanopyLAI.text                = "0.0"
        child_CanopyDirectThroughfall.text  = "1.0"
        child_SoilAlb.text                  = "0.2"
        child_ErosionLevel.text             = "0.0"

    child_locRef                = xroot.find(caaml_ns + "locRef")
    
    """Add name to location reference if missing"""
    child_name                  = child_locRef.find(caaml_ns + "name")
    if child_name.text is None:
        child_name.text = "nameless"

    """Add slope angle to location reference if missing"""
    child_validSlopeAngle       = child_locRef.find(caaml_ns + "validSlopeAngle")
    if child_validSlopeAngle is None:
        child_validSlopeAngle    = et.SubElement(child_locRef,caaml_ns + "validSlopeAngle")
        child_SlopeAnglePosition = et.SubElement(child_validSlopeAngle,caaml_ns + "SlopeAnglePosition")
        child_position_angle     = et.SubElement(child_SlopeAnglePosition,caaml_ns + "position")
        child_position_angle.text = str(0)

    """Add aspect to location reference if missing"""
    child_validAspect    = child_locRef.find(caaml_ns + "validAspect")
    if child_validAspect is None:
        child_validAspect    = et.SubElement(child_locRef,caaml_ns + "validAspect")
        child_AspectPosition = et.SubElement(child_validAspect,caaml_ns + "AspectPosition")
        child_positionA      = et.SubElement(child_AspectPosition,caaml_ns + "position")
        child_positionA.text = "N"
    

    return tree


def validate_and_prepare_for_snp(file, new_file):
    """
    Adds density profile based on grain shape and hardness parameterization alla Monti (2014)

    Arguments:
        file / new_file (str):          Path to xml file of snowprofile
        remove_observed_density (bool): Relevant when density is already provided from the observations (already included in .caaml file)

    Returns nothing, but saves new profile

    Comments:
        - If density exists from observation or simulation it is removed
    """

    try:
        tree = et.parse(file)
    except:
        print(f"[i]  Reading of observed snow profile {file.split('/')[-1]} failed. File will be removed.")
        os.remove(file)
        return

    tree = add_snp_metadata(tree)
    xroot = tree.getroot()
    
    """Add density section"""
    child_snowProfileResultsOf = xroot.find(caaml_ns + 'snowProfileResultsOf')
    # if child_snowProfileResultsOf is None:
    #     return # No profile data
    child_SnowProfileMeasurements = child_snowProfileResultsOf.find(caaml_ns + 'SnowProfileMeasurements')
    try:
        child_densityProfile = child_SnowProfileMeasurements.find(caaml_ns + 'densityProfile')
        child_SnowProfileMeasurements.remove(child_densityProfile)
    except:
        pass
    child_densityProfile = et.SubElement(child_SnowProfileMeasurements, caaml_ns + 'densityProfile')
    child_densityMetaData = et.SubElement(child_densityProfile, caaml_ns + 'densityMetaData')
    child_methodOfMeas = et.SubElement(child_densityMetaData, caaml_ns + 'methodOfMeas')
    child_methodOfMeas.text = 'other'
    
    """Check temperature profile (needs to be ascending order for SNP)"""
    child_tempProfile = child_SnowProfileMeasurements.find(caaml_ns + 'tempProfile')
    if child_tempProfile is None:
        print(f"[i]  Temperature profile missing for observed snow profile {file.split('/')[-1]}. File will be removed.")
        os.remove(file)
        return

    # depths    = []
    # temp_vals = []
    # for obs in child_tempProfile.iter(caaml_ns + 'Obs'):
    #     depths.append(obs.find(caaml_ns + 'depth').text)
    #     temp_vals.append(obs.find(caaml_ns + 'snowTemp').text)
    # print(depths)
    # print(temp_vals)


    """Add density layer for each stratigraphy layer (and check for valid grain size)"""
    child_stratProfile = child_SnowProfileMeasurements.find(caaml_ns + 'stratProfile')
    for layer in child_stratProfile.iter(caaml_ns + 'Layer'):
        

        child_densityLayer = et.SubElement(child_densityProfile, caaml_ns + 'Layer')
        child_depthTop  = et.SubElement(child_densityLayer,caaml_ns + 'depthTop')
        child_thickness = et.SubElement(child_densityLayer,caaml_ns + 'thickness')

        child_density   = et.SubElement(child_densityLayer,caaml_ns + 'density')

        child_depthTop.set('uom', 'cm')
        child_thickness.set('uom', 'cm')
        child_density.set('uom', 'kgm-3')
    
        for val_depth in layer.iter(caaml_ns + 'depthTop'):
            child_depthTop.text = val_depth.text
        for val_thickness in layer.iter(caaml_ns + 'thickness'):
            child_thickness.text = val_thickness.text
        
        """Get grain forms and hardness to estimate density"""
        for form in layer.iter(caaml_ns + 'grainFormPrimary'):
            form.text = check_grainType(form.text)
            grainFormPrimary = form.text
        grainFormSecondary = ""
        for form in layer.iter(caaml_ns + 'grainFormSecondary'):
            form.text = check_grainType(form.text)
            grainFormSecondary = form.text
        for val_hardness in layer.iter(caaml_ns + 'hardness'):
            val_hardness.text = check_handHardness(val_hardness.text)
            handhardness      = val_hardness.text
        
        child_density.text = str(monti_density(grainFormPrimary, grainFormSecondary, handhardness))

        """Check LWC/wetness (only full values are allowed)"""
        for val_wetness in layer.iter(caaml_ns + 'wetness'):
            val_wetness.text = check_wetness(val_wetness.text)

        """Check if grain size is available"""
        if grainFormPrimary in ["PP","DF","FC","RG"]:
            child_grainType = layer.find(caaml_ns + 'grainSize')
            if child_grainType is None:
                print(f"[i]  Grain size is missing for observed snow profile {file.split('/')[-1]}. File will be removed.")
                os.remove(file)
                return
                
    tree.write(new_file,encoding='UTF-8', method="xml", doctype='<?xml version="1.0" encoding="UTF-8"?>')


def monti_density(form, form2, hardness_str):
    """Density parameterisation from Monti 2014"""

    """
    #### - SNOW CRYSTAL SHAPES - ####
    Precip particles = PP
    Decomp/fragm     = df
    Rounded Grains   = RG
    Faceted cystals  = FC
    Depth hoar       = DH
    Surface hoar     = SH
    Meltforms        = MF
    Faceted rounded  = FCxr
    Graupel          = PPgp
    Melt freez crust = MFcr
    
    #### - STR OF HARDNESS TO FLOAT - ####
    1   = F- = Faust
    1-2 = F+
    2   = 4F = 4Finger
    3   = 1F = 1Finger
    4-4 = 1F+
    4   = P  = Bleistift
    4-5 = P+
    5   = K  = Knife 
    5-6 = K+ 
    6   = I  = Ice
    """
    
    dict_hardness = {'F' : 1,
                     'F+' : 1.5,
                     '4F' : 2,
                     '4F+' : 2.5,
                     '1F' : 3,
                     '1F+': 3.5,
                     'P'  : 4,
                     'P+' : 4.5,
                     'K'  : 5,
                     'K+' : 5.5,
                     'I'  : 6}
    
    dict_hardness_i =  {'F' : 0,
                        'F+' : 1,
                        '4F' : 2,
                        '4F+': 3,
                        '1F' : 4,
                        '1F+': 5,
                        'P'  : 6,
                        'P+' : 7,
                        'K'  : 8,
                        'K+' : 9,
                        'I'  : 10}
    

    ice = 910
    sho = 50
    mfc = 520
    dict_shapes = {'PP'  : [48,96,144,144,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                   'PPgp': [48,96,144,144,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                   'DF'  : [71,142,232,250,269,269,np.nan,np.nan,np.nan,np.nan,np.nan],
                   'SH'  : [sho,sho,sho,sho,sho,sho,sho,sho,sho,sho,sho],
                   'DH'  : [96,192,307,325,345,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                   'FC'  : [82,164,272,295,347,373,440,479,519,np.nan,np.nan],
                   'FCxr': [86,172,282,304,350,373,426,455,485,np.nan,np.nan],
                   'RG'  : [63,126,219,248,308,338,394,418,444,np.nan,np.nan],
                   'MF'  : [73,146,248,283,348,377,518,629,740,np.nan,np.nan],
                   'MFcr': [mfc,mfc,mfc,mfc,mfc,mfc,mfc,mfc,mfc,mfc,mfc],
                   'IF'  : [ice,ice,ice,ice,ice,ice,ice,ice,ice,ice,ice]}
    
    """
    1/3 of range used for full hardness, 2/3 used for .5 hardness
    
    [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6] -> nx = 11

    matrix = [[48,96,144,144,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
          [71,142,232,250,269,269,np.nan,np.nan,np.nan,np.nan,np.nan],
          [63,126,219,248,308,338,394,418,444,np.nan,np.nan],
          [82,164,272,295,347,373,440,479,519,np.nan,np.nan],
          [96,192,307,325,345,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
          [73,146,248,283,348,377,518,629,740,np.nan,np.nan],
          [86,172,282,304,350,373,426,455,485,np.nan,np.nan]]
    """
        
    hardness   = dict_hardness[hardness_str]
    hardness_i = dict_hardness_i[hardness_str]
    density    = dict_shapes[form][hardness_i]
    
    return density


def check_wetness(wet_input):
    if wet_input == "D-M":
        wet_output = "M" # LWC: 0.01
    elif wet_input == "M-W":
        wet_output = "W" # LWC: 0.03
    elif wet_input == "W-V":
        wet_output = "V" # LWC: 0.08
    elif wet_input == "V-S":
        wet_output = "S" # LWC: 0.15
    else:
        wet_output = wet_input
    return wet_output


def check_grainType(gt_input):
    if gt_input == "PPdc":
        gt_output = "PP"
    elif gt_input == "DFdc":
        gt_output = "DF"
    elif gt_input == "MFcl":
        gt_output = "MF"
    elif gt_input == "MFpc":
        gt_output = "MF"
    elif gt_input == "MFsl":
        gt_output = "MF"
    elif gt_input == "RGsr":
        gt_output = "RG"
    elif gt_input == "RGlr":
        gt_output = "RG"
    elif gt_input == "RGwp":
        gt_output = "RG"
    elif gt_input == "RGxf":
        gt_output = "RG"
    elif gt_input == "DHxr":
        gt_output = "DH"
    elif gt_input == "SHxr":
        gt_output = "SH"
    else:
        gt_output = gt_input
    return gt_output


def check_handHardness(hh_input):
    if hh_input == "F-":
        hh_output = "F"
    elif hh_input == "F-4F":
        hh_output = "F+"
    elif hh_input == "4F-":
        hh_output = "F+"
    elif hh_input == "4F-1F":
        hh_output = "4F+"
    elif hh_input == "1F-":
        hh_output = "4F+"
    elif hh_input == "1F-P":
        hh_output = "1F+"
    elif hh_input == "P-":
        hh_output = "1F+"
    elif hh_input == "P-K":
        hh_output = "P+"
    elif hh_input == "K-I":
        hh_output = "K+"
    else:
        hh_output = hh_input

    return hh_output


if __name__ == "__main__":
    print("This script provides useful functions to process CAAMLv6 snow profiles. Call them directly.")

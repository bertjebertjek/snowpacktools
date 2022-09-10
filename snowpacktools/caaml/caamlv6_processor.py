################################################################################
# Copyright 2022 Avalanche Warning Service Tyrol                               #
################################################################################
# This is free software you can redistribute/modify under the terms of the     #
# GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
################################################################################

import sys
import os
import numpy as np
from lxml import etree as et

def get_prof_metadata(file):
    '''Returns dictionary of relevant metadata for provided snowprofile.'''

    tree = et.parse(file)
    xroot = tree.getroot()
    
    """Altitude, location, aspect"""
    child_locRef = xroot.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}locRef')
    child_validElevation = child_locRef.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}validElevation')
    child_ElevationPosition = child_validElevation.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}ElevationPosition')
    child_position = child_ElevationPosition.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}position')
    alt = float(child_position.text)

    child_pointLocation = child_locRef.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}pointLocation')
    child_Point = child_pointLocation.find('{http://www.opengis.net/gml}Point') # be careful with namespace (gml)
    child_pos = child_Point.find('{http://www.opengis.net/gml}pos')
    lon_lat = child_pos.text

    lon_lat = lon_lat.split(' ')
    lon= float(lon_lat[0])
    lat= float(lon_lat[1])

    child_validAspect = child_locRef.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}validAspect')
    child_AspectPosition = child_validAspect.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}AspectPosition')
    child_positionA = child_AspectPosition.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}position')
    aspect = child_positionA.text

    """Date"""
    child_timeRef = xroot.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}timeRef')
    child_recordTime =  child_timeRef.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}recordTime')
    child_TimeInstant = child_recordTime.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}TimeInstant')
    child_timePosition = child_TimeInstant.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}timePosition')
    date = child_timePosition.text # string

    """Name"""
    child_name = child_locRef.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}name')
    name = child_name.text

    prof_meta = {'lon'     : lon,
                 'lat'     : lat,
                 'alt'     : alt,
                 'datetime': date,
                 'name'    : name,
                 'aspect'  : aspect}

    return prof_meta


def add_monti_density(file, new_file):
    '''
    Adds density profile based on grain shape and hardness parameterization alla Monti (2014)

    Arguments:
        file / new_file (str):          Path to xml file of snowprofile
        remove_observed_density (bool): Relevant when density is already provided from the observations (already included in .caaml file)

    Returns nothing, but saves new profile

    Comments:
        - If density exists from observation or simulation it is removed
    '''
    
    tree = et.parse(file)
    xroot = tree.getroot()
    
    """Add density section"""
    child_snowProfileResultsOf = xroot.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}snowProfileResultsOf')
    child_SnowProfileMeasurements = child_snowProfileResultsOf.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}SnowProfileMeasurements')
    try:
        child_densityProfile = child_SnowProfileMeasurements.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}densityProfile')
        child_SnowProfileMeasurements.remove(child_densityProfile)
    except:
        pass
    child_densityProfile = et.SubElement(child_SnowProfileMeasurements, '{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}densityProfile')
    child_densityMetaData = et.SubElement(child_densityProfile,'{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}densityMetaData')
    child_methodOfMeas = et.SubElement(child_densityMetaData,'{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}methodOfMeas')
    child_methodOfMeas.text = 'other'
    
    """Add density layer for each steatigraphy layer"""
    child_stratProfile = child_SnowProfileMeasurements.find('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}stratProfile')
    for layer in child_stratProfile.iter('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}Layer'):
        child_densityLayer = et.SubElement(child_densityProfile,'{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}Layer')
        child_depthTop  = et.SubElement(child_densityLayer,'{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}depthTop')
        child_thickness = et.SubElement(child_densityLayer,'{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}thickness')
        child_density   = et.SubElement(child_densityLayer,'{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}density')

        child_depthTop.set('uom', 'cm')
        child_thickness.set('uom', 'cm')
        child_density.set('uom', 'kgm-3')
    
        for val_depth in layer.iter('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}depthTop'):
            child_depthTop.text = val_depth.text
        for val_thickness in layer.iter('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}thickness'):
            child_thickness.text = val_thickness.text
        
        """Get grain forms and hardness to estimate density"""
        for form in layer.iter('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}grainFormPrimary'):
            grainFormPrimary = form.text
        for form in layer.iter('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}grainFormSecondary'):
            grainFormSecondary = form.text
        for val_hardness in layer.iter('{http://caaml.org/Schemas/SnowProfileIACS/v6.0.3}hardness'):
            hardness = val_hardness.text
        
        child_density.text = str(_monti_density(grainFormPrimary, grainFormSecondary, hardness))
    
    tree.write(new_file,encoding='UTF-8', method="xml", doctype='<?xml version="1.0" encoding="UTF-8"?>')


def _monti_density(form, form2, hardness_str):
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
    
    dict_hardness = {'F-' : 1,
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
    
    dict_hardness_i =  {'F-' : 0,
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
                   'RG'  : [63,126,219,248,308,338,394,418,444,np.nan,np.nan],
                   'FC'  : [82,164,272,295,347,373,440,479,519,np.nan,np.nan],
                   'DH'  : [96,192,307,325,345,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                   'MF'  : [73,146,248,283,348,377,518,629,740,np.nan,np.nan],
                   'FCxr': [86,172,282,304,350,373,426,455,485,np.nan,np.nan],
                   'IF'  : [ice,ice,ice,ice,ice,ice,ice,ice,ice,ice,ice],
                   'SH'  : [sho,sho,sho,sho,sho,sho,sho,sho,sho,sho,sho],
                   'MFcr': [mfc,mfc,mfc,mfc,mfc,mfc,mfc,mfc,mfc,mfc,mfc]}
    
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


if __name__ == "__main__":
    print("This script provides useful functions to process CAAMLv6 snow profiles")

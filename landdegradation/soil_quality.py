from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import ee

from landdegradation import stats, GEEIOError
from landdegradation.util import TEImage
from landdegradation.schemas.schemas import BandInfo

def soil_quality(depth, geometry, EXECUTION_ID,logger):
    """
    ===========================================================================================
                 SOIL QUALITY INDEX (CQI)
    ===========================================================================================
    The impact of the soil factor to the process of desertification is determined by the strength of cohesion between soil particles, water retention ability of the soil, soil texture, and structure. The soil quality index (SQI), developed by OSS, that will be used is based on four parameters:

    - parent material
    - soil depth
    - soil texture
    - slope

    The formula used to compute the SQI from the above-mentioned parameters is as shown below:

    SQI =(Parent material* Depth* Texture *slope)^1/4
    """

    logger.debug("Entering soil quality function.")

    # define slope 
    srtm = ee.Image("USGS/SRTMGL1_003")

    slope = ee.Terrain.slope(srtm).clip(geometry)

    # define soil texture 
    texture = ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02").clip(geometry)

    # select texture based on depth 
    if depth < 10:
        texture = texture.select(['b0'])
    elif depth >= 10 or depth < 30:
        texture = texture.select(['b10'])
    elif depth >= 30 or depth < 60:
        texture = texture.select(['b30'])
    elif depth >= 60 or depth < 100:
        texture = texture.select(['b60'])
    elif depth >= 100 or depth < 200:
        texture = texture.select(['b100'])
    else:
        texture = texture.select(['b200'])
    

    # define soil parent material 
    lithology = ee.Image("users/miswagrace/global_lithology").clip(geometry).select('b1')

    # calculate SQI
    image = ee.Image()
    sqi = image.expression('(slope * parent_material * soil_texture * soil_depth) ** (1/4)', {
    'slope':slope,
    'parent_material':lithology,
    'soil_texture':texture,
    'soil_depth':depth
    }).rename('sqi')

    return TEImage(sqi,
        [BandInfo("Soil Quality Index", add_to_map=True, metadata={'depth':depth})])




    

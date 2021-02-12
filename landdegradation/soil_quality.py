from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import ee

from landdegradation import stats, GEEIOError
from landdegradation.util import TEImage
from landdegradation.schemas.schemas import BandInfo


def soil_quality(depth, texture_matrix, geometry, EXECUTION_ID, logger):
    """
    ===========================================================================================
                 SOIL QUALITY INDEX (SQI)
    ===========================================================================================
    The impact of the soil factor to the process of desertification is determined by the strength of cohesion between soil particles, water retention ability of the soil, soil texture, and structure. The soil quality index (SQI), developed by OSS, that will be used is based on four parameters:

    - parent material
    - soil depth (0 to 200cm)
    - soil texture
    - slope
    - rock fragment
    - drainage

    The formula used to compute the SQI from the above-mentioned parameters is as shown below:

    SQI = (parent material*soil depth*soil texture* slope* rock fragment*drainage)^1/6
    """

    logger.debug("Entering soil quality function.")
    srtm = ee.Image("USGS/SRTMGL1_003")

    # ==========================
    # PARENT MATERIAL
    # ==========================

    # parent_material_map = [
    #   [1,2,3,4,5,6,7,8,9, 10, 11, 12, 13, 14, 15],pmaterial_matrix]

    parent_material = ee.Image("users/miswagrace/parent_material_northafrica").clip(geometry)
    parent_material = parent_material.remap([1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 2.0],[1.0, 1.7, 1.7, 1.7, 1.7, 1.7, 2.0])

    # resample parent material to srtm projection
    parent_material = parent_material.reproject(crs=srtm.projection())
    # remap parent material 
    # parent_material = parent_material.remap(parent_material_map[0], parent_material_map[1])

    # ==========================
    # SLOPE
    # ==========================
    slope = ee.Image("users/miswagrace/slope_north_africa").clip(geometry)
    slope = slope.reproject(crs=srtm.projection())
   
    # ==========================
    # TEXTURE
    # ==========================    
    # soil_texture = ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02").clip(geometry)
    soil_texture = ee.Image("users/miswagrace/texture_north_africa").clip(geometry)
    remap_texture = [
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],texture_matrix]

    # remap texture 
    soil_texture_remap = soil_texture.remap(remap_texture[0], remap_texture[1])
    # soil_texture_remap = soil_texture_remap \
    #     .where(soil_texture_remap.eq(1), 1) \
    #     .where(soil_texture_remap.eq(2), 1.2) \
    #     .where(soil_texture_remap.eq(3), 1.6) \
    #     .where(soil_texture_remap.eq(4), 2) \
    #     .rename("Soil Texture")

    # soil_texture_remap = soil_texture_remap.updateMask(soil_texture_remap.neq(-32768))

    # ==========================
    # DEPTH
    # ==========================
    if depth<15 :
        depthIndex = 4
    elif depth>=15 and depth<30:
        depthIndex = 3
    elif depth>=30 and depth<75:
        depthIndex = 2
    elif depth>=75:
        depthIndex = 1
    else:
        logger.debug("Unexpected depth value")

    # ==========================
    # ROCK FRAGMENT
    # ==========================
    # classify rock fragment based on 3 classes
    rock_fragment = ee.Image("users/miswagrace/rock_fragment_north_africa")
    fragmentClass = rock_fragment \
        .where(rock_fragment.lt(20), 2) \
        .where(rock_fragment.gte(20).And(rock_fragment.lte(60)), 1.3) \
        .where(rock_fragment.gt(60), 1) \
        .rename("Rock Fragments")

    # ==========================
    # SOIL DRAINAGE
    # ==========================
 
    drainage = ee.Image("users/miswagrace/drainage_northafrica")
    soil_drainage = drainage.remap([1,2,3,4,5,6,7],[2.0, 2.0, 1.4, 1.2, 1.0, 1.7, 2.0])
    # ==========================
    # SOIL QUALITY INDEX
    # ========================== 
    # compute soil quality index based on datasetes generated above
    img = ee.Image()
    sqi = img.expression('(slope * parent_material * soil_texture * soil_depth * rock_fragment *drainage) ** (1/6)', {
        'slope':slope,
        'parent_material':parent_material,
        'soil_texture':soil_texture_remap,
        'soil_depth':depthIndex,
        'rock_fragment':fragmentClass.clip(geometry),
        'drainage':soil_drainage.clip(geometry)
        }).rename('sqi')  

    # classify output sqi into 3 classes 
    sqi = sqi \
        .where(sqi.lt(1.13), 1) \
        .where(sqi.gte(1.13).And(sqi.lte(1.45)), 2) \
        .where(sqi.gt(1.45), 3)
    

    return TEImage(sqi,
        [BandInfo("Soil Quality Index (cm deep)", add_to_map=True, metadata={'depth':depth})])


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import ee

from landdegradation import stats, GEEIOError
from landdegradation.util import TEImage
from landdegradation.schemas.schemas import BandInfo


def soil_quality(depth, texture_matrix, pmaterial_matrix, geometry, EXECUTION_ID, logger):
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

    # ==========================
    # PARENT MATERIAL
    # ==========================

    parent_material_map = [
      [1,2,3,4,5,6,7,8,9, 10, 11, 12, 13, 14, 15],pmaterial_matrix]

    parent_material = ee.Image("users/miswagrace/tunisia_parent_material")

    # remap parent material 
    parent_material = parent_material.remap(parent_material_map[0], parent_material_map[1])

    # ==========================
    # SLOPE
    # ==========================

    strm = ee.Image("USGS/SRTMGL1_003")
    slope = ee.Terrain.slope(srtm).clip(geometry)

    # convert slope to radians
    slopeRad = slope.multiply(0.0174533)

    # convert slope radians to percentage
    slopePercent = slopeRad.atan().multiply(100)

    # reclassify slope according to ranges
    slopeClass = ee.Image(-32768) \
        .where(slopePercent.lt(6), 1) \
        .where(slopePercent.gte(6).And(slopePercent.lte(18)), 1.2) \
        .where(slopePercent.gte(18).And(slopePercent.lte(35)), 1.5) \
        .where(slopePercent.gt(35), 2) \
        .rename("Slope Class")

    # remove no data values
    slopeClass = slopeClass.where(slopeClass.eq(9999), -32768)
    slopeClass = slopeClass.updateMask(slopeClass.neq(-32768))

    # ==========================
    # TEXTURE
    # ==========================    
    soil_texture = ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02").clip(geometry)
    remap_texture = [
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],texture_matrix]

    # remap texture 
    soil_texture_remap = soil_texture.remap(remap_texture[0], remap_texture[1])
    soil_texture_remap = soil_texture_remap \
        .where(soil_texture_remap.eq(1), 1) \
        .where(soil_texture_remap.eq(2), 1.2) \
        .where(soil_texture_remap.eq(3), 1.6) \
        .where(soil_texture_remap.eq(4), 2) \
        .rename("Soil Texture")

    soil_texture_remap = soil_texture_remap.updateMask(soil_texture_remap.neq(-32768))

    # ==========================
    # DEPTH
    # ==========================
    if depth<15 :
        depthIndex = 1
    elif depth>=15 and depth<30:
        depthIndex = 2
    elif depth>=30 and depth<75:
        depthIndex = 3
    elif depth>=75:
        depthIndex = 4
    else:
        logger.debug("Unexpected depth value")

    # ==========================
    # ROCK FRAGMENT
    # ==========================
    # classify rock fragment based on 3 classes
    rock_fragment = ee.Image("users/miswagrace/rock_fragments")
    fragmentClass = rock_fragment \ 
        .where(rock_fragment.lt(20), 2) \
        .where(rock_fragment.gte(20).And(rock_fragment.lte(60)), 1.3) \
        .where(rock_fragment.gt(60), 1) \
        .rename("Rock Fragments")

    fragmentClass = fragmentClass.unmask(-32768)
    fragmentClass = fragmentClass.where(fragmentClass.eq(-32768), 2)
    
    # ==========================
    # SOIL DRAINAGE
    # ==========================
    # classify rock fragment based on 3 classes
    drainage = ee.Image("users/miswagrace/soil_drainage")
    soil_drainage = ee.Image(-32768) \
        .where(drainage.lte(2), 2) \
        .where(drainage.eq(3), 1.2) \
        .where(drainage.gt(3), 1) \
        .rename("Soil Drainage")

    # ==========================
    # SOIL QUALITY INDEX
    # ========================== 
    # compute soil quality index based on datasetes generated above
    img = ee.Image()
    sqi = img.expression('(slope * parent_material * soil_texture * soil_depth * rock_fragment *drainage) ** (1/6)', {
        'slope':slopeClass,
        'parent_material':parent_material,
        'soil_texture':soil_texture_remap,
        'soil_depth':depthIndex,
        'rock_fragment':fragmentClass.clip(geometry),
        'drainage':soil_drainage.clip(geometry)
        }).rename('sqi')  

    # classify output sqi into 3 classes 
    sqi = sqi
        .where(sqi.lt(1.13), 1) \
        .where(sqi.gte(1.13).And(sqi.lte(1.45)), 2) \
        .where(sqi.gt(1.46), 3)
    
    return TEImage(sqi,
        [BandInfo("Soil Quality Index", add_to_map=True, metadata={'depth':depth})])


# def soil_quality(depth, geometry, EXECUTION_ID,logger):
    # """
    # ===========================================================================================
    #              SOIL QUALITY INDEX (CQI)
    # ===========================================================================================
    # The impact of the soil factor to the process of desertification is determined by the strength of cohesion between soil particles, water retention ability of the soil, soil texture, and structure. The soil quality index (SQI), developed by OSS, that will be used is based on four parameters:

    # - parent material
    # - soil depth (0 to 200cm)
    # - soil texture
    # - slope
    # - rock fragment
    # - drainage

    # The formula used to compute the SQI from the above-mentioned parameters is as shown below:

    # SQI = (parent material*soil depth*soil texture* slope* rock fragment*drainage)^1/6
    # """

    # logger.debug("Entering soil quality function.")

    # # define slope 
    # srtm = ee.Image("USGS/SRTMGL1_003")

    # slope = ee.Terrain.slope(srtm).clip(geometry)

    # # define soil texture 
    # texture = ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02").clip(geometry)

    # # select texture based on depth 
    # if depth < 10:
    #     texture = texture.select(['b0'])
    # elif depth >= 10 or depth < 30:
    #     texture = texture.select(['b10'])
    # elif depth >= 30 or depth < 60:
    #     texture = texture.select(['b30'])
    # elif depth >= 60 or depth < 100:
    #     texture = texture.select(['b60'])
    # elif depth >= 100 or depth < 200:
    #     texture = texture.select(['b100'])
    # else:
    #     texture = texture.select(['b200'])
    

    # # define soil parent material 
    # lithology = ee.Image("users/miswagrace/global_lithology").clip(geometry).select('b1')

    # # calculate SQI
    # image = ee.Image()
    # sqi = image.expression('(slope * parent_material * soil_texture * soil_depth) ** (1/4)', {
    # 'slope':slope,
    # 'parent_material':lithology,
    # 'soil_texture':texture,
    # 'soil_depth':depth
    # }).rename('sqi')

    # return TEImage(sqi,
    #     [BandInfo("Soil Quality Index", add_to_map=True, metadata={'depth':depth})])




    

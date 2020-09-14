from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import ee

from landdegradation import stats, GEEIOError
from landdegradation.util import TEImage
from landdegradation.schemas.schemas import BandInfo

def forest_fire(prefire_start,prefire_end,postfire_start,postfire_end, geometry, platform, logger):
    """
    ===========================================================================================
                 BURN SEVERITY MAPPING USING THE NORMALIZED BURN RATIO (NBR)
    ===========================================================================================
     Normalized Burn Ratio will be applied to imagery from before and after a wild fire. By
     calculating the difference afterwards (dNBR) Burn Severity is derived, showing the spatial
     impact of the disturbance. Imagery used in this process comes from either Sentinel-2 or 
     Landsat 8.
    """

    logger.debug("Entering forest_fire function.")

    # SELECT one of the following:   'L8'  or 'S2' 

    if platform == 'S2' or platform == 's2':
        ImCol = 'COPERNICUS/S2'
        pl = 'Sentinel-2'
    else:
        ImCol = 'LANDSAT/LC08/C01/T1_SR'
        pl = 'Landsat 8'

    logger.debug(ee.String('Data selected for analysis: ').cat(pl))
    logger.debug(ee.String('Fire incident occurred between ').cat(prefire_end).cat(' and ').cat(postfire_start))

    # Location
    area = ee.FeatureCollection(geometry)

    # add image collection 
    imagery = ee.ImageCollection(ImCol)

    prefireImCol = ee.ImageCollection(imagery
        # Filter by dates.
        .filterDate(prefire_start, prefire_end)
        # Filter by location.
        .filterBounds(area))

    postfireImCol = ee.ImageCollection(imagery
        # Filter by dates.
        .filterDate(postfire_start, postfire_end)
        # Filter by location.
        .filterBounds(area))

    logger.debug("Pre-fire Image Collection: "+prefireImCol)
    logger.debug("Post-fire Image Collection: "+postfireImCol)

    def maskS2sr(image):
        # Bits 10 and 11 are clouds and cirrus, respectively.
        cloudBitMask = ee.Number(2).pow(10).int()
        cirrusBitMask = ee.Number(2).pow(11).int()
        # Get the pixel QA band.
        qa = image.select('QA60')
        # All flags should be set to zero, indicating clear conditions.
        mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))
        # Return the masked image, scaled to TOA reflectance, without the QA bands.
        return image.updateMask(mask).copyProperties(image, ["system:time_start"])    

    def maskL8sr(image):
        # Bits 3 and 5 are cloud shadow and cloud, respectively.
        cloudShadowBitMask = 1 << 3
        cloudsBitMask = 1 << 5
        snowBitMask = 1 << 4
        # Get the pixel QA band.
        qa = image.select('pixel_qa')
        # All flags should be set to zero, indicating clear conditions.
        mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(qa.bitwiseAnd(cloudsBitMask).eq(0)).And(qa.bitwiseAnd(snowBitMask).eq(0))
        # Return the masked image, scaled to TOA reflectance, without the QA bands.
        return image.updateMask(mask).select("B[0-9]*").copyProperties(image, ["system:time_start"])

    # Apply platform-specific cloud mask
    if platform == 'S2' or platform == 's2':
        prefire_CM_ImCol = prefireImCol.map(maskS2sr)
        postfire_CM_ImCol = postfireImCol.map(maskS2sr)
    else:
        prefire_CM_ImCol = prefireImCol.map(maskL8sr)
        postfire_CM_ImCol = postfireImCol.map(maskL8sr)

    pre_mos = prefireImCol.mosaic().clip(area)
    post_mos = postfireImCol.mosaic().clip(area)

    pre_cm_mos = prefire_CM_ImCol.mosaic().clip(area)
    post_cm_mos = postfire_CM_ImCol.mosaic().clip(area)

    if platform == 'S2' or platform == 's2':
        preNBR = pre_cm_mos.normalizedDifference(['B8', 'B12'])
        postNBR = post_cm_mos.normalizedDifference(['B8', 'B12'])
    else:
        preNBR = pre_cm_mos.normalizedDifference(['B5', 'B7'])
        postNBR = post_cm_mos.normalizedDifference(['B5', 'B7'])

    # The result is called delta NBR or dNBR
    dNBR_unscaled = preNBR.subtract(postNBR)

    # Scale product to USGS standards
    dNBR = dNBR_unscaled.multiply(1000)

    return TEImage(pre_mos.addBands(post_mos).addBands(dNBR),
                    [BandInfo("Pre-fire image", add_to_map=True, metadata={'prefire_start':prefire_start,'prefire_end':prefire_end, 'postfire_start':postfire_start, 'postfire_end':postfire_end})]
                    [BandInfo("Post-fire image", add_to_map=True, metadata={'prefire_start':prefire_start,'prefire_end':prefire_end, 'postfire_start':postfire_start, 'postfire_end':postfire_end})]
                    [BandInfo("dNBR classified",  add_to_map=True,metadata={'prefire_start':prefire_start,'prefire_end':prefire_end, 'postfire_start':postfire_start, 'postfire_end':postfire_end})])

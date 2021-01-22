from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import ee

from landdegradation import stats, GEEIOError
from landdegradation.util import TEImage
from landdegradation.schemas.schemas import BandInfo

def vegetation_quality(year, ndvi_start, ndvi_end, drought_matrix, fire_matrix, erosion_matrix, geometry, EXECUTION_ID,logger):

    drought_remap_matrix =  [
        [
            10,11,12,20,30,40,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,
            130,140,150,151,152,153,160,170,180,200,201,202
        ],
        drought_matrix
    ]

    fire_remap_matrix =  [
        [
            10,11,12,20,30,40,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,
            130,140,150,151,152,153,160,170,180,200,201,202
        ],
        fire_matrix
    ]

    erosion_remap_matrix =  [
        [
            10,11,12,20,30,40,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,
            130,140,150,151,152,153,160,170,180,200,201,202
        ],
        erosion_matrix
    ]

    ## land cover
    lc = ee.Image("users/geflanddegradation/toolbox_datasets/lcov_esacc_1992_2018").clip(geometry)
    lc = lc.where(lc.eq(9999), -32768)
    lc = lc.updateMask(lc.neq(-32768))

    # DEFINE DROUGHT RESISTANCE 
    lc_remapped_drought = lc.select('y'+year).remap(drought_remap_matrix[0], drought_remap_matrix[1]).divide(10)

    # DEFINE FIRE RISK
    lc_remapped_fire = lc.select('y'+year).remap(fire_remap_matrix[0], fire_remap_matrix[1]).divide(10)

    # DEFINE EROSION PROTECTION
    lc_remapped_erosion = lc.select('y'+year).remap(erosion_remap_matrix[0], erosion_remap_matrix[1]).divide(10)

    # DEFINE PLANT COVER
    max_ndvi = ee.ImageCollection('VITO/PROBAV/C1/S1_TOC_100M') \
                  .filter(ee.Filter.date(ndvi_start, ndvi_end)) \
                  .max() \
                  .select('NDVI')

    plant_cover_range = max_ndvi.divide(255).multiply(100)

    plant_cover_class = plant_cover_range \
        .where(plant_cover_range.gte(80), 1.0) \
        .where(plant_cover_range.lt(80).And(plant_cover_range.gte(72)), 1.1) \
        .where(plant_cover_range.lt(72).And(plant_cover_range.gte(62)), 1.2) \
        .where(plant_cover_range.lt(62).And(plant_cover_range.gte(50)), 1.3) \
        .where(plant_cover_range.lt(50).And(plant_cover_range.gte(38)), 1.4) \
        .where(plant_cover_range.lt(38).And(plant_cover_range.gte(16)), 1.5) \
        .where(plant_cover_range.lt(26).And(plant_cover_range.gte(18)), 1.6) \
        .where(plant_cover_range.lt(18).And(plant_cover_range.gte(13)), 1.7) \
        .where(plant_cover_range.lt(13).And(plant_cover_range.gte(11)), 1.8) \
        .where(plant_cover_range.lt(11).And(plant_cover_range.gte(10)), 1.9) \
        .where(plant_cover_range.lt(10), 2) \
        .clip(geometry) \
        .rename("Plant Cover")


    # CALCULATE VEGETATION QUALITY INDEX
    vqi = ee.Image()

    vqi = vqi.expression('(fire_risk * erosion_protection * drought_resistance * plant_cover)**(1/4)',{
        'fire_risk':lc_remapped_fire,
        'erosion_protection':lc_remapped_erosion,
        'drought_resistance':lc_remapped_drought,
        'plant_cover':plant_cover_class,
        }
    ).rename("Vegetation Quality Index")

    vqi_range = vqi \
        .where(vqi.lte(1.13), 1) \
        .where(vqi.lte(1.38).And(vqi.gt(1.13)), 2) \
        .where(vqi.gt(1.38), 3)

    return TEImage(vqi_range.clip(geometry),
        [BandInfo("Vegetation Quality Index", add_to_map=True, metadata={'year':year})])
  



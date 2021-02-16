from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import ee

from landdegradation import stats, GEEIOError
from landdegradation.util import TEImage
from landdegradation.schemas.schemas import BandInfo

def climate_quality(year, geometry, EXECUTION_ID,logger):
    """
    ===========================================================================================
                 CLIMATE QUALITY INDEX (CQI)
    ===========================================================================================

    """
    logger.debug("Entering climate quality function.")

    terra_climate = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") \
                .filter(ee.Filter.date('{}-01-01'.format(year), '{}-12-31'.format(year))) \
                .mean() \
                .clip(geometry)

    precipitation = terra_climate.select('pr')

    precipitation_class = precipitation \
        .where(precipitation.gt(650), 1) \
        .where(precipitation.gte(570).And(precipitation.lt(650)), 1.05) \
        .where(precipitation.gte(490).And(precipitation.lt(570)), 1.15) \
        .where(precipitation.gte(440).And(precipitation.lt(490)), 1.25) \
        .where(precipitation.gte(390).And(precipitation.lt(440)), 1.35) \
        .where(precipitation.gte(345).And(precipitation.lt(390)), 1.50) \
        .where(precipitation.gte(310).And(precipitation.lt(345)), 1.65) \
        .where(precipitation.gte(280).And(precipitation.lt(310)), 1.80) \
        .where(precipitation.lt(280), 2) \
        .rename("Precipitation")

    evapotrans = terra_climate.select('pet')

    aridityIndex = ee.Image().expression('(precipitation / evapotranspiration)', {
    'precipitation':precipitation,
    'evapotranspiration':evapotrans,
    }).rename('aridityIndex')

    aridityIndex = aridityIndex \
        .where(aridityIndex.gte(1.0), 1) \
        .where(aridityIndex.gte(0.75).And(aridityIndex.lt(1.0)), 1.05) \
        .where(aridityIndex.gte(0.65).And(aridityIndex.lt(0.75)), 1.15) \
        .where(aridityIndex.gte(0.50).And(aridityIndex.lt(0.65)), 1.25) \
        .where(aridityIndex.gte(0.35).And(aridityIndex.lt(0.50)), 1.35) \
        .where(aridityIndex.gte(0.20).And(aridityIndex.lt(0.35)), 1.45) \
        .where(aridityIndex.gte(0.10).And(aridityIndex.lt(0.20)), 1.55) \
        .where(aridityIndex.gte(0.03).And(aridityIndex.lt(0.10)), 1.75) \
        .where(aridityIndex.lt(0.03), 2) \
        .rename("Aridity Index")


    cqi = ee.Image().expression('(precipitation * aridity_index) ** 1/2', {
        'precipitation':precipitation_class,
        'aridity_index':aridityIndex,
        })

    cqi_class = cqi \
        .where(cqi.lt(1.15), 1) \
        .where(cqi.gte(1.15).And(cqi.lte(1.81)), 2) \
        .where(cqi.gt(1.81), 3) \
        .rename('Climate Quality Reclass')

    srtm_proj = ee.Image("USGS/SRTMGL1_003").projection()

    # resample parent material to srtm projection
    cqi_class = cqi_class.reproject(crs=srtm_proj)
    
    return TEImage(cqi_class.clip(geometry),
        [BandInfo("Climate Quality Index (year)", add_to_map=True, metadata={'year':year})]
    )

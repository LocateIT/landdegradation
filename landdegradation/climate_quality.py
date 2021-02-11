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



# def climate_quality(month,next_month, geometry, EXECUTION_ID,logger):
#     """
#     ===========================================================================================
#                  CLIMATE QUALITY INDEX (CQI)
#     ===========================================================================================
#     Adverse climate conditions such as recurrent and prolonged drought, increase the susceptibility of land to desertification. 
#     The Climate Quality Index(CQI) is analyzed using the following parameters: rainfall, aridity index(1), and aspect, using the formula: CQI = (rainfall*aridity*aspect)^1/3
#     """

#     logger.debug("Entering climate quality function.")

#     climate = ee.Image('WORLDCLIM/V1/BIO').clip(geometry)
#     annualMeanRainfall = climate.select('bio12')

#     # reclassify rainfall values into 3 classes
#     rainfallClass =ee.Image(-32768) \
#         .where(annualMeanRainfall.gt(650), 1) \
#         .where(annualMeanRainfall.gte(280).And(annualMeanRainfall.lte(650)), 2) \
#         .where(annualMeanRainfall.lt(280), 4) \
#         .rename("Annual Rainfall")

#     # remove no data values
#     rainfallClass = rainfallClass.where(rainfallClass.eq(9999), -32768)
#     rainfallClass = rainfallClass.updateMask(rainfallClass.neq(-32768))

#     # // --------------------------------------------------------------
#     # // DEFINE ASPECT in 2 classes
#     # // --------------------------------------------------------------

#     # aspect(exposure)
#     srtm = ee.Image("USGS/SRTMGL1_003").clip(geometry)

#     logger.debug("Calculating field orientation")
#     # // converted to millimeteres
#     aspect = ee.Terrain.aspect(srtm)
#     fieldOrientation = ee.Image(-32768) \
#         .where(aspect.gte(0).And(aspect.lt(112.5)), 1) \
#         .where(aspect.gte(247.5).And(aspect.lt(360)), 2) \
#         .where(aspect.gte(112.5).And(aspect.lt(247.5)), 3)

#     # // combine fieldOrientation class 1 and 2
#     fieldOrientation = fieldOrientation \
#         .where(fieldOrientation.eq(1).And(fieldOrientation.eq(2)), 1) \
#         .where(fieldOrientation.eq(3), 2) \
#         .rename("Field Orientation")

#     # // remove no data values
#     # // fieldOrientation = fieldOrientation.where(fieldOrientation.eq(9999), -32768)
#     fieldOrientation = fieldOrientation.updateMask(fieldOrientation.neq(-32768))

#     ecmwf = ee.ImageCollection("ECMWF/ERA5/MONTHLY") \
#         .filter(ee.Filter.date('{}-01'.format(month), '{}-01'.format(next_month))) \
#         .first() \
#         .clip(geometry)

#     aridity_index =  ecmwf.expression('(2*mean_2m_air_temperature - total_precipitation)* 1',{
#         'mean_2m_air_temperature':ecmwf.select(['mean_2m_air_temperature']).subtract(273.15),
#         'total_precipitation':ecmwf.select(['total_precipitation']).multiply(1000)
#         }
#     )

#     # reclassify rainfall values into 3 classes
#     aridityClass =ee.Image(-32768) \
#         .where(aridity_index.lt(50), 1) \
#         .where(aridity_index.gte(50).And(aridity_index.lt(75)), 1.1) \
#         .where(aridity_index.gte(75).And(aridity_index.lt(100)), 1.2) \
#         .where(aridity_index.gte(100).And(aridity_index.lt(125)), 1.4) \
#         .where(aridity_index.gte(125).And(aridity_index.lt(150)), 1.8) \
#         .where(aridity_index.gte(150), 2) \
#         .rename("Aridity Index")
  
#     # remove no data values
#     aridityClass = aridityClass.where(aridityClass.eq(9999), -32768)
#     aridityClass = aridityClass.updateMask(aridityClass.neq(-32768))

#     # // -------------------------------------------------------------------
#     # // DEFINE CLIMATE QUALITY INDEX in 3 classes
#     # // -------------------------------------------------------------------
#     logger.debug("Calculating climate quality")

#     cqi = ee.Image()
#     cqi = cqi.expression('(rainfall * aspect * aridity_index) ** (1/3)', {
#         'rainfall':rainfallClass,
#         'aspect':fieldOrientation,
#         'aridity_index':aridityClass,
#         }).rename('cqi')

#     # classify climate quality values based on ranges
#     cqi = cqi \
#         .where(cqi.lt(1.15), 1) \
#         .where(cqi.gte(1.15).And(cqi.lte(1.81)), 2) \
#         .where(cqi.gt(1.81), 3) \
#         .rename('Climate Quality')

#     # remove no data values
#     cqi = cqi.where(cqi.eq(9999), -32768)
#     cqi = cqi.updateMask(cqi.neq(-32768))
    
#     return TEImage(cqi.addBands(rainfallClass).addBands(fieldOrientation).addBands(aridityClass).clip(geometry),
#         [BandInfo("Climate Quality Index (month)", add_to_map=True, metadata={'month':month}),
#         BandInfo("Rainfall", add_to_map=True, metadata={'month':month}),
#         BandInfo("Field orientation", add_to_map=True, metadata={'month':month}),
#         BandInfo("Aridity", add_to_map=True, metadata={'month':month})])


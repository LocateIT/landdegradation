from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import ee

from landdegradation import stats, GEEIOError
from landdegradation.util import TEImage
from landdegradation.schemas.schemas import BandInfo

def management_quality(year, lu_matrix, geometry, EXECUTION_ID,logger):
    logger.debug("Entering management quality function.")

    lu_remap_matrix =  [
        [
            10,11,12,20,30,40,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,
            130,140,150,151,152,153,160,170,180,200,201,202
        ],
        lu_matrix
    ]

    ## land cover
    lc = ee.Image("users/geflanddegradation/toolbox_datasets/lcov_esacc_1992_2018").clip(geometry)
    lc = lc.where(lc.eq(9999), -32768)
    lc = lc.updateMask(lc.neq(-32768))

    # DEFINE LAND USE INTENSITY
    # Remap LC according to input matrix
    lc_remapped_lu = lc.select('y'+year).remap(lu_remap_matrix[0], lu_remap_matrix[1]).divide(10)

    # DEFINE POPULATION INTENSITY
    population_density = ee.ImageCollection("CIESIN/GPWv411/GPW_Population_Density").toList(5)

    if(year <= 2000):
        population_density = ee.Image(population_density.get(0))
    elif(year > 2000 and year<= 2005):
        population_density = ee.Image(population_density.get(1))
    elif(year > 2005 and year<= 2010):
        population_density = ee.Image(population_density.get(2))
    elif(year > 2010 and year<= 2015):
        population_density = ee.Image(population_density.get(3))
    elif(year > 2015):
        population_density = ee.Image(population_density.get(4))
    else:
        logger.debug("Invalid date range")

    population_density_class = population_density \
        .select("population_density") \
        .where(population_density.lt(4), 1.0) \
        .where(population_density.gte(4).And(population_density.lt(30)), 1.1) \
        .where(population_density.gte(30).And(population_density.lt(80)), 1.2) \
        .where(population_density.gte(80).And(population_density.lt(170)), 1.3) \
        .where(population_density.gte(170).And(population_density.lt(300)), 1.4) \
        .where(population_density.gte(300).And(population_density.lt(500)), 1.5) \
        .where(population_density.gte(500).And(population_density.lt(850)), 1.6) \
        .where(population_density.gte(850).And(population_density.lt(1400)), 1.7) \
        .where(population_density.gte(1400).And(population_density.lt(2000)), 1.8) \
        .where(population_density.gte(2000).And(population_density.lt(2700)), 1.9) \
        .where(population_density.gte(2700), 2.0) \
        .clip(geometry)

    mqi = ee.Image()
    mqi = mqi.expression('(landuse_intensity * population_density)**(1/2)',{
        'landuse_intensity':lc_remapped_lu,
        'population_density':population_density_class
        }
    ).rename("Management Quality Index")

    mqi_range = mqi \
        .where(mqi.lte(1.25), 1) \
        .where(mqi.lte(1.50).And(mqi.gt(1.25)), 2) \
        .where(mqi.gt(1.50), 3)

    return TEImage(mqi_range.clip(geometry),
        [BandInfo("Management Quality Index", add_to_map=True, metadata={'year':year})])


    





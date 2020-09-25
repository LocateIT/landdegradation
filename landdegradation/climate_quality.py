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
    Adverse climate conditions such as recurrent and prolonged drought, increase the susceptibility of land to desertification. 
    The Climate Quality Index(CQI) is analyzed using the aridity index (AI), using the formula: AI = P/PET
    Where AI is the Aridity index; P is the yearly mean precipitation; and PET is the mean potential evapotranspiration
    """

    logger.debug("Entering climate quality function.")
    # params include AOI, year, 
    terra_climate = (ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
                  .filter(ee.Filter.date('{year}-01-01', '{year}-12-31'))
                 )

    # calculate yearly mean
    terra_climate_mean = terra_climate.mean()

    # retrieve yearly mean precipitation and mean potential evapotranspiration
    prec = terra_climate_mean.select('pr')
    pet = terra_climate_mean.select('pet')

    logger.debug("Calculating Aridity Index.")
    # calculate Aridity index (p/pet) results range from 0 to 1
    aridityIndex = prec.divide(pet)

    return TEImage(aridityIndex.clip(geometry),
        [BandInfo("Aridity Index", add_to_map=True, metadata={'year':year)])




    

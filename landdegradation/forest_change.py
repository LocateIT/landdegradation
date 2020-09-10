from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import ee

from landdegradation import stats, GEEIOError
from landdegradation.util import TEImage
from landdegradation.schemas.schemas import BandInfo

# This dataset is updated yearly, so we get the latest version.
gfc2019 = ee.Image("UMD/hansen/global_forest_change_2019_v1_7")

"""
Calulate forest loss for a given period over a given region
"""
def forest_loss( year, geojson,EXECUTION_ID, logger):
    logger.debug("Entering Forest Loss function.")
    
    # Make sure the bounding box of the poly is used, and not the geodesic 
    # version, for the clipping
    poly = ee.Geometry(geojson, opt_geodesic=False)
    
    lossYear = gfc2019.select(['lossyear']).eq(int(str(year)[2:]))

    lossImageAOI = lossYear.clip(poly)

    # calculate area in metres 
    lossAreaImage = lossImageAOI.multiply(ee.Image.pixelArea())

    return TEImage(lossAreaImage.updateMask(lossAreaImage), 
                    [BandInfo("Forest Loss in {0}".format(year), add_to_map=True, metadata={year:year})])

"""
Calulate forest gain for a given period over a given region
"""
def forest_gain(year_start, year_end, ndvi_gee_dataset, geojson,EXECUTION_ID, logger):
    logger.debug("Entering Forest Gain function.")
    gainImage = gfc2019.select(['gain'])



"""
Calulate forest cover for a given period over a given region
"""
def forest_cover(year_start, year_end, ndvi_gee_dataset, geojson,EXECUTION_ID, logger):
    logger.debug("Entering Forest Cover function.")
    treeCover = gfc2019.select(['treecover2000'])


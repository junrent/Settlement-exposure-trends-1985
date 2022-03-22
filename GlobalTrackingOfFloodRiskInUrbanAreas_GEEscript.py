# Given:
# i) the WSF evolution dataset (either outlining the yearly or the 5-year step global settlement extent) stored as Google Earth Engine (GEE) Image Collection asset,
# ii) a reference global flood hazard map (e.g., Fathom Pluvial/Fluvial Undefended, JRC Coastal) stored as GEE Image Collection asset,
# iii) a reference shapefile outlining the study regions of interests (e.g., ADM0, ADM1, ...) stored as GEE Feature Collection asset,
# this script allows to compute the total settlement extent for different classes of interests defined on the basis of expected water depth.

import ee
import sys

ee.Initialize()

WSFevolution = ee.ImageCollection('XXX').mosaic(); # <- Replace XXX with the path to the WSF evolution Image Collection asset
HazardMap = ee.ImageCollection('XXX').mosaic();    # <- Replace XXX with the path to the flood hazard Image Collection asse
HazardMap = HazardMap.unmask(0);                   # <- depending on the specific flood hazard map, it might be necessary to unmask nodata and/or set them to 0
ReferenceSHP = ee.FeatureCollection('XXX');        # <- Replace XXX with the path to the reference shapefile Feature Collection asset



# ***********************************************
# Simple routine to associate each feature in the reference shapefile with a specific ID number; this will identify the corresponding output CSV file (see further below);
# please note that any attribute in the original reference shapefile will also be reported in the output CSV (hence including any name or code that allows to univocally identify a given feature)
indexes = ee.List(ReferenceSHP.aggregate_array('system:index'));
ids = ee.List.sequence(1, indexes.size());
idByIndex = ee.Dictionary.fromLists(indexes, ids);

def setId (feature): \
  return feature.set('id', idByIndex.get(feature.get('system:index')));

ReferenceSHP_withID = ReferenceSHP.map(setId);
# ***********************************************



# ***********************************************
# Defining the spatial masks outlining the different classes of interest on the basis of the expected water depth  
mask_0 = HazardMap.eq(0); # {0}
mask_015 = HazardMap.lte(0.15).And(HazardMap.gt(0)); # (0÷0.15m]
mask_050 = HazardMap.lte(0.5).And(HazardMap.gt(0.15)); # (0.15÷0.50m]
mask_150 = HazardMap.lte(1.5).And(HazardMap.gt(0.5)); # (0.5÷1.50m]
mask_150p = HazardMap.gt(1.5); # (1.50m÷+∞)
# ***********************************************



# ***********************************************
# Function computing - for the given reference year and feature under analysis  - the total settlement extent (SE) for the different classes of interests defined above
def computeSettlementExtent (year):
  def computeSettlementExtentFeature (selectedFeature): \
    WSFevolution_subset = WSFevolution.clip(selectedFeature.geometry()); \
    WSFevolution_subset = WSFevolution_subset.updateMask(WSFevolution_subset.lte(year)); \
    WSFevolution_subset = WSFevolution_subset.gt(0); \
    WSFevolution_subset_0 = WSFevolution_subset.updateMask(mask_0.eq(1)); \
    WSFevolution_subset_015 = WSFevolution_subset.updateMask(mask_015.eq(1)); \
    WSFevolution_subset_050 = WSFevolution_subset.updateMask(mask_050.eq(1)); \
    WSFevolution_subset_150 = WSFevolution_subset.updateMask(mask_150.eq(1)); \
    WSFevolution_subset_150p = WSFevolution_subset.updateMask(mask_150p.eq(1)); \
     \
    settlementExtent = ee.Image.pixelArea().multiply(0.000001).multiply(WSFevolution_subset.unmask(0)) \
    .reduceRegion(reducer=ee.Reducer.sum(), geometry=selectedFeature.geometry(), bestEffort=True, maxPixels=1e18, tileScale=1, scale=30).get('area'); \
     \
    settlementExtent_0 = ee.Image.pixelArea().multiply(0.000001).multiply(WSFevolution_subset_0.unmask(0)) \
    .reduceRegion(reducer=ee.Reducer.sum(), geometry=selectedFeature.geometry(), bestEffort=True, maxPixels=1e18, tileScale=1, scale=30).get('area'); \
     \
    settlementExtent_015 = ee.Image.pixelArea().multiply(0.000001).multiply(WSFevolution_subset_015.unmask(0)) \
    .reduceRegion(reducer=ee.Reducer.sum(), geometry=selectedFeature.geometry(), bestEffort=True, maxPixels=1e18, tileScale=1, scale=30).get('area'); \
     \
    settlementExtent_050 = ee.Image.pixelArea().multiply(0.000001).multiply(WSFevolution_subset_050.unmask(0)) \
    .reduceRegion(reducer=ee.Reducer.sum(), geometry=selectedFeature.geometry(), bestEffort=True, maxPixels=1e18, tileScale=1, scale=30).get('area'); \
     \
    settlementExtent_150 = ee.Image.pixelArea().multiply(0.000001).multiply(WSFevolution_subset_150.unmask(0)) \
    .reduceRegion(reducer=ee.Reducer.sum(), geometry=selectedFeature.geometry(), bestEffort=True, maxPixels=1e18, tileScale=1, scale=30).get('area'); \
     \
    settlementExtent_150p = ee.Image.pixelArea().multiply(0.000001).multiply(WSFevolution_subset_150p.unmask(0)) \
    .reduceRegion(reducer=ee.Reducer.sum(), geometry=selectedFeature.geometry(), bestEffort=True, maxPixels=1e18, tileScale=1, scale=30).get('area'); \
     \
    selectedFeature = selectedFeature.set('SE_' + str(year), settlementExtent); \
    selectedFeature = selectedFeature.set('SE_0_' + str(year), settlementExtent_0); \
    selectedFeature = selectedFeature.set('SE_015_' + str(year), settlementExtent_015); \
    selectedFeature = selectedFeature.set('SE_050_' + str(year), settlementExtent_050); \
    selectedFeature = selectedFeature.set('SE_150_' + str(year), settlementExtent_150); \
    selectedFeature = selectedFeature.set('SE_150p_' + str(year), settlementExtent_150p); \
    return selectedFeature;
  return computeSettlementExtentFeature;
# ***********************************************



# ***********************************************
# The for-loop below allows to export for each feature in the reference shapefile a CSV file to the 'Results' folder of the Google Drive containing - in addition to the original attributes - the total settlement extent (SE)
# for the different classes of interests on a yerly basis or at 5-year steps depending on the specific version of the WSF evolution used (this proved to be far more effective than processing the whole reference shapefile at once).
#
for j in range (0,ReferenceSHP_withID.size().getInfo(),1): # <- GEE allows at maximum 3000 tasks in the queue. In case the number of features is greater, one shall manually set the range of features to be processed [e.g., (1,3001,1); (3001,6001,1), etc.)
  AOI = ReferenceSHP_withID.filterMetadata('id','equals', j+1);
  countriesEdited = AOI.map(computeSettlementExtent(1985));
  for i in range (1986,2016,1): # <- use this in the case of the yearly WSF evolution
  #for i in range (1990,2016,5): # <- use this in the case of the 5-year WSF evolution
    countriesEdited = countriesEdited.map(computeSettlementExtent(i));
  task = ee.batch.Export.table.toDrive(collection = countriesEdited.select([".*"], None, False),
                                       description = 'WSFevolution_ZonalStatistics_' + str(j+1),
                                       folder = 'Results',
                                      )
  task.start()
# ***********************************************


 


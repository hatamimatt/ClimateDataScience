// ============================================================
// 1. ASSET LIST & PROJECT SETTINGS
// ============================================================
var us_se = ee.Geometry.Rectangle([-94, 24, -75, 38]);
var hist_range = ['1981-01-01', '2010-12-31'];
var fut_range = ['2071-01-01', '2100-12-31'];

// Load all 4 thresholds as an Image Stack
var t10  = ee.Image("projects/ptscalingproject/assets/se10yr24ha");
var t25  = ee.Image("projects/ptscalingproject/assets/se25yr24ha");
var t50  = ee.Image("projects/ptscalingproject/assets/se50yr24ha");
var t100 = ee.Image("projects/ptscalingproject/assets/se100yr24ha");

var atlasStack = ee.Image.cat([t10, t25, t50, t100])
                   .rename(['t10', 't25', 't50', 't100'])
                   .clip(us_se)
                   .divide(1000).multiply(25.4);

// ============================================================
// 2. MULTI-THRESHOLD ANALYSIS FUNCTION
// ============================================================
var runMultiThreshold = function(modelName, scenario) {
  var cmip6 = ee.ImageCollection("NASA/GDDP-CMIP6")
                .filterBounds(us_se)
                .filter(ee.Filter.eq('model', modelName))
                .select('pr');

  var toDaily = function(img) { return img.multiply(86400).copyProperties(img, ['system:time_start']); };
  var hist = cmip6.filter(ee.Filter.eq('scenario', 'historical')).filterDate(hist_range[0], hist_range[1]).map(toDaily);
  var fut  = cmip6.filter(ee.Filter.eq('scenario', scenario)).filterDate(fut_range[0], fut_range[1]).map(toDaily);

  // Bias Correction (Using the 10yr threshold as the baseline anchor)
  var annMax = ee.ImageCollection.fromImages(ee.List.sequence(1981, 2010).map(function(y){
    return hist.filterDate(ee.Date.fromYMD(y,1,1).advance(0, 'year'), ee.Date.fromYMD(y,1,1).advance(1, 'year')).max();
  }));
  var q90 = annMax.reduce(ee.Reducer.percentile([90]));
  
  // Note: For multi-threshold, you apply the specific CF for each threshold level
  var cfStack = atlasStack.divide(q90).clamp(0.5, 5.0);

  // Calculate Exceedances for all 4 bands at once
  var hEx = hist.map(function(img){ return img.multiply(cfStack).gt(atlasStack); }).sum().divide(30);
  var fEx = fut.map(function(img){ return img.multiply(cfStack).gt(atlasStack); }).sum().divide(30);

  return fEx.subtract(hEx).set({
    'model': modelName,
    'scenario': scenario,
    'desc': '4-Band Exceedance Change'
  });
};

// ============================================================
// 3. EXECUTION
// ============================================================
var scenario = 'ssp585';
var models = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CESM2', 'CESM2-WACCM', 
              'CMCC-CM2-SR5', 'CMCC-ESM2', 'CNRM-CM6-1', 'CNRM-ESM2-1', 'CanESM5', 
              'EC-Earth3', 'EC-Earth3-Veg-LR', 'FGOALS-g3', 'GFDL-CM4', 'GFDL-ESM4', 
              'GISS-E2-1-G', 'HadGEM3-GC31-LL', 'IITM-ESM', 'INM-CM4-8', 'INM-CM5-0', 
              'IPSL-CM6A-LR', 'KACE-1-0-G', 'KIOST-ESM', 'MIROC-ES2L', 'MIROC6', 
              'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NESM3', 'NorESM2-LM', 
              'NorESM2-MM', 'TaiESM1', 'UKESM1-0-LL'];

models.forEach(function(m) {
  var result = runMultiThreshold(m, scenario);
  Export.image.toDrive({
    image: result,
    description: m + '_' + scenario + '_MultiThreshold',
    folder: 'CMIP6_Multi_Band_Results',
    scale: 25000,
    region: us_se,
    crs: 'EPSG:4326'
  });
});
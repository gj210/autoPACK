#include as follow : execfile('pathto/.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import MultiCylindersIngr
LPO128= MultiCylindersIngr( 
packingMode = 'random',
color = [1, 0, 0],
radii = [[3.1682549560422122, 3.1682549560422122, 2.9334277157949957]],
cutoff_boundary = 0,
positions2 = [[[4.6484614885770359, 3.3700000231082621, -9.1453846153846161], [-2.9276923216306248, -2.285384613734025, 1.3961538461538463], [3.1066666692495346, 8.5324999888738002, -14.699999999999998]]],
Type = 'MultiCylinder',
useLength = False,
gradient = '',
jitterMax = (0.20000000000000001, 0.10000000000000001, 0.20000000000000001),
packingPriority = 0,
nbJitter = 5,
molarity = 1.0,
meshFile = 'http://autofill.googlecode.com/svn/data/Lipids/recipe/LPO128.c4d',#google
perturbAxisAmplitude = 0.1,
principalVector = [-0.29690364003181458, 0.76161420345306396, -0.57601392269134521],
name = 'LPO128',
positions = [[[1.7900000174840291, 1.2399999935179948, -2.0233333333333334], [1.7900000174840291, 1.2399999935179948, -2.0233333333333334], [4.6484614885770359, 3.3700000231082621, -9.1453846153846161]]],
placeType = 'jitter',
cutoff_surface = 0,
nbMol = 0,
)
recipe.addIngredient(LPO128)

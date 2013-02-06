#include as follow : execfile('pathto/LPO118.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO118= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO118.sph',
radii = [[1.97, 2.8399999999999999, 1.22, 3.1400000000000001, 3.3900000000000001, 4.0599999999999996, 2.52, 3.4399999999999999]],
cutoff_boundary = 0,
Type = 'MultiSphere',
cutoff_surface = 0,
gradient = '',
jitterMax = [0.5, 0.5, 0.10000000000000001],
packingPriority = 0,
rotAxis = [0.0, 2.0, 1.0],
nbJitter = 5,
molarity = 1.0,
rotRange = 6.2831,
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO118.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO118',
positions = [[(-4.4500000000000002, -1.6799999999999999, -6.5199999999999996), (-4.9500000000000002, -0.84999999999999998, -11.970000000000001), (-3.2400000000000002, -3.2999999999999998, -2.7200000000000002), (5.0999999999999996, 4.2400000000000002, -8.8499999999999996), (2.3300000000000001, 0.91000000000000003, -20.010000000000002), (1.1799999999999999, -0.97999999999999998, -24.579999999999998), (3.1299999999999999, 4.1100000000000003, -14.9), (-3.2999999999999998, -2.29, -18.170000000000002)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO118)

#include as follow : execfile('pathto/POP58.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP58= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP58.sph',
radii = [[1.55, 3.0099999999999998, 1.8799999999999999, 2.8300000000000001, 3.2200000000000002, 4.2999999999999998, 3.1800000000000002, 2.9300000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP58.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP58',
positions = [[(-3.0, -1.4299999999999999, 24.440000000000001), (0.050000000000000003, -4.8799999999999999, 12.26), (-1.0900000000000001, 1.3200000000000001, 23.34), (0.76000000000000001, -3.4399999999999999, 18.690000000000001), (0.94999999999999996, 3.5899999999999999, 11.83), (0.84999999999999998, 4.1799999999999997, 3.7200000000000002), (0.40000000000000002, -2.27, 5.46), (0.82999999999999996, 1.1100000000000001, 18.710000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP58)

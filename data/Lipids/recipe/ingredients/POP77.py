#include as follow : execfile('pathto/POP77.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP77= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP77.sph',
radii = [[3.0499999999999998, 4.7000000000000002, 1.9199999999999999, 1.9399999999999999, 2.3199999999999998, 0.77000000000000002, 4.2999999999999998, 2.6699999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP77.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP77',
positions = [[(-0.45000000000000001, 1.04, -18.100000000000001), (-4.8499999999999996, -0.95999999999999996, -5.1699999999999999), (10.52, 0.23999999999999999, -7.0800000000000001), (3.0299999999999998, 1.55, -12.970000000000001), (6.0700000000000003, 1.3100000000000001, -8.8100000000000005), (11.460000000000001, -1.1000000000000001, -4.1500000000000004), (-5.2699999999999996, -0.22, -14.949999999999999), (-1.1599999999999999, -0.82999999999999996, -23.600000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP77)

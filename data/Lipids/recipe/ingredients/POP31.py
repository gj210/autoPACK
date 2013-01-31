#include as follow : execfile('pathto/POP31.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP31= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP31.sph',
radii = [[4.5999999999999996, 4.8700000000000001, 3.7200000000000002, 6.0800000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP31.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP31',
positions = [[(-8.4399999999999995, -0.01, -8.5600000000000005), (-2.1699999999999999, -2.0499999999999998, -16.789999999999999), (-2.6000000000000001, -2.6600000000000001, -9.0399999999999991), (-1.3700000000000001, 1.1699999999999999, -23.760000000000002)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP31)

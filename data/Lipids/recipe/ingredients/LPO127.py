#include as follow : execfile('pathto/LPO127.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO127= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO127.sph',
radii = [[2.46, 1.7, 3.6299999999999999, 3.7599999999999998, 2.3799999999999999, 2.52, 2.5600000000000001, 4.0899999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO127.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO127',
positions = [[(3.52, -3.46, 3.2200000000000002), (0.45000000000000001, -5.7000000000000002, 24.199999999999999), (-1.5600000000000001, -1.6599999999999999, 13.1), (1.6200000000000001, 5.9299999999999997, 10.49), (-0.14000000000000001, -2.4100000000000001, 22.579999999999998), (-2.54, 0.40000000000000002, 19.0), (0.11, 4.9100000000000001, 16.960000000000001), (0.01, 0.68000000000000005, 6.2400000000000002)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO127)

#include as follow : execfile('pathto/POP17.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP17= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP17.sph',
radii = [[2.6800000000000002, 3.7599999999999998, 3.3999999999999999, 1.5, 4.0300000000000002, 1.51, 3.2799999999999998, 1.3]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP17.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP17',
positions = [[(3.21, 1.4099999999999999, 13.57), (3.8100000000000001, 3.3500000000000001, 6.21), (-5.9100000000000001, 2.3900000000000001, 4.5199999999999996), (1.04, -5.2300000000000004, 23.469999999999999), (-3.0800000000000001, 0.059999999999999998, 11.94), (0.75, -1.6200000000000001, 21.760000000000002), (0.89000000000000001, 0.089999999999999997, 19.02), (0.28999999999999998, -2.3399999999999999, 25.030000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP17)

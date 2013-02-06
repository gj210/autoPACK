#include as follow : execfile('pathto/POP83.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP83= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP83.sph',
radii = [[1.29, 2.3900000000000001, 3.4300000000000002, 4.21, 0.0, 2.3599999999999999, 2.8900000000000001, 3.1600000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP83.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP83',
positions = [[(3.8100000000000001, -4.5099999999999998, -4.2400000000000002), (1.8200000000000001, -3.23, -7.7599999999999998), (-2.2599999999999998, 4.2199999999999998, -16.91), (2.27, -0.25, -18.359999999999999), (5.9199999999999999, -5.0800000000000001, -2.8700000000000001), (-1.22, 2.8999999999999999, -10.6), (-2.9500000000000002, -1.1299999999999999, -8.0800000000000001), (-0.91000000000000003, -2.0299999999999998, -13.01)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP83)

#include as follow : execfile('pathto/POP63.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP63= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP63.sph',
radii = [[2.46, 3.4100000000000001, 1.3, 2.46, 3.9300000000000002, 3.1299999999999999, 2.75, 2.3500000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP63.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP63',
positions = [[(-4.29, -4.4900000000000002, 9.3200000000000003), (0.31, -5.75, 14.27), (-0.16, 8.7200000000000006, 6.6299999999999999), (-1.46, 6.6100000000000003, 10.210000000000001), (0.73999999999999999, -4.4100000000000001, 20.41), (-1.1899999999999999, 1.3500000000000001, 20.59), (2.73, 0.42999999999999999, 24.93), (-0.66000000000000003, 4.8499999999999996, 15.67)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP63)

#include as follow : execfile('pathto/POP72.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP72= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP72.sph',
radii = [[3.1299999999999999, 1.9299999999999999, 4.8700000000000001, 1.6200000000000001, 1.29, 2.4100000000000001, 4.71, 2.73]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP72.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP72',
positions = [[(0.63, -2.5299999999999998, -20.079999999999998), (2.0899999999999999, -3.1600000000000001, -9.7100000000000009), (-5.0999999999999996, 3.7799999999999998, -7.6200000000000001), (5.3099999999999996, -0.72999999999999998, -21.420000000000002), (2.3900000000000001, -3.3100000000000001, -5.3499999999999996), (2.1800000000000002, -3.3999999999999999, -14.779999999999999), (-2.7799999999999998, 1.52, -17.199999999999999), (2.1899999999999999, 1.1899999999999999, -22.510000000000002)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP72)

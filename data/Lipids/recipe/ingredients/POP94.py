#include as follow : execfile('pathto/POP94.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP94= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP94.sph',
radii = [[3.23, 3.73, 2.0800000000000001, 2.4700000000000002, 3.8999999999999999, 1.6100000000000001, 2.48, 3.2799999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP94.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP94',
positions = [[(3.2999999999999998, 2.2999999999999998, 10.800000000000001), (3.2200000000000002, 3.0, 2.8199999999999998), (4.3600000000000003, -2.7999999999999998, 21.84), (3.1400000000000001, 0.22, 17.27), (-10.48, 2.1499999999999999, 10.0), (6.6100000000000003, -4.0800000000000001, 19.73), (-0.66000000000000003, -1.6200000000000001, 19.329999999999998), (-5.5999999999999996, -0.33000000000000002, 15.710000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP94)

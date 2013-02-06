#include as follow : execfile('pathto/POP34.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP34= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP34.sph',
radii = [[1.3100000000000001, 4.9199999999999999, 3.6000000000000001, 3.6200000000000001, 1.96, 3.52, 1.02, 4.1600000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP34.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP34',
positions = [[(1.6000000000000001, -2.73, 22.52), (0.10000000000000001, 6.6799999999999997, 6.0999999999999996), (-1.3899999999999999, -0.95999999999999996, 16.699999999999999), (-1.49, -2.7400000000000002, 8.7899999999999991), (-0.80000000000000004, 0.48999999999999999, 21.699999999999999), (-1.52, -7.4699999999999998, 2.8100000000000001), (0.46999999999999997, -1.3400000000000001, 22.739999999999998), (3.0600000000000001, 3.5, 15.33)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP34)

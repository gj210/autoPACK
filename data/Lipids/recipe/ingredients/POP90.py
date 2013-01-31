#include as follow : execfile('pathto/POP90.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP90= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP90.sph',
radii = [[4.4500000000000002, 3.6499999999999999, 7.0700000000000003, 6.0700000000000003]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP90.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP90',
positions = [[(-0.42999999999999999, -3.7400000000000002, 12.69), (-4.8700000000000001, -3.8799999999999999, 6.0599999999999996), (1.1899999999999999, 2.5499999999999998, 19.809999999999999), (-0.88, 4.21, 7.2599999999999998)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP90)

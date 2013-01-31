#include as follow : execfile('pathto/POP71.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP71= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP71.sph',
radii = [[2.9900000000000002, 5.8600000000000003, 4.8899999999999997, 6.0999999999999996]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP71.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP71',
positions = [[(-0.34999999999999998, 2.79, -0.59999999999999998), (0.56999999999999995, 0.96999999999999997, -19.41), (5.8200000000000003, -1.98, -12.220000000000001), (0.72999999999999998, 0.27000000000000002, -7.5999999999999996)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP71)

#include as follow : execfile('pathto/POP8.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP8= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP8.sph',
radii = [[2.3599999999999999, 3.4900000000000002, 2.2200000000000002, 2.9199999999999999, 2.2000000000000002, 4.0099999999999998, 2.98, 3.0099999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP8.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP8',
positions = [[(-2.29, 1.3500000000000001, -15.16), (-0.27000000000000002, -2.48, -15.92), (-1.6499999999999999, 1.4399999999999999, -19.84), (5.5999999999999996, -0.26000000000000001, -0.90000000000000002), (-2.0899999999999999, -1.6299999999999999, -20.73), (3.48, -1.0900000000000001, -8.5999999999999996), (-2.5600000000000001, 1.74, -2.96), (-1.3300000000000001, 2.6200000000000001, -9.1500000000000004)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP8)

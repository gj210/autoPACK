#include as follow : execfile('pathto/POP21.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP21= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP21.sph',
radii = [[5.6200000000000001, 6.1200000000000001, 6.1399999999999997, 3.1699999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP21.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP21',
positions = [[(-0.080000000000000002, 0.96999999999999997, 19.010000000000002), (-0.56000000000000005, -5.9500000000000002, 9.0199999999999996), (-4.1699999999999999, 7.0800000000000001, 8.9800000000000004), (-1.47, -0.56999999999999995, 23.82)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP21)

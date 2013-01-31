#include as follow : execfile('pathto/POP74.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP74= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP74.sph',
radii = [[4.1699999999999999, 4.2699999999999996, 8.3000000000000007, 5.3399999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP74.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP74',
positions = [[(-0.28000000000000003, 2.4199999999999999, -21.510000000000002), (-3.1499999999999999, 1.49, -13.050000000000001), (-2.6800000000000002, -4.0800000000000001, -6.8099999999999996), (2.1800000000000002, -4.0, -16.66)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP74)

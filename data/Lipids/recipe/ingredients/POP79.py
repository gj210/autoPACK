#include as follow : execfile('pathto/POP79.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP79= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP79.sph',
radii = [[4.5099999999999998, 5.6600000000000001, 4.3899999999999997, 6.7000000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP79.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP79',
positions = [[(-0.55000000000000004, 4.5300000000000002, -11.859999999999999), (0.39000000000000001, -0.96999999999999997, -16.98), (0.17000000000000001, 5.8600000000000003, -2.1800000000000002), (5.9199999999999999, -4.2800000000000002, -5.04)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP79)

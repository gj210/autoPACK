#include as follow : execfile('pathto/LPO118.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO118= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO118.sph',
radii = [[7.4800000000000004, 5.4100000000000001, 5.6600000000000001, 3.3500000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO118.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO118',
positions = [[(-3.4500000000000002, 0.72999999999999998, -8.6600000000000001), (5.2999999999999998, 6.7000000000000002, -11.029999999999999), (1.05, 1.9199999999999999, -19.48), (2.3999999999999999, 1.6599999999999999, -25.109999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO118)

#include as follow : execfile('pathto/POP60.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP60= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP60.sph',
radii = [[5.9699999999999998, 5.0099999999999998, 4.9800000000000004, 5.25]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP60.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP60',
positions = [[(0.029999999999999999, 0.070000000000000007, 19.100000000000001), (-1.05, -2.1899999999999999, 9.4000000000000004), (3.1000000000000001, 4.3399999999999999, 8.0800000000000001), (0.14999999999999999, 2.4300000000000002, 0.38)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP60)

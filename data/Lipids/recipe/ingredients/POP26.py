#include as follow : execfile('pathto/POP26.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP26= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP26.sph',
radii = [[3.8999999999999999, 6.2199999999999998, 6.8300000000000001, 4.8600000000000003]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP26.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP26',
positions = [[(-0.98999999999999999, 3.7400000000000002, -12.199999999999999), (2.0, 0.23000000000000001, -18.809999999999999), (2.0800000000000001, -4.9699999999999998, -5.9100000000000001), (-1.0800000000000001, 4.9500000000000002, -2.6899999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP26)

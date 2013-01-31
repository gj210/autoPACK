#include as follow : execfile('pathto/POP87.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP87= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP87.sph',
radii = [[4.4400000000000004, 4.9800000000000004, 5.5300000000000002, 6.8600000000000003]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP87.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP87',
positions = [[(-2.6800000000000002, -3.4900000000000002, 21.350000000000001), (-5.75, -6.1699999999999999, 11.859999999999999), (-0.76000000000000001, 5.0599999999999996, 8.9800000000000004), (-0.40999999999999998, 2.98, 22.329999999999998)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP87)

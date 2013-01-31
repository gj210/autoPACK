#include as follow : execfile('pathto/POP41.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP41= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP41.sph',
radii = [[4.0499999999999998, 3.0800000000000001, 6.4199999999999999, 3.9700000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP41.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP41',
positions = [[(-2.04, -2.2999999999999998, 13.369999999999999), (-0.57999999999999996, 2.8999999999999999, 24.59), (-2.1000000000000001, -4.4400000000000004, 6.2300000000000004), (-0.44, -1.0600000000000001, 19.75)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP41)

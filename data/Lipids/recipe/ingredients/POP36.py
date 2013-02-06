#include as follow : execfile('pathto/POP36.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP36= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP36.sph',
radii = [[5.0599999999999996, 1.9099999999999999, 4.3499999999999996, 2.3900000000000001, 2.3199999999999998, 2.5099999999999998, 3.6299999999999999, 1.3100000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP36.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP36',
positions = [[(-3.7400000000000002, 8.4800000000000004, 9.7400000000000002), (-3.8500000000000001, -3.7000000000000002, 9.3699999999999992), (0.10000000000000001, 4.5999999999999996, 17.780000000000001), (3.4300000000000002, -3.9700000000000002, 22.390000000000001), (0.23000000000000001, -3.8100000000000001, 12.18), (2.2200000000000002, -4.75, 16.109999999999999), (2.2000000000000002, -0.75, 21.780000000000001), (-4.5, -5.4000000000000004, 6.2199999999999998)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP36)

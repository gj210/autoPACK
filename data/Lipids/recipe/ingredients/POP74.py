#include as follow : execfile('pathto/POP74.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP74= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP74.sph',
radii = [[1.73, 1.98, 2.9500000000000002, 4.2199999999999998, 1.5700000000000001, 2.7599999999999998, 4.5899999999999999, 3.4199999999999999]],
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
positions = [[(0.050000000000000003, 4.8600000000000003, -22.550000000000001), (-2.3199999999999998, 2.21, -17.91), (-2.1200000000000001, 2.3199999999999998, -11.85), (3.1899999999999999, -3.29, -16.899999999999999), (3.02, 4.0599999999999996, -21.780000000000001), (0.55000000000000004, 1.0600000000000001, -20.850000000000001), (2.0800000000000001, -6.2800000000000002, -7.4500000000000002), (-5.9199999999999999, 0.46000000000000002, -6.6500000000000004)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP74)

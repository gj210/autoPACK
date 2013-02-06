#include as follow : execfile('pathto/POP80.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP80= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP80.sph',
radii = [[2.96, 1.9099999999999999, 2.5699999999999998, 3.2000000000000002, 2.52, 4.04, 2.2799999999999998, 1.9399999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP80.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP80',
positions = [[(-4.9900000000000002, -0.40999999999999998, -8.5800000000000001), (0.60999999999999999, -2.2799999999999998, -6.1900000000000004), (-2.0299999999999998, -0.23000000000000001, -14.210000000000001), (3.1699999999999999, -2.2999999999999998, -17.129999999999999), (-5.3600000000000003, 2.9700000000000002, -2.7799999999999998), (2.48, 0.93999999999999995, -23.809999999999999), (0.01, 0.84999999999999998, -19.199999999999999), (1.4399999999999999, -2.1200000000000001, -11.109999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP80)

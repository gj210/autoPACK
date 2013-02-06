#include as follow : execfile('pathto/POP32.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP32= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP32.sph',
radii = [[2.1499999999999999, 3.71, 2.4199999999999999, 1.3500000000000001, 2.3300000000000001, 3.4500000000000002, 4.0899999999999999, 3.2400000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP32.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP32',
positions = [[(0.41999999999999998, 2.1899999999999999, -20.920000000000002), (-1.8500000000000001, 5.7000000000000002, -4.1900000000000004), (-3.9399999999999999, 0.59999999999999998, -17.300000000000001), (3.1099999999999999, 0.58999999999999997, -22.010000000000002), (0.01, -2.3900000000000001, -19.109999999999999), (1.48, -3.8399999999999999, -13.01), (3.1800000000000002, -4.6699999999999999, -5.0499999999999998), (-2.9500000000000002, 2.79, -11.01)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP32)

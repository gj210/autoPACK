#include as follow : execfile('pathto/POP44.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP44= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP44.sph',
radii = [[4.96, 3.6600000000000001, 4.5099999999999998, 1.95, 0.77000000000000002, 2.9399999999999999, 2.27, 2.9199999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP44.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP44',
positions = [[(2.71, -8.1099999999999994, -10.33), (-0.47999999999999998, 2.6099999999999999, -17.420000000000002), (-2.4199999999999999, -3.9500000000000002, -18.129999999999999), (1.5600000000000001, 7.8899999999999997, -6.21), (3.23, 10.550000000000001, -3.9100000000000001), (0.23000000000000001, 3.77, -10.050000000000001), (-1.4399999999999999, 3.3700000000000001, -23.100000000000001), (-0.93000000000000005, 0.34999999999999998, -23.440000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP44)

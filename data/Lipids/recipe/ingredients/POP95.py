#include as follow : execfile('pathto/POP95.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP95= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP95.sph',
radii = [[0.76000000000000001, 2.8599999999999999, 3.8799999999999999, 4.9000000000000004, 1.9199999999999999, 2.7400000000000002, 1.8999999999999999, 2.9500000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP95.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP95',
positions = [[(-4.5700000000000003, 3.4399999999999999, 8.4800000000000004), (-3.9399999999999999, 2.8700000000000001, 17.940000000000001), (2.6800000000000002, -5.5099999999999998, 13.91), (3.6499999999999999, -7.0, 5.1799999999999997), (-3.4900000000000002, 4.3200000000000003, 5.4199999999999999), (0.95999999999999996, 4.2400000000000002, 21.350000000000001), (-4.9100000000000001, 3.0899999999999999, 12.18), (1.02, -0.20999999999999999, 19.100000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP95)

#include as follow : execfile('pathto/POP70.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP70= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP70.sph',
radii = [[7.1900000000000004, 4.0700000000000003, 6.6900000000000004, 4.8899999999999997]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP70.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP70',
positions = [[(0.93000000000000005, -0.40999999999999998, -22.870000000000001), (-2.8199999999999998, 2.1600000000000001, -17.940000000000001), (4.9900000000000002, -3.0, -9.9399999999999995), (-3.54, 3.5699999999999998, -7.8499999999999996)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP70)

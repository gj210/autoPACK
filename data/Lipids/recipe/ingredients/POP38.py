#include as follow : execfile('pathto/POP38.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP38= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP38.sph',
radii = [[1.9399999999999999, 1.3600000000000001, 5.54, 3.8300000000000001, 3.4500000000000002, 3.4300000000000002, 1.9399999999999999, 2.79]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP38.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP38',
positions = [[(-0.029999999999999999, -1.9399999999999999, 23.030000000000001), (-1.71, -4.3600000000000003, 22.440000000000001), (5.6100000000000003, -3.7400000000000002, 6.6900000000000004), (-7.25, 5.3799999999999999, 5.2800000000000002), (4.9199999999999999, -1.5700000000000001, 15.869999999999999), (-5.1500000000000004, 5.0599999999999996, 13.16), (1.6499999999999999, -1.01, 19.800000000000001), (-0.55000000000000004, 1.99, 17.489999999999998)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP38)

#include as follow : execfile('pathto/POP3.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP3= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP3.sph',
radii = [[4.8899999999999997, 5.4900000000000002, 1.7, 2.1200000000000001, 3.27, 4.6799999999999997, 2.3900000000000001, 0.0]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP3.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP3',
positions = [[(-1.3400000000000001, -0.31, -13.359999999999999), (3.3900000000000001, 3.02, -5.0300000000000002), (0.48999999999999999, 0.81999999999999995, -26.140000000000001), (2.1800000000000002, -0.48999999999999999, -23.760000000000002), (1.03, 2.0499999999999998, -18.219999999999999), (-4.7400000000000002, -4.1399999999999997, -7.2699999999999996), (-0.92000000000000004, -1.21, -19.140000000000001), (3.96, 0.01, -24.25)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP3)

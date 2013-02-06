#include as follow : execfile('pathto/POP70.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP70= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP70.sph',
radii = [[4.1100000000000003, 2.3700000000000001, 4.0099999999999998, 2.9100000000000001, 1.3, 1.6200000000000001, 1.98, 3.9100000000000001]],
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
positions = [[(-0.38, -0.94999999999999996, -19.949999999999999), (4.7300000000000004, -3.0499999999999998, -8.9600000000000009), (-4.3300000000000001, 3.4100000000000001, -6.7699999999999996), (2.79, -4.4000000000000004, -14.220000000000001), (6.8399999999999999, -1.25, -5.3499999999999996), (1.26, 0.13, -25.969999999999999), (-1.21, 1.03, -24.010000000000002), (-3.21, 3.1000000000000001, -15.25)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP70)

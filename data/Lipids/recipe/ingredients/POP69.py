#include as follow : execfile('pathto/POP69.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP69= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP69.sph',
radii = [[3.46, 2.3500000000000001, 2.6000000000000001, 1.3, 2.7400000000000002, 3.04, 1.9199999999999999, 2.3399999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP69.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP69',
positions = [[(2.02, -0.66000000000000003, -19.5), (-4.3300000000000001, 4.3399999999999999, -12.539999999999999), (-3.4100000000000001, 1.5600000000000001, -17.140000000000001), (-3.1000000000000001, 5.9000000000000004, -8.4700000000000006), (2.6699999999999999, -3.0499999999999998, -13.02), (2.0899999999999999, -2.9700000000000002, -6.6200000000000001), (-1.0600000000000001, 7.04, -5.4000000000000004), (-0.35999999999999999, -6.3499999999999996, -3.21)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP69)

#include as follow : execfile('pathto/POP46.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP46= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP46.sph',
radii = [[0.73999999999999999, 2.71, 4.3099999999999996, 5.6600000000000001, 1.6899999999999999, 1.6699999999999999, 3.71, 3.5499999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP46.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP46',
positions = [[(1.1200000000000001, 0.28000000000000003, -21.309999999999999), (-0.88, -1.6000000000000001, -17.030000000000001), (-3.1800000000000002, 1.45, -14.380000000000001), (-4.4000000000000004, 3.4300000000000002, -3.8399999999999999), (0.01, 1.25, -20.370000000000001), (-0.93000000000000005, 2.3900000000000001, -23.25), (8.2899999999999991, -4.8200000000000003, -5.7000000000000002), (3.9500000000000002, -3.9300000000000002, -12.99)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP46)

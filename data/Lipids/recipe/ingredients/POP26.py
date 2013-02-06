#include as follow : execfile('pathto/POP26.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP26= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP26.sph',
radii = [[3.7999999999999998, 3.1699999999999999, 2.75, 2.9900000000000002, 2.52, 2.7000000000000002, 2.4500000000000002, 2.54]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP26.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP26',
positions = [[(1.3600000000000001, -4.3499999999999996, -8.8000000000000007), (2.4399999999999999, 1.55, -20.600000000000001), (0.59999999999999998, -3.54, -16.41), (-2.0499999999999998, 3.2200000000000002, -13.35), (0.62, -6.7699999999999996, -1.8600000000000001), (-2.0600000000000001, 4.0099999999999998, -6.6399999999999997), (-0.80000000000000004, 0.29999999999999999, -18.329999999999998), (-2.0899999999999999, 4.9000000000000004, -0.35999999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP26)

#include as follow : execfile('pathto/POP96.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP96= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP96.sph',
radii = [[2.5800000000000001, 1.27, 3.2999999999999998, 3.6800000000000002, 3.25, 2.5499999999999998, 2.54, 2.0]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP96.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP96',
positions = [[(3.27, 1.46, 16.640000000000001), (3.1800000000000002, 1.72, 11.66), (0.65000000000000002, 1.95, 22.789999999999999), (-3.2200000000000002, -1.9299999999999999, 18.75), (-4.4400000000000004, -3.8500000000000001, 11.02), (6.2599999999999998, 2.5800000000000001, 3.3199999999999998), (-6.7999999999999998, -6.0, 5.1799999999999997), (4.7400000000000002, 2.9900000000000002, 8.1199999999999992)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP96)

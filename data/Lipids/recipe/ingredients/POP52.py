#include as follow : execfile('pathto/POP52.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP52= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP52.sph',
radii = [[2.8700000000000001, 3.3399999999999999, 3.04, 2.2400000000000002, 1.8999999999999999, 3.8199999999999998, 1.5900000000000001, 3.6899999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP52.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP52',
positions = [[(-6.6799999999999997, -1.3999999999999999, -12.75), (6.8899999999999997, -3.8199999999999998, -4.6799999999999997), (1.03, 2.1600000000000001, -15.859999999999999), (-4.2599999999999998, 2.8199999999999998, -16.84), (-0.93999999999999995, 0.69999999999999996, -20.329999999999998), (-5.21, -3.0600000000000001, -5.4500000000000002), (1.8400000000000001, 1.8999999999999999, -20.149999999999999), (5.4900000000000002, 1.1100000000000001, -10.0)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP52)

#include as follow : execfile('pathto/POP18.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP18= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP18.sph',
radii = [[2.7400000000000002, 1.7, 1.97, 2.9900000000000002, 3.96, 2.8100000000000001, 3.0099999999999998, 3.6499999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP18.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP18',
positions = [[(0.87, -1.1399999999999999, 15.52), (0.33000000000000002, -2.1800000000000002, 19.239999999999998), (-1.05, 0.44, 19.420000000000002), (1.71, -2.5800000000000001, 9.8599999999999994), (-0.48999999999999999, -3.4300000000000002, 2.3700000000000001), (-0.65000000000000002, 2.54, 13.390000000000001), (-0.71999999999999997, 1.9299999999999999, 6.8799999999999999), (0.25, 4.4900000000000002, -0.46000000000000002)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP18)

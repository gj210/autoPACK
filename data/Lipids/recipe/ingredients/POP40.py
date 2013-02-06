#include as follow : execfile('pathto/POP40.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP40= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP40.sph',
radii = [[3.75, 1.9299999999999999, 3.3399999999999999, 4.04, 3.96, 3.0, 2.4199999999999999, 1.6799999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP40.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP40',
positions = [[(1.72, 0.33000000000000002, 12.25), (2.8799999999999999, -0.79000000000000004, 23.329999999999998), (-1.72, 0.75, 18.120000000000001), (-2.9399999999999999, -0.46999999999999997, 11.34), (-3.3500000000000001, -0.62, 4.1500000000000004), (2.8700000000000001, 1.47, 19.390000000000001), (-5.5300000000000002, -0.52000000000000002, 7.6600000000000001), (5.8399999999999999, -0.42999999999999999, 22.030000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP40)

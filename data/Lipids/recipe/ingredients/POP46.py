#include as follow : execfile('pathto/POP46.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP46= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP46.sph',
radii = [[6.3300000000000001, 3.1000000000000001, 5.6600000000000001, 6.0499999999999998]],
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
positions = [[(9.3000000000000007, -4.8899999999999997, -7.9500000000000002), (1.8400000000000001, 1.3100000000000001, -21.969999999999999), (-2.29, 3.0800000000000001, -3.8399999999999999), (0.58999999999999997, -0.78000000000000003, -15.27)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP46)

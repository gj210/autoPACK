#include as follow : execfile('pathto/POP1.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP1= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP1.sph',
radii = [[3.8100000000000001, 6.2199999999999998, 5.0700000000000003, 3.0299999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP1.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP1',
positions = [[(-2.1000000000000001, -1.51, -22.41), (-0.48999999999999999, -4.6799999999999997, -8.7200000000000006), (-3.8199999999999998, -4.4100000000000001, -16.579999999999998), (3.6000000000000001, -0.97999999999999998, -24.75)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP1)

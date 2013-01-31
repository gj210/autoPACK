#include as follow : execfile('pathto/POP5.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP5= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP5.sph',
radii = [[6.5499999999999998, 4.2999999999999998, 6.79, 4.7599999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP5.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP5',
positions = [[(-1.3100000000000001, -0.63, -18.559999999999999), (2.3100000000000001, -6.7599999999999998, -2.8300000000000001), (-6.5499999999999998, 2.3500000000000001, -6.1600000000000001), (3.5099999999999998, -2.1699999999999999, -11.0)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP5)

#include as follow : execfile('pathto/POP51.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP51= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP51.sph',
radii = [[1.3, 1.3500000000000001, 2.46, 2.6699999999999999, 3.1800000000000002, 2.9300000000000002, 4.6500000000000004, 3.5699999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP51.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP51',
positions = [[(-4.7300000000000004, -2.7400000000000002, -15.48), (-3.0899999999999999, -0.60999999999999999, -17.359999999999999), (-5.2199999999999998, -4.1600000000000001, -10.84), (1.28, 0.48999999999999999, -16.969999999999999), (2.5699999999999998, 2.8300000000000001, -10.539999999999999), (1.3100000000000001, -0.71999999999999997, -20.84), (4.1500000000000004, 4.7599999999999998, -2.5099999999999998), (-3.3500000000000001, -3.4199999999999999, -4.1600000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP51)

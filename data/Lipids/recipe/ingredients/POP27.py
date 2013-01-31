#include as follow : execfile('pathto/POP27.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP27= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP27.sph',
radii = [[5.9900000000000002, 4.7300000000000004, 3.6800000000000002, 5.9000000000000004]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP27.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP27',
positions = [[(2.2000000000000002, 0.83999999999999997, -19.309999999999999), (-4.9299999999999997, -2.1299999999999999, -11.01), (-5.8300000000000001, -5.1200000000000001, -2.6899999999999999), (7.0999999999999996, -0.75, -8.9900000000000002)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP27)

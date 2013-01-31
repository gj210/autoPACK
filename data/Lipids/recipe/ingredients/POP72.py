#include as follow : execfile('pathto/POP72.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP72= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP72.sph',
radii = [[5.54, 4.8700000000000001, 6.0700000000000003, 4.71]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP72.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP72',
positions = [[(2.3799999999999999, -4.0800000000000001, -9.6099999999999994), (-4.8899999999999997, 2.8500000000000001, -7.6200000000000001), (2.9199999999999999, -1.99, -20.829999999999998), (-2.5699999999999998, 0.59999999999999998, -17.199999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP72)

#include as follow : execfile('pathto/LPO114.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO114= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO114.sph',
radii = [[1.5700000000000001, 3.6499999999999999, 3.54, 1.3, 3.0699999999999998, 3.6099999999999999, 1.9099999999999999, 2.5800000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO114.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO114',
positions = [[(-1.1100000000000001, 4.0899999999999999, -18.260000000000002), (-2.9199999999999999, -2.8999999999999999, -7.71), (1.6599999999999999, 0.14000000000000001, -11.84), (-1.6399999999999999, 5.2400000000000002, -14.65), (-0.46999999999999997, 3.0699999999999998, -5.1200000000000001), (1.95, -3.8199999999999998, -18.359999999999999), (0.14000000000000001, 4.75, -11.1), (-0.82999999999999996, -0.080000000000000002, -18.859999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO114)

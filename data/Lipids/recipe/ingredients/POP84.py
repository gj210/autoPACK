#include as follow : execfile('pathto/POP84.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP84= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP84.sph',
radii = [[6.46, 2.6600000000000001, 5.9400000000000004, 5.9400000000000004]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP84.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP84',
positions = [[(-0.070000000000000007, -0.17999999999999999, 13.949999999999999), (-2.7400000000000002, -0.34000000000000002, 19.460000000000001), (1.55, -8.2100000000000009, 5.9199999999999999), (-4.9900000000000002, 4.8700000000000001, 3.3599999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP84)

#include as follow : execfile('pathto/POP21.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP21= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP21.sph',
radii = [[3.1200000000000001, 1.3200000000000001, 3.4500000000000002, 2.5, 2.6899999999999999, 1.9299999999999999, 3.1400000000000001, 3.4199999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP21.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP21',
positions = [[(-1.0600000000000001, 5.5800000000000001, 13.970000000000001), (-1.5600000000000001, -0.54000000000000004, 23.109999999999999), (-3.4500000000000002, 7.4699999999999998, 6.7400000000000002), (1.4299999999999999, 3.23, 19.829999999999998), (1.5, -1.75, 20.710000000000001), (0.78000000000000003, -0.65000000000000002, 25.010000000000002), (1.53, -3.98, 14.109999999999999), (0.59999999999999998, -6.9400000000000004, 7.29)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP21)

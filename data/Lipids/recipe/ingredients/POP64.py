#include as follow : execfile('pathto/POP64.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP64= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP64.sph',
radii = [[0.0, 4.0700000000000003, 0.5, 4.3399999999999999, 2.04, 3.8900000000000001, 1.3300000000000001, 2.7400000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP64.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP64',
positions = [[(5.4400000000000004, -2.3799999999999999, 21.0), (-1.6599999999999999, -0.48999999999999999, 0.68000000000000005), (4.7199999999999998, -3.3999999999999999, 20.600000000000001), (0.050000000000000003, 2.2799999999999998, 13.34), (1.5600000000000001, -0.95999999999999996, 20.640000000000001), (-1.6200000000000001, 1.8500000000000001, 7.0), (3.54, -2.7599999999999998, 21.600000000000001), (-0.12, -1.1599999999999999, 15.56)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP64)

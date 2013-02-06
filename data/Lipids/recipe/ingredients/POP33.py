#include as follow : execfile('pathto/POP33.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP33= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP33.sph',
radii = [[4.0999999999999996, 2.7999999999999998, 3.6099999999999999, 1.29, 2.4199999999999999, 1.3100000000000001, 2.71, 2.3500000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP33.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP33',
positions = [[(3.23, 4.5, 7.9000000000000004), (0.88, -1.28, 22.870000000000001), (-0.37, 4.21, 15.27), (-3.54, -3.1600000000000001, 14.92), (-1.5600000000000001, -0.68999999999999995, 17.940000000000001), (-2.2400000000000002, -4.0199999999999996, 11.93), (-0.28000000000000003, -3.2200000000000002, 1.6000000000000001), (-1.1799999999999999, -2.8900000000000001, 7.79)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP33)

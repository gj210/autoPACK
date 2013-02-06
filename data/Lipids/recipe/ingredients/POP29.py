#include as follow : execfile('pathto/POP29.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP29= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP29.sph',
radii = [[2.9300000000000002, 4.0499999999999998, 1.6000000000000001, 1.9099999999999999, 2.3500000000000001, 3.0699999999999998, 1.9099999999999999, 4.7599999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP29.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP29',
positions = [[(-2.8999999999999999, 1.1899999999999999, -20.710000000000001), (1.21, 4.96, -16.300000000000001), (-3.04, -2.3999999999999999, -22.0), (0.81999999999999995, -9.6500000000000004, -3.6499999999999999), (0.73999999999999999, -5.4299999999999997, -5.9199999999999999), (-0.31, -0.84999999999999998, -16.109999999999999), (0.78000000000000003, -3.5, -10.630000000000001), (2.3900000000000001, 5.4500000000000002, -6.7199999999999998)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP29)

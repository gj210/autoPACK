#include as follow : execfile('pathto/POP30.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP30= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP30.sph',
radii = [[0.76000000000000001, 3.25, 4.3600000000000003, 1.26, 2.96, 1.9299999999999999, 1.9399999999999999, 5.5700000000000003]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP30.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP30',
positions = [[(-1.97, 5.1500000000000004, -2.1400000000000001), (-1.03, 3.2999999999999998, -18.059999999999999), (1.71, -1.74, -16.010000000000002), (-1.4199999999999999, 2.8199999999999998, -4.1399999999999997), (-2.5699999999999998, 0.93999999999999995, -22.949999999999999), (-1.5, 2.0800000000000001, -7.79), (-1.78, 3.7400000000000002, -12.57), (4.2400000000000002, -5.9800000000000004, -6.9299999999999997)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP30)

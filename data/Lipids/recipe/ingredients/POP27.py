#include as follow : execfile('pathto/POP27.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP27= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP27.sph',
radii = [[1.95, 2.5600000000000001, 3.1400000000000001, 2.5699999999999998, 1.29, 3.8999999999999999, 2.9199999999999999, 3.5299999999999998]],
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
positions = [[(-6.5199999999999996, -1.1599999999999999, -11.58), (-6.6600000000000001, -4.9299999999999997, -1.76), (-0.28999999999999998, 0.25, -18.23), (-6.7400000000000002, -2.3700000000000001, -6.79), (-3.6600000000000001, -0.080000000000000002, -14.75), (6.2000000000000002, 0.96999999999999997, -6.9000000000000004), (1.9099999999999999, 3.6899999999999999, -21.260000000000002), (5.1399999999999997, -1.3200000000000001, -14.68)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP27)

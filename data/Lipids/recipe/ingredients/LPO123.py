#include as follow : execfile('pathto/LPO123.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO123= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO123.sph',
radii = [[3.0899999999999999, 2.7999999999999998, 4.3099999999999996, 3.29, 3.9199999999999999, 1.48, 2.0099999999999998, 2.29]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO123.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO123',
positions = [[(0.20000000000000001, -0.83999999999999997, -4.1399999999999997), (1.6599999999999999, -0.070000000000000007, -17.039999999999999), (-1.3799999999999999, -5.04, -15.34), (0.26000000000000001, -3.0600000000000001, -22.780000000000001), (0.20000000000000001, 0.02, -10.289999999999999), (-2.3399999999999999, 7.3200000000000003, -21.030000000000001), (0.90000000000000002, 1.1100000000000001, -21.199999999999999), (0.14999999999999999, 4.5999999999999996, -21.82)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO123)

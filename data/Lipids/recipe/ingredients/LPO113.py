#include as follow : execfile('pathto/LPO113.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO113= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO113.sph',
radii = [[2.7599999999999998, 3.7999999999999998, 2.2999999999999998, 1.95, 3.7400000000000002, 2.0899999999999999, 2.6499999999999999, 3.2599999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO113.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO113',
positions = [[(2.0600000000000001, -0.90000000000000002, -23.449999999999999), (-0.60999999999999999, 6.1399999999999997, -3.6400000000000001), (-0.56000000000000005, -2.8500000000000001, -18.690000000000001), (-0.16, 3.8100000000000001, -16.739999999999998), (-1.3799999999999999, -5.6500000000000004, -4.4699999999999998), (0.81000000000000005, 1.3100000000000001, -20.0), (-0.41999999999999998, 5.6299999999999999, -11.19), (-1.3799999999999999, -5.1900000000000004, -12.49)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO113)

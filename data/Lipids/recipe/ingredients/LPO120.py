#include as follow : execfile('pathto/LPO120.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO120= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO120.sph',
radii = [[1.95, 2.7000000000000002, 2.6800000000000002, 3.46, 2.5800000000000001, 3.3500000000000001, 4.0199999999999996, 1.3100000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO120.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO120',
positions = [[(1.76, 1.8100000000000001, 2.6600000000000001), (-1.1100000000000001, 3.3900000000000001, 7.2400000000000002), (0.68999999999999995, -0.71999999999999997, 17.460000000000001), (2.2599999999999998, -6.3200000000000003, 5.5599999999999996), (2.2799999999999998, -4.0099999999999998, 12.369999999999999), (-2.54, 3.6499999999999999, 14.279999999999999), (-2.0299999999999998, 0.98999999999999999, 21.379999999999999), (4.7300000000000004, -0.16, 0.089999999999999997)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO120)

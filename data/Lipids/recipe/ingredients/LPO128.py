#include as follow : execfile('pathto/LPO128.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO128= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO128.sph',
radii = [[3.4700000000000002, 2.75, 2.73, 2.7400000000000002, 3.1400000000000001, 3.3599999999999999, 1.9299999999999999, 1.9199999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO128.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'LPO128',
positions = [[(-4.5599999999999996, -4.8300000000000001, 22.870000000000001), (-1.3600000000000001, 0.34999999999999998, 19.949999999999999), (4.2400000000000002, -0.31, 10.609999999999999), (0.76000000000000001, -3.0699999999999998, 20.23), (1.1100000000000001, 3.1400000000000001, 14.02), (2.3399999999999999, 4.6399999999999997, 7.7000000000000002), (-0.25, 8.6699999999999999, 4.9199999999999999), (3.2799999999999998, -2.3999999999999999, 15.75)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO128)

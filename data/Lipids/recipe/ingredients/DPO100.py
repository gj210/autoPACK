#include as follow : execfile('pathto/DPO100.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO100= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO100.sph',
radii = [[2.9700000000000002, 3.6299999999999999, 2.98, 2.4199999999999999, 3.7999999999999998, 2.52, 2.9100000000000001, 2.0]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/DPO100.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'DPO100',
positions = [[(2.4900000000000002, -0.70999999999999996, 15.58), (0.71999999999999997, 1.3899999999999999, 22.530000000000001), (-2.02, 1.8300000000000001, 2.8300000000000001), (-4.1500000000000004, -0.42999999999999999, 8.3200000000000003), (2.4100000000000001, 0.13, 2.6800000000000002), (-0.93999999999999995, -1.51, 17.75), (1.1200000000000001, -0.80000000000000004, 9.25), (-2.7799999999999998, -1.8799999999999999, 13.449999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO100)

#include as follow : execfile('pathto/DPO110.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO110= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO110.sph',
radii = [[3.46, 0.76000000000000001, 1.98, 1.9299999999999999, 2.9199999999999999, 2.98, 1.96, 4.1100000000000003]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/DPO110.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'DPO110',
positions = [[(-0.28999999999999998, -0.23000000000000001, 13.880000000000001), (3.8199999999999998, -2.0699999999999998, 14.76), (2.8300000000000001, -3.4100000000000001, 11.890000000000001), (0.92000000000000004, -4.4699999999999998, 8.0700000000000003), (-1.3700000000000001, -3.1499999999999999, 2.9300000000000002), (0.28000000000000003, 6.2300000000000004, 9.6099999999999994), (3.0299999999999998, -0.26000000000000001, 17.559999999999999), (-1.6100000000000001, 1.1299999999999999, 19.489999999999998)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO110)

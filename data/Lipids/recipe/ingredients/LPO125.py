#include as follow : execfile('pathto/LPO125.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO125= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO125.sph',
radii = [[5.8099999999999996, 6.0999999999999996, 3.2200000000000002, 5.7300000000000004]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO125.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO125',
positions = [[(-0.93999999999999995, 0.14999999999999999, 17.239999999999998), (-5.0800000000000001, -6.9000000000000004, 7.2400000000000002), (-0.059999999999999998, -0.22, 23.440000000000001), (-1.5, 10.449999999999999, 13.039999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO125)

#include as follow : execfile('pathto/POP61.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP61= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP61.sph',
radii = [[5.4800000000000004, 4.1500000000000004, 5.9299999999999997, 3.96]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP61.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP61',
positions = [[(1.54, -6.9000000000000004, 7.3799999999999999), (-0.79000000000000004, 2.4100000000000001, 6.9199999999999999), (-0.40000000000000002, -0.31, 16.059999999999999), (1.5900000000000001, 2.3199999999999998, 22.699999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP61)

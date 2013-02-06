#include as follow : execfile('pathto/POP4.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP4= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP4.sph',
radii = [[3.9900000000000002, 0.76000000000000001, 0.76000000000000001, 4.9299999999999997, 3.6299999999999999, 2.3799999999999999, 1.3, 4.4100000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP4.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP4',
positions = [[(1.8500000000000001, -0.93000000000000005, -24.190000000000001), (6.8799999999999999, -3.8199999999999998, -3.3999999999999999), (6.9500000000000002, -1.74, -4.8099999999999996), (-8.6400000000000006, 1.49, -10.51), (2.4100000000000001, 1.75, -17.41), (4.4900000000000002, 1.8899999999999999, -10.470000000000001), (6.3099999999999996, 0.17999999999999999, -6.7699999999999996), (-3.98, -1.1799999999999999, -17.93)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP4)

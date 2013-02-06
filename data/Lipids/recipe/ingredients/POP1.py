#include as follow : execfile('pathto/POP1.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP1= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP1.sph',
radii = [[2.75, 3.1600000000000001, 5.5099999999999998, 0.70999999999999996, 2.5, 1.3400000000000001, 3.3399999999999999, 4.0099999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP1.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP1',
positions = [[(2.1899999999999999, 2.54, -24.52), (-4.8600000000000003, -2.5299999999999998, -12.0), (1.52, -1.5700000000000001, -8.3100000000000005), (4.6500000000000004, 1.78, -25.350000000000001), (-5.1100000000000003, -1.71, -18.07), (6.1799999999999997, 1.6299999999999999, -24.280000000000001), (-2.02, 0.68000000000000005, -23.620000000000001), (0.080000000000000002, 0.96999999999999997, -18.239999999999998)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP1)

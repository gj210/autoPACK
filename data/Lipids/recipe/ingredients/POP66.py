#include as follow : execfile('pathto/POP66.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP66= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP66.sph',
radii = [[3.23, 1.9299999999999999, 3.3900000000000001, 3.1400000000000001, 2.8100000000000001, 2.6200000000000001, 3.1299999999999999, 3.0499999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP66.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP66',
positions = [[(-3.48, 1.02, 16.82), (-6.6500000000000004, -2.5499999999999998, 6.0199999999999996), (4.3499999999999996, 0.38, 11.119999999999999), (1.97, -0.17000000000000001, 22.120000000000001), (5.21, -0.070000000000000007, 17.73), (-6.5999999999999996, -0.68000000000000005, 10.81), (1.0600000000000001, 2.0899999999999999, 5.29), (-0.93999999999999995, -0.66000000000000003, 22.98)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP66)

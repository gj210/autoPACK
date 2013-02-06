#include as follow : execfile('pathto/POP79.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP79= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP79.sph',
radii = [[1.52, 3.2999999999999998, 2.9399999999999999, 2.5899999999999999, 3.1400000000000001, 3.8100000000000001, 1.8700000000000001, 2.48]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP79.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP79',
positions = [[(-0.14000000000000001, -1.45, -19.329999999999998), (4.8700000000000001, -4.9800000000000004, -0.93999999999999995), (-3.0, 4.9100000000000001, -7.9000000000000004), (-1.9299999999999999, 3.6400000000000001, -13.93), (-1.1200000000000001, 5.9000000000000004, -1.04), (3.52, -4.0199999999999996, -8.6300000000000008), (-3.2599999999999998, -0.75, -18.100000000000001), (-1.1100000000000001, -1.5700000000000001, -13.59)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP79)

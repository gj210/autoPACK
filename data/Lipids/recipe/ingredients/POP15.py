#include as follow : execfile('pathto/POP15.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP15= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP15.sph',
radii = [[4.1399999999999997, 1.97, 2.5099999999999998, 2.73, 1.9099999999999999, 2.3799999999999999, 3.29, 3.5800000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP15.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP15',
positions = [[(-1.0, 3.5, 7.8200000000000003), (0.92000000000000004, -6.0300000000000002, 7.46), (-1.51, -7.8099999999999996, 2.7799999999999998), (-0.56000000000000005, 2.3100000000000001, 20.120000000000001), (0.23000000000000001, -3.1600000000000001, 16.359999999999999), (0.94999999999999996, -5.4100000000000001, 12.35), (0.44, 6.1900000000000004, 0.97999999999999998), (1.1799999999999999, 1.9399999999999999, 16.059999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP15)

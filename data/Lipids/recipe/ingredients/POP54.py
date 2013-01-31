#include as follow : execfile('pathto/POP54.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP54= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP54.sph',
radii = [[6.2400000000000002, 6.0599999999999996, 4.3600000000000003, 4.8499999999999996]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP54.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP54',
positions = [[(1.28, -5.9900000000000002, -7.7199999999999998), (1.1399999999999999, -2.1000000000000001, -19.719999999999999), (6.3200000000000003, 7.4800000000000004, -4.8300000000000001), (3.6800000000000002, 4.5599999999999996, -13.789999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP54)

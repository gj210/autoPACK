#include as follow : execfile('pathto/POP53.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP53= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP53.sph',
radii = [[3.0299999999999998, 2.8399999999999999, 3.4399999999999999, 2.5499999999999998, 3.21, 1.6699999999999999, 2.0099999999999998, 3.71]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP53.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP53',
positions = [[(-1.25, -2.98, -20.710000000000001), (-1.8200000000000001, 6.8899999999999997, -12.49), (1.2, -5.2000000000000002, -7.3600000000000003), (0.34999999999999998, 7.2300000000000004, -6.3700000000000001), (1.3600000000000001, -4.9100000000000001, -14.960000000000001), (3.6499999999999999, -1.3200000000000001, -21.829999999999998), (0.87, -0.54000000000000004, -23.629999999999999), (-3.0600000000000001, 2.9199999999999999, -18.699999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP53)

#include as follow : execfile('pathto/POP49.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP49= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP49.sph',
radii = [[7.1299999999999999, 6.8399999999999999, 3.6000000000000001, 4.4900000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP49.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP49',
positions = [[(-5.4699999999999998, -1.4099999999999999, -8.9299999999999997), (-0.02, 1.5900000000000001, -21.079999999999998), (6.4800000000000004, 3.3900000000000001, -4.4699999999999998), (4.4900000000000002, 2.48, -12.789999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP49)

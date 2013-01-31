#include as follow : execfile('pathto/POP23.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP23= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP23.sph',
radii = [[4.2800000000000002, 5.6399999999999997, 5.0099999999999998, 6.4699999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP23.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP23',
positions = [[(-2.3199999999999998, 1.45, 17.949999999999999), (6.3200000000000003, -2.8799999999999999, 10.210000000000001), (-5.4100000000000001, 2.7000000000000002, 8.1600000000000001), (0.60999999999999999, -2.8599999999999999, 21.77)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP23)

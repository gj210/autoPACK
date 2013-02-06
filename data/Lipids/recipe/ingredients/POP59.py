#include as follow : execfile('pathto/POP59.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP59= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP59.sph',
radii = [[2.4399999999999999, 2.0600000000000001, 3.1800000000000002, 3.7599999999999998, 3.6899999999999999, 3.0, 1.9199999999999999, 3.6699999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP59.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP59',
positions = [[(3.1400000000000001, 0.33000000000000002, 18.809999999999999), (0.69999999999999996, -0.58999999999999997, 21.899999999999999), (-1.22, -3.3799999999999999, 13.220000000000001), (-3.5699999999999998, -2.71, 5.6900000000000004), (1.3799999999999999, 4.3899999999999997, 4.2599999999999998), (-0.23000000000000001, 1.8799999999999999, 25.399999999999999), (-0.85999999999999999, -3.2000000000000002, 19.190000000000001), (1.4299999999999999, 1.0, 11.74)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP59)

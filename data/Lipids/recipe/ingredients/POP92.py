#include as follow : execfile('pathto/POP92.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP92= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP92.sph',
radii = [[1.6000000000000001, 4.5599999999999996, 4.2300000000000004, 2.73, 2.5299999999999998, 1.95, 2.48, 2.8599999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP92.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP92',
positions = [[(-1.9299999999999999, -5.0899999999999999, 2.02), (-3.4300000000000002, 6.1100000000000003, 6.2599999999999998), (-1.95, 3.7200000000000002, 15.539999999999999), (-0.029999999999999999, -0.62, 19.699999999999999), (2.8300000000000001, -2.9100000000000001, 10.33), (1.05, -3.8700000000000001, 5.5599999999999996), (1.6699999999999999, -2.9199999999999999, 16.34), (2.77, -1.54, 22.129999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP92)

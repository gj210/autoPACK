#include as follow : execfile('pathto/POP68.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP68= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP68.sph',
radii = [[2.7400000000000002, 3.3799999999999999, 3.21, 2.02, 3.1299999999999999, 1.6599999999999999, 3.6200000000000001, 2.2599999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP68.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP68',
positions = [[(1.26, 1.5600000000000001, 14.289999999999999), (1.3899999999999999, 6.2199999999999998, 10.050000000000001), (-1.5800000000000001, -6.6200000000000001, 10.99), (-1.5600000000000001, 0.28000000000000003, 19.120000000000001), (1.53, -3.5699999999999998, 5.1299999999999999), (-3.0699999999999998, -2.3300000000000001, 18.280000000000001), (1.6100000000000001, 6.6399999999999997, 2.4199999999999999), (-0.029999999999999999, -3.48, 16.100000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP68)

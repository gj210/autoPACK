#include as follow : execfile('pathto/POP50.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP50= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP50.sph',
radii = [[5.9299999999999997, 6.1699999999999999, 5.7000000000000002, 2.75]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP50.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP50',
positions = [[(-0.02, 2.8900000000000001, -5.0899999999999999), (-1.1599999999999999, 0.56999999999999995, -16.920000000000002), (-7.75, -1.4299999999999999, -8.3499999999999996), (0.97999999999999998, -2.7000000000000002, -22.079999999999998)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP50)

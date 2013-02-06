#include as follow : execfile('pathto/POP82.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP82= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP82.sph',
radii = [[3.77, 3.3900000000000001, 1.29, 2.3500000000000001, 1.25, 4.2199999999999998, 2.77, 3.2400000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP82.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP82',
positions = [[(6.9699999999999998, -3.4700000000000002, -4.0599999999999996), (-3.0699999999999998, -0.20000000000000001, -13.199999999999999), (-0.20000000000000001, -0.29999999999999999, 0.41999999999999998), (-0.31, -0.54000000000000004, -7.3399999999999999), (-1.22, -0.62, -2.7000000000000002), (3.6800000000000002, -1.03, -11.23), (-3.8300000000000001, 2.1499999999999999, -21.420000000000002), (-1.1000000000000001, 1.97, -16.09)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP82)

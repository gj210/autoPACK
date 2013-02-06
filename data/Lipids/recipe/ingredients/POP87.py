#include as follow : execfile('pathto/POP87.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP87= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP87.sph',
radii = [[1.28, 3.8300000000000001, 4.8300000000000001, 3.3100000000000001, 2.3399999999999999, 1.8899999999999999, 1.6200000000000001, 0.77000000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP87.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP87',
positions = [[(-3.5499999999999998, -6.7999999999999998, 8.0899999999999999), (0.38, 2.9500000000000002, 17.09), (1.3500000000000001, 4.9000000000000004, 7.9199999999999999), (1.8899999999999999, 2.2400000000000002, 23.879999999999999), (-0.14999999999999999, -3.3100000000000001, 22.629999999999999), (-2.3799999999999999, -5.5700000000000003, 18.260000000000002), (-4.1100000000000003, -6.9100000000000001, 14.039999999999999), (-4.0099999999999998, -6.6500000000000004, 10.73)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP87)

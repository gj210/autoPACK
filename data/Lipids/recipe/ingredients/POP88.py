#include as follow : execfile('pathto/POP88.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP88= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP88.sph',
radii = [[4.96, 1.5800000000000001, 2.5, 3.6899999999999999, 1.95, 2.6000000000000001, 1.8999999999999999, 3.21]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP88.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP88',
positions = [[(-6.4900000000000002, -5.1500000000000004, 9.3699999999999992), (-1.54, 0.64000000000000001, 20.739999999999998), (5.8499999999999996, 2.54, 3.5), (-3.7999999999999998, -0.60999999999999999, 16.379999999999999), (8.5299999999999994, 0.040000000000000001, -0.27000000000000002), (5.3799999999999999, 2.6099999999999999, 9.1199999999999992), (0.41999999999999998, 3.0299999999999998, 19.579999999999998), (2.5099999999999998, 1.6399999999999999, 14.66)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP88)

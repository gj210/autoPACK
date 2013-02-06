#include as follow : execfile('pathto/POP85.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP85= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP85.sph',
radii = [[2.6400000000000001, 2.5, 2.96, 3.1699999999999999, 2.98, 3.0, 1.9399999999999999, 3.5]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP85.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP85',
positions = [[(-0.66000000000000003, 0.070000000000000007, 20.149999999999999), (8.1099999999999994, 1.71, 5.4800000000000004), (-4.2999999999999998, -0.12, 16.579999999999998), (3.0899999999999999, 0.19, 15.720000000000001), (-2.8199999999999998, 3.21, 3.96), (-0.39000000000000001, -2.9500000000000002, 23.07), (4.8799999999999999, 0.68000000000000005, 9.8200000000000003), (-3.2999999999999998, -0.26000000000000001, 9.0700000000000003)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP85)

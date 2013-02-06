#include as follow : execfile('pathto/POP71.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP71= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP71.sph',
radii = [[4.7300000000000004, 1.1599999999999999, 0.93999999999999995, 4.8799999999999999, 0.76000000000000001, 3.6800000000000002, 2.9900000000000002, 1.8200000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP71.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP71',
positions = [[(-0.85999999999999999, -1.45, -16.07), (-0.71999999999999997, 2.52, -20.359999999999999), (-1.6599999999999999, 2.6499999999999999, -23.010000000000002), (-0.42999999999999999, -0.17000000000000001, -6.75), (-0.38, 1.2, -22.98), (4.71, -2.3700000000000001, -11.26), (-1.8500000000000001, 2.3300000000000001, -0.59999999999999998), (0.13, 1.3700000000000001, -19.399999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP71)

#include as follow : execfile('pathto/POP19.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP19= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP19.sph',
radii = [[4.04, 2.98, 1.9199999999999999, 2.8700000000000001, 3.1600000000000001, 2.5299999999999998, 3.1699999999999999, 1.3300000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP19.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP19',
positions = [[(0.26000000000000001, -1.3999999999999999, 19.399999999999999), (-2.1299999999999999, -3.8300000000000001, 21.809999999999999), (0.040000000000000001, 7.9299999999999997, 3.2000000000000002), (-0.85999999999999999, 5.1600000000000001, 8.0899999999999999), (0.65000000000000002, -3.7200000000000002, 13.18), (2.5099999999999998, -1.6899999999999999, 6.8499999999999996), (-0.42999999999999999, 4.0199999999999996, 14.85), (3.1000000000000001, -0.34000000000000002, 2.02)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP19)

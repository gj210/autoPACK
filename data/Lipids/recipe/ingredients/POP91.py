#include as follow : execfile('pathto/POP91.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP91= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP91.sph',
radii = [[2.3799999999999999, 4.21, 1.79, 3.52, 2.0099999999999998, 1.9099999999999999, 2.79, 2.54]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP91.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP91',
positions = [[(-2.25, 3.9300000000000002, 12.789999999999999), (-6.4400000000000004, 6.7599999999999998, 10.68), (-0.28999999999999998, -1.48, 21.77), (6.2000000000000002, -6.04, 4.6299999999999999), (-1.5, 0.059999999999999998, 19.399999999999999), (0.17000000000000001, 1.1100000000000001, 15.73), (5.7400000000000002, -4.7199999999999998, 11.449999999999999), (2.3399999999999999, -3.3999999999999999, 16.059999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP91)

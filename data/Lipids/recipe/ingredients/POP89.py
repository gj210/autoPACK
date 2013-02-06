#include as follow : execfile('pathto/POP89.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP89= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP89.sph',
radii = [[4.5899999999999999, 3.5800000000000001, 2.0, 4.2000000000000002, 4.6200000000000001, 0.81999999999999995, 1.24, 0.0]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP89.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP89',
positions = [[(0.54000000000000004, 0.23999999999999999, 12.09), (2.3399999999999999, -1.48, 17.510000000000002), (-0.25, -1.3799999999999999, 21.379999999999999), (-2.8700000000000001, -3.4900000000000002, 6.0899999999999999), (-0.48999999999999999, 4.8200000000000003, 4.0899999999999999), (0.16, 1.6599999999999999, 22.73), (-0.64000000000000001, 0.51000000000000001, 23.34), (1.1899999999999999, 2.04, 23.460000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP89)

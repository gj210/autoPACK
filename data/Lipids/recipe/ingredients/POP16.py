#include as follow : execfile('pathto/POP16.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP16= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP16.sph',
radii = [[3.8700000000000001, 2.3300000000000001, 3.0099999999999998, 2.3999999999999999, 3.8599999999999999, 2.8700000000000001, 1.6000000000000001, 1.8999999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP16.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP16',
positions = [[(1.1499999999999999, -2.9500000000000002, 1.96), (0.59999999999999998, 2.6299999999999999, 18.77), (-2.8599999999999999, 0.34999999999999998, 4.2800000000000002), (2.21, -0.29999999999999999, 15.289999999999999), (-0.12, 0.25, 9.8300000000000001), (0.64000000000000001, 4.3799999999999999, 14.69), (-0.23999999999999999, -3.0099999999999998, 18.460000000000001), (-1.3100000000000001, -0.68999999999999995, 20.719999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP16)

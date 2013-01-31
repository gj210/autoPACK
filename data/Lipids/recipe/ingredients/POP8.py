#include as follow : execfile('pathto/POP8.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP8= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP8.sph',
radii = [[5.1399999999999997, 4.2400000000000002, 3.8700000000000001, 5.0300000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP8.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP8',
positions = [[(2.9300000000000002, 2.1400000000000001, -10.66), (7.1100000000000003, 1.5700000000000001, -2.0600000000000001), (-0.76000000000000001, 4.1600000000000001, -4.0099999999999998), (0.059999999999999998, 1.48, -18.57)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP8)

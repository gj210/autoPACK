#include as follow : execfile('pathto/POP12.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP12= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP12.sph',
radii = [[2.4300000000000002, 3.9199999999999999, 3.0800000000000001, 1.3200000000000001, 0.76000000000000001, 1.99, 5.8799999999999999, 1.5800000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP12.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP12',
positions = [[(3.1699999999999999, 0.95999999999999996, -18.710000000000001), (-2.6400000000000001, -1.9299999999999999, -16.07), (-2.1600000000000001, 1.46, -21.539999999999999), (4.2400000000000002, 2.0699999999999998, -10.32), (4.54, 2.25, -7.3499999999999996), (4.7800000000000002, 3.4300000000000002, -4.2800000000000002), (-2.8100000000000001, -3.9300000000000002, -6.0199999999999996), (5.5700000000000003, 2.0299999999999998, -14.109999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP12)

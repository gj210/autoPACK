#include as follow : execfile('pathto/POP13.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP13= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP13.sph',
radii = [[3.0, 3.1800000000000002, 3.1200000000000001, 2.6099999999999999, 2.5600000000000001, 1.8500000000000001, 2.8199999999999998, 3.1899999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP13.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP13',
positions = [[(-1.96, 9.2699999999999996, 5.5300000000000002), (-1.8300000000000001, 7.1900000000000004, 12.17), (-0.76000000000000001, 4.6399999999999997, 18.920000000000002), (-0.57999999999999996, -0.68000000000000005, 18.539999999999999), (1.3999999999999999, -8.9900000000000002, 9.2400000000000002), (-1.77, -9.25, 5.9000000000000004), (1.98, -1.5700000000000001, 21.109999999999999), (1.73, -4.5300000000000002, 13.23)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP13)

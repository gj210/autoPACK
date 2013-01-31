#include as follow : execfile('pathto/POP10.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP10= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP10.sph',
radii = [[6.7000000000000002, 4.0999999999999996, 4.8099999999999996, 4.2000000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP10.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP10',
positions = [[(0.35999999999999999, -7.1600000000000001, -9.8399999999999999), (-3.1499999999999999, -3.7799999999999998, -15.380000000000001), (1.3700000000000001, -0.19, -21.09), (4.8200000000000003, -1.3300000000000001, -13.050000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP10)

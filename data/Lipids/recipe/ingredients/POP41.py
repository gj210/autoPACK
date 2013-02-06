#include as follow : execfile('pathto/POP41.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP41= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP41.sph',
radii = [[3.2400000000000002, 3.1099999999999999, 3.0499999999999998, 2.9199999999999999, 1.96, 2.8999999999999999, 2.5, 2.29]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP41.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP41',
positions = [[(-5.0199999999999996, -3.23, 5.96), (-0.79000000000000004, -0.81000000000000005, 17.649999999999999), (1.1599999999999999, 0.83999999999999997, 20.800000000000001), (0.62, 4.54, 24.760000000000002), (1.8, 0.28999999999999998, 14.279999999999999), (3.0600000000000001, -3.4500000000000002, 4.75), (-2.9100000000000001, -1.9099999999999999, 12.380000000000001), (0.93999999999999995, -1.1399999999999999, 9.3300000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP41)

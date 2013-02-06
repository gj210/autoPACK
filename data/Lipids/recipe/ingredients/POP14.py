#include as follow : execfile('pathto/POP14.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP14= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP14.sph',
radii = [[2.4100000000000001, 4.0800000000000001, 2.5600000000000001, 1.3500000000000001, 2.5, 2.6699999999999999, 2.2599999999999998, 3.1099999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP14.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP14',
positions = [[(2.0099999999999998, 6.1399999999999997, 1.73), (-2.77, -0.84999999999999998, 21.539999999999999), (1.3600000000000001, 3.1899999999999999, 12.130000000000001), (2.9199999999999999, -5.6799999999999997, 3.5600000000000001), (1.0900000000000001, -5.5099999999999998, 8.2100000000000009), (0.48999999999999999, 2.0299999999999998, 18.370000000000001), (2.5800000000000001, 4.6500000000000004, 6.9500000000000002), (-0.51000000000000001, -4.3099999999999996, 14.41)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP14)

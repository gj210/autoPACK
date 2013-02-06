#include as follow : execfile('pathto/POP76.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP76= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP76.sph',
radii = [[0.73999999999999999, 3.2000000000000002, 1.3200000000000001, 2.9100000000000001, 3.6600000000000001, 5.1100000000000003, 4.0999999999999996, 2.6699999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP76.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP76',
positions = [[(-1.4299999999999999, -0.81000000000000005, -23.0), (6.5099999999999998, 8.3000000000000007, -7.5800000000000001), (-3.21, 2.6699999999999999, -23.850000000000001), (5.2599999999999998, 3.8599999999999999, -12.32), (0.23999999999999999, 1.6699999999999999, -15.82), (-2.1099999999999999, -6.04, -7.2199999999999998), (-2.5299999999999998, -4.4900000000000002, -16.789999999999999), (-1.6100000000000001, 0.45000000000000001, -21.280000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP76)

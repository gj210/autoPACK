#include as follow : execfile('pathto/LPO122.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO122= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO122.sph',
radii = [[3.1000000000000001, 2.9300000000000002, 2.5899999999999999, 3.5800000000000001, 2.5499999999999998, 2.0, 2.25, 2.6000000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO122.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO122',
positions = [[(-3.77, -2.3900000000000001, 18.420000000000002), (4.4800000000000004, -1.72, 16.530000000000001), (-3.9199999999999999, -1.1399999999999999, 11.789999999999999), (1.02, -0.070000000000000007, 22.190000000000001), (-2.02, 1.46, 6.3600000000000003), (1.8999999999999999, -1.6699999999999999, 25.23), (2.1000000000000001, 1.3899999999999999, 13.44), (-0.34000000000000002, 5.5599999999999996, 15.220000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO122)

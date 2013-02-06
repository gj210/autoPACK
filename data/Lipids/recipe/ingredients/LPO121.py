#include as follow : execfile('pathto/LPO121.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO121= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO121.sph',
radii = [[2.5699999999999998, 1.6299999999999999, 4.1699999999999999, 3.2400000000000002, 1.54, 2.4100000000000001, 2.7000000000000002, 5.9900000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO121.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'LPO121',
positions = [[(-4.8300000000000001, -1.95, 12.289999999999999), (-5.4100000000000001, -2.8399999999999999, 17.440000000000001), (2.29, 0.23999999999999999, 15.73), (-4.5199999999999996, -0.45000000000000001, 5.6699999999999999), (3.3500000000000001, -4.6399999999999997, 25.66), (-1.8400000000000001, -2.0499999999999998, 20.010000000000002), (1.02, -1.6899999999999999, 24.329999999999998), (3.9100000000000001, 6.1799999999999997, 7.46)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO121)

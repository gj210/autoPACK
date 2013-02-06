#include as follow : execfile('pathto/DPO102.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO102= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO102.sph',
radii = [[3.0499999999999998, 2.0, 3.6400000000000001, 2.29, 3.6499999999999999, 3.0299999999999998, 1.95, 2.79]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/DPO102.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'DPO102',
positions = [[(-1.8500000000000001, 3.3900000000000001, -16.260000000000002), (2.71, -4.8600000000000003, -3.6699999999999999), (-2.8900000000000001, -0.48999999999999999, -3.7000000000000002), (1.99, -2.6200000000000001, -12.1), (1.71, 0.93999999999999995, -23.23), (2.0600000000000001, -0.25, -17.75), (1.3400000000000001, -2.3599999999999999, -7.1399999999999997), (-4.9299999999999997, 1.8700000000000001, -10.07)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO102)

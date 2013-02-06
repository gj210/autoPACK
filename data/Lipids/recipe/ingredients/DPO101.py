#include as follow : execfile('pathto/DPO101.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO101= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO101.sph',
radii = [[3.0, 1.9099999999999999, 3.1600000000000001, 3.3300000000000001, 2.4300000000000002, 2.9500000000000002, 2.5499999999999998, 2.2999999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/DPO101.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'DPO101',
positions = [[(-5.0, -2.3500000000000001, -5.9900000000000002), (-0.01, -0.73999999999999999, -19.440000000000001), (3.7400000000000002, 0.28999999999999998, -9.2400000000000002), (2.8599999999999999, 2.96, -16.079999999999998), (-4.4199999999999999, -2.8500000000000001, -11.92), (0.81000000000000005, 1.3899999999999999, -22.149999999999999), (3.1099999999999999, -0.78000000000000003, -3.1299999999999999), (-3.8500000000000001, -1.5, -17.210000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO101)

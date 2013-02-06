#include as follow : execfile('pathto/DPO97.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO97= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO97.sph',
radii = [[2.3300000000000001, 2.5299999999999998, 3.3900000000000001, 2.5800000000000001, 2.52, 2.3399999999999999, 2.9900000000000002, 3.0499999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/DPO97.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'DPO97',
positions = [[(-1.26, 0.19, -17.600000000000001), (-2.8999999999999999, -5.4100000000000001, -1.5), (2.7599999999999998, 3.9700000000000002, -11.1), (-1.27, -0.33000000000000002, -11.539999999999999), (-1.5700000000000001, -2.29, -6.2199999999999998), (2.5, 1.29, -17.41), (1.0900000000000001, 2.4399999999999999, -4.9000000000000004), (-0.35999999999999999, -0.85999999999999999, -23.41)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO97)

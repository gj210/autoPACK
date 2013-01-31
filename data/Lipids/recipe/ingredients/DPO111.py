#include as follow : execfile('pathto/DPO111.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO111= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO111.sph',
radii = [[7.0700000000000003, 4.6500000000000004, 3.8399999999999999, 6.4800000000000004]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/DPO111.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'DPO111',
positions = [[(8.8399999999999999, 6.2400000000000002, -8.3200000000000003), (-4.5800000000000001, -0.46999999999999997, -11.779999999999999), (-9.0299999999999994, 0.42999999999999999, -3.3399999999999999), (0.93000000000000005, 2.3599999999999999, -20.120000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO111)

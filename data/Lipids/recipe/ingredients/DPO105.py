#include as follow : execfile('pathto/DPO105.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO105= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO105.sph',
radii = [[2.3399999999999999, 2.5600000000000001, 2.9700000000000002, 2.1499999999999999, 3.46, 2.29, 3.0099999999999998, 3.8700000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/DPO105.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'DPO105',
positions = [[(3.5099999999999998, 1.51, -8.8200000000000003), (2.9700000000000002, 1.47, -14.68), (-4.7300000000000004, -2.6699999999999999, -13.15), (-2.21, 4.7400000000000002, -20.66), (0.38, 1.6000000000000001, -21.23), (-2.1000000000000001, -2.2200000000000002, -18.82), (6.5599999999999996, -0.98999999999999999, -4.0700000000000003), (-4.04, -3.52, -5.6900000000000004)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO105)

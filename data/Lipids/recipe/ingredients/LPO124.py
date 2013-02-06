#include as follow : execfile('pathto/LPO124.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO124= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/LPO124.sph',
radii = [[2.1000000000000001, 1.26, 3.54, 2.0800000000000001, 3.0699999999999998, 4.0800000000000001, 2.3900000000000001, 3.23]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/LPO124.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO124',
positions = [[(-2.3599999999999999, -3.4500000000000002, -23.710000000000001), (-0.94999999999999996, 4.4000000000000004, -3.9399999999999999), (2.8300000000000001, 1.23, -12.16), (-1.1000000000000001, -2.8599999999999999, -27.039999999999999), (4.75, -1.29, -5.2800000000000002), (-0.5, 0.60999999999999999, -20.300000000000001), (-0.54000000000000004, 2.3399999999999999, -7.7999999999999998), (-2.4100000000000001, 1.1599999999999999, -14.16)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO124)

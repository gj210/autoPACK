#include as follow : execfile('pathto/DPO112.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO112= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO112.sph',
radii = [[2.8599999999999999, 1.6000000000000001, 2.0099999999999998, 3.3700000000000001, 3.8799999999999999, 3.6000000000000001, 2.3999999999999999, 1.95]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/DPO112.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'DPO112',
positions = [[(0.70999999999999996, 0.050000000000000003, 16.530000000000001), (-2.5600000000000001, -1.1000000000000001, 10.25), (-0.40999999999999998, -6.9000000000000004, 2.27), (5.4800000000000004, 2.1899999999999999, 1.6499999999999999), (-2.6400000000000001, 2.1099999999999999, 20.920000000000002), (3.0499999999999998, 1.8, 9.2799999999999994), (-1.8300000000000001, -4.4100000000000001, 6.9699999999999998), (-3.1600000000000001, -0.28999999999999998, 15.01)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO112)

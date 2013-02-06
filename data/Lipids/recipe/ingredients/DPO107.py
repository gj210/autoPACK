#include as follow : execfile('pathto/DPO107.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO107= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO107.sph',
radii = [[3.2799999999999998, 1.98, 3.54, 3.1400000000000001, 2.0299999999999998, 1.9299999999999999, 3.4199999999999999, 3.73]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/DPO107.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'DPO107',
positions = [[(2.5299999999999998, -1.45, 18.559999999999999), (3.5, 0.54000000000000004, 24.559999999999999), (-1.76, -0.55000000000000004, 14.130000000000001), (2.1000000000000001, 1.8899999999999999, 8.8399999999999999), (0.52000000000000002, 1.1100000000000001, 27.34), (1.53, 3.1699999999999999, 2.75), (-0.81000000000000005, 0.35999999999999999, 20.899999999999999), (-5.79, -3.3399999999999999, 7.8799999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO107)

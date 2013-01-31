#include as follow : execfile('pathto/DPO101.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO101= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO101.sph',
radii = [[4.8899999999999997, 6.4400000000000004, 4.2800000000000002, 4.8899999999999997]],
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
positions = [[(1.77, 2.3199999999999998, -14.26), (-1.1899999999999999, -1.6299999999999999, -5.0800000000000001), (-0.94999999999999996, 0.46000000000000002, -21.129999999999999), (-5.7699999999999996, -2.6899999999999999, -13.52)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO101)

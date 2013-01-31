#include as follow : execfile('pathto/DPO97.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO97= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO97.sph',
radii = [[5.3799999999999999, 4.6399999999999997, 3.0499999999999998, 5.21]],
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
positions = [[(-2.0600000000000001, -2.3900000000000001, -3.3300000000000001), (1.3100000000000001, 3.48, -8.4299999999999997), (-0.53000000000000003, -0.040000000000000001, -23.41), (0.28999999999999998, 1.6599999999999999, -16.219999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO97)

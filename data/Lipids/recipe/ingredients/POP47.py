#include as follow : execfile('pathto/POP47.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP47= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP47.sph',
radii = [[5.4199999999999999, 6.5599999999999996, 4.0700000000000003, 5.46]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP47.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP47',
positions = [[(-5.1900000000000004, 0.53000000000000003, -5.6500000000000004), (2.8599999999999999, -0.76000000000000001, -20.170000000000002), (-2.73, 0.050000000000000003, -15.57), (6.0499999999999998, -0.81000000000000005, -8.1300000000000008)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP47)

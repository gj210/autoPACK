#include as follow : execfile('pathto/POP18.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP18= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP18.sph',
radii = [[5.8700000000000001, 2.5800000000000001, 3.52, 4.9100000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP18.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP18',
positions = [[(-0.17000000000000001, 0.20000000000000001, 0.67000000000000004), (-0.34999999999999998, -2.0600000000000001, 19.329999999999998), (0.40999999999999998, -0.84999999999999998, 14.609999999999999), (0.17000000000000001, -2.1000000000000001, 8.0)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP18)

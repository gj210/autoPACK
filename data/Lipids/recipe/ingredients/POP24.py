#include as follow : execfile('pathto/POP24.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP24= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP24.sph',
radii = [[2.4100000000000001, 3.2200000000000002, 3.6000000000000001, 2.02, 3.9399999999999999, 1.29, 3.5499999999999998, 1.6599999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP24.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP24',
positions = [[(2.6600000000000001, 1.9399999999999999, 8.0), (-6.71, -3.3599999999999999, 5.8200000000000003), (-3.6699999999999999, -1.27, 12.75), (-1.0600000000000001, 0.87, 21.120000000000001), (1.8200000000000001, -1.0700000000000001, 18.890000000000001), (0.52000000000000002, 0.42999999999999999, 3.8599999999999999), (6.0999999999999996, 0.46000000000000002, 13.140000000000001), (-1.1799999999999999, 3.7400000000000002, 22.629999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP24)

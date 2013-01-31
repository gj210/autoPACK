#include as follow : execfile('pathto/POP19.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP19= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP19.sph',
radii = [[3.7599999999999998, 4.8899999999999997, 6.3499999999999996, 6.3899999999999997]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP19.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP19',
positions = [[(-1.26, 7.8399999999999999, 4.6500000000000004), (-1.0900000000000001, 5.2000000000000002, 13.1), (1.6599999999999999, -0.80000000000000004, 6.7999999999999998), (-1.4199999999999999, -1.78, 19.82)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP19)

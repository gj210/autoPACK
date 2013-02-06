#include as follow : execfile('pathto/POP2.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP2= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP2.sph',
radii = [[1.6899999999999999, 3.4300000000000002, 3.23, 1.95, 2.73, 1.97, 1.8899999999999999, 4.0099999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP2.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP2',
positions = [[(-5.25, 3.3999999999999999, -2.04), (3.1000000000000001, -1.97, -12.75), (-4.6100000000000003, 0.84999999999999998, -11.56), (-5.9199999999999999, 3.1499999999999999, -6.7699999999999996), (1.97, 0.80000000000000004, -22.329999999999998), (3.1699999999999999, -2.29, -6.1500000000000004), (6.1600000000000001, -1.45, -2.8199999999999998), (-0.34000000000000002, -1.1100000000000001, -18.0)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP2)

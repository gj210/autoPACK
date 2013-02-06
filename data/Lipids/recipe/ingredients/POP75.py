#include as follow : execfile('pathto/POP75.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP75= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP75.sph',
radii = [[3.71, 1.8999999999999999, 4.0, 1.72, 2.9300000000000002, 2.8500000000000001, 1.9399999999999999, 2.5899999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP75.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP75',
positions = [[(0.23999999999999999, 0.91000000000000003, -17.48), (-1.8500000000000001, -3.0899999999999999, -13.880000000000001), (4.2300000000000004, 5.46, -6.6699999999999999), (-1.72, -4.5899999999999999, -5.3099999999999996), (-2.0899999999999999, -0.01, -21.84), (4.46, 4.9500000000000002, -13.91), (-2.0499999999999998, -3.9399999999999999, -9.5700000000000003), (-3.52, -7.1799999999999997, -0.97999999999999998)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP75)

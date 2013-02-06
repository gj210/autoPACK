#include as follow : execfile('pathto/POP25.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP25= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP25.sph',
radii = [[1.95, 2.8199999999999998, 2.7200000000000002, 1.1599999999999999, 3.8700000000000001, 3.8500000000000001, 2.98, 3.6600000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP25.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP25',
positions = [[(-2.9300000000000002, -0.39000000000000001, -15.31), (-1.75, 4.5099999999999998, -20.539999999999999), (2.8199999999999998, 2.4399999999999999, -17.280000000000001), (-0.64000000000000001, 0.51000000000000001, -17.960000000000001), (-4.5199999999999996, -5.8099999999999996, -4.3099999999999996), (5.3399999999999999, -0.12, -3.3599999999999999), (-5.0700000000000003, -3.7599999999999998, -11.529999999999999), (5.4000000000000004, -0.37, -11.42)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP25)

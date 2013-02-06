#include as follow : execfile('pathto/POP93.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP93= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP93.sph',
radii = [[2.9300000000000002, 1.6200000000000001, 1.6399999999999999, 4.5199999999999996, 3.2799999999999998, 1.3100000000000001, 3.8700000000000001, 2.9500000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP93.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP93',
positions = [[(2.6400000000000001, 2.8599999999999999, 18.57), (3.2599999999999998, 6.4400000000000004, 14.0), (3.27, 2.3100000000000001, 8.3800000000000008), (-4.2300000000000004, -2.54, 3.8799999999999999), (-1.0900000000000001, -0.81999999999999995, 19.059999999999999), (4.0599999999999996, 5.7699999999999996, 10.359999999999999), (-5.6100000000000003, -2.1000000000000001, 12.460000000000001), (3.2999999999999998, -2.4500000000000002, 21.399999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP93)

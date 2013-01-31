#include as follow : execfile('pathto/POP42.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP42= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP42.sph',
radii = [[5.6799999999999997, 7.1799999999999997, 4.5899999999999999, 5.6200000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP42.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP42',
positions = [[(-0.34000000000000002, -0.66000000000000003, 23.140000000000001), (-1.1699999999999999, -8.5199999999999996, 12.119999999999999), (-2.4399999999999999, 4.9400000000000004, 18.719999999999999), (-7.3799999999999999, 6.54, 10.02)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP42)

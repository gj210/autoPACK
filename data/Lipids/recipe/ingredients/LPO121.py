#include as follow : execfile('pathto/LPO121.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO121= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/LPO121.sph',
radii = [[6.1200000000000001, 7.0199999999999996, 3.8599999999999999, 5.04]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/LPO121.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'LPO121',
positions = [[(-1.76, 0.47999999999999998, 17.710000000000001), (3.77, 7.5700000000000003, 8.2400000000000002), (1.72, -0.51000000000000001, 24.739999999999998), (-4.4900000000000002, 1.1499999999999999, 7.4800000000000004)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO121)

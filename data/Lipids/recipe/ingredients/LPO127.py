#include as follow : execfile('pathto/LPO127.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO127= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/LPO127.sph',
radii = [[5.0800000000000001, 6.3899999999999997, 6.1799999999999997, 5.0999999999999996]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/LPO127.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO127',
positions = [[(2.1200000000000001, -1.76, 22.379999999999999), (4.2300000000000004, 0.48999999999999999, 5.2000000000000002), (2.96, 5.8600000000000003, 15.039999999999999), (0.52000000000000002, -0.46999999999999997, 14.48)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO127)

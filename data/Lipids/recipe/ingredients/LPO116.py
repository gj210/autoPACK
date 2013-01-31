#include as follow : execfile('pathto/LPO116.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO116= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/LPO116.sph',
radii = [[4.7400000000000002, 5.6600000000000001, 6.7300000000000004, 4.75]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/LPO116.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO116',
positions = [[(-0.59999999999999998, -4.1200000000000001, 19.789999999999999), (2.6899999999999999, 3.0600000000000001, 7.4199999999999999), (2.1499999999999999, 3.3100000000000001, 19.960000000000001), (0.12, -5.4299999999999997, 9.9299999999999997)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO116)

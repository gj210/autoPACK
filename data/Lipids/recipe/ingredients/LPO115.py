#include as follow : execfile('pathto/LPO115.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO115= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/LPO115.sph',
radii = [[5.7400000000000002, 2.9399999999999999, 3.4100000000000001, 4.5499999999999998]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/LPO115.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'LPO115',
positions = [[(-1.1100000000000001, 0.94999999999999996, -3.0800000000000001), (0.39000000000000001, -2.46, -22.73), (-0.52000000000000002, 0.72999999999999998, -17.460000000000001), (-2.1499999999999999, 0.73999999999999999, -10.9)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO115)

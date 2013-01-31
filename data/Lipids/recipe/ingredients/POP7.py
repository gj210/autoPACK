#include as follow : execfile('pathto/POP7.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP7= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/POP7.sph',
radii = [[6.2999999999999998, 5.8600000000000003, 6.1500000000000004, 5.3099999999999996]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/POP7.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP7',
positions = [[(1.3700000000000001, -3.8300000000000001, -8.3000000000000007), (-1.1100000000000001, 0.56000000000000005, -17.239999999999998), (-5.5099999999999998, 2.2999999999999998, -6.6399999999999997), (2.1899999999999999, -0.33000000000000002, -23.920000000000002)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP7)

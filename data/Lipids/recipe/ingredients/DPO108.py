#include as follow : execfile('pathto/DPO108.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO108= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/DPO108.sph',
radii = [[5.7400000000000002, 7.1299999999999999, 2.8399999999999999, 2.7799999999999998]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/DPO108.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'DPO108',
positions = [[(0.20000000000000001, -0.040000000000000001, -16.350000000000001), (0.80000000000000004, -2.2999999999999998, -5.5099999999999998), (-2.1499999999999999, 0.26000000000000001, -22.18), (-4.5499999999999998, 2.0699999999999998, -20.18)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO108)

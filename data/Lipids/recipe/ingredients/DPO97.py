#include as follow : execfile('pathto/DPO97.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO97= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/DPO97.sph',
radii = [[5.7300000000000004, 4.2599999999999998, 3.0499999999999998, 4.2699999999999996]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/DPO97.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'DPO97',
positions = [[(0.85999999999999999, 2.8700000000000001, -7.8200000000000003), (0.40999999999999998, 1.75, -16.510000000000002), (-0.53000000000000003, -0.040000000000000001, -23.41), (-2.79, -3.77, -3.0299999999999998)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO97)

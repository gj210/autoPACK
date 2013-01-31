#include as follow : execfile('pathto/POP77.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP77= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/POP77.sph',
radii = [[5.9800000000000004, 5.2800000000000002, 6.6299999999999999, 3.5]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/POP77.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP77',
positions = [[(9.9800000000000004, -0.029999999999999999, -7.6500000000000004), (-3.1000000000000001, -1.6699999999999999, -5.7400000000000002), (-0.26000000000000001, 0.050000000000000003, -16.059999999999999), (0.48999999999999999, -1.4199999999999999, -23.289999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP77)

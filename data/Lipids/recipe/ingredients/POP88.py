#include as follow : execfile('pathto/POP88.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP88= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/POP88.sph',
radii = [[5.4199999999999999, 6.4400000000000004, 5.6399999999999997, 5.4299999999999997]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/POP88.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP88',
positions = [[(-6.9299999999999997, -6.4699999999999998, 9.8599999999999994), (6.0800000000000001, 0.34000000000000002, 3.8799999999999999), (-3.0499999999999998, -1.23, 18.809999999999999), (1.3799999999999999, 0.48999999999999999, 15.98)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP88)

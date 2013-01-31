#include as follow : execfile('pathto/POP85.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP85= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/POP85.sph',
radii = [[6.7000000000000002, 4.2699999999999996, 5.8600000000000003, 4.9800000000000004]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/POP85.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP85',
positions = [[(8.1099999999999994, 1.1299999999999999, 9.0299999999999994), (-1.8799999999999999, 0.040000000000000001, 15.529999999999999), (2.2000000000000002, -1.3400000000000001, 21.219999999999999), (-0.78000000000000003, 1.6699999999999999, 5.7800000000000002)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP85)

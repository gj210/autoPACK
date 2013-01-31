#include as follow : execfile('pathto/POP21.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP21= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/POP21.sph',
radii = [[6.9299999999999997, 4.8300000000000001, 4.7400000000000002, 6.1100000000000003]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/POP21.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP21',
positions = [[(-0.52000000000000002, -5.75, 9.5099999999999998), (-4.5199999999999996, 7.3899999999999997, 7.8700000000000001), (-0.93999999999999995, 4.4299999999999997, 17.550000000000001), (-0.84999999999999998, -0.91000000000000003, 22.449999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP21)

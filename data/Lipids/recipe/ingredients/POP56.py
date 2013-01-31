#include as follow : execfile('pathto/POP56.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP56= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/POP56.sph',
radii = [[3.1499999999999999, 4.8200000000000003, 5.4400000000000004, 4.5599999999999996]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/POP56.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP56',
positions = [[(3.5600000000000001, 2.7000000000000002, -22.140000000000001), (-0.29999999999999999, -1.9099999999999999, -19.289999999999999), (-0.28999999999999998, 0.13, -9.8800000000000008), (7.3200000000000003, -1.6399999999999999, -12.109999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP56)

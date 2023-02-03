include("../src/NARDS.jl")

using .NARDS, Test

println("Running NARDS tests...")

# Gamma tubulin test with real IDR.
z_scores = @time nards(collect("AAEQDSYLDDVLVDDENMVGELEEDLDADGDHKLV"), num_null_models=100_000) 
display(z_scores)

# Simple tests to make sure it isn't doing anything super dumb. 
z_scores = @time nards(collect("AAAAAAAAAAAAAAAAAAAAAAAPPPPPPPPPPPPPPPPPPPPPPPPP"), num_null_models=100_000)
display(z_scores)
@test z_scores[6, 6] ≈ 15.2 atol=0.2
@test z_scores[7, 7] ≈ 15.2 atol=0.2
@test z_scores[6, 7] ≈ 15.2 atol=0.2

z_scores = @time nards(collect("APAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAP"), num_null_models=100_000)
display(z_scores)
@test z_scores[6, 6] ≈ -1.46 atol=0.2
@test z_scores[7, 7] ≈ -1.46 atol=0.2
@test z_scores[6, 7] ≈ -1.46 atol=0.2
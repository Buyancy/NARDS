include("../src/NARDS.jl")

using Documenter
using .NARDS

makedocs(
    sitename = "NARDS",
    format = Documenter.HTML(),
    modules = [NARDS]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

using Documenter
using CropGrowth

makedocs(sitename = "CropGrowth", format = Documenter.HTML(), modules = [CropGrowth])

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

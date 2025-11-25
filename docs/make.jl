using CropGrowth
using Documenter

DocMeta.setdocmeta!(CropGrowth, :DocTestSetup, :(using CropGrowth); recursive = true)

makedocs(;
    modules = [CropGrowth],
    authors = "jeffersonparil@gmail.com",
    sitename = "CropGrowth.jl",
    format = Documenter.HTML(;
        canonical = "https://jeffersonfparil.github.io/CropGrowth.jl",
        edit_link = "main",
        assets = String[],
        size_threshold = 1000000,
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/jeffersonfparil/CropGrowth.jl", devbranch = "main")

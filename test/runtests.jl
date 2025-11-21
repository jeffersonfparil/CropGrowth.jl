using CropGrowth
using Test
using Documenter, DocumenterTools

try
    Documenter.doctest(CropGrowth)
catch
    DocumenterTools.generate()
    Pkg.activate("docs/")
    Documenter.doctest(CropGrowth)
end

@testset "CropGrowth.jl" begin
    @test 1 == 1
end

using Pkg
# Develop the local package first
Pkg.develop(PackageSpec(; path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter, Literate
using CircleMedianFilter

DocMeta.setdocmeta!(
    CircleMedianFilter, :DocTestSetup, :(using CircleMedianFilter); recursive=true
)

# Generate tutorials from example files
@info "current dir = $(@__DIR__)"
tutorials_folder = (@__DIR__) * "/../examples"
docs_tutorials_folder = (@__DIR__) * "/src/examples"
@info "Processing examples from: $tutorials_folder"

# Process example files with Literate.jl if examples folder exists
if isdir(tutorials_folder)
    for f in readdir(tutorials_folder; join=true)
        if endswith(f, ".jl")
            @info "Processing example: $f"
            Literate.markdown(f, docs_tutorials_folder)
        end
    end
end

makedocs(;
    sitename="CircleMedianFilter.jl",
    modules=[CircleMedianFilter],
    authors="Oleg Soloviev",
    repo="https://github.com/olejorik/CircleMedianFilter.jl/blob/{commit}{path}#L{line}",
    checkdocs=:exports,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://olejorik.github.io/CircleMedianFilter.jl/stable/",
        assets=["assets/favicon.ico"],
        highlights=["yaml"],
    ),
    clean=false,
    pages=[
        "Home" => "index.md",
        "About" => "about.md",
        "API Reference" => "api.md",
        "Examples" => ["examples/circular_median_demo.md"],
    ],
)

# Uncomment for deployment to GitHub Pages
# deploydocs(
#     repo = "github.com/olejorik/CircleMedianFilter.jl.git",
# )

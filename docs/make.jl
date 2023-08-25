using Documenter
using OrnsteinZernike

ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"

push!(LOAD_PATH,"../src/")
makedocs(sitename="OrnsteinZernike.jl Documentation",
         pages = [
            "Index" => "index.md",
            "Tutorials" => ["SingleCompLJ.md",
                            "HardSphereMixture.md",
                            "HighDensities.md",
                            "OtherDimensions.md",
                            "Accuracy.md",
                            "ThermodynamicConsistency.md"
                            ],
            "Basics"=> ["GeneralWorkflow.md",
                        "Potentials.md",
                        "Systems.md",
                        "Closures.md",
                        "Solvers.md"],
            "Theory"=>"Theory.md",
            "Extending" => ["ExtendingClosures.md", 
                            "ExtendingPotentials.md"
                            ],
            "API" => "API.md",
            "From Python" => "FromPython.md"
         ],
         format = Documenter.HTML(prettyurls = false)
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/IlianPihlajamaa/OrnsteinZernike.jl.git",
    devbranch = "main"
)

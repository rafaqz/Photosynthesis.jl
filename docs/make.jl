using Documenter, Photosynthesis

makedocs(
    modules = [Photosynthesis],
    doctest = false,
    clean = false,
    sitename = "Photosynthesis.jl",
    format = :html,
    pages = Any[
        "Introduction" => "index.md",
    ]
)

# deploydocs(
#     repo = "github.com/rafaqz/Photosynthesis.jl.git",
#     osname = "linux",
#     julia = "0.6",
#     target = "build",
#     deps = nothing,
#     make = nothing
# )

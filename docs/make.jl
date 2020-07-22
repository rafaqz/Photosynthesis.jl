using Documenter, Photosynthesis

makedocs(
    modules = [Photosynthesis],
    sitename = "Photosynthesis.jl",
    pages = Any[
        "Home" => "index.md",
    ],
    clean = false,
)

deploydocs(
    repo = "github.com/rafaqz/Photosynthesis.jl.git",
)

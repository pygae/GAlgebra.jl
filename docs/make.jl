push!(LOAD_PATH,"../")
using Documenter, GAlgebra

root = normpath(joinpath(dirname(pathof(GAlgebra)), ".."))
src = joinpath(root, "README.md")
dst = joinpath(root, "docs", "src",  "index.md")
cp(src, dst; force = true)

makedocs(
    # root = joinpath(dirname(pathof(GAlgebra)), ".."),
    # source = "./",
    # build = "./docs/build",
    modules = [GAlgebra],
    clean = false,
    sitename = "GAlgebra.jl",
    pages = Any[
        "Home" => "index.md",
        "API" => "api.md"
        # "Manual" => Any[
        #     "Plotting" => "man/plotting.md",
        #     "Compositing" => "man/compositing.md",
        #     "Backends" => "man/backends.md",
        #     "Themes" => "man/themes.md",
        # ]
    ]
)

# deploydocs(
#     repo   = "github.com/xxx/yyy.git",
# )
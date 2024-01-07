using PhysicsEducation
using Documenter

DocMeta.setdocmeta!(PhysicsEducation, :DocTestSetup, :(using PhysicsEducation); recursive=true)

makedocs(;
    modules=[PhysicsEducation],
    authors="wwangnju <wwangnju@163.com> and contributors",
    repo="https://github.com/wwangnju/PhysicsEducation.jl/blob/{commit}{path}#{line}",
    sitename="PhysicsEducation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wwangnju.github.io/PhysicsEducation.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/wwangnju/PhysicsEducation.jl",
    devbranch="master",
)

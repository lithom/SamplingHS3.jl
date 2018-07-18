using Documenter, SamplingHS3

makedocs(modules=[SamplingHS3],
        doctest=true)

deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/lithom/SamplingHS3.git",
    julia  = "0.6.2",
    osname = "Windows")

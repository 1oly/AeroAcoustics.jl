using Documenter, AeroAcoustics

makedocs(
    modules = [AeroAcoustics],
    #format = Documenter.HTML(),
    checkdocs = :none,
    sitename = "AeroAcoustics.jl",
    authors = "Oliver Lylloff",
    pages = ["index.md","API Reference" => "reference.md"]
)

deploydocs(
  repo   = "github.com/1oly/AeroAcoustics.jl.git",
  #target = "build",
)


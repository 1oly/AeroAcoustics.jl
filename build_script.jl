Pkg.update()
Pkg.clone(pwd(), \"AeroAcoustics\")
#Pkg.clone("https://gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl","AeroAcoustics")
Pkg.build("AeroAcoustics")

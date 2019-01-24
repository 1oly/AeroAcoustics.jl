using Test
using Dates
using BinaryProvider

const testdatadir = @__DIR__

REPO_URL = "https://gitlab.windenergy.dtu.dk/ollyl/aeroacousticsdata/raw/master/"

remotefiles = [
    ("data/test1_csm.h5", "8f2ba2948f852cf7e406d4c8ad1988beecc5e437f38c1d9369d502891e8e947a"),
    ("data/test1_timeseries.h5", "6b940853919c07d3db2e10a9c5c39495698309e0d60ed77dedac9e4a3406a8cc")
    ]

for (f, sha) in remotefiles
    # create the directories if they don't exist
    currdir = dirname(f)
    isdir(currdir) || mkpath(currdir)
    # download the file if it is not there or if it has a different checksum
    currfile = normpath(joinpath(testdatadir, f))
    url = REPO_URL * f * "?inline=false"
    download_verify(url, sha, currfile; force=true)
end

@testset "AeroAcoustics" begin
    cd(dirname(@__FILE__)) do
        @time include("monopole_noflow.jl")
    end
end

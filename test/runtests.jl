using Test
using Dates
using BinaryProvider

const testdatadir = @__DIR__

REPO_URL = "https://gitlab.windenergy.dtu.dk/ollyl/aeroacousticsdata/raw/master/"

remotefiles = [
    ("data/test1_csm.h5", "8f2ba2948f852cf7e406d4c8ad1988beecc5e437f38c1d9369d502891e8e947a"),
    ("data/test1_timeseries.h5", "1f897db4039138f0c31f9cb340cce3304b85949f9f745200cc161f47426e45ec")
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

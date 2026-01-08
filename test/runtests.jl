using Test
using Dates
using Downloads
using SHA

const testdatadir = @__DIR__

REPO_URL = "https://gitlab.windenergy.dtu.dk/ollyl/aeroacousticsdata/raw/master/"

remotefiles = [
    (
        "data/test1_csm.h5",
        "8f2ba2948f852cf7e406d4c8ad1988beecc5e437f38c1d9369d502891e8e947a",
        REPO_URL * "data/test1_csm.h5?inline=false",
    ),
    (
        "data/task6_DTU_AeroAcousticsjl_CBF_srcint.csv",
        "c75ab472c9767435e801523e157526b936f31d53858119c58833c2a410699ab5",
        "https://raw.githubusercontent.com/MicrophoneArrayBenchmarking/airfoil-in-kevlar-walled-windtunnel/main/data/task6_DTU_AeroAcousticsjl_CBF_srcint.csv",
    ),
    #(
    #    "data/test1_timeseries.h5",
    #    "1f897db4039138f0c31f9cb340cce3304b85949f9f745200cc161f47426e45ec",
    #    REPO_URL * "data/test1_timeseries.h5?inline=false",
    #)
]

# Fix instead of using BinaryProvider 
# https://github.com/yeesian/ArchGDAL.jl/blob/master/test/remotefiles.jl
function verify(path::AbstractString, hash::AbstractString)
    @assert occursin(r"^[0-9a-f]{64}$", hash)
    hash = lowercase(hash)
    if isfile(path)
        calc_hash = open(path) do file
            return bytes2hex(sha256(file))
        end
        @assert occursin(r"^[0-9a-f]{64}$", calc_hash)
        if calc_hash != hash
            @error "Hash Mismatch! Expected: $hash, Calculated: $calc_hash\n"
            return false
        else
            return true
        end
    else
        error("File read error: $path")
    end
end

function download_verify(
    url::AbstractString,
    hash::Union{AbstractString,Nothing},
    dest::AbstractString,
)
    file_existed = false
    # verify if file exists
    if isfile(dest)
        file_existed = true
        if hash !== nothing && verify(dest, hash)
            # hash verified
            return true
        else
            # either hash is nothing or couldn't pass the SHA test
            @error(
                "Failed to verify file: $dest with hash: $hash. Re-downloading file..."
            )
        end
    end
    # if the file exists but some problem exists, we delete it to start from scratch
    file_existed && Base.rm(dest; force = true)
    # Make sure the containing folder exists
    mkpath(dirname(dest))
    # downloads the file at dest
    Downloads.download(url, dest)
    # hash exists and verification fails
    if hash !== nothing && !verify(dest, hash)
        if file_existed
            # the file might be corrupted so we start from scracth
            Base.rm(dest; force = true)
            Downloads.download(url, dest)
            if hash !== nothing && !verify(dest, hash)
                error("Verification failed")
            end
        else
            error("Verification failed. File not created after download.")
        end
    end
    return !file_existed
end

for (f, sha, url) in remotefiles
    # create the directories if they don't exist
    currdir = dirname(f)
    isdir(currdir) || mkpath(currdir)
    # download the file if it is not there or if it has a different checksum
    currfile = normpath(joinpath(testdatadir, f))
    download_verify(url, sha, currfile)
end

@testset "AeroAcoustics" begin
    cd(dirname(@__FILE__)) do
        include("utility_tests.jl")
        include("csm_test.jl")
        include("monopole_noflow.jl")
        include("benchmark_airfoil.jl")
    end
end

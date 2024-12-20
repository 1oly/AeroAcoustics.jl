using Test
using AeroAcoustics
import DSP

# Generate a test signal, a 2 Vrms sine wave at 1234 Hz, corrupted by
# 0.001 V**2/Hz of white noise sampled at 10 kHz.
# From scipy.signal.welch docs:
# https://github.com/scipy/scipy/blob/f2ec91c4908f9d67b5445fbfacce7f47518b35d1/scipy/signal/spectral.py#L223

@testset "CSM test:" begin
    for N in (1e4,1e5,1e6)
        for fs in (10e3,10e4,10e5)
            for n = (2^10,2^11,2^12,2^13)
                amp = 2*sqrt(2)
                fr = 1234.0
                noise_power = 0.001*fs/2
                time = (1:N)/fs
                x = amp*sin.(2*pi*fr.*time)
                x .+= sqrt.(noise_power).*rand(length(time))
                @time Pxx_d = csm(x;n=n,noverlap=div(n,2),win=DSP.hanning(n),fs=fs,scaling="density")
                fc = Pxx_d.fc
                Pxx_d1 = real.(Pxx_d[1,1,:])
                Pxx_w = DSP.welch_pgram(x,n;onesided=true,fs=fs,window=DSP.hanning(n))
                @test Pxx_d1 ≈ Pxx_w.power
                @test Pxx_d ≈ AeroAcoustics.csm_slow(x;n=n,noverlap=div(n,2),win=DSP.hanning(n),fs=fs,scaling="density")
                # spectrum scaling
                Pxx_s = csm(x;n=n,noverlap=div(n,2),win=DSP.hanning(n),fs=fs,scaling="spectrum")
                fc = Pxx_s.fc
                Pxx_s1 = real.(Pxx_s[1,1,:])
                ENBW = enbw(fs,DSP.hanning(n))
                @test Pxx_s1 ≈ Pxx_w.power*ENBW
                @test Pxx_s ≈ AeroAcoustics.csm_slow(x;n=n,noverlap=div(n,2),win=DSP.hanning(n),fs=fs,scaling="spectrum")
            end
        end
    end
end

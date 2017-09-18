# Test shear layer computation
using AeroAcoustics

function f(xi,fvec,xm,M,h)
    a = sqrt(xi[1]^2+(1-M^2)*(xi[2]^2+xi[3]^2))
    b = sqrt((xm[1]-xi[1])^2+(xm[2]-xi[2])^2+(xm[3]-xi[3])^2)
    fvec[1] = (1/a)*xi[1]-(1/b)*(1-M^2)*(xm[1]-xi[1]) - M
    fvec[2] = (1/a)*xi[2]-(1/b)*(1-M^2)*(xm[2]-xi[2])
    fvec[3] = xi[3] - h
end

Ma = 0.3
h = 0.25

res = nlsolve((x,fvec)->f(x,fvec,xm,M,h), [ 1.; 1.; 1.])

xi = res.zero
c0 = 343

r1 = sqrt(xi[1]^2+xi[2]^2+xi[3]^2)
d = xi[1]*M*c0/r1
c1 = d + sqrt(d^2+c0^2-(M*c0)^2)
r2 = sqrt((xm[1]-xi[1])^2+(xm[2]-xi[2])^2+(xm[3]-xi[3])^2)
ta = r1/c1 + r2/c0

dx = 0.025
dy = dx
x0,x1 = 0.0,0.2
y0,y1 = -0.1,0.4
rx = x0:dx:x1
ry = y0:dy:y1
Nx = length(rx)
Ny = length(ry)
z0 = 0.5

X,Y = (Float64[i for i in rx, j in ry],Float64[j for i in rx, j in ry])

M = size(rn,1)    # Number of microphones
const omega = 2pi*f         # Angular frequency
const c = 343.0       # Speed of sound
Ma = 0.2
const h = 0.25
# CSM[eye(Bool,M)] = 0;        # Naive diagonal removal

# Allocation of arrays
gj = Array{C}(M)
gjs = Array{C}(Nx,Ny,M)
b = Array{T}(Nx,Ny)
rn = [0. 0. ; 1 1]
xm = zeros(3)
# Compute transfer functions
i = 1
j = 1
r0 = sqrt(X[i,j]^2 + Y[i,j]^2 + z0^2)
m = 1
xm = [X[i,j]-rn[m,1], Y[i,j]-rn[m,2], z0]
res = nlsolve((x,fvec)->shear(x,fvec,xm,Ma,h), [ 1.; 1.; 1.])
xi = res.zero
r1 = sqrt(xi[1]^2+xi[2]^2+xi[3]^2)
d = xi[1]*M*c/r1
c^2-(Ma*c)^2
c1 = d + sqrt(d^2+c^2-(M*c)^2)
r2 = sqrt((xm[1]-xi[1])^2+(xm[2]-xi[2])^2+(xm[3]-xi[3])^2)
ta = r1/c1 + r2/c
gj[m] = (1/M)*(ta*c/r0)*exp(-im*omega*ta) # TYPE II? Steering vector:

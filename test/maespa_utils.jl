using Photosynthesis: quad, Lower, Upper

# Quadratic solvers
# quadm: quad
quadm_fortran = Libdl.dlsym(photosynlib, :quadm_)
a = 0.5
b = -0.5
c = 0.05
quad_ref = ccall(quadm_fortran, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Int32}), a, b, c, 1)/4.0
quad_test = quad(Lower(), a,b,c)/4
@test quad_ref ≈ quad_test atol=1e-5

# quadm: quap
quadp_fortran = Libdl.dlsym(photosynlib, :quadp_)
a = 0.5
b = -0.5
c = 0.05
quad_ref = ccall(quadp_fortran, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Int32}), a, b, c, 1)/4.0
quad_test = quad(Upper(), a,b,c)/4
@test quad_ref ≈ quad_test atol=1e-5


# Arrenius equations
# arrhfn: arrhenius
arrhfn_fortran = Libdl.dlsym(photosynlib, :arrh_)
arrh_ref = ccall(arrhfn_fortran, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 42.75, 37830.0, 30.0, 25.0)
arrh = arrhenius(42.75u"μmol*mol^-1", 37830.0u"J*mol^-1", 30.0u"°C", 25.0u"°C")
@test arrh.val ≈ arrh_ref

# lattice constant, in angstrom
l = 1.42
E0 = [3/2 3/2; -sqrt(3)/2 sqrt(3)/2]
P0 = [1; 0]

# Graphene elastic constants (eV/unit cell)
K = 2.0*34.759  # Bulk modulus
G = 2.0*23.676  # Shear modulus

# From Koshino 2017
# K = 40.8626
# G = 45.5776

fac = 1 # enhance the GSFE

# Graphene GSFE functional
function GSFE(X::Array{Float64,2})
    s = view(X, 1, :)
    t = view(X, 2, :)

    return 4.063543396e-3 .* (cos.(s) .+ cos.(t) .+ cos.(s.+t)) .-
                0.374e-3 .* (cos.(s.+2. .*t) .+ cos.(s-t) .+ cos.(2. .*s .+ t)) .-
                0.094514083e-3 .* (cos.(2. .*s) .+ cos.(2. .*t) .+ cos.(2. .*(s .+ t)))
end

# In-place variant - uses storage2 and outputs results in storage array
function GSFE!(storage::SubArray{Float64,1}, X::SubArray{Float64,2})
    s = view(X, 1, :)
    t = view(X, 2, :)

    for i=1:length(storage)
      @inbounds storage[i] = fac*(4.063543396e-3 * (cos(s[i]) + cos(t[i]) + cos(s[i]+t[i])) -
                0.374e-3 * (cos(s[i]+2. *t[i]) + cos(s[i]-t[i]) + cos(2. *s[i] + t[i])) -
                0.094514083e-3 * (cos(2. *s[i]) + cos(2. *t[i]) + cos(2. *(s[i] + t[i]))))
    end
    nothing
end


# Graphene GSFE functional gradient
function gradient_GSFE(X::Array{Float64,2})
    s = view(X, 1, :)
    t = view(X, 2, :)

    G = similar(X)
    for i=1:size(X,2)
      @inbounds G[1,i] = fac*(-4.063543396e-3 * (sin(s[i]) + sin(s[i]+t[i])) +
                 0.374e-3 * (sin(s[i]+2. *t[i]) + sin(s[i]-t[i]) + 2. *sin(2. *s[i] + t[i])) +
                 0.189028166e-3 * (sin(2. *s[i]) + sin(2. *(s[i]+t[i]))))
      @inbounds G[2,i] = fac*(-4.063543396e-3 * (sin(t[i]) + sin(s[i]+t[i])) +
                  0.374e-3 * (2. *sin(s[i]+2. *t[i]) + sin(t[i]-s[i]) + sin(2. *s[i]+t[i])) +
                  0.189028166e-3 * (sin(2. *t[i]) + sin(2. *(s[i]+t[i]))))
    end

    return G
end

# In-place variant - outputs results in storage array
function gradient_GSFE!(storage::SubArray{Float64,2}, X::SubArray{Float64,2})
    s = view(X, 1, :)
    t = view(X, 2, :)

    for i=1:size(storage,2)
      @inbounds storage[1,i] = fac*(-4.063543396e-3 * (sin(s[i]) + sin(s[i]+t[i])) +
                 0.374e-3 * (sin(s[i]+2. *t[i]) + sin(s[i]-t[i]) + 2. *sin(2. *s[i] + t[i])) +
                 0.189028166e-3 * (sin(2. *s[i]) + sin(2. *(s[i]+t[i]))))
      @inbounds storage[2,i] = fac*(-4.063543396e-3 * (sin(t[i]) + sin(s[i]+t[i])) +
                  0.374e-3 * (2. *sin(s[i]+2. *t[i]) + sin(t[i]-s[i]) + sin(2. *s[i]+t[i])) +
                  0.189028166e-3 * (sin(2. *t[i]) + sin(2. *(s[i]+t[i]))))
    end

    nothing
end

nothing

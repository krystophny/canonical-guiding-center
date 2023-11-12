const B0 = 1.0

function Bmod(x)
    return B0
end

function dBmod(x)
    return [0.0, 0.0, 0.0]
end

function A(x)
    return B0 .* [0.0, x[1], 0.0]
end

function dA(x)
    dA = zeros(3, 3)
    dA[1, 2] = B0
    return dA
end

function B(x)
    return [0.0, 0.0, B0]
end

function dB(x)
    return zeros(3, 3)
end

function e1(x)
    return [0.0, 1.0, 0.0]
end

function e2(x)
    return [1.0, 0.0, 0.0]
end

function de1(x)
    return zeros(3, 3)
end

function de2(x)
    return zeros(3, 3)
end

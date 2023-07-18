local matrice = {}

local function matmul(vec, b)
    local a1, a2, a3, a4 = vec[1],vec[2],vec[3],vec[4]
    return {
          a1 * b[1]  + a2 * b[2]  + a3 * b[3]  + a4 * b[4],
          a1 * b[5]  + a2 * b[6]  + a3 * b[7]  + a4 * b[8],
          a1 * b[9]  + a2 * b[10] + a3 * b[11] + a4 * b[12],
          a1 * b[13] + a2 * b[14] + a3 * b[15] + a4 * b[16]
    }
end

local function attacher(self,matrice)
    return setmetatable(matmul(self,matrice),{
        __mul=attacher,
    })
end

function matrice.new(N, M, ...)
    local m = { ... }
    m.matrix_width = N
    m.matrix_height = M
    return setmetatable(m,{__mul = attacher})
end

function matrice.vector(...)
    local m = { ... }
    return matrice.new(#m, 1, ...)
end

return matrice
## B-splines, recursive definition
B0(u::Vector{Float64}, i::Integer) = x::Real -> (u[i] <= x < u[i+1]) ? float(1) : float(0) ## non-continuous
B(p::Integer,u::Vector{Float64},i::Integer) = (p==0) ? B0(u,i) : x -> float(((x-u[i])/(u[i+p]-u[i]))*B(p-1,u,i)(x) + ((u[i+p+1]-x)/(u[i+p+1]-u[i+1]))*B(p-1,u,i+1)(x))

totalNumNodes(numNodes::Integer, p::Integer) = numNodes + 2p ## adding p nodes on each side
numOfBasisElements(numNodes::Integer,p::Integer) = numNodes + p - 1 ## number of non-zero spline functions

## returns an index-function of basis-functions, which takes as input i in 1:k
## creates evenly-spaced nodes
function splineBasisFuns(p::Integer, range::(Real,Real), numNodes::Integer)
    intervalLen = (range[2]-range[1])/(numNodes-1)
    u = linspace(range[1] - p*intervalLen, range[2]+ p*intervalLen, totalNumNodes(numNodes,p))
    i -> B(p,u,i)
end

## B-splines with pre-specific nodes
## WARNING: all functions in this basis will be 0 at the leftmost and rightmost nodes, and the same is true of their first p-1 derivatives.
function splineBasisFuns(p::Integer, nodes::Vector{Float64})
    i::Integer -> B(p,nodes,i)
end


cubicSplineBasisFuns(range::(Real,Real), numNodes::Integer) = splineBasisFuns(3, range, numNodes)
quadraticSplineBasisFuns(range::(Real,Real), numNodes::Integer) = splineBasisFuns(2, range, numNodes)
linearSplineBasisFuns(range::(Real,Real), numNodes::Integer) = splineBasisFuns(1, range, numNodes)


## BFE "Basis Function Expansion"

## evenly-spaced nodes
function splineBFE(p::Integer, x::Vector{Float64}, range::(Real,Real), numNodes::Integer)
    funs = splineBasisFuns(p, range, numNodes)
    num = numOfBasisElements(numNodes,p)
    Float64[funs(i)(xj) for xj in x, i in 1:num]
end

## pre-specified nodes
function splineBFE(p::Integer, x::Vector{Float64}, nodes::Vector{Float64})
    funs = splineBasisFuns(p, nodes)
    n = length(nodes) - p - 1
    Float64[funs(i)(xj) for xj in x, i in 1:n]
end


cubicSplineBFE(x::Vector{Float64}, range::(Real,Real), numNodes::Integer) =
    splineBFE(3, x, range, numNodes)
quadraticSplineBFE(x::Vector{Float64}, range::(Real,Real), numNodes::Integer) =
    splineBFE(2, x, range, numNodes)
linearSplineBFE(x::Vector{Float64}, range::(Real,Real), numNodes::Integer) =
    splineBFE(1, x, range, numNodes)

cubicSplineBFE(x::Vector{Float64}, nodes::Vector{Float64}) = splineBFE(3, x, nodes)
quadraticSplineBFE(x::Vector{Float64}, nodes::Vector{Float64}) = splineBFE(2, x, nodes)
linearSplineBFE(x::Vector{Float64}, nodes::Vector{Float64}) = splineBFE(1, x, nodes)

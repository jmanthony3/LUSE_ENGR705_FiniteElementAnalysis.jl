using LinearAlgebra
using LinearSolve
using Symbolics

struct LineElements
    nodes   ::Tuple{Integer, Integer}
    A       ::Float64
    E       ::Float64
    L       ::Float64
    k       ::Float64
    θ       ::Float64
    K       ::Matrix{Float64}
end

struct BeamElements
    nodes   ::Tuple{Integer, Integer}
    E       ::Float64
    I       ::Float64
    L       ::Float64
    k       ::Float64
    θ       ::Float64
    K       ::Matrix{Float64}
end

struct FiniteElementAnalysisSolution
    F::Vector{Float64}
    K::Matrix{Float64}
    U::Vector{Float64}
end

function transformationmatrix(θ::Real, isdegrees::Bool=true)::Matrix{Float64}
    C, S = isdegrees ? (cosd(θ), sind(θ)) : (cos(θ), sin(θ))
    return float.([
        [C S 0 0];
        [-S C 0 0];
        [0 0 C S];
        [0 0 -S C]
    ])
end

function transformationmatrixstar(θ::Real, isdegrees::Bool=true)::Matrix{Float64}
    C, S = isdegrees ? (cosd(θ), sind(θ)) : (cos(θ), sin(θ))
    return float.([[C S 0 0]; [0 0 C S]])
end

"""
    prepareelements_line(elements, elementareas, elementmoduli, elementlengths, nodeangles=0, isdegrees=true; dims=1, stiffnesscoefficients=nothing)::Dict{String, LineElements}

Return the dictionary with properties specific to those **line** elements between node pairs.

See also: `globalstiffnessmatrix`, `solve`, and `axialstresses`.

# Arguments
- `elements::Dict{String, Tuple{Int, Int}}`: the dictionary of "element"=>(node 1, node 2) for the element between node pairs.
- `elementareas::Union{<:Real, Vector{<:Real}}`: the cross-sectional area(s) of elements.
- `elementmoduli::Union{<:Real, Vector{<:Real}}`: the elastic modul(us/i) of elements.
- `elementlengths::Union{<:Real, Vector{<:Real}}`: the length(s) of elements.
- `nodeangles::Union{<:Real, Vector{<:Real}}=0`: the angle between the global and local coordinate systems of elements.
- `isdegrees::Bool=true`: control value whether elements of `nodeangles` is degrees (`true`) or radians (`false`).
- `dims::Integer=1`: the number of dimensions being analyzed.
- `stiffnesscoefficients::Union{<:Real, Vector{<:Real}}`: the stiffness coefficient, ``k`` associated with each element.

# Examples
```jldoctest; output=false
A = [0.75, 0.5, 1]   # [in²]
E = 10.6e6           # [psi]
L = [12, 9, 8]       # [in]
elements = Dict("1"=>(1, 2), "2"=>(2, 3), "3"=>(3, 4))
prepareelements_line(elements, A, E, L)

# output

Dict{String, FiniteElementAnalysis.LineElements} with 3 entries:
  "1" => LineElements((1, 2), 0.75, 1.06e7, 12.0, 662500.0, 0.0, [662500.0 -662500.0; -662500.0 662500.0])
  "2" => LineElements((2, 3), 0.5, 1.06e7, 9.0, 5.88889e5, 0.0, [5.88889e5 -5.88889e5; -5.88889e5 5.88889e5])
  "3" => LineElements((3, 4), 1.0, 1.06e7, 8.0, 1.325e6, 0.0, [1.325e6 -1.325e6; -1.325e6 1.325e6])
```
"""
function prepareelements_line(
    elements                ::Dict{String, Tuple{Int, Int}},
    elementareas            ::Union{<:Real, Vector{<:Real}},
    elementmoduli           ::Union{<:Real, Vector{<:Real}},
    elementlengths          ::Union{<:Real, Vector{<:Real}},
    nodeangles              ::Union{<:Real, Vector{<:Real}}=0,
    isdegrees               ::Bool=true;
    dims                    ::Integer=1,
    stiffnesscoefficients   ::Union{Nothing, <:Real, Vector{<:Real}}=nothing,
)::Dict{String, LineElements}
    preparedelementsdictionary = Dict{String, LineElements}()
    for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
        i = parse(Int64, key)
        A = float(length(elementareas) == 1 ? elementareas[1] : elementareas[i])
        E = float(length(elementmoduli) == 1 ? elementmoduli[1] : elementmoduli[i])
        L = float(length(elementlengths) == 1 ? elementlengths[1] : elementlengths[i])
        θ = if isdegrees
            float(length(nodeangles) == 1 ? nodeangles[1] : nodeangles[i])
        else
            rad2deg(length(nodeangles) == 1 ? nodeangles[1] : nodeangles[i])
        end
        k = isnothing(stiffnesscoefficients) ? A*E/L : float(length(stiffnesscoefficients) == 1 ? stiffnesscoefficients[1] : stiffnesscoefficients[i])
        if dims == 2
            T = transformationmatrix(θ, isdegrees)
            K = T'*[[1 0 -1 0]; [0 0 0 0]; [-1 0 1 0]; [0 0 0 0]]*T
        else
            K = [[1 -1]; [-1 1]]
        end
        K *= k
        preparedelementsdictionary[key] = LineElements((value[1], value[2]), A, E, L, k, θ, K)
    end
    return preparedelementsdictionary
end

"""
    prepareelements_beam(elements, elementmoduli, elementareasofinertia, elementlengths, nodeangles=0, isdegrees=true; dims=2, stiffnesscoefficients=nothing)::Dict{String, BeamElements}

Return the dictionary with properties specific to those **beam** elements between node pairs.

See also: `globalstiffnessmatrix`, `solve`, and `axialstresses`.

# Arguments
- `elements::Dict{String, Tuple{Int, Int}}`: the dictionary of "element"=>(node 1, node 2) for the element between node pairs.
- `elementmoduli::Union{<:Real, Vector{<:Real}}`: the elastic modulus of each element.
- `elementareasofinertia::Union{<:Real, Vector{<:Real}}`: the cross-sectional area of each element.
- `elementlengths::Union{<:Real, Vector{<:Real}}`: the length of each element.
- `nodeangles::Union{<:Real, Vector{<:Real}}=0`: the angle between the global and local coordinate systems of each element.
- `isdegrees::Bool=true`: control value whether elements of `nodeangles` is degrees (true) or radians (false).
- `dims::Integer=2`: the number of dimensions being analyzed.
- `stiffnesscoefficients::Union{<:Real, Vector{<:Real}}`: the stiffness coefficient, ``k`` associated with each element.

# Examples
```jldoctest; output=false
E, I, L = 30e6, 200, 20     # [psi, in⁴, ft]
elements = Dict("1"=>(1, 2),"2"=>(2, 3))
prepareelements_beam(elements, E, I, L)

# output

Dict{String, FiniteElementAnalysis.BeamElements} with 2 entries:
  "1" => BeamElements((1, 2), 3.0e7, 200.0, 20.0, 750000.0, 0.0, [9.0e6 9.0e7 -9.0e6 9.0e7; 9.0e7 1.2e9 -9.0e7 6.0e8; -9.0e6 -9.0e7 9.0e6 -9.0e7; 9.0e7 6.0e8 -9.0e7 1.2e9])
  "2" => BeamElements((2, 3), 3.0e7, 200.0, 20.0, 750000.0, 0.0, [9.0e6 9.0e7 -9.0e6 9.0e7; 9.0e7 1.2e9 -9.0e7 6.0e8; -9.0e6 -9.0e7 9.0e6 -9.0e7; 9.0e7 6.0e8 -9.0e7 1.2e9])
```
"""
function prepareelements_beam(
    elements                ::Dict{String, Tuple{Int, Int}},
    elementmoduli           ::Union{<:Real, Vector{<:Real}},
    elementareasofinertia   ::Union{<:Real, Vector{<:Real}},
    elementlengths          ::Union{<:Real, Vector{<:Real}},
    nodeangles              ::Union{<:Real, Vector{<:Real}}=0,
    isdegrees               ::Bool=true;
    dims                    ::Integer=2,
    stiffnesscoefficients   ::Union{Nothing, <:Real, Vector{<:Real}}=nothing,
)::Dict{String, BeamElements}
    preparedelementsdictionary = Dict{String, BeamElements}()
    for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
        i = parse(Int64, key)
        E = float(length(elementmoduli) == 1 ? elementmoduli[1] : elementmoduli[i])
        I = float(length(elementareasofinertia) == 1 ? elementareasofinertia[1] : elementareasofinertia[i])
        L = float(length(elementlengths) == 1 ? elementlengths[1] : elementlengths[i])
        θ = if isdegrees
            float(length(nodeangles) == 1 ? nodeangles[1] : nodeangles[i])
        else
            rad2deg(length(nodeangles) == 1 ? nodeangles[1] : nodeangles[i])
        end
        k = isnothing(stiffnesscoefficients) ? E*I/L^3 : float(length(stiffnesscoefficients) == 1 ? stiffnesscoefficients[1] : stiffnesscoefficients[i])
        if dims == 2
            # T = transformationmatrix(θ, isdegrees)
            # K = T'*[[12 6L -12 6L]; [6L 4L^2 -6L 2L^2]; [-12 -6L 12 -6L]; [6L 2L^2 -6L 4L^2]]*T
            K = [[12 6L -12 6L]; [6L 4L^2 -6L 2L^2]; [-12 -6L 12 -6L]; [6L 2L^2 -6L 4L^2]]
        else
            K = [[12 6L -12 6L]; [6L 4L^2 -6L 2L^2]; [-12 -6L 12 -6L]; [6L 2L^2 -6L 4L^2]]
        end
        K *= k
        preparedelementsdictionary[key] = BeamElements((value[1], value[2]), E, I, L, k, θ, K)
    end
    return preparedelementsdictionary
end

"""
    globalstiffnessmatrix(elements::Dict{String, LineElements}, dims::Integer=1)::Matrix{Float64}
    globalstiffnessmatrix(elements::Dict{String, BeamElements}, dims::Integer=2)::Matrix{Float64}

Assembles the global stiffness matrix from the dictionary, `elements`.

See also: `prepareelements_line`, `prepareelements_beam`, and `transformationmatrix`.

# Arguments
- `elements`: entries infer `localstiffnessmatrices` and `nodeconnections` for each element in local coordinate system transformed by angle.
- `dims`: the number of dimensions being analyzed.

# Examples
```jldoctest; output=false
A = [0.75, 0.5, 1]   # [in²]
E = 10.6e6           # [psi]
L = [12, 9, 8]       # [in]
elements = Dict("1"=>(1, 2), "2"=>(2, 3), "3"=>(3, 4))
globalstiffnessmatrix(prepareelements_line(elements, A, E, L))

# output

4×4 Matrix{Float64}:
  662500.0  -662500.0         0.0         0.0
 -662500.0        1.25139e6  -5.88889e5   0.0
       0.0       -5.88889e5   1.91389e6  -1.325e6
       0.0        0.0        -1.325e6     1.325e6
```

```jldoctest; output=false
E, I, L = 30e6, 510, 60 # [psi, in⁴, in]
W(x) = -12\\1e3(x)      # [lb/in]
F_applied = ([(W, 1, 2, 3)])
elements = Dict("1"=>(1, 2), "2"=>(2, 3))
globalstiffnessmatrix(prepareelements_beam(elements, E, I, L))

# output

6×6 Matrix{Float64}:
  850000.0      2.55e7  -850000.0      2.55e7        0.0      0.0
       2.55e7   1.02e9       -2.55e7   5.1e8         0.0      0.0
 -850000.0     -2.55e7        1.7e6    0.0     -850000.0      2.55e7
       2.55e7   5.1e8         0.0      2.04e9       -2.55e7   5.1e8
       0.0      0.0     -850000.0     -2.55e7   850000.0     -2.55e7
       0.0      0.0           2.55e7   5.1e8        -2.55e7   1.02e9
```
"""
function globalstiffnessmatrix(
    elements    ::Dict{String, LineElements};
    dims        ::Integer=1,
)::Matrix{Float64}
    nodeconnections = Vector{Tuple{Integer, Integer}}([])
    for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
        push!(nodeconnections, value.nodes)
    end
    n = maximum(maximum(nodeconnections))
    K = zeros((n, n) .* dims)
    for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
        if dims == 2
            p = 0
            for node ∈ value.nodes
                for i ∈ 1:1:2
                    q = 0
                    for n ∈ value.nodes
                        for j ∈ 1:1:2
                            K[2(node - 1) + i, 2(n - 1) + j] += value.K[i + p, j + q]
                        end
                        q += 2
                    end
                end
                p += 2
            end
        else
            matrixdims = size(value.K)
            for i ∈ 1:1:matrixdims[1], j ∈ 1:1:matrixdims[2]
                K[value.nodes[i], value.nodes[j]] += value.K[i, j]
            end
        end
    end
    return K
end

function globalstiffnessmatrix(
    elements    ::Dict{String, BeamElements};
    dims        ::Integer=2,
)::Matrix{Float64}
    nodeconnections = Vector{Tuple{Integer, Integer}}([])
    for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
        push!(nodeconnections, value.nodes)
    end
    n = maximum(maximum(nodeconnections))
    # K = zeros((n, n) .* dims)
    K = (dims == 2 ? zeros((2n, 2n)) : zeros((4n, 4n)))
    for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
        k = parse(Int64, key)
        if dims == 2
            p = 0
            for node ∈ value.nodes
                for i ∈ 1:1:2
                    q = 0
                    for n ∈ value.nodes
                        for j ∈ 1:1:2
                            # println("node=$node, n=$n: [2(node - 1) + i, 2(n - 1) + j] = $([2(node - 1) + i, 2(n - 1) + j]) <- [i + p, j + q] = [$(i + p), $(j + q)]")
                            K[2(node - 1) + i, 2(n - 1) + j] += value.K[i + p, j + q]
                        end
                        q += 2
                    end
                end
                p += 2
            end
        else
            p = 0
            for node ∈ value.nodes
                for i ∈ 1:1:2
                    q = 0
                    for n ∈ value.nodes
                        for j ∈ 1:1:2
                            # println("node=$node, n=$n: [2(node - 1) + i, 2(n - 1) + j] = $([2(node - 1) + i, 2(n - 1) + j]) <- [i + p, j + q] = [$(i + p), $(j + q)]")
                            K[2(node - 1) + i, 2(n - 1) + j] += value.K[i + p, j + q]
                        end
                        q += 2
                    end
                end
                p += 2
            end
        end
    end
    return K
end

"""
    reducedglobalstiffnessmatrix(K::Matrix{Float64}, nodes::Union{Vector{Integer}, Vector{Tuple{Integer, Integer}}}, dims::Integer=1)::Matrix{Float64}

Eliminate the rows and columns of `K` according to `nodes` which implies the displacement of the appropriate node is known: e.g. is fixed or a prescribed displacement.

See also: `globalstiffnessmatrix`.

# Examples
```jldoctest; output=false
A = [0.75, 0.5, 1]   # [in²]
E = 10.6e6           # [psi]
L = [12, 9, 8]       # [in]
elements = Dict("1"=>(1, 2), "2"=>(2, 3), "3"=>(3, 4))
K = globalstiffnessmatrix(prepareelements_line(elements, A, E, L))
skipnodes = Vector{Integer}([1, 4])
reducedglobalstiffnessmatrix(K, skipnodes)

# output

2×2 Matrix{Float64}:
  1.25139e6  -5.88889e5
 -5.88889e5   1.91389e6
```
"""
function reducedglobalstiffnessmatrix(
    K       ::Matrix{Float64},
    nodes   ::Union{Vector{Integer}, Vector{Tuple{Integer, Integer}}};
    dims    ::Integer=1
)::Matrix{Float64}
    k = 0
    for node ∈ nodes
        if dims == 2
            if node[2] == 0
                K = K[begin:end .!= 2node[1] - 1 - k, begin:end .!= 2node[1] - 1 - k]
                k += 1
                K = K[begin:end .!= 2node[1] - k, begin:end .!= 2node[1] - k]
            else
                K = K[begin:end .!= 2node[1] - 1 + node[2] - 1 - k, begin:end .!= 2node[1] - 1 + node[2] - 1 - k]
            end
        else
            K = K[begin:end .!= node - k, begin:end .!= node - k]
        end
        k += 1
    end
    return K
end

"""
    solve(elements::Dict{String, LineElements}, nodeboundaryconditions, nodeforces::Vector{<:Tuple{Integer, Vararg{<:Real}}}, dims::Integer=1)
    solve(elements::Dict{String, BeamElements}, nodeboundaryconditions, nodeforces::Vector{<:Tuple{Integer, Vararg{<:Real}}}, dims::Integer=2)

Perform Finite Element Analysis (FEA) with entries of `elements` and applied boundary conditions by the Direct Stiffness Method according to Hooke's Law of linear-elastic deformation.

If **line** elements, stiffness coefficient defined by: ``k ≡ AE/L`` (typical of trusses with two-force members).
If **beam** elements, stiffness coefficient defined by: ``k ≡ EI/L³``.

See also: `prepareelements_line`, `prepareelements_beam`, `globalstiffnessmatrix`, and `reducedglobalstiffnessmatrix`.

# Arguments
- `elements`: the prepared dictionary according to element type.
- `nodeboundaryconditions`: the tuple of nodes that either experience zero displacement (tuple item of integer for node) or prescribed displacement constrained by the stiffness of a certain element (tuple item of node integer and tuple of real displacement and integer element).
- `nodeforces::Vector{<:Tuple{Integer, Vararg{<:Real}}}`: the vector of tuple pairs for nodes experiencing some applied force.
- `dims::Integer=1`: the number of dimensions being analyzed.
- `alg=nothing`: algorithm of choice to solve linear system of equations according to [LinearSolve.jl](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/).

# Examples
## 1D
### `nodeboundaryconditions`:= tuple of nodes with zero displacement.
```jldoctest; output=false
A, E, L = 4, 30e6, 30               # [in², psi, in]
elements = Dict("1"=>(1, 2), "2"=>(2, 3), "3"=>(3, 4))
preparedelements = prepareelements_line(elements, A, E, L)
nodeboundaryconditions = (1, 4)
F_applied = [(2, 5e3), (3, -10e3)]  # (node, force [lb])
solve(preparedelements, nodeboundaryconditions, F_applied).U

# output

4-element Vector{Float64}:
  0.0
  0.0
 -0.00125
  0.0
```

### `nodeboundaryconditions`:= tuple of nodes with zero and prescribed displacements.
The tuple items are identified in the following way:
1. Node (1) with zero displacement.
2. Node (3) has a prescribed displacement of `25e-3` constrained by the deformation of element `2`.
```jldoctest; output=false
A, E, L = 4e-4, 210e9, 2                    # [m², Pa, m]
elements = Dict("1"=>(1, 2), "2"=>(2, 3))
preparedelements = prepareelements_line(elements, A, E, L)
nodeboundaryconditions = [1, (3, 25e-3, 2)] # (node, displacement [m], element)
F_applied = [(2, -5e3)]                     # (node, force [N])
solve(preparedelements, nodeboundaryconditions, F_applied).U

# output

3-element Vector{Float64}:
 0.0
 0.012440476190476191
 0.025
```

## 2D
```jldoctest; output=false
A, E, L = 1, 10e6, 100          # [in², psi, in]
angles = [120, 0, 210]          # [°]
L = L ./ abs.(cosd.(angles))    # [in]
elements = Dict("1"=>(1, 2), "2"=>(1, 3), "3"=>(1, 4))
preparedelements = prepareelements_line(elements, A, E, L, angles, dims=2)
nodeboundaryconditions = [2, 3, 4]
F_applied = [(1, (1e3, 1e3))]   # (node, (force_x [lb], force_y [lb]))
solve(preparedelements, nodeboundaryconditions, F_applied, dims=2).U

# output

8-element Vector{Float64}:
 0.004226497308103743
 0.01577350269189626
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
```
"""
function solve(
    elements    ::Dict{String, LineElements},
    nodeboundaryconditions,
    nodeforces  ::Union{Vector{<:Tuple{Integer, <:Real}}, Vector{<:Tuple{Integer, Tuple{<:Real, <:Real}}}};
    dims        ::Integer=1,
    alg=nothing
)::FiniteElementAnalysisSolution
    nodeconnections, stiffnesscoefficients = Vector{Tuple{Integer, Integer}}([]), Vector{Float64}([])
    for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
        push!(nodeconnections, value.nodes)
        push!(stiffnesscoefficients, value.k)
    end
    K = globalstiffnessmatrix(elements, dims=dims)
    n = maximum(maximum(nodeconnections))
    F = zeros(n * dims)
    for (node, force) ∈ nodeforces
        if dims == 1
            F[node] = force
        elseif dims == 2
            F[2node - 1]    = force[1]
            F[2node]        = force[2]
        end
    end
    if dims == 1
        fixednodes = Vector{Integer}([])
    elseif dims == 2
        fixednodes = Vector{Tuple{Integer, Integer}}([])
    end
    for nbc ∈ nodeboundaryconditions
        if dims == 1
            if length(nbc) == 1 || nbc[2] == 0
                push!(fixednodes, nbc[1])
            end
        elseif dims == 2
            if length(nbc) == 1 || nbc[2] == (0, 0)
                push!(fixednodes, (nbc[1], 0))
            end
        end
    end
    skipnodes = deepcopy(fixednodes)
    if length(fixednodes) != length(nodeboundaryconditions)
        for nbc ∈ nodeboundaryconditions
            if nbc[1] ∉ fixednodes
                if dims == 2
                    if length(nbc) != 1
                        if nbc[2][1] == Inf
                            push!(skipnodes, (nbc[1], 2))
                        elseif nbc[2][2] == Inf
                            push!(skipnodes, (nbc[1], 1))
                        end
                    end
                end
                inspectnodes = findall(x->x!=nbc[1], nodeconnections)
                for node ∈ inspectnodes
                    if node ∉ fixednodes
                        if dims == 1
                            push!(skipnodes, nbc[1])
                            if length(stiffnesscoefficients) == 1
                                F[node] += stiffnesscoefficients[1]*nbc[2]
                            else
                                F[node] += stiffnesscoefficients[nbc[3]]*nbc[2]
                            end
                        elseif dims == 2
                            if length(nbc) != 1
                                if nbc[2][1] != 0 && nbc[2][1] != Inf
                                    if length(stiffnesscoefficients) == 1
                                        F[node] += [(stiffnesscoefficients[1]*nbc[2][1])...]
                                    else
                                        F[node] += [(stiffnesscoefficients[nbc[3]]*nbc[2][1])...]
                                    end
                                elseif nbc[2][2] != 0 && nbc[2][2] != Inf
                                    if length(stiffnesscoefficients) == 1
                                        F[node] += [(stiffnesscoefficients[1]*nbc[2][2])...]
                                    else
                                        F[node] += [(stiffnesscoefficients[nbc[3]]*nbc[2][2])...]
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    skipnodes = sort(collect(skipnodes), by=x->x[1])
    K_reduced = (skipnodes != [] ? reducedglobalstiffnessmatrix(K, skipnodes, dims=dims) : K)
    k = 0
    for node ∈ skipnodes
        if dims == 1
            F = F[begin:end .!= node - k]
        elseif dims == 2
            if node[2] == 0
                F = F[begin:end .!= 2node[1] - 1 - k]
                k += 1
                F = F[begin:end .!= 2node[1] - k]
            else
                F = F[begin:end .!= 2node[1] - 1 + node[2] - 1 - k]
            end
        end
        k += 1
    end
    U_reduced = try
        LinearSolve.solve(LinearProblem(K_reduced, F), alg).u
    catch error
        if isa(error, SingularException)
            throw(ErrorException("Error occurred in 'LinearSolve.jl': $error\nThis typically happens with a badly conditioned matrix (cond(K) = $(cond(K_reduced)))."))
        end
    end
    U = U_reduced
    for nbc ∈ nodeboundaryconditions
        if length(nbc) == 1
            if dims == 1
                insert!(U, nbc, 0)
            elseif dims == 2
                insert!(U, 2nbc-1, 0)
                insert!(U, 2nbc, 0)
            end
        else
            if dims == 1
                insert!(U, nbc[1], nbc[2])
            elseif dims == 2
                if length(nbc) == 1
                    insert!(U, 2nbc[1]-1, 0)
                    insert!(U, 2nbc[1], 0)
                else
                    if nbc[2][1] != Inf
                        insert!(U, 2nbc[1]-1, nbc[2][1])
                    end
                    if nbc[2][2] != Inf
                        insert!(U, 2nbc[1], nbc[2][2])
                    end
                end
            end
        end
    end
    F = K*U
    return FiniteElementAnalysisSolution(F, K, U)
end

function solve(
    elements    ::Dict{String, BeamElements},
    nodeboundaryconditions,
    nodeforces  ::Vector{<:Union{Tuple{Function, Vararg{Integer}}, Tuple{Integer, Tuple{<:Real, <:Real}}, Tuple{Integer, Tuple{<:Real, <:Real, <:Real}}}};
    dims        ::Integer=2,
    alg=nothing
)::FiniteElementAnalysisSolution
    nodeconnections, stiffnesscoefficients = Vector{Tuple{Integer, Integer}}([]), Vector{Float64}([])
    for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
        push!(nodeconnections, value.nodes)
        push!(stiffnesscoefficients, value.k)
    end
    K = globalstiffnessmatrix(elements, dims=dims)
    n = maximum(maximum(nodeconnections))
    # K = zeros((n, n) .* dims)
    F = zeros(dims == 2 ? 2n : 4n)
    for nodeforce ∈ nodeforces
        if isa(nodeforce[1], Function)
            @variables x
            d_nodeforces = Differential(x)(nodeforce[1](x))
            d_nodeforces_expanded = expand_derivatives(d_nodeforces)
            remainingvariables = Symbolics.get_variables(d_nodeforces_expanded)
            for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
                if in(value.nodes[1], nodeforce[2:end]) && in(value.nodes[2], nodeforce[2:end])
                    (i, j), L = value.nodes, value.L
                    if length(remainingvariables) == 0
                        F[2i - 1]    += nodeforce[1](L)/2
                        F[2i]        += nodeforce[1](L)*L/12
                        F[2j - 1]    += nodeforce[1](L)/2
                        F[2j]        += -nodeforce[1](L)*L/12
                    elseif length(remainingvariables) == 1
                        F[2i - 1]    += 3nodeforce[1](L)/20
                        F[2i]        += nodeforce[1](L)*L/30
                        F[2j - 1]    += 7nodeforce[1](L)/20
                        F[2j]        += -nodeforce[1](L)*L/30
                    end
                end
            end
        else
            if dims == 1
                i = 0
                for force ∈ nodeforce[2]
                    F[2nodeforce[1] - 1 + i] += force
                    i += 1
                end
            elseif dims == 2
                i = 0
                for force ∈ nodeforce[2]#[2:3]
                    F[2nodeforce[1] - 1 + i] += force
                    i += 1
                end
            end
        end
    end
    if dims == 1
        fixednodes = Vector{Tuple{Integer, Integer}}([])
    elseif dims == 2
        fixednodes = Vector{Tuple{Integer, Integer}}([])
    end
    for nbc ∈ nodeboundaryconditions
        if dims == 1
            if length(nbc) == 1 || nbc[2] == (0, 0)
                push!(fixednodes, (nbc[1], 0))
            end
        elseif dims == 2
            if length(nbc) == 1 || nbc[2] == (0, 0)
                push!(fixednodes, (nbc[1], 0))
            end
        end
    end
    skipnodes = deepcopy(fixednodes)
    if length(fixednodes) != length(nodeboundaryconditions)
        for nbc ∈ nodeboundaryconditions
            if nbc[1] ∉ fixednodes
                if dims == 2
                    if length(nbc) != 1
                        if nbc[2][1] == Inf ||  nbc[2][2] == Inf
                            push!(skipnodes, (nbc[1], 1))
                        end
                    end
                end
                inspectnodes = findall(x->x!=nbc[1], nodeconnections)
                for node ∈ inspectnodes
                    if node ∉ fixednodes
                        if dims == 1
                            if length(nbc) != 1
                                if nbc[2][1] != 0 && nbc[2][1] != Inf
                                    if length(stiffnesscoefficients) == 1
                                        F[node] += [(stiffnesscoefficients[1]*nbc[2][1])...]
                                    else
                                        F[node] += [(stiffnesscoefficients[nbc[3]]*nbc[2][1])...]
                                    end
                                elseif nbc[2][2] != 0 && nbc[2][2] != Inf
                                    if length(stiffnesscoefficients) == 1
                                        F[node] += [(stiffnesscoefficients[1]*nbc[2][2])...]
                                    else
                                        F[node] += [(stiffnesscoefficients[nbc[3]]*nbc[2][2])...]
                                    end
                                end
                            end
                        elseif dims == 2
                            if length(nbc) != 1
                                if nbc[2][1] != 0 && nbc[2][1] != Inf
                                    if length(stiffnesscoefficients) == 1
                                        F[node] += [(stiffnesscoefficients[1]*nbc[2][1])...]
                                    else
                                        F[node] += [(stiffnesscoefficients[nbc[3]]*nbc[2][1])...]
                                    end
                                elseif nbc[2][2] != 0 && nbc[2][2] != Inf
                                    if length(stiffnesscoefficients) == 1
                                        F[node] += [(stiffnesscoefficients[1]*nbc[2][2])...]
                                    else
                                        F[node] += [(stiffnesscoefficients[nbc[3]]*nbc[2][2])...]
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    skipnodes = sort(collect(skipnodes), by=x->x[1])
    K_reduced = (skipnodes != [] ? reducedglobalstiffnessmatrix(K, skipnodes, dims=dims) : K)
    k = 0
    for node ∈ skipnodes
        if dims == 1
            if node[2] == 0
                F = F[begin:end .!= 2node[1] - 1 - k]
                k += 1
                F = F[begin:end .!= 2node[1] - k]
            else
                F = F[begin:end .!= 2node[1] - 1 + node[2] - 1 - k]
            end
        elseif dims == 2
            if node[2] == 0
                F = F[begin:end .!= 2node[1] - 1 - k]
                k += 1
                F = F[begin:end .!= 2node[1] - k]
            else
                F = F[begin:end .!= 2node[1] - 1 + node[2] - 1 - k]
            end
        end
        k += 1
    end
    U_reduced = try
        LinearSolve.solve(LinearProblem(K_reduced, F), alg).u
    catch error
        if isa(error, SingularException)
            throw(ErrorException("Error occurred in 'LinearSolve.jl': $error\nThis typically happens with a badly conditioned matrix (cond(K) = $(cond(K_reduced)))."))
        end
    end
    U = U_reduced
    for skipnode ∈ skipnodes
        if skipnode[2] == 0
            insert!(U, 2skipnode[1] - 1, 0)
            insert!(U, 2skipnode[1], 0)
        else
            insert!(U, 2skipnode[1] - 1 + skipnode[2] - 1, 0)
        end
        # if length(nbc) == 1
        # else
        #     if dims == 1
        #         if length(nbc) == 1
        #             insert!(U, 2nbc[1]-1, 0)
        #             insert!(U, 2nbc[1], 0)
        #         else
        #             if nbc[2][1] != Inf
        #                 insert!(U, 2nbc[1]-1, nbc[2][1])
        #             end
        #             if nbc[2][2] != Inf
        #                 insert!(U, 2nbc[1], nbc[2][2])
        #             end
        #         end
        #     elseif dims == 2
        #         if length(nbc) == 1
        #             insert!(U, 2nbc[1]-1, 0)
        #             insert!(U, 2nbc[1], 0)
        #         else
        #             if nbc[2][1] != Inf
        #                 insert!(U, 2nbc[1]-1, nbc[2][1])
        #             end
        #             if nbc[2][2] != Inf
        #                 insert!(U, 2nbc[1], nbc[2][2])
        #             end
        #         end
        #     end
        # end
    end
    F = K*U
    return FiniteElementAnalysisSolution(F, K, U)
end

"""
    axialstresses(elements::Dict{String, LineElements}, nodedisplacements::Vector{Float64}, dims=1)
    axialstresses(elements::Dict{String, BeamElements}, bendingmoments::Vector{Float64}, distancesfromneutralaxis::Union{<:Real, Vector{<:Real}}, dims=2)

Determine the axial stresses present in element-type system.

If **line** elements, of two-force members in their local coordinate systems.
If **beam** elements, of beam elements which assumes no direct, axial forces.

See also: `prepareelements_line`, `prepareelements_beam`, `solve`, and `transformationmatrixstar`.

# Arguments
- `elements`: the prepared dictionary according to element type.
- `nodedisplacements`: the vector for displacement of each node in the global coordinate system.
- `bendingmoments`: the vector of reaction moments at each node from `solve.F`.
- `distancesfromneutralaxis`: the maximum distance(s) of nodes from neutral axis.
- `dims=1`: the number of dimensions being analyzed.

# Examples
```jldoctest; output=false
A, E, L = 1, 10e6, 100          # [in², psi, in]
angles = [120, 0, 210]          # [°]
L = L ./ abs.(cosd.(angles))    # [in]
elements = Dict("1"=>(1, 2), "2"=>(1, 3), "3"=>(1, 4))
preparedelements = prepareelements_line(elements, A, E, L, angles, dims=2)
nodeboundaryconditions = [2, 3, 4]
F_applied = [(1, (1e3, 1e3))]   # (node, (force_x [lb], force_y [lb]))
U = solve(preparedelements, nodeboundaryconditions, F_applied, dims=2).U
axialstresses(preparedelements, U, dims=2)

# output

3-element Vector{Float64}:
 -577.3502691896259
 -422.6497308103743
 1000.0000000000001
```

```jldoctest; output=false
A, E, I, L = 0.004887968, 200e9, 102e6 * (1/1e3)^4, 6               # [m², Pa, m⁴, m]
Q, c, t_w = 323473.474 * (1/1e3)^3, 176.5e-3, 6.48e-3               # [m³, m, m]
W(x) = -10e3(x)                                                     # [N/m]
F_applied = ([(W, 1, 2, 3, 4)])
numberofelements = 10
L /= numberofelements
x = 0:L:18
elements, nodeboundaryconditions = Dict{String, Tuple{Int, Int}}(), []
for node ∈ 1:1:numberofelements+1
    if node == 1
        push!(nodeboundaryconditions, node)
    else
        elements["\$(node-1)"] = (node-1, node)
        if x[node]%6 == 0
            push!(nodeboundaryconditions, (node, (0, Inf)))
        end
    end
end
preparedelements = prepareelements_beam(elements, E, I, L)
sol = solve(preparedelements, nodeboundaryconditions, F_applied)    # [m]
V = [sol.F[2i - 1] for i ∈ 1:1:length(sol.F)÷2]                     # [N]
M = [sol.F[2i] for i ∈ 1:1:length(sol.F)÷2]                         # [N-m]
maximum(axialstresses(preparedelements, M, c))                      # [Pa]

# output

519117.64705877303
```
"""
function axialstresses(
    elements            ::Dict{String, LineElements},
    nodedisplacements   ::Vector{Float64};
    dims                ::Integer=1
)::Vector{Float64}
    axialstresses = Vector{Float64}([])
    for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
        i, (m, n), k, A = parse(Int64, key), value.nodes, value.k, value.A
        if dims == 1
            push!(axialstresses, k*(nodedisplacements[n] - nodedisplacements[m])/A)
        elseif dims == 2
            C = [-1 1]*transformationmatrixstar(value.θ)
            d = [nodedisplacements[2m - 1:2m]... nodedisplacements[2n - 1:2n]...]'
            push!(axialstresses, k*(C*d)[1]/A)
        end
    end
    return axialstresses
end

function axialstresses(
    elements                    ::Dict{String, BeamElements},
    bendingmoments              ::Vector{Float64},
    distancesfromneutralaxis    ::Union{<:Real, Vector{<:Real}};
    dims                        ::Integer=2
)::Vector{Float64}
    axialstresses = Vector{Float64}([])
    for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
        i, (m, n), k, I = parse(Int64, key), value.nodes, value.k, value.I
        M = bendingmoments[n]
        c = float(length(distancesfromneutralaxis) == 1 ? distancesfromneutralaxis[1] : distancesfromneutralaxis[i])
        if dims == 2
            push!(axialstresses, M*c/I)
        # elseif dims == 2
        #     C = [-1 1]*transformationmatrixstar(value.θ)
        #     d = [nodedisplacements[2m - 1:2m]... nodedisplacements[2n - 1:2n]...]'
        #     push!(axialstresses, k*(C*d)[1]/A)
        end
    end
    return axialstresses
end

"""
    shearstresses(elements, shearforces, firstmomentofareas, bendingthicknesses, dims=2)

Determine the axial stresses present in **beam** element system.

See also: `prepareelements_beam`, `solve`, and `transformationmatrixstar`.

# Arguments
- `elements::Dict{String, BeamElements}`: the prepared dictionary according to element type.
- `shearforces::Vector{Float64}`: the vector for reaction, shear forces at each node in the global coordinate system.
- `firstmomentofarea::Union{<:Real, Vector{<:Real}}`: the first moment of area(s), ``Q`` above the neutral axis.
- `bendingthicknesses::Union{<:Real, Vector{<:Real}}`: the web thickness of I-beam.
- `dims::Integer=2`: the number of dimensions being analyzed.

# Examples
```jldoctest; output=false
A, E, I, L = 0.004887968, 200e9, 102e6 * (1/1e3)^4, 6               # [m², Pa, m⁴, m]
Q, c, t_w = 323473.474 * (1/1e3)^3, 176.5e-3, 6.48e-3               # [m³, m, m]
W(x) = -10e3(x)                                                     # [N/m]
F_applied = ([(W, 1, 2, 3, 4)])
numberofelements = 10
L /= numberofelements
x = 0:L:18
elements, nodeboundaryconditions = Dict{String, Tuple{Int, Int}}(), []
for node ∈ 1:1:numberofelements+1
    if node == 1
        push!(nodeboundaryconditions, node)
    else
        elements["\$(node-1)"] = (node-1, node)
        if x[node]%6 == 0
            push!(nodeboundaryconditions, (node, (0, Inf)))
        end
    end
end
preparedelements = prepareelements_beam(elements, E, I, L)
sol = solve(preparedelements, nodeboundaryconditions, F_applied)    # [m]
V = [sol.F[2i - 1] for i ∈ 1:1:length(sol.F)÷2]                     # [N]
M = [sol.F[2i] for i ∈ 1:1:length(sol.F)÷2]                         # [N-m]
maximum(shearstresses(preparedelements, V, Q, t_w))                 # [Pa]

# output

366682.5532475419
```
"""
function shearstresses(
    elements            ::Dict{String, BeamElements},
    shearforces         ::Vector{Float64},
    firstmomentofareas  ::Union{<:Real, Vector{<:Real}},
    bendingthicknesses  ::Union{<:Real, Vector{<:Real}};
    dims                ::Integer=2
)::Vector{Float64}
    shearstresses = Vector{Float64}([])
    for (key, value) ∈ sort(collect(pairs(elements)), by=x->x[1])
        i, (m, n), k, I = parse(Int64, key), value.nodes, value.k, value.I
        V = shearforces[n]
        Q = float(length(firstmomentofareas) == 1 ? firstmomentofareas[1] : firstmomentofareas[i])
        t = float(length(bendingthicknesses) == 1 ? bendingthicknesses[1] : bendingthicknesses[i])
        if dims == 2
            push!(shearstresses, V*Q/I/t)
        end
    end
    return shearstresses
end
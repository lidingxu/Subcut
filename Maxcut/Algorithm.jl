#=========================================================
 Algorithm
=========================================================#
using LinearAlgebra, JuMP, MosekTools

mutable struct Stat
    sol_val::Float64
    sol_bd::Float64
    time::Float64

    function Stat()
        stat = new()
        return stat
    end   
end

mutable struct Problem
    # parameters
    nv::Int
    ne::Int
    weights

 
    # construtor
    function Problem(nv::Int, ne::Int, weights)
        problem = new(nv, ne,  weights)
        return problem
    end
end

function solve!(problem::Problem, time_limit)
    # prepare the problem
 
    msk = Model(MosekTools.Optimizer)
    set_optimizer_attribute(msk, "INTPNT_MULTI_THREAD", 0)
    set_optimizer_attribute(msk, "OPTIMIZER_MAX_TIME", time_limit) 
    

    stat = solveSDP!(problem, msk)


    #verification(problem, sol)

    #println("Optimal value: ", objective_value(cflg))
    println(stat)

    return stat
end




function solveSDP!(problem::Problem, model)
    nv = problem.nv
    ne = problem.ne
    weights  = problem.weights
    w = zeros(nv, nv)
    ns = collect(1:nv)
    ms = collect(1:ne)
    for i  in ms
        i_ = weights[i][1]
        j_ = weights[i][2]
        w_ = weights[i][3]
        w[i_,j_] = w_
        w[j_,i_] = w_
    end
    laplacian = LinearAlgebra.diagm(0 => w * ones(nv)) - w
    @variable(
        model,
        X[i = 1:nv, j = 1:nv],
        PSD
    )
    @objective(model, Max, 1 / 4 * LinearAlgebra.dot(laplacian, X))
    @constraint(model, LinearAlgebra.diag(X) .== 1)
    optimize!(model)
    stat = Stat()
    stat.sol_val = objective_value(model)
    stat.sol_bd = objective_bound(model)
    stat.time = solve_time(model)
    return stat
end

function readMaxcut(absolute_path)
    print(absolute_path)
    lines = readlines(absolute_path)
    line1 = split.(lines[1], " ")
    nv = parse(Int, line1[1])
    ne = parse(Int, line1[2])
    weights = []
    row = 1
    for line in lines[2:end]
        line = Vector(split.(line, " "))
        v1 = parse(Int, line[1])
        v2 = parse(Int, line[2])
        w = parse(Float64, line[3])
        push!(weights,(v1, v2,w))
    end
    #print(weights)
    return Problem(nv, ne,  weights)
end

function main(args)
    @assert(length(args) == 2)
    time_limit =  parse(Float64, args[1])
    instance = args[2]
    problem = readMaxcut( string("benchmark/", instance ))
    stat = solve!(problem, time_limit)
    istat_info = string("instance: ", instance, "\n", "bound: ", string(stat.sol_val), "\n", "time: ", string(stat.time))
    #print(stat)
end


main(ARGS)

using DynamicPolynomials
using EigenvalueSolver

function bench()
    @polyvar t[1:3]
    β = [-13 -1 -1 24 -1; -13 -1 -1 24 -1; -13 -1 -1 24 -1]
    mons1 = [1 t[2]^2 t[3]^2 t[2]*t[3] t[2]^2*t[3]^2]
    mons2 = [1 t[3]^2 t[1]^2 t[3]*t[1] t[3]^2*t[1]^2]
    mons3 = [1 t[1]^2 t[2]^2 t[1]*t[2] t[1]^2*t[2]^2]
    f = [β[1,:]'*mons1';β[2,:]'*mons2';β[3,:]'*mons3'][:]
    @time sol, A₀, E, D = EigenvalueSolver.solve_CI_mixed(f,t;verbose=true)
    BWEs = EigenvalueSolver.get_residual(f,sol,t)
    BWE = maximum(BWEs)
    println("Maximal backward error $(BWE)")
end
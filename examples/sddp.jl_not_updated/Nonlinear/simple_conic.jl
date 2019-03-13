#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Ipopt, Base.Test

"""
A trivial two-stage stochastic programming problem to work out the
least-squares estimator for the set of `data`. The first stage chooses `x`,
and the second-stage evaluates the decision.
"""
function least_squares(data::Vector{Float64},radius=0.0)
    m = SDDPModel(
        stages          = 2,
        sense           = :Min,
        objective_bound = 0.0,
        solver          = IpoptSolver(print_level=0),
        risk_measure    = DRO(radius)
            ) do sp, t

        @state(sp, x′>=0, x==0)
        if t == 2
            @stageobjective(sp, ω=data, (x - ω)^2)
        else
            @stageobjective(sp, 0.0)
        end
    end
end

"""
Solve the model `m` and return the estimator for `data`.
"""
function solvemodel!(m)
    status = solve(m, iteration_limit=20, print_level=0)
    @test status == :iteration_limit
    d = simulate(m, 1, [:x′])
    return d[1][:x′][1]
end

# =============================================
#   Tests
# =============================================
srand(1234)
N = 50
data = rand(N)
m = least_squares(data)
solvemodel!(m)
@test isapprox(getbound(m), mean((data - mean(data)).^2), atol=1e-6)

# =============================================
#   SDDP.jl chapter
# =============================================
function getfittedvalue(N, radius)
    data = rand(N)
    m = least_squares(data, radius)
    solvemodel!(m)
end

fits = Vector{Float64}[]
for r in 0.0:0.1:0.5
    srand(1234)
    push!(fits, [getfittedvalue(N, r) for i in 1:100])
end
fitted_estimates = hcat(fits...)

using Plots, StatPlots
fntsm = Plots.font("times", 10.0Plots.pt)
fntlg = Plots.font("times", 12.0Plots.pt)
default(
    titlefont=fntlg, guidefont=fntlg, tickfont=fntsm,legendfont=fntsm,
    left_margin=6Plots.mm,bottom_margin=6Plots.mm,
    legend=false,
    size=(500,300),
    c="#00467F",
    xticks=(1:6, 0.0:0.1:0.5),
    xlabel="DRO Radius"
)
p1 = boxplot(fitted_estimates,
    title="(a)",
    ylabel="Estimator",
    ylims=(0.35, 0.65)
)
p2 = plot(
    (std(fitted_estimates, 1)').^2,
    w=3,
    title="(b)",
    ylabel="VAR(Estimator)",
    ylims=(0, 0.002)
)
plot(p1, p2, layout=(1,2), size=(1000, 300))
savefig("conic.pdf")

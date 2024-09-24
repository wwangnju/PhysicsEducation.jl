using DifferentialEquations, Plots
plotlyjs() 
function dampedoscillator!(d²x, dx, x, p, t)
    ξ, ω₀, f, ω = p
    d²x .= -2*ξ*ω₀*dx .- ω₀^2 * x .+ f*cos(ω*t)
end
x₀ = [3.0]
dx₀ = [0.10]
ξ, ω₀, f= 0.3, 0.8, 1.0
ω =  ω₀*sqrt(1-2*ξ^2)
p = [ξ, ω₀, f, ω]
tspan = (0.0, 40.0)

prob = SecondOrderODEProblem(dampedoscillator!, dx₀, x₀, tspan, p)
sol = solve(prob)
p₂ = [ξ, ω₀, 1.0, 1.9*√(1-2*ξ^2)]
prob₂ = SecondOrderODEProblem(dampedoscillator!, dx₀, x₀, tspan, p₂)
sol₂ = solve(prob₂)
plt=plot(sol, idxs=[2], linewidth = 2, ylims=(-4,4), lc=:black, label = "共振")
plot!(plt, sol₂, idxs=[2], ls=:dash, lw=2, label="非共振", xlabel = "时间/s", ylabel = "x/m", lc=:black, legend=:topright)
savefig(plt,"受迫振动数值解.png")



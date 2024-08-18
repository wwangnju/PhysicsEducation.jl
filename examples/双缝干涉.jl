using PhysicsEducation
using GLMakie
using Makie.Colors
function doubleslit_figure()
    A, λ = 1.0, 600.0e-9 #unit m
    d = 0.08e-3 #m
    t0 = 0.0
    p = [A, λ, d, 0.0]
    xy = Vector{Float64}[collect(0.01:1e-3:1.0), collect(-0.05:1e-5:0.05) ]
 
    doubleslit = WaveFunction(doubleslit_intensity, t0, xy, p)
    dso = WaveObservable([doubleslit]; current_step=0, Δt=0.01)

    ps = [Dict(2=>4e-7:1e-8:8e-6, 3=>1e-5:2e-6:0.001, 4=>0.0:0.1:2pi)]
    pnames = [Dict(2=>"波长λ", 3=>"缝距d", 4=>"初相位差θ")]
    colors = [:grays]#[:viridis]

    fig = Figure(;resolution = (1000, 1000))

    intensity_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        colors=colors, 
        axis=(; limits=((minimum(xy[1]), maximum(xy[1])), (minimum(xy[2]),maximum(xy[2]))),     xminorgridvisible=true, yminorgridvisible=true), 
        plotkwargs=(xminorgridvisible=true, yminorgridvisible=true,colorrange=(0,1e-3)),
        scale=(x->x),
        sliderkwargs=(;linestyle=:dash, color=:white)
    )
    Label(fig[0, :], "Double Slit", fontsize = 30)
    # Label(fig[0, :], "杨氏双缝干涉", fontsize = 30)
    return fig, dso
end
fig, = doubleslit_figure()
fig

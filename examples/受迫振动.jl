using PhysicsEducation
using GLMakie
using Makie.Colors
function damped_figure()
    Δt, tail = 0.01, 5000

    ω₀, ξ = 0.5, 0.5
    v₀, x₀ = 0.0, 3.0
    f, ω = 1.0, 0.4
    t0, xmax = 0.0, 5
    p = [ω₀, ξ, x₀, v₀, f, ω]
    dv = DynamicalFunction(forcevibration, t0, p)
    dso = TrajectoryObservable([dv]; tail=tail, current_step=0, Δt=Δt)
    
    fig = Figure(;resolution = (1200, 800))
    
    ps = [Dict(1=>0.1:0.1:5, 2=>0.0:0.1:2.0, 3=> 0.0:0.1:5, 4=>0.0:0.1:5, 5=>0.0:0.5:4, 6=>0.1:0.1:10)]
    pnames = [Dict(1=>"固有角频率ω₀", 2=>"阻尼比"*"ξ", 3=>"初始位置"*"x₀", 4=>"初始速度"*"v₀", 5=>"驱动力"*"f", 6=>"驱动角频率"*"ω")]
    colors = [:black, ]
    trajectory_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        colors=colors, 
        axis=(; limits=((-xmax, xmax), (-xmax, xmax)), xminorgridvisible=true, yminorgridvisible=true), 
        plotkwargs=(xminorgridvisible=true, yminorgridvisible=true)
    )
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    timeserieslayout = fig[:, 2] = GridLayout()
    fs = [1, ]
    labels = L"x(t)"
    axs=timeseries_plots!(
        timeserieslayout, dso, fs; 
        colors=colors, 
        timeseries_names=["x/m", ], 
        linekwargs=(linewidth=3.0, label=labels, fontsize=30),
        timeseries_ylims=[(-xmax, xmax), (-xmax, xmax)], 
        axis=(;xminorgridvisible=true, yminorgridvisible=true),
        timeunit = 1,
        timelabel="时间/s"
    )
    
    axislegend(axs[1])
    Label(fig[0,1], L"\frac{d^2x}{dt^2} + 2ξω_0\frac{dx}{dt} + ω_0^2x = f\cos(ωt)", tellwidth=false, tellheight=true, fontsize = 20)
    # Label(fig[0,2], labels2, tellwidth=false, tellheight=true, fontsize = 20)
    Label(fig[-1, :], "受迫振动", fontsize = 20)
    return fig, dso
end
fig, = damped_figure()
fig
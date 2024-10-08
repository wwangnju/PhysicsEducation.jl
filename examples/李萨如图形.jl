using PhysicsEducation
using GLMakie
using Makie.Colors
function lissajous_figure()
    Δt, tail = 0.01, 5000

    A₁,  θ₁ = 2., 0.0
    A₂, θ₂ = 2., 0.0
    ω₁, ω₂ = 0.2, 0.2

    t0, x0 = 0.0, 5
    
    p = [ω₁, A₁, θ₁, ω₂, A₂, θ₂]
    
    lsa = DynamicalFunction(lissajous, t0, p)
    
    fig = Figure(;resolution = (1000, 1000))
    
    dso = TrajectoryObservable([lsa]; tail=tail, current_step=0, Δt=Δt)
    
    ps = [Dict(1=>0.1:0.01:1, 2=>1.0:1:4.0, 3=> 0.0:0.01:2, 4=>0.1:0.01:1, 5=>1.0:1:4.0, 6=> 0.0:0.01:2)]
    pnames = [Dict(1=>(L"f_x"), 2=>(L"A_x"), 3=>L"θ_x", 4=>L"f_y", 5=>L"A_y", 6=>L"θ_y")]
    formats = [Dict(1=>"{:.2f}Hz", 2=>"{:}m", 3=>"{:.2f}rad/π", 4=>"{:.2f}Hz", 5=>"{:}m", 6=>"{:.2f}rad/π")]
    colors = [:blue, ]
    trajectory_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        formats=formats,
        colors=colors, 
        axis=(;aspect=1.0, limits=((-x0, x0), (-x0, x0)), xminorgridvisible=true, yminorgridvisible=true), 
        plotkwargs=(xminorgridvisible=true, yminorgridvisible=true)
    )
    timeserieslayout = fig[:, 2] = GridLayout()
    fs = [1, 2]

    timeseries_plots!(
        timeserieslayout, dso, fs; 
        colors=colors, 
        timeseries_names=["x/m", "y/m"], 
        linekwargs=(linewidth=3.0, label=nothing),
        timeseries_ylims=[(-x0, x0), (-x0, x0)], 
        axis=(;xminorgridvisible=true, yminorgridvisible=true),
        timeunit = 1,
        timelabel="时间/s"
    )
    Label(fig[0, :], "李萨如图形", fontsize = 30)
    return fig, dso
end
fig, =lissajous_figure()
fig
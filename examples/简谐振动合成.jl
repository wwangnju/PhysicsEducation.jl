using PhysicsEducation
using GLMakie
using Makie.Colors
function harmonicmotion()
    Δt = 0.01
    A₁,  θ₁ = 2., 0.0
    A₂, θ₂ = 2., 0.0
    ω₁, ω₂ = 1/20, 1/20
    t, x0 = 0.0, 5
    
    p₁ = [ω₁, A₁, θ₁]
    p₂ = [ω₂, A₂, θ₂]
    
    shm₁ = DynamicalFunction(simpleharmonic, t, p₁)
    shm₂ = DynamicalFunction(simpleharmonic, t, p₂)
    
    shms = [shm₁, shm₂]
    
    tail = 5000
    fig = Figure(; resolution= (1000, 1000), fontsize=20, px_per_unit=4)
    
    dso = TrajectoryObservable(shms; tail=tail, current_step=0, Δt=Δt, mapdf=Base.:+)
    
    ps = [Dict(1=>0.01:0.01:2, 2=>0.0:1:4.0, 3=> 0.0:0.05:2pi), Dict(1=>0.01:0.01:2, 2=>0.0:1:4.0, 3=> 0.0:0.05:2pi)]
    pnames = [Dict(1=>"频率f₁", 2=>"振幅A₁", 3=>"初相位θ₁"), Dict(1=>"频率f₂", 2=>"振幅A₂", 3=>"初相位θ₂")]
    colors = [:black, :black, :black]#[:blue, :purple, :green]
    ax = trajectory_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        colors=colors, 
        axis=(;aspect=1.0,limits=((-x0, x0), (-x0, x0)), xminorgridvisible=true, yminorgridvisible=true), 
        plotkwargs=[(xminorgridvisible=true, yminorgridvisible=true, linestyle=:dashdot), (xminorgridvisible=true, yminorgridvisible=true, linewidth=3.0, linestyle=:dot), (xminorgridvisible=true, yminorgridvisible=true, linestyle=:solid)]
    )
    
    superposition_plot!(ax, dso.state_observable; colors=colors)
    timeserieslayout = fig[:, 2] = GridLayout()
    fs = [1, ]
    label3 = L"x_3=x_1+x_2"
    axtime=timeseries_plots!(
        timeserieslayout, dso, fs; 
        colors=colors, 
        timeseries_names=["位移/m"], 
        linekwargs=[[
            (linewidth=3.0, linestyle=:dashdot, label=L"x_1=A_1\cos(ω_1t+θ_1)"),
            ((linewidth=3.0, linestyle=:dot, label=L"x_2=A_2\cos(ω_2t+θ_2)")), 
            ((linewidth=3.0, linestyle=:solid, label=label3))
        ],],
        timeseries_ylims=[(-x0, x0)], 
        axis=(;xminorgridvisible=true, yminorgridvisible=true),
        timeunit = 1+0*Δt*200,
        timelabel="时间/s"
    )
    axislegend(axtime[end])
    Label(fig[0, :], "简谐振动的合成", fontsize = 30)
    return fig, dso
end
fig, = harmonicmotion()
fig

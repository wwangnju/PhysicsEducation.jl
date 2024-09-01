
using PhysicsEducation
using GLMakie
using Makie.Colors
function standingwave_figure()
    Δt = 0.01
    A₁,  θ₁ = 2., 0.0
    A₂, θ₂ = 2., 0.0
    ω₁, ω₂ = 0.1, 0.1
    k₁, k₂ = -2pi/2.0, 2pi/2
    t0, xmax = 0.0, 4
    
    p₁ = [ω₁, k₁, A₁, θ₁, 4.0]
    p₂ = [ω₂, k₂, A₂, θ₂, -4.0]
    xy = Vector{Float64}[collect(-8:0.01:8), ]
    sw₁ = WaveFunction(harmonicwave1D, t0, xy, p₁)
    sw₂ = WaveFunction(harmonicwave1D, t0, xy, p₂)
    dso = WaveObservable([sw₁, sw₂]; mapdf=Base.:+, current_step=0, Δt=Δt)

    fig = Figure(;resolution = (1000, 1000))

    ps = [Dict(1=>0.01:0.01:pi, 2=>-2pi:0.01:2pi, 3=>1.0:1:4.0, 4=> 0.0:0.05:2pi), Dict(1=>0.01:0.01:pi, 2=>-2pi:0.01:2pi, 3=>1.0:1:4.0, 4=> 0.0:0.05:2pi)]
    pnames = [Dict(1=>"角频率ω₁",2=>"波矢k₁", 3=>"振幅A₁", 4=>"初相位θ₁"), Dict(1=>"角频率ω₂", 2=>"波矢k₂",3=>"振幅A₂", 4=>"初相位θ₂")]
    formats=[Dict(1=>"{:.2f}rad/s",2=>"{:.2f}rad/m", 3=>"{:}m", 4=>"{:.2f}rad"), Dict(1=>"{:.2f}rad/s",2=>"{:.2f}rad/m", 3=>"{:}m", 4=>"{:.2f}rad")]
    colors = [:blue, :purple, :green]
    ax=trajectory_plot!(fig, dso; 
        parameter_sliders=ps, 
        parameter_names=pnames, 
        formats=formats,
        colors=colors, 
        axis=(; limits=((-xmax, xmax), (-xmax, xmax)), xminorgridvisible=true, yminorgridvisible=true), 
        plotkwargs=[
            (xminorgridvisible=true, yminorgridvisible=true, linestyle=:dashdot,linewidth=3.0,label=L"y_1=A_1\cos(\omega_1t-k_1x+\theta_1)"),
            (xminorgridvisible=true, yminorgridvisible=true, linestyle=:dot,linewidth=3.0,label=L"y_2=A_2\cos(\omega_2t-k_2x+\theta_2)"),
            (xminorgridvisible=true, yminorgridvisible=true, linestyle=:solid,linewidth=3.0,label=L"y_3=y_2+y_1")
        ]
    )
    axislegend(ax)
    Label(fig[0, :], "驻波", tellwidth=false,fontsize = 30)
    return fig, dso
end
fig, = standingwave_figure()
fig

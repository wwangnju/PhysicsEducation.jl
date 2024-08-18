module Wave
using Makie
using SpecialFunctions: besselj1
using ..Vibration: _trajectory_plot_controls!, _add_ds_param_controls!
# using ..PhysicsEducation: AbstractStateObservable

import ..PhysicsEducation: reinit!, current_state, current_parameters, step!, set_parameter!,  set_state!, trajectory_plot!

export harmonicwave1D, harmonicwave2D, dopplereffect
export WaveObservable, WaveFunction, intensity_plot!, doubleslit_intensity,michelsoninterferometer_intensity
export doubleslit_fraunhofer, singleslit_fraunhofer, rectangular_fraunhofer, circular_fraunhofer, slit_grating, newton_ring_intensity


mutable struct WaveFunction{T<:Real, P}
    f::Function
    t::T
    xycoord::Vector{Vector{Float64}}
    p::Vector{P}
    function WaveFunction(f::Function, t::Real, xycoord::Vector{Vector{Float64}}, p::Vector)
        new{typeof(t), eltype(p)}(f, t, xycoord, p)  
    end
end
current_state(df::WaveFunction) = df.f(df.t, df.xycoord, df.p)
current_parameters(df::WaveFunction) = df.p
step!(ds::WaveFunction, Δt::Real) = (ds.t += Δt;)
reinit!(df::WaveFunction, t0::Real) = (df.t = t0;)
set_parameter!(df::WaveFunction, index::Int, value) = (df.p[index] = value;)
function set_parameter!(df::WaveFunction, value::Vector)  
    for i in eachindex(value)
        set_parameter!(df, i, value[i])
    end
    return nothing
end


struct WaveObservable{T<:WaveFunction} 
    wfs::Vector{T} # reference to the dynamical system
    state_observables::Vector{Observable} # [Observable(Matrix),...]
    param_observables::Vector{Observable}
    current_step::Observable{Int}
    Δt::Real # a default value for `step!`
end
function WaveObservable(wfs::Vector{S}; current_step::Int=0, Δt::Real=0.001, mapdf::Union{Function, Nothing}=nothing) where S<:WaveFunction
    stateobs = Observable[]
    for wf in wfs
        u = current_state(wf)
        push!(stateobs, Observable(u))
    end   
    !isa(mapdf, Nothing) && push!(stateobs, map(mapdf, stateobs...))
    param_observables = Observable[Observable(deepcopy(current_parameters(wf))) for wf in wfs]
    return WaveObservable{S}(wfs, stateobs, param_observables, Observable(current_step), Δt)
end
"""
    step!(dso::WaveObservable, n::Int = 1) -> Nothing

To update the states which is stored in `dso`.
"""
function step!(dso::WaveObservable, n::Int = 1)
    Δt = dso.Δt
    for _ in 1:n
        for i in eachindex(dso.wfs)
            step!(dso.wfs[i], Δt)
            last_state = current_state(dso.wfs[i])
            dso.state_observables[i][] =  copy(last_state)
        end
    end
    dso.current_step.val = dso.current_step[] + n
    notify.(dso.state_observables)
    return nothing
end
function set_state!(dso::WaveObservable, t0::Real, i::Int = 1)
    dso.current_step.val = 0
    reinit!(dso.wfs[i], t0)
    u = current_state(dso.wfs[i])
    dso.state_observables[i][] = u
    # for j in eachindex(val); val[j] = u; end
    notify(dso.state_observables[i])
    return nothing
end
function set_parameter!(dso::WaveObservable, i::Int, index::Int, value)
    dso.param_observables[i][][index] = value
    set_parameter!(dso.wfs[i], index, value)
    notify(dso.param_observables[i])
    return nothing
end

"""
    trajectory_plot!(fig::Figure, dso::WaveObservable;
        parameter_sliders=nothing, 
        parameter_names=isnothing(parameter_sliders) ? nothing : [Dict(keys(ps) .=> string.(keys(ps)).*"(j)") for (j, ps) in enumerate(parameter_sliders)], 
        colors=[:blue for _ in 1:length(dso.state_observables)], 
        axis=NamedTuple(), 
        plotkwargs=NamedTuple(),
        is3D::Bool=false
    )-> Axis

Plot the trajectories.
"""
function trajectory_plot!(fig::Figure, dso::WaveObservable;
    parameter_sliders=nothing, 
    parameter_names=isnothing(parameter_sliders) ? nothing : [Dict(keys(ps) .=> string.(keys(ps)).*"($j)") for (j, ps) in enumerate(parameter_sliders)], 
    colors=[:blue for _ in 1:length(dso.state_observables)], 
    axis=NamedTuple(), 
    plotkwargs=NamedTuple(),
    is3D::Bool=false
)
    t0 = [deepcopy(ds.t) for ds in dso.wfs]
    p0 = [deepcopy(current_parameters(ds)) for ds in dso.wfs]
    dims = dso.wfs[1].xycoord|>length 
    xycoords = Observable[Observable(ds.xycoord) for ds in dso.wfs]

    statespacelayout = fig[1,1] = GridLayout() 
    statespaceax = _init_statespace_plot!(statespacelayout, xycoords, dso.state_observables, dims, colors, axis, plotkwargs, is3D)
    reset, run, step, stepslider = _trajectory_plot_controls!(statespacelayout)
    isrunning = Observable(false)
    on(run) do c; isrunning[] = !isrunning[]; end
    on(run) do c
        @async while isrunning[]
            step[] = step[] + 1
            isopen(fig.scene) || break # ensures computations stop if closed window
            yield() 
        end
    end
    on(step) do clicks
        n = stepslider[]
        step!(dso, n)
    end   
    on(reset) do clicks
        for j in eachindex(t0)
            set_state!(dso, copy(t0[j]), j)
        end
    end
    if !isnothing(parameter_sliders)
        paramlayout = fig[2, :] = GridLayout(tellheight = true, tellwidth = false)
        _add_param_controls!(fig, paramlayout, dso, p0, parameter_sliders, parameter_names)
    end
    return statespaceax
end
function _init_statespace_plot!(layout::GridLayout, xycoords::Vector{Observable}, state_observables::Vector{Observable}, dims::Int, colors, axis=NamedTuple(), plotkwargs=NamedTuple(), is3D::Bool=false)
    statespaceax = !is3D ? Axis(layout[1, 1]; xlabel = "x/m", ylabel ="y/m", axis... ) : Axis3(layout[1, 1]; xlabel = "x/m", ylabel ="y/m", zlabel="z/m", axis... )
    xss = [map(x->first(x), xycoord) for xycoord in xycoords ]
    yss = dims==2 ? [map(x->x[2], xycoord) for xycoord in xycoords ] : nothing

    for (i, ob) in enumerate(state_observables)
        pk = plotkwargs isa Vector ? plotkwargs[i] : plotkwargs 
        if dims == 1
            x = to_color(colors[i])
            ys = map(x->x[:, 1], ob)
            Makie.lines!(statespaceax, xss[1], ys; color = x, linewidth = 3.0, transparency = true, pk...)
        elseif dims == 2
            if is3D
                i == length(state_observables) && Makie.surface!(statespaceax, xss[1], yss[1], ob; colormap=colors[i], transparency=true, pk...)
            else
                # Makie.contour!(statespaceax, xss[1], yss[1], ob; levels=0.1:0.1:0.1, color=x, pk...)
                i == length(state_observables) && (hm=Makie.heatmap!(statespaceax, xss[1], yss[1], ob; colormap=colors[i], pk...); Colorbar(layout[0,:], hm, vertical=false))
            end
        else
            error("_init_statespace_plot error: not support 3-D.")
        end
    end
    is3D && (statespaceax.protrusions = 50)
    return statespaceax
end

function  _add_param_controls!(fig::Figure, paramlayout::GridLayout, dso::WaveObservable, p0::Vector{<:Vector}, parameter_sliders::Vector{<:Dict}, parameter_names=[Dict(keys(ps) .=> string.(keys(ps)).*"($j)") for (j, ps) in enumerate(parameter_sliders)])
    
    slidervals, sliders = _add_ds_param_controls!(
        paramlayout, parameter_sliders, parameter_names, [current_parameters(ds) for ds in dso.wfs]
    )
    update = Button(fig, label = "update", tellwidth = false, tellheight = true)
    resetp = Button(fig, label = "reset p", tellwidth = false, tellheight = true)
    gl = paramlayout[2, :] = GridLayout()
    gl[1,1] = update
    gl[1,2] = resetp
    # [update, resetp]
    on(update.clicks) do clicks
        for (i, sdvals) in enumerate(slidervals)
            for l in keys(sdvals)
                v = sdvals[l][]
                set_parameter!(dso, i, l, v)
            end
        end
    end
    on(resetp.clicks) do clicks
        for i in 1:length(dso.wfs)
            set_parameter!(dso.wfs[i], p0[i])
        # Also **visually** reset sliders to initial parameters
            for k in keys(p0[i])
                haskey(sliders[i], k) || continue
                set_close_to!(sliders[i][k], p0[i][k])
            end
        end
    end
    return nothing
end


function Θ(x::Real)
    x < 0.0 && return 1.0
    return 0.0
end

function harmonicwave1D(t::Real, xy::Vector{Vector{Float64}}, p::Vector)
    ω₁, k₁, A₁, ϕ₁ = p[1], p[2], p[3], p[4]
    u = ω₁/k₁
    x = xy[1]
    phase = -k₁*x .+ (ϕ₁ + ω₁*t )
    z = zeros(length(x), 1)
    z[:, 1] = A₁*cos.(phase) .* Θ.(sign(u)*x .- (abs(u)*t + sign(u)*p[5]) )
    return z
end
function harmonicwave2D(t::Real, xy::Vector{Vector{Float64}}, p::Vector)
    ω₁, kx, ky, A₁, ϕ₁ = p[1], p[2], p[3], p[4], p[5]
    xs, ys = xy[1], xy[2]
    phase = Float64[kx*x + ky*y + (ϕ₁ + ω₁*t ) for x in 1:length(xs), y in 1:length(ys)]
    z = zeros(length(xs), length(ys))
    z[:, :] = A₁*cos.(phase)
    return z
end

function dopplereffect(t::Real, xy::Vector{Vector{Float64}}, p::Vector)
    ν₀, u, νₛ, ϕ₁ = p[1], p[2], p[3], p[4]
    ω₁ = 2pi*ν₀
    x0 = νₛ*t
    # u = ω₁/k₁
    z = zeros(length(xy[1]), 1)
    for (i, x) in enumerate(xy[1])
        z[i, 1] = cos(ω₁*(t-abs(x-x0)/u)+ϕ₁)*(t-abs(x-x0)>=0.0)
    end
end

#optic intensity

function intensity_plot!(fig::Figure, dso::WaveObservable;
    parameter_sliders=nothing, 
    parameter_names=isnothing(parameter_sliders) ? nothing : [Dict(keys(ps) .=> string.(keys(ps)).*"($j)") for (j, ps) in enumerate(parameter_sliders)], 
    colors=[:blue for _ in 1:length(dso.state_observables)], 
    axis=NamedTuple(), 
    plotkwargs=NamedTuple(),
    scale=identity,
    slider::Symbol=:x,
    sliderkwargs=NamedTuple()
)
    p0 = [deepcopy(current_parameters(ds)) for ds in dso.wfs]
    @assert  2 == dso.wfs[end].xycoord|>length "intensity_plot! error: only support two dimenison."
    xs, ys = dso.wfs[end].xycoord[1], dso.wfs[end].xycoord[2]
    nx, ny = length(xs), length(ys)
    xycoords = Observable[Observable(ds.xycoord) for ds in dso.wfs]
    statespacelayout = fig[1, 1] = GridLayout() 
    statespaceax = _init_statespace_plot!(statespacelayout, xycoords, dso.state_observables, 2, colors, axis, plotkwargs, false)
    if slider == :x
        layout2 = fig[:, 2] = GridLayout()
        axint = Axis(layout2[:, 1]; ylabel = "相对强度", xlabel="y/m")
        sl_x = Slider(fig.layout[1, 1][2,1], range=1:1:nx, startvalue=nx,color_active=:black, color_active_dimmed=:gray)
        int_y = lift(sl_x.value, dso.state_observables[end]) do x, y 
            scale.(y[x, :])
        end
        x0 = lift(sl_x.value, xycoords[end]) do i, x
            x[1][i]
        end 
        vlines!(statespaceax, x0; linewidth=3.0, sliderkwargs...)
        lines!(axint, ys, int_y, linewidth=3.0, color=:black)
        ylims!(axint, 0, maximum(int_y[])+1e-9)
    elseif slider == :y
        sl_y = Slider(fig.layout[1, 1][1,2], range=1:1:ny, horizontal = false, startvalue=ny, color_active=:black, color_active_dimmed=:gray)
        int_y = lift(sl_y.value, dso.state_observables[end]) do x, y 
            scale.(y[:, x])
        end
        y0 = lift(sl_y.value, xycoords[end]) do i, x
            x[2][i]
        end 
        hlines!(statespaceax, y0; linewidth=3.0, sliderkwargs...)
        layout2 = fig[:, 2] = GridLayout()
        axint = Axis(layout2[:, 1]; ylabel = "相对强度", xlabel="x/m")
        lines!(axint, xs, int_y, linewidth=3.0, color=:black)
        ylims!(axint, 0, maximum(int_y[])+1e-9)
    end
    
    on(int_y) do y
        ylims!(axint, 0, maximum(y)+1e-9)
    end    
    if !isnothing(parameter_sliders)
        paramlayout = fig[2, :] = GridLayout(tellheight = true, tellwidth = false)
        slidervals, sliders = _add_ds_param_controls!(
            paramlayout, parameter_sliders, parameter_names, [current_parameters(ds) for ds in dso.wfs]
        )
        update = Button(fig, label = "update", tellwidth = false, tellheight = true)
        resetp = Button(fig, label = "reset p", tellwidth = false, tellheight = true)
        gl = paramlayout[2, :] = GridLayout()
        gl[1,1] = update
        gl[1,2] = resetp
        # [update, resetp]
        on(update.clicks) do clicks
            for (i, sdvals) in enumerate(slidervals)
                for l in keys(sdvals)
                    v = sdvals[l][]
                    set_parameter!(dso, i, l, v)
                end
            end
            step!(dso, 1)
            notify(int_y)    
        end
        on(resetp.clicks) do clicks
            for i in 1:length(dso.wfs)
                set_parameter!(dso.wfs[i], p0[i])
                for k in keys(p0[i])
                    haskey(sliders[i], k) || continue
                    set_close_to!(sliders[i][k], p0[i][k])
                end
            end
            step!(dso, 1)
            notify(int_y)
        end
    end
    return statespaceax, axint
end


"""

Young's double slit interference.
"""
function doubleslit_intensity(t::Float64, xy::Vector{Vector{Float64}}, p::Vector)
    A₁, λ₁ = p[1], p[2]
    d, θ = p[3], p[4]
    xs, ys = xy[1], xy[2]
    k₁ = 2pi/λ₁
    r₁ = [sqrt(x^2 + (y - d/2)^2) for x in xs, y in ys]
    r₂ = [sqrt(x^2 + (y + d/2)^2) for x in xs, y in ys]
    AA₁ = [r ≈ 0.0 ? A₁ : A₁/r for r in r₁]
    AA₂ = [r ≈ 0.0 ? A₁ : A₁/r for r in r₂]
    wf₁ = AA₁.*exp.(-im*(k₁*r₁ .+ θ)) + AA₂.*exp.(-im*k₁*r₂)
    int = abs2.(wf₁)
    return (int/maximum(int))
end
function newton_ring_intensity(t::Float64, xy::Vector{Vector{Float64}}, p::Vector)
    λ, R = p[1], p[2]
    k = 2pi/λ
    xs, ys = xy[1], xy[2]
    #e^2 -2R e+r^2=0
    e = [R - sqrt(R^2-x^2-y^2) for x in xs, y in ys]
    e = [(x^2+y^2)/R for x in xs, y in ys]
    cos.(k*(e .+ λ/4)).^2
end
function michelsoninterferometer_intensity(t::Float64, xy::Vector{Vector{Float64}}, p::Vector)
    A₁, λ₁ = p[1], p[2]
    d = p[3]
    Lz = 1.0
    xs, ys = xy[1], xy[2]
    k₁ = 2pi/λ₁
    r₁ = [sqrt(x^2 + (y)^2 + (Lz)^2) for x in xs, y in ys]
    r₂ = [sqrt(x^2 + (y)^2+(Lz + 2*d)^2) for x in xs, y in ys]
    AA₁ = [r ≈ 0.0 ? A₁ : A₁/r for r in r₁]
    AA₂ = [r ≈ 0.0 ? A₁ : A₁/r for r in r₂]
    wf₁ = AA₁.*exp.(-im*(k₁*r₁ )) + AA₂.*exp.(-im*k₁*r₂)
    int = abs2.(wf₁)
    return (int)
end
function doubleslit_fraunhofer(t::Float64, xy::Vector{Vector{Float64}}, p::Vector)
    λ, d, W, L = p[1], p[2], p[3], p[4]#d is the distance between two slits, W is the width of a slit.
    θs = sin.(xy[1]/L) # x/L and x << L
    ys = xy[2]
    [(cos(pi/λ* d*θ)*sinc(W*θ/λ))^2 for θ in θs, _ in ys]
end
function singleslit_fraunhofer(t::Float64, xy::Vector{Vector{Float64}}, p::Vector)
    λ, W, L = p[1], p[2], p[3]
    θs = sin.(xy[1]/L)
    ys = xy[2]
    int = [(sinc(W*θ/λ))^2 for θ in θs, _ in ys]
    return (int)
end
function rectangular_fraunhofer(t::Float64, xy::Vector{Vector{Float64}}, p::Vector)
    λ, W, H, L = p[1], p[2], p[3], p[4]
    xs = xy[1]
    ys = xy[2]
    int = [(sinc(W*x/λ/L)*sinc(H*y/(λ*L)))^2 for x in xs, y in ys] #
    return (int)
end
"""
    circular_aperture(t::Float64, xy::Vector{Vector{Float64}}, p::Vector)

[see](https://en.wikipedia.org/wiki/Fraunhofer_diffraction_equation)
"""
function circular_fraunhofer(t::Float64, xy::Vector{Vector{Float64}}, p::Vector)
    λ, W, L = p[1], p[2], p[3]
    xs = xy[1]
    ys = xy[2]

    int = [sin(atan(sqrt(x^2+y^2), L)) ≈ 0.0 ? 1/4 : (besselj1(pi*W/λ*sin(atan(sqrt(x^2+y^2), L)))/(pi*W/λ*sin(atan(sqrt(x^2+y^2), L))))^2 for x in xs, y in ys]
    return (int)
end
function slit_grating(t::Float64, xy::Vector{Vector{Float64}}, p::Vector)
    λ, d, W, L, N = p[1], p[2], p[3], p[4], p[5]
    xs = xy[1]
    ys = xy[2]
    [isapprox(d/λ*sin(atan(x,L)), round(d/λ*sin(atan(x,L))),atol=1e-14) ? N^2*sinc(W/λ*sin(atan(x,L)))^2 : sinc(W/λ*sin(atan(x,L)))^2*sin(pi*N*d/λ*sin(atan(x,L)))^2/(sin(pi*d/λ*sin(atan(x,L))))^2 for x in xs, y in ys]
end
end# module
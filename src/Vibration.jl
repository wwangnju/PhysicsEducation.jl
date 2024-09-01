module Vibration
using Makie
using DataStructures: CircularBuffer

import ..PhysicsEducation: set_parameter!, trajectory_plot!, current_state, current_parameters,  step!, set_state!, reinit!

export DynamicalFunction, TrajectoryObservable
export  simpleharmonic, timeseries_plots!, superposition_plot!, lissajous
export dampedvibration, forcevibration, _trajectory_plot_controls!, _add_ds_param_controls!

"""
    simpleharmonic(t::Real, p::Vector{Float64}) -> Vector{Float64}

Define the simple harmonic osscillator
"""
function simpleharmonic(t::Real, p::Vector{Float64})
    ω, A, ϕ = p[1], p[2], p[3]
    x = A*cos(2pi*ω*t+ϕ)
    y = A*sin(2pi*ω*t+ϕ)
    return Float64[x, y]
end
"""
    lissajous(t::Real, p::Vector{Float64}) -> Vector

Define the lissajous picture.
"""
function lissajous(t::Real, p::Vector{Float64})
    ω, A, ϕ = p[1], p[2], p[3]
    ω₂, A₂, ϕ₂ = p[4], p[5], p[6]
    x = A*cos(2pi*ω*t+ϕ)
    y = A₂*cos(2pi*ω₂*t+ϕ₂)
    return Float64[x, y]
end
"""
    dampedvibration(t::Real, p::Vector{<:Real}) -> Vector

``\frac{d^2x}{dt^2} + γ\frac{dx}{dt} + ω_0^2x = 0``
"""
function dampedvibration(t::Real, p::Vector{<:Real})
    ω₀, γ, x₀, v₀ = p[1], p[2], p[3], p[4]
    ω₁ = sqrt(complex(ω₀^2 - γ^2/4))
    temp = abs(ω₁) ≈ 0.0 ? t : sin(ω₁*t)/ω₁
    x = exp(-γ/2*t)*(x₀*cos(ω₁*t) + (v₀+γ/2*x₀)*temp)
    return Float64[real(x), 0.0] 
end
"""
    forcevibration(t::Real, p::Vector{<:Real}) -> Vector

```math
\\frac{d^2x}{dt^2} + 2ξω_0\\frac{dx}{dt} + ω_0^2x = f\\cos(ωt)
```
"""
function forcevibration(t::Real, p::Vector{<:Real})
    ω₀, ξ, x₀, v₀ = p[1], p[2], p[3], p[4]
    f, ω = p[5], p[6]
    γ = 2*ω₀*ξ
    ω₁ = sqrt(complex(ω₀^2 - γ^2/4))
    A = f*γ*ω/((ω₀^2-ω^2)^2+(γ*ω)^2)
    B = f*(ω₀^2-ω^2)/((ω₀^2-ω^2)^2+(γ*ω)^2)
    ω₁ = sqrt(complex(ω₀^2 - γ^2/4))
    temp = abs(ω₁) ≈ 0.0 ? t : sin(ω₁*t)/ω₁
    x = A*sin(ω*t)+B*cos(ω*t) + exp(-γ/2*t)*((x₀-B)*cos(ω₁*t) + (v₀-A*ω+γ/2*(x₀-B))*temp)
    return Float64[real(x), 0.0]
end
"""
    DynamicalFunction{T<:Real, P}

Define dynamical function. 
"""
mutable struct DynamicalFunction{T<:Real, P}
    f::Function
    t::T
    p::Vector{P}
    function DynamicalFunction(f::Function, t::Real, p::Vector)
        new{typeof(t), eltype(p)}(f, t, p)  
    end
end
current_state(df::DynamicalFunction) = df.f(df.t, df.p)
current_parameters(df::DynamicalFunction) = df.p
"""
    step!(ds::DynamicalFunction, Δt::Real)

Update the times.
"""
step!(ds::DynamicalFunction, Δt::Real) = (ds.t += Δt;)
reinit!(df::DynamicalFunction, t0::Real) = (df.t = t0;)
set_parameter!(df::DynamicalFunction, index::Int, value) = (df.p[index] = value;)
function set_parameter!(df::DynamicalFunction, value::Vector)  
    for i in eachindex(value)
        set_parameter!(df, i, value[i])
    end
    return nothing
end

"""
    TrajectoryObservable{T<:DynamicalFunction}

Construct the trajectory.
"""
struct TrajectoryObservable{T<:DynamicalFunction} 
    dss::Vector{T} # reference to the dynamical system
    state_observable::Observable
    tail_observables::Vector{Observable}
    param_observables::Vector{Observable}
    current_step::Observable{Int}
    Δt::Real # a default value for `step!`
end
function TrajectoryObservable(dss::Vector{S}; tail::Int=1000, current_step::Int=0, Δt::Real=0.001, mapdf::Union{Function, Nothing}=nothing) where S<:DynamicalFunction
    tailobs = Observable[]
    for ds in dss
        u = current_state(ds)
        T = typeof(u)
        cb = CircularBuffer{T}(tail)
        fill!(cb, u)
        push!(tailobs, Observable(cb))
    end
    
    !isa(mapdf, Nothing) && push!(tailobs, map(mapdf, tailobs...))
    finalpoints = Observable([x[][end] for x in tailobs])
    param_observables = Observable[Observable(deepcopy(current_parameters(ds))) for ds in dss]
    return TrajectoryObservable{S}(dss, finalpoints, tailobs, param_observables, Observable(current_step), Δt)
end
"""
    step!(dso::TrajectoryObservable, n::Int = 1) -> Nothing

To update the states which is stored in `dso`.
"""
function step!(dso::TrajectoryObservable, n::Int = 1)
    Δt = dso.Δt
    for _ in 1:n
        for i in eachindex(dso.dss)
            step!(dso.dss[i], Δt)
            ob = dso.tail_observables[i]
            last_state = current_state(dso.dss[i])
            push!(ob[], copy(last_state))
        end
    end
    dso.current_step.val = dso.current_step[] + n
    notify.(dso.tail_observables)
    dso.state_observable[] = [x[][end] for x in dso.tail_observables]
    return nothing
end
"""
    set_state!(dso::DynamicalSystemObservable, u, i::Int = 1) -> Nothing

To set state with `u`. `i` is the index of Vector `dso.dss`.
"""
function set_state!(dso::TrajectoryObservable, t0::Real, i::Int = 1)
    dso.current_step.val = 0
    reinit!(dso.dss[i], t0)
    u = current_state(dso.dss[i])
    val = dso.tail_observables[i][]
    for j in eachindex(val); val[j] = u; end
    notify(dso.tail_observables[i])
    dso.state_observable.val[i] = u
    notify(dso.state_observable)
    return nothing
end
function set_parameter!(dso::TrajectoryObservable, i::Int, index::Int, value)
    dso.param_observables[i][][index] = value
    set_parameter!(dso.dss[i], index, value)
    notify(dso.param_observables[i])
    return nothing
end

"""
    trajectory_plot!(fig::Figure, dso::TrajectoryObservable;
        parameter_sliders=nothing, 
        parameter_names=isnothing(parameter_sliders) ? nothing : [Dict(keys(ps) .=> string.(keys(ps)).*"(j)") for (j, ps) in enumerate(parameter_sliders)], 
        formats=isnothing(parameter_sliders) ? nothing : [Dict(keys(ps) .=> "{:.3f}单位") for (j, ps) in enumerate(parameter_sliders)],
        colors=[:blue for _ in 1:length(dso.tail_observables)], 
        axis=NamedTuple(), 
        plotkwargs=NamedTuple()
    ) -> Axis

Plot the trajectories.
"""
function trajectory_plot!(fig::Figure, dso::TrajectoryObservable;
    parameter_sliders=nothing, 
    parameter_names=isnothing(parameter_sliders) ? nothing : [Dict(keys(ps) .=> string.(keys(ps)).*"($j)") for (j, ps) in enumerate(parameter_sliders)], 
    formats=isnothing(parameter_sliders) ? nothing : [Dict(keys(ps) .=> "{:.3f}单位") for (j, ps) in enumerate(parameter_sliders)],
    colors=[:blue for _ in 1:length(dso.tail_observables)], 
    axis=NamedTuple(), 
    plotkwargs=NamedTuple()
)
    t0 = [deepcopy(ds.t) for ds in dso.dss]
    p0 = [deepcopy(current_parameters(ds)) for ds in dso.dss]
    statespacelayout = fig[1,1] = GridLayout() 
    statespaceax = _init_statespace_plot!(statespacelayout, dso.tail_observables, dso.state_observable, colors, axis, plotkwargs)
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
        dso.state_observable[] = [x[][end] for x in dso.tail_observables]
    end
    if !isnothing(parameter_sliders)
        paramlayout = fig[2, :] = GridLayout(tellheight = true, tellwidth = false)
        slidervals, sliders = _add_ds_param_controls!(
            paramlayout, parameter_sliders, parameter_names, [current_parameters(ds) for ds in dso.dss], formats
        )
        update = Button(fig, label = "update", tellwidth = false, tellheight = true)
        resetp = Button(fig, label = "reset p", tellwidth = false, tellheight = true)
        # paramlayout[2, 1] = update
        # paramlayout[2, 2] = resetp
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
            for i in 1:length(dso.dss)
                set_parameter!(dso.dss[i], p0[i])
            # Also **visually** reset sliders to initial parameters
                for k in keys(p0[i])
                    haskey(sliders[i], k) || continue
                    set_close_to!(sliders[i][k], p0[i][k])
                end
            end
        end
    end
    return statespaceax
end
"""
    superposition_plot!(statespaceax::Axis, plotted_finalpoints::Observable; colors=[:blue for _ in 1:length(plotted_finalpoints[])])

Plot some lines for a simple harmonic vibration.
"""
function superposition_plot!(statespaceax::Axis, plotted_finalpoints::Observable; colors=[:blue for _ in 1:length(plotted_finalpoints[])])
    xcoords = map(x->Point2f[[y[1], 0.0] for y in x], plotted_finalpoints)
    linx = lift(plotted_finalpoints) do ob
        res = [(Point2f(ob[i][1], 0.0),  Point2f(ob[i])) for i in eachindex(ob)]
        append!(res, [(Point2f(ob[i]), Point2f(ob[end])) for i in 1:length(ob)-1])
    end
    scatter!(statespaceax, xcoords, markersize=20.0, marker=:diamond, color=colors)
    arrows!(statespaceax, [Point2f(0,0) for _ in 1:length(plotted_finalpoints[])], plotted_finalpoints, linewidth=4.0, arrowsize = 20, color=colors)
    linesegments!(statespaceax, linx, linestyle=:dash, linewidth=2.0, color=:black)
end

function _init_statespace_plot!(layout::GridLayout, tailobs::Vector{Observable}, finalpoints::Observable, colors, axis=NamedTuple(), plotkwargs= NamedTuple())
    statespaceax = Axis(layout[1, 1]; xlabel = "x/m", ylabel ="y/m", axis... )
    plotted_tailobs = [map(x -> Point2f[y[[1, 2]] for y in x], ob) for ob in tailobs]
    plotted_finalpoints = map(x -> Point2f[y[[1, 2]] for y in x], finalpoints)
    tail = first(plotted_tailobs)[]|>length
    for (i, ob) in enumerate(plotted_tailobs)
        pk = plotkwargs isa Vector ? plotkwargs[i] : plotkwargs
        x = to_color(colors[i])
        x = [RGBAf(x.r, x.g, x.b, i/tail) for i in 1:tail]
        Makie.lines!(statespaceax, ob; color = x, linewidth = 3.0, transparency = true, pk...)
    end
    Makie.scatter!(statespaceax, plotted_finalpoints; color = colors, markersize = 20, marker = :circle)
    return statespaceax
end
function _trajectory_plot_controls!(layout)
    layout[2, 1] = controllayout = GridLayout(tellwidth = false)
    reset = Button(controllayout[1, 0]; label = "reset")
    run = Button(controllayout[1, 1]; label = "run")
    step = Button(controllayout[1, 2]; label = "step")
    slider_vals = vcat(1:10, 100:100:1000)
    sg = SliderGrid(controllayout[1,3],
        (label = "steps =", range = slider_vals, startvalue = 1,color_active=:gray, color_active_dimmed=:gray),
    )   
    return reset.clicks, run.clicks, step.clicks, sg.sliders[1].value
end
function _add_ds_param_controls!(paramlayout, parameter_sliders, pnames, p0, formats)
    slidervals = [Dict{keytype(ps), Observable}() for ps in parameter_sliders]# directly has the slider observables
    sliders = [Dict{keytype(ps), Any}() for ps in parameter_sliders] # for updating via reset parameters
    tuples_for_slidergrid = []
    for (j, ps) in enumerate(parameter_sliders)
        for (i, (l, vals)) in enumerate(ps)
            startvalue = p0[j][l]
            label = string(pnames[j][l])
            format=string(formats[j][l])
            push!(tuples_for_slidergrid, (;label=label, range = vals, format=format, startvalue=startvalue,color_active=:gray, color_active_dimmed=:gray))
        end
    end
    sg = SliderGrid(paramlayout[1,1], tuples_for_slidergrid...; tellheight = true)
    count = 0
    for (j, ps) in enumerate(parameter_sliders)
        for (i, (l, vals)) in enumerate(ps)
            count += 1
            slidervals[j][l] = sg.sliders[count].value
            sliders[j][l] = sg.sliders[count]
        end
    end
    return slidervals, sliders
end

# timeseries
"""
    timeseries_plots!(
        layout::GridLayout, dso::DynamicalSystemObservable, fs::Vector;  
        linekwargs = (linewidth = 3,),
        timeseries_names = [_timeseries_name(f) for f in fs],
        colors = [:blue for _ in 1:length(dso.tail_observables)],
        timeseries_ylims = [(0, 1) for f in fs],
        timelabel = "time", timeunit = 1,
    )

Plot timeseries.
"""
function timeseries_plots!(
    layout::GridLayout, dso::TrajectoryObservable, fs::Vector;  
    linekwargs = (linewidth = 3,),
    timeseries_names = [_timeseries_name(f) for f in fs],
    colors = [:blue for _ in 1:length(dso.tail_observables)],
    timeseries_ylims = [(0, 1) for f in fs],
    timelabel = "time",
    timeunit = 1,
    axis = NamedTuple()
)
    tsnames, tslims = timeseries_names, timeseries_ylims
# First, create axis
    axs = [Axis(layout[i, 1]; ylabel = tsnames[i], axis...) for i in 1:length(fs)]
    for i in 1:length(fs); ylims!(axs[i], tslims[i]); end
    linkxaxes!(axs...)
    for i in 1:length(fs)-1; hidexdecorations!(axs[i]; grid = false); end
    axs[end].xlabel = timelabel
# Create and plot the observables of the timeseries
    T = length(dso.tail_observables[1][])
    for (j, f) in enumerate(fs)
        for (i, tail) in enumerate(dso.tail_observables)
            observed_data = map(tail, dso.current_step) do x, n
                [Point2f(max(0, dso.Δt*(n - T + k)/timeunit), _obtain_data(x[k], f)) for k in eachindex(x)]
            end
            # plot them
            lk = linekwargs isa AbstractVector ? linekwargs[j][i] : linekwargs
            lines!(axs[j], observed_data; color = colors[i], lk...)
        end
    # Add a last observable trigger that changes the axis xspan
        on(dso.tail_observables[end]) do x
            n = dso.current_step[]
            xlims!(axs[end],
                max(0, dso.Δt*(n - T)/timeunit),
                max(T*dso.Δt/timeunit, dso.Δt*n/timeunit)
            )
        end
    end
    
    return axs
end

_obtain_data(x::AbstractVector, f::Int) = x[f]
_obtain_data(x::AbstractVector, f::Function) = f(x)
_timeseries_name(f::Int) = "x"*subscript(f)
_timeseries_name(f) = string(f)

end # module

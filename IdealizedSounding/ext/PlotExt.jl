# plotting for IdealizedSounding; loaded when PythonPlot is loaded too
module PlotExt

using IdealizedSounding, PythonPlot
import IdealizedSounding: plot_idealized_sounding

"""
plot_idealized_sounding(; p=exp.(range(log(1010), log(10); length=400)), fontsize=14, kw...) -> fig
Profiles of T and Td [°C], θ [K], q [g/kg] and rh [%], and u, v [m/s] against log p, lines at pressures p [hPa],
dots at the vertices. kw override sounding_params.
"""
function plot_idealized_sounding(; p=exp.(range(log(1010), log(10); length=400)), fontsize=14, kw...)
    s, v = idealized_sounding(p; kw...), idealized_vertices(; kw...)
    fig, axs = subplots(1, 4; figsize=(14, 6), sharey=true)
    panels = (((:T, -273.15, 1, "C3", "T"), (:Td, -273.15, 1, "C2", "Td")),
              ((:θ, 0, 1, "C0", "θ"),),
              ((:q, 0, 1e3, "C4", "q [g/kg]"), (:rh, 0, 100, "C2", "rh [%]")),
              ((:u, 0, 1, "C0", "u"), (:v, 0, 1, "C3", "v")))
    for (j, (pan, xl)) in enumerate(zip(panels, ("T, Td [°C]", "θ [K]", "q [g/kg], rh [%]", "wind [m/s]")))
        ax = axs[j-1]
        for (f, off, sc, c, lab) in pan
            ax.plot(sc .* (s[f] .+ off), s.p; color=c, linewidth=1.5, label=lab)
            ax.plot(sc .* (v[f] .+ off), v.p; color=c, marker=".", linestyle="none")
        end
        ax.set_xlabel(xl; fontsize=fontsize)
        ax.tick_params(axis="both"; labelsize=fontsize)
        length(pan) > 1 && ax.legend(; fontsize=fontsize-4, frameon=false)
    end
    axs[3].axvline(0; color="k", linewidth=0.5)
    axs[1].set_xscale("log")                                         # θ spans 300-800 K
    axs[0].set_yscale("log"); axs[0].set_ylim(1010, 10)
    ticks = [1000, 850, 700, 500, 300, 200, 100, 50, 20, 10]
    axs[0].set_yticks(ticks); axs[0].set_yticklabels(string.(ticks))
    axs[0].set_ylabel("p [hPa]"; fontsize=fontsize)
    fig.tight_layout()
    fig
end

end # module PlotExt

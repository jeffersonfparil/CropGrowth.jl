using WGLMakie
using Bonito, Markdown
using DataFrames, CSV, StatsBase

df = Da


WGLMakie.activate!(; use_html_widgets = true)

fig = Figure(size=(1000, 500))
gl = GridLayout(fig[2, 1], tellwidth = false)
subgl = GridLayout(gl[1, 1])

cb1 = Makie.Checkbox(subgl[1, 1], checked = false)
cb2 = Makie.Checkbox(subgl[2, 1], checked = true)
cb3 = Makie.Checkbox(subgl[3, 1], checked = true)

Label(subgl[1, 2], "Dataset A", halign = :left)
Label(subgl[2, 2], "Dataset B", halign = :left)
Label(subgl[3, 2], "Dataset C", halign = :left)
rowgap!(subgl, 8)
colgap!(subgl, 8)

ax1 = Axis(fig[1, 1])

for cb in [cb1, cb2, cb3]
    lines!(ax1, cumsum(randn(1000)), alpha = @lift($(cb.checked) ? 1.0 : 0.1))
end

ax2 = Axis(fig[1, 2])
x = rand(100)
y = x .* rand(100)
scatter!(ax2, x, y)


fig




# Page() # for Franklin, you still need to configure
# WGLMakie.activate!()
# Makie.inline!(true) # Make sure to inline plots into Documenter output!
# scatter(1:4, color=1:4)
# N = 60
# function xy_data(x, y)
#     r = sqrt(x^2 + y^2)
#     r == 0.0 ? 1f0 : (sin(r)/r)
# end
# l = range(-10, stop = 10, length = N)
# z = Float32[xy_data(x, y) for x in l, y in l]
# surface(
#     -1..1, -1..1, z,
#     colormap = :Spectral
# )
# using Observables

# App() do session::Session
#     n = 10
#     index_slider = Slider(1:n)
#     volume = rand(n, n, n)
#     slice = map(index_slider) do idx
#         return volume[:, :, idx]
#     end
#     fig = Figure()
#     ax, cplot = contour(fig[1, 1], volume)
#     rectplot = linesegments!(ax, Rect(-1, -1, 12, 12), linewidth=2, color=:red)
#     on(index_slider) do idx
#         translate!(rectplot, 0,0,idx)
#     end
#     heatmap(fig[1, 2], slice)
#     slider = DOM.div("z-index: ", index_slider, index_slider.value)
#     return Bonito.record_states(session, DOM.div(slider, fig))
# end

# App() do session
#     f, ax, pl = scatter(1:4, markersize=100, color=Float32[0.3, 0.4, 0.5, 0.6])
#     custom_info = ["a", "b", "c", "d"]
#     on_click_callback = js"""(plot, index) => {
#         // the plot object is currently just the raw THREEJS mesh
#         console.log(plot)
#         // Which can be used to extract e.g. position or color:
#         const {pos, color} = plot.geometry.attributes
#         console.log(pos)
#         console.log(color)
#         const x = pos.array[index*2] // everything is a flat array in JS
#         const y = pos.array[index*2+1]
#         const c = Math.round(color.array[index] * 10) / 10 // rounding to a digit in JS
#         const custom = $(custom_info)[index]
#         // return either a string, or an HTMLNode:
#         return "Point: <" + x + ", " + y + ">, value: " + c + " custom: " + custom
#     }
#     """

#     # ToolTip(figurelike, js_callback; plots=plots_you_want_to_hover)
#     tooltip = WGLMakie.ToolTip(f, on_click_callback; plots=pl)
#     return DOM.div(f, tooltip)
# end

# App() do session::Session
#     # We can now use this wherever we want:
#     fig = Figure(size=(300, 300))
#     contour(fig[1,1], rand(4,4))
#     card = Card(Grid(
#         Centered(DOM.h1("Hello"); style=Styles("grid-column" => "1 / 3")),
#         StylableSlider(1:100; style=Styles("grid-column" => "1 / 3")),
#         DOM.img(src="https://julialang.org/assets/infra/logo.svg"),
#         fig; columns="1fr 1fr", justify_items="stretch"
#     ))
#     # Markdown creates a DOM as well, and you can interpolate
#     # arbitrary jsrender'able elements in there:
#     return DOM.div(card)
# end

# spinner = WGLMakie.CircleSpinner(;
#     size=100,                                    # diameter in pixels
#     stroke=10,                                   # border width in pixels
#     color="red",                       # color of the spinning part
#     background_color="rgba(1, 0, 0, 0.9)",      # color of the background circle
#     duration=1                                  # rotation speed in seconds
# )
# WGLMakie.activate!(; spinner=spinner)

# WGLMakie.activate!(; use_html_widgets = true)

# fig = Figure(size=(1000, 500))
# ax = Axis(fig[1, 1])
# ax2 = Axis(fig[1, 2])
# sl = Makie.Slider(fig[2, 1], range = 0:0.1:10, startvalue = 5, tellwidth=false)
# btn = Makie.Button(fig[3, 1], label = "Click me", tellwidth=false)
# x = 0:10
# lines!(ax, x, map(y -> sin.(y .* x), sl.value))

# scatter!(ax2, 1:4, color=1:4)


# fig # Will get rendered with HTML widgets


# ########################################################
# using GLMakie
# f = Figure()

# gl = GridLayout(f[2, 1], tellwidth = false)
# subgl = GridLayout(gl[1, 1])

# cb1 = Checkbox(subgl[1, 1], checked = false)
# cb2 = Checkbox(subgl[2, 1], checked = true)
# cb3 = Checkbox(subgl[3, 1], checked = true)
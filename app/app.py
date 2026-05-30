# fitgrowthmodels(
#   df::DataFrame;
#   A = Dict(
#       :init=>minimum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
#       :lower=>0.0,
#       :upper=>maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
#   ),
#   K = Dict(
#       :init=>maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
#       :lower=>0.0,
#       :upper=>2*maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
#   ),
#   C = Dict(:init=>1.0, :lower=>1.0, :upper=>1.0),
#   Q = Dict(:init=>1.0, :lower=>1.0, :upper=>1.0),
#   B = Dict(:init=>1.0, :lower=>0.0, :upper=>10.0),
#   v = Dict(:init=>1.0, :lower=>1e-5, :upper=>10.0),
#   min_t::Int64 = 3,
#   frac_of_final::Vector{Float64} = [0.5, 0.9],
#   fit_statistic::String = "R²",
#   maxiters::Int64 = 10_000,
#   seed::Int64 = 42,
#   show_plots::Bool = false,
#   verbose::Bool = false,
# )::Tuple{DataFrame, Vector{String}}
# ```

# ### Arguments

# - `df::DataFrame`: Input data containing the required columns specified in `REQUIRED_COLUMNS = ["entries", "sites", "replications", "growing_periods", "time_points"]` and at least one trait column.
# - `A::Dict`: Search space for the parameter `A` (lower asymptote). Contains `:init`, `:lower`, and `:upper` keys. Defaults to the minimum and maximum of the trait column with `init=minimum`.
# - `K::Dict`: Search space for the parameter `K` (upper asymptote). Contains `:init`, `:lower`, and `:upper` keys. Defaults to the minimum and 2×maximum of the trait column with `init=maximum`.
# - `C::Dict`: Search space for the parameter `C`. Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=1.0`, `upper=1.0`.
# - `Q::Dict`: Search space for the parameter `Q`. Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=1.0`, `upper=1.0`.
# - `B::Dict`: Search space for the parameter `B` (growth rate). Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=0.0`, `upper=10.0`.
# - `v::Dict`: Search space for the parameter `v` (asymmetry parameter). Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=1e-5`, `upper=10.0`.
# - `min_t::Int64`: Minimum number of time points required to fit the growth model for a specific combination of entry, site, replication, and growing period. Defaults to `3`.
# - `frac_of_final::Vector{Float64}`: Percentages of the final value for which the time to reach these fraction will be calculated. Defaults to `[0.5, 0.9]`.
# - `fit_statistic::String`: The fit statistic to be used for evaluating the model. Must be one of `["R²", "RMSE", "MSE", "MAE", "ρ"]`. Defaults to `"R²"`.
# - `maxiters::Int64`: Maximum number of iterations allowed for the optimization process. Defaults to `10_000`.
# - `seed::Int64`: Random seed for reproducibility. Defaults to `42`.
# - `show_plots::Bool`: Whether to show fitted growth curve plots. Defaults to `false`.
# - `verbose::Bool`: Whether to display progress and additional information during the fitting process. Defaults to `false`.

# ### Returns

# `Tuple{DataFrame, Vector{String}}`: 
#   + The first element is a `DataFrame` containing the fitted parameters (`A`, `K`, `C`, `Q`, `B`, `v`), fit statistics, value of the growth models at ``t=0`` (`y_t0`), maximum value of the growth model (`y_max`), and time to reach specified fraction of the final value for each combination of `entries`, `sites`, `replications`, and `growing_periods`.
#   + The second element is a `Vector{String}` containing the combinations that were skipped due to insufficient data points.

# ### Notes

# - The input `DataFrame` must contain the required columns specified in the global variable `REQUIRED_COLUMNS = ["entries", "sites", "replications", "growing_periods", "time_points"]`, as well as at least one additional trait column.
# - If the `DataFrame` contains more than one trait column, only the first trait column will be used.
# - Combinations with fewer than `min_t` time points will be skipped.
# - The function uses a progress bar to indicate the fitting process if `verbose=true`.
# - The optimisation is performed using the `BBO_adaptive_de_rand_1_bin_radiuslimited()` algorithm ([details of the optimisation algorithm](https://docs.sciml.ai/Optimization/stable/optimisation_packages/blackboxoptim/)).
# - The optimisation algorithm minimises the mean squared error between the observed data `y` and the generalised logistic model. 
import dash
from dash import dcc, html, callback, Input, Output
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

app = dash.Dash(__name__)

# Load your results dataframe
df_results = pd.read_csv('sim_data.csv')

app.layout = html.Div([
	html.Div([
		html.Label("X-axis:"),
		dcc.Dropdown(id='x-axis', style={'width': '45%', 'display': 'inline-block'}),
		html.Label("Y-axis:", style={'marginLeft': '5%'}),
		dcc.Dropdown(id='y-axis', style={'width': '45%', 'display': 'inline-block', 'marginLeft': '5%'}),
	], style={'marginBottom': '20px'}),
	dcc.Graph(id='scatter-plot'),
	dcc.Store(id='selected-point'),
	html.Div(id='growth-curve-container')
])

@callback(
	[Output('x-axis', 'options'), Output('y-axis', 'options')],
	Input('x-axis', 'id')
)
def update_dropdowns(_):
	params = ['A', 'K', 'B', 'v', 'C', 'Q', 'R²', 'y_t0', 'y_max']
	options = [{'label': p, 'value': p} for p in params]
	return options, options

@callback(
	Output('scatter-plot', 'figure'),
	[Input('x-axis', 'value'), Input('y-axis', 'value')]
)
def update_scatter(x_col, y_col):
	if not x_col or not y_col:
		return go.Figure()
	
	hover_text = df_results.apply(
		lambda row: f"Entry: {row['entries']}<br>Site: {row['sites']}<br>Rep: {row['replications']}<br>Period: {row['growing_periods']}",
		axis=1
	)
	
	fig = px.scatter(df_results, x=x_col, y=y_col, hover_data=hover_text)
	fig.update_traces(hovertemplate='%{customdata}<extra></extra>')
	return fig

@callback(
	Output('growth-curve-container', 'children'),
	Input('scatter-plot', 'clickData')
)
def show_growth_curve(clickData):
	if not clickData:
		return html.Div()
	
	# Extract selected point index and create growth curve plot
	# You'll need to join with raw data to get observed values
	point_num = clickData['points'][0]['pointNumber']
	selected = df_results.iloc[point_num]
	
	return dcc.Graph(figure=go.Figure())  # Replace with actual growth curve plot

if __name__ == '__main__':
	app.run(debug=True)
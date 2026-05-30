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
from dash import dcc, html, callback, Input, Output, State
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

# ==========================================
# 1. MOCK DATA GENERATION
# ==========================================
np.random.seed(42)
n_entries = 100

df_results = pd.DataFrame({
    'entries': [f'Entry_{i}' for i in range(n_entries)],
    'sites': np.random.choice(['Site_A', 'Site_B', 'Site_C'], n_entries),
    'replications': np.random.choice([1, 2, 3], n_entries),
    'growing_periods': np.random.choice(['2022', '2023'], n_entries),
    'A': np.random.uniform(0, 5, n_entries),
    'K': np.random.uniform(50, 100, n_entries),
    'C': [1.0] * n_entries,
    'Q': [1.0] * n_entries,
    'B': np.random.uniform(0.1, 0.5, n_entries),
    'v': np.random.uniform(0.8, 1.2, n_entries),
    'R²': np.random.uniform(0.85, 0.99, n_entries),
    'y_t0': np.random.uniform(0, 5, n_entries),
    'y_max': np.random.uniform(50, 100, n_entries)
})

raw_data_list = []
for _, row in df_results.iterrows():
    t_points = np.linspace(0, 20, 10)
    # y(t) = A + (K - A) / (C + Q * exp(-B * t))^(1/v)
    y_vals = row['A'] + (row['K'] - row['A']) / ((row['C'] + row['Q'] * np.exp(-row['B'] * t_points))**(1/row['v']))
    y_vals += np.random.normal(0, 2, len(t_points)) 
    
    # Intentionally add a massive outlier to Entry_0 to test exclusions
    if row['entries'] == 'Entry_0':
        y_vals[4] += 50 
        
    for t, y in zip(t_points, y_vals):
        raw_data_list.append({
            'entries': row['entries'],
            'sites': row['sites'],
            'replications': row['replications'],
            'growing_periods': row['growing_periods'],
            'time_points': float(t),
            'trait_value': max(0, float(y)) 
        })
        
df_raw = pd.DataFrame(raw_data_list)


# ==========================================
# 2. MATHEMATICAL MODEL
# ==========================================
def generalized_logistic(t, A, K, C, Q, B, v):
    """
    Generalized Logistic Function
    y(t) = A + (K - A) / (C + Q * exp(-B * t))^(1/v)
    """
    exponent = np.clip(-B * t, -100, 100)
    base = np.maximum(C + Q * np.exp(exponent), 1e-10) # Prevent domain errors
    return A + (K - A) / (base**(1/v))


# ==========================================
# 3. APP LAYOUT & CONFIGURATION
# ==========================================
app = dash.Dash(__name__)

PARAMS = ['A', 'K', 'B', 'v', 'C', 'Q', 'R²', 'y_t0', 'y_max']
xy_options = [{'label': p, 'value': p} for p in PARAMS]

def get_filter_options(column):
    return [{'label': str(val), 'value': val} for val in sorted(df_results[column].unique())]

app.layout = html.Div(style={'fontFamily': 'Arial, sans-serif', 'padding': '20px'}, children=[
    
    # State Manager for the Active Entry
    dcc.Store(id='active-curve-state', data={'entry': None, 'site': None, 'rep': None, 'period': None}),
    
    html.H2("Growth Parameter Dashboard"),
    
    html.Div(style={'backgroundColor': '#f9f9f9', 'padding': '20px', 'borderRadius': '8px', 'marginBottom': '20px'}, children=[
        html.Div(style={'display': 'flex', 'gap': '15px', 'marginBottom': '20px'}, children=[
            html.Div(style={'flex': 1}, children=[
                html.Label("Filter Entries:", style={'fontWeight': 'bold'}),
                dcc.Dropdown(id='filter-entries', options=get_filter_options('entries'), multi=True)
            ]),
            html.Div(style={'flex': 1}, children=[
                html.Label("Filter Sites:", style={'fontWeight': 'bold'}),
                dcc.Dropdown(id='filter-sites', options=get_filter_options('sites'), multi=True)
            ]),
            html.Div(style={'flex': 1}, children=[
                html.Label("Filter Replications:", style={'fontWeight': 'bold'}),
                dcc.Dropdown(id='filter-reps', options=get_filter_options('replications'), multi=True)
            ]),
            html.Div(style={'flex': 1}, children=[
                html.Label("Filter Periods:", style={'fontWeight': 'bold'}),
                dcc.Dropdown(id='filter-periods', options=get_filter_options('growing_periods'), multi=True)
            ]),
        ]),
        html.Hr(style={'borderColor': '#eaeaea'}),
        html.Div(style={'display': 'flex', 'gap': '20px', 'marginTop': '15px'}, children=[
            html.Div(style={'flex': 1}, children=[
                html.Label("X-axis Parameter:", style={'fontWeight': 'bold'}),
                dcc.Dropdown(id='x-axis', options=xy_options, value='A', clearable=False)
            ]),
            html.Div(style={'flex': 1}, children=[
                html.Label("Y-axis Parameter:", style={'fontWeight': 'bold'}),
                dcc.Dropdown(id='y-axis', options=xy_options, value='K', clearable=False)
            ]),
        ]),
    ]),
    
    dcc.Graph(id='scatter-plot', style={'height': '50vh'}),
    
    html.Hr(),
    
    html.H3("Interactive Growth Curve", id='curve-title'),
    
    # NEW: Dynamic Time Point Control Panel
    html.Div(id='timepoint-control-panel', style={'display': 'none'}, children=[
        html.Label("Active Time Points (uncheck to exclude & refit):", style={'fontWeight': 'bold', 'color': '#2c3e50'}),
        dcc.Checklist(
            id='timepoint-selector',
            inline=True,
            style={'marginTop': '10px', 'display': 'flex', 'gap': '20px', 'flexWrap': 'wrap'}
        )
    ]),
    
    dcc.Graph(id='growth-curve-plot', style={'height': '50vh', 'display': 'none'})
])


# ==========================================
# 4. CALLBACKS
# ==========================================

@app.callback(
    Output('scatter-plot', 'figure'),
    [Input('x-axis', 'value'), Input('y-axis', 'value'),
     Input('filter-entries', 'value'), Input('filter-sites', 'value'),
     Input('filter-reps', 'value'), Input('filter-periods', 'value')]
)
def update_scatter(x_col, y_col, f_entries, f_sites, f_reps, f_periods):
    if not x_col or not y_col: raise PreventUpdate
    
    filtered_df = df_results.copy()
    if f_entries: filtered_df = filtered_df[filtered_df['entries'].isin(f_entries)]
    if f_sites: filtered_df = filtered_df[filtered_df['sites'].isin(f_sites)]
    if f_reps: filtered_df = filtered_df[filtered_df['replications'].isin(f_reps)]
    if f_periods: filtered_df = filtered_df[filtered_df['growing_periods'].isin(f_periods)]
        
    if filtered_df.empty:
        return go.Figure().update_layout(title="No data matches filters.", xaxis_visible=False, yaxis_visible=False)
    
    fig = px.scatter(
        filtered_df, x=x_col, y=y_col, 
        custom_data=['entries', 'sites', 'replications', 'growing_periods'],
        title=f"Scatterplot: {x_col} vs {y_col}"
    )
    fig.update_traces(
        hovertemplate=(
            "<b>Entry:</b> %{customdata[0]}<br><b>Site:</b> %{customdata[1]}<br>"
            "<b>Rep:</b> %{customdata[2]}<br><b>Period:</b> %{customdata[3]}<br>"
            "<b>X:</b> %{x:.3f}<br><b>Y:</b> %{y:.3f}<extra></extra>"
        ),
        marker=dict(size=10, opacity=0.7, line=dict(width=1, color='DarkSlateGrey'))
    )
    fig.update_layout(margin={'l': 40, 'b': 40, 't': 40, 'r': 0}, hovermode='closest')
    return fig


# --- STATE UPDATER: Triggered when Scatter Plot is clicked ---
@app.callback(
    [Output('active-curve-state', 'data'),
     Output('timepoint-selector', 'options'),
     Output('timepoint-selector', 'value'),
     Output('timepoint-control-panel', 'style')],
    Input('scatter-plot', 'clickData'),
    prevent_initial_call=True
)
def load_new_curve_data(clickData):
    if not clickData or 'customdata' not in clickData['points'][0]: 
        raise PreventUpdate

    custom_data = clickData['points'][0]['customdata']
    entry, site, rep, period = custom_data[0], custom_data[1], custom_data[2], custom_data[3]

    # Find available time points for this combination
    raw_condition = (
        (df_raw['entries'] == entry) & (df_raw['sites'] == site) &
        (df_raw['replications'] == rep) & (df_raw['growing_periods'] == period)
    )
    t_points = sorted(df_raw[raw_condition]['time_points'].unique())

    # Build options for the checklist (rounded for cleaner UI)
    options = [{'label': f" t={t:.1f} ", 'value': t} for t in t_points]
    
    # State payload
    state = {'entry': entry, 'site': site, 'rep': rep, 'period': period}
    
    # Reveal the control panel with styling
    panel_style = {'display': 'block', 'marginTop': '15px', 'padding': '15px', 'backgroundColor': '#e8f4f8', 'borderRadius': '5px'}

    # Return state, checklist options, default values (all checked), and panel style
    return state, options, t_points, panel_style


# --- CURVE RENDERER: Triggered by State Change OR Checklist Change ---
@app.callback(
    [Output('growth-curve-plot', 'figure'),
     Output('growth-curve-plot', 'style'),
     Output('curve-title', 'children')],
    [Input('active-curve-state', 'data'),
     Input('timepoint-selector', 'value')],
    prevent_initial_call=True
)
def draw_and_fit_curve(state, active_t_points):
    if not state.get('entry') or active_t_points is None:
        raise PreventUpdate
        
    entry, site, rep, period = state['entry'], state['site'], state['rep'], state['period']
    
    orig_params = df_results[
        (df_results['entries'] == entry) & (df_results['sites'] == site) & 
        (df_results['replications'] == rep) & (df_results['growing_periods'] == period)
    ].iloc[0]
    
    df_point_raw = df_raw[
        (df_raw['entries'] == entry) & (df_raw['sites'] == site) & 
        (df_raw['replications'] == rep) & (df_raw['growing_periods'] == period)
    ].copy()
    
    if df_point_raw.empty:
        return go.Figure(), {'height': '50vh', 'display': 'none'}, "No raw data found."

    # Separate data based on user's checklist selection
    mask = df_point_raw['time_points'].isin(active_t_points)
    df_active = df_point_raw[mask]
    df_excluded = df_point_raw[~mask]
    
    has_excluded_points = len(df_excluded) > 0
    
    # Refit the Curve using SciPy (Requires at least 6 points for 6 parameters)
    fit_success = False
    popt = [orig_params['A'], orig_params['K'], orig_params['C'], orig_params['Q'], orig_params['B'], orig_params['v']]
    
    if len(df_active) >= 3: # Minimum points to attempt fitting (ideally should be >= 6 for 6 parameters, but we can try with fewer)
        try:
            max_y = df_active['trait_value'].max() if df_active['trait_value'].max() > 0 else 1.0
            
            # Bounds: [A, K, C, Q, B, v]
            bounds = (
                [0.0, 0.0, 1e-5, 1e-5, 0.0, 1e-5],              
                [max_y, 2 * max_y, 10.0, 10.0, 10.0, 10.0]     
            )
            
            popt, _ = curve_fit(
                generalized_logistic, 
                df_active['time_points'], df_active['trait_value'],
                p0=popt, bounds=bounds, maxfev=10000
            )
            fit_success = True
        except Exception:
            pass # Fails safely, uses original parameters 
    
    # Build Figure
    fig = go.Figure()
    
    # Active Points
    if not df_active.empty:
        fig.add_trace(go.Scatter(
            x=df_active['time_points'], y=df_active['trait_value'],
            mode='markers', name='Included Data',
            marker=dict(color='black', size=10, line=dict(color='white', width=1))
        ))
    
    # Excluded Points (Shown faintly for reference)
    if not df_excluded.empty:
        fig.add_trace(go.Scatter(
            x=df_excluded['time_points'], y=df_excluded['trait_value'],
            mode='markers', name='Excluded Data',
            marker=dict(color='red', size=8, symbol='x', opacity=0.5)
        ))
        
    # Fitted Line
    max_t = df_point_raw['time_points'].max() * 1.1 
    t_smooth = np.linspace(0, max_t, 200)
    y_smooth = generalized_logistic(t_smooth, *popt)
    
    line_name = 'Refitted Model' if fit_success and has_excluded_points else 'Original Model'
    line_color = '#2ecc71' if fit_success and has_excluded_points else 'red'
    
    fig.add_trace(go.Scatter(
        x=t_smooth, y=y_smooth, mode='lines', name=line_name, hoverinfo='skip',
        line=dict(color=line_color, width=3, dash='solid' if fit_success else 'dot')
    ))
    
    # Layout with explicitly modified X-Axis ticks
    fig.update_layout(
        xaxis=dict(
            title="Time Points",
            tickmode='array',
            tickvals=active_t_points # Forces X-axis to only show ticks for active points
        ),
        yaxis_title="Biomass (y)",
        margin={'l': 40, 'b': 40, 't': 40, 'r': 0},
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01)
    )
    
    title_text = f"Growth Curve - Entry: {entry} | Site: {site} | Rep: {rep}"
    if not fit_success and has_excluded_points:
        title_text += " (Warning: Insufficient points or optimization failed)"
        
    return fig, {'height': '50vh', 'display': 'block'}, title_text

if __name__ == '__main__':
    app.run(debug=True)
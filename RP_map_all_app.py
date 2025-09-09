#suppressing a warning - UserWarning: Wswiglal-redir-stdio
import warnings
warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

import sys
sys.path.insert(0, "../")

from regular_precession_no_match import *
from systems_lib import *

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams["figure.dpi"] = 200

import math

from astropy.cosmology import FlatLambdaCDM

def redshift_to_luminosity_distance(z)->float:
    """
    Convert redshift to luminosity distance in standard Lambda CDM cosmology.
    _______________________________________________________________________________________    
    Parameters used:
    z : float : redshift
    _______________________________________________________________________________________
    Returns:    
    luminosity_distance : float : luminosity distance in Gpc
    """
    # Define the cosmology
    cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315) # Planck 2018 results
    luminosity_distance = cosmo.luminosity_distance(z).to('Gpc').value
    return luminosity_distance


def redshifted_new_params(z, rp_params)->dict:
    """
    Get redshifted params
    _______________________________________________________________________________________    
    Parameters used:
    z : float : redshift
    rp_params : dict : Regular precession parameters
    _______________________________________________________________________________________
    Returns:    
    rp_params : dict : Regular precession parameters
    """
    rp_params["dist"] = redshift_to_luminosity_distance(z)*giga_parsec
    old_chirp = rp_params["mcz"]
    #new chirp mass
    rp_params["mcz"] = old_chirp*(1+z)
    return rp_params


# Define the default parameters for the system 1
default_precession_params_sys1_RP = redshifted_new_params(0.3, default_precession_params_sys1_RP)
default_precession_params_sys1_NP = redshifted_new_params(0.3, default_precession_params_sys1_NP)

# Define the default parameters for the system 2
default_precession_params_sys2_RP = redshifted_new_params(0.3, default_precession_params_sys2_RP)
default_precession_params_sys2_NP = redshifted_new_params(0.3, default_precession_params_sys2_NP)

# Define the default parameters for the system 3
default_precession_params_sys3_RP = redshifted_new_params(0.3, default_precession_params_sys3_RP)
default_precession_params_sys3_NP = redshifted_new_params(0.3, default_precession_params_sys3_NP)



# Default parameters
default_precession_params_sys2_RP['omega_tilde'] = 2
default_precession_params_sys2_RP['theta_tilde'] = 5
default_precession_params_sys2_RP['gamma_P'] = 0

import plotly.graph_objects as go
f_min = 20
delta_f = 0.05

def plot_cos_L_phi_L(chirp_mass, mass_ratio, redshift, theta_S, phi_S, theta_J, phi_J, theta_tilde, omega_tilde, gamma_P):
    params_NP = default_precession_params_sys2_NP.copy()
    params_NP['theta_S'] = theta_S
    params_NP['phi_S'] = phi_S
    params_NP['theta_J'] = theta_J
    params_NP['phi_J'] = phi_J
    params_NP['mcz'] = chirp_mass*solar_mass
    params_NP['eta'] = mass_ratio/((1+mass_ratio)**2)
    params_NP['dist'] = redshift_to_luminosity_distance(redshift)*giga_parsec
    
    params_NP = redshifted_new_params(redshift, params_NP)
    
    params_RP = params_NP.copy()
    params_RP['theta_tilde'] = theta_tilde
    params_RP['omega_tilde'] = omega_tilde
    params_RP['gamma_P'] = gamma_P
    
    
    f_cut = Regular_precession(params_RP).get_f_cut()
    f_range = np.arange(f_min, f_cut, delta_f)
    
    # Create instance
    precession_initial_RP = Regular_precession(params_RP)

    # Calculating modulation of L around the direction of the J
    cos_L_RP = precession_initial_RP.cos_theta_L(f_range)
    
    phi_L_RP = precession_initial_RP.phi_L(f_range)
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=phi_L_RP, y=cos_L_RP, mode='lines', name='RP', line=dict(color='blue')))
    fig.add_trace(go.Scatter(x=[phi_J], y=[np.cos(theta_J)], mode='markers', marker=dict(color='red'), name='J'))
    
    fig.update_xaxes(tickvals=[np.pi/4, 3*np.pi/8, np.pi/2, 5*np.pi/8, 3*np.pi/4], ticktext=['π/4', '3π/8', 'π/2', '5π/8', '3π/4'])
    fig.update_yaxes(range=[-0.8, 0.8])

    fig.update_layout(
    xaxis_title=r'$φ_L$',
    yaxis_title=r'$cos θ_L$',
    title='Modulation of L',
    showlegend=False,  # Suppress the legend
    grid=dict()
    )
    
    return fig

def plot_amplitude(chirp_mass, mass_ratio, redshift, theta_S, phi_S, theta_J, phi_J, theta_tilde, omega_tilde, gamma_P):
    params_NP = default_precession_params_sys2_NP.copy()
    params_NP['theta_S'] = theta_S
    params_NP['phi_S'] = phi_S
    params_NP['theta_J'] = theta_J
    params_NP['phi_J'] = phi_J
    params_NP['mcz'] = chirp_mass*solar_mass
    params_NP['eta'] = mass_ratio/(1+mass_ratio)**2
    params_NP['dist'] = redshift_to_luminosity_distance(redshift)*giga_parsec
    
    params_NP = redshifted_new_params(redshift, params_NP)
    
    params_RP = params_NP.copy()
    params_RP['theta_tilde'] = theta_tilde
    params_RP['omega_tilde'] = omega_tilde
    params_RP['gamma_P'] = gamma_P

    f_cut = Regular_precession(params_RP).get_f_cut()
    f_range = np.arange(f_min, f_cut, delta_f)
    
    # Create instances
    precession_initial_NP = Regular_precession(params_NP)
    precession_initial_RP = Regular_precession(params_RP)

    # Get amplitudes
    Amp_NP = precession_initial_NP.amplitude(f_range)
    Amp_RP = precession_initial_RP.amplitude(f_range)

    
    frac_t_2_RP = (Amp_RP / Amp_NP) - 1
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(x=f_range, y=frac_t_2_RP, mode='lines', name='RP', line=dict(color='blue')))
    fig.add_hline(y=0, line=dict(color='black', dash='dash'))
    
    fig.update_xaxes(type='log', tickvals=[20, 40, 60, 80, 100, 140], ticktext=['20', '40', '60', '80', '100', '140'])
    fig.update_yaxes(range=[-1.05, 0.95])
    
    fig.update_layout(
        xaxis_title=r'f (Hz)',
        yaxis_title=r'$B_GP/B_NP - 1$',
        title='Amplitude Difference',
        showlegend=True,
        grid=dict()
    )
    
    return fig


def plot_phase(chirp_mass, mass_ratio, redshift, theta_S, phi_S, theta_J, phi_J, theta_tilde, omega_tilde, gamma_P):
    params_NP = default_precession_params_sys2_NP.copy()
    params_NP['theta_S'] = theta_S
    params_NP['phi_S'] = phi_S
    params_NP['theta_J'] = theta_J
    params_NP['phi_J'] = phi_J
    params_NP['mcz'] = chirp_mass*solar_mass
    params_NP['eta'] = mass_ratio/(1+mass_ratio)**2
    params_NP['dist'] = redshift_to_luminosity_distance(redshift)*giga_parsec
    
    params_NP = redshifted_new_params(redshift, params_NP)
    
    params_RP = params_NP.copy()
    params_RP['theta_tilde'] = theta_tilde
    params_RP['omega_tilde'] = omega_tilde
    params_RP['gamma_P'] = gamma_P

    f_cut = Regular_precession(params_RP).get_f_cut()
    f_range = np.arange(f_min, f_cut, delta_f)
    
    # Create instances
    precession_initial_RP = Regular_precession(params_RP)
    
    # Get amplitudes
    Phase_RP = precession_initial_RP.phase_phi_P(f_range) + 2*precession_initial_RP.phase_delta_phi(f_range)

    phase_diff_t_2_RP = (Phase_RP)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=f_range, y=-phase_diff_t_2_RP, mode='lines', name='RP', line=dict(color='blue')))
    fig.add_hline(y=0, line=dict(color='black', dash='dash'))
    
    fig.update_xaxes(type='log', tickvals=[20, 40, 60, 80, 100, 140], ticktext=['20', '40', '60', '80', '100', '140'])
    
    fig.update_layout(
        xaxis_title='f (Hz)',
        yaxis_title=r'$Φ_{GP} - Φ_{NP}$',
        title='Phase Difference',
        showlegend=False,
        grid=dict()
    )
    
    return fig

from plotly.subplots import make_subplots
from dash import Dash, dcc, html
from dash.dependencies import Input, Output


def plot_combined(chirp_mass, mass_ratio, redshift, theta_S, phi_S, theta_J, phi_J, theta_tilde, omega_tilde, gamma_P):
    
    cos_i_JN = (np.sin(theta_J) * np.sin(theta_S) * np.cos(phi_J - phi_S)) + (np.cos(theta_J) * np.cos(theta_S))
    
    label = (
        f'The changes in the gravitational wave signal from an inspiralling binary black hole system with regular precession shown in blue as compared to a signal with no precession.'
        f'<br>Modulation</b> of the orbital angular momentum (L) around the total angular momentum (J) of the binary in the left panel, along with <b>the fractional change in amplitude</b> due to regular precession in the middle panel,'
        f'<br>and <b>phase differences</b> between a regularly precessing and non-precessing binary in the right panel. The red dot in the left panel shows the direction of J.')

    fig = make_subplots(rows=1, cols=3, subplot_titles=(r"Modulation", "Amplitude", "Phase"), vertical_spacing=0.25)

    fig_cos_L = plot_cos_L_phi_L(chirp_mass, mass_ratio, redshift, theta_S, phi_S, theta_J, phi_J, theta_tilde, omega_tilde, gamma_P)
    fig_amp = plot_amplitude(chirp_mass, mass_ratio, redshift, theta_S, phi_S, theta_J, phi_J, theta_tilde, omega_tilde, gamma_P)
    fig_phase = plot_phase(chirp_mass, mass_ratio, redshift, theta_S, phi_S, theta_J, phi_J, theta_tilde, omega_tilde, gamma_P)
    

    # Add traces and update axes for fig_cos_L
    for trace in fig_cos_L['data']:
        fig.add_trace(trace, row=1, col=1)

        fig.update_xaxes(title_text=r'φ_L', row=1, col=1)
        fig.update_layout(showlegend=False)
        fig.update_yaxes(title_text=r'cos θ_L', title_standoff=1, row=1, col=1)
        fig.update_yaxes(range=[-0.8, 0.8], row=1, col=1)
        fig.update_xaxes(range=[np.pi/4, 3*np.pi/4], tickvals=[np.pi/4, 3*np.pi/8, np.pi/2, 5*np.pi/8, 3*np.pi/4], ticktext=['π/4', '3π/8', 'π/2', '5π/8', '3π/4'], row=1, col=1)

    # Add traces and update axes for fig_amp
    for trace in fig_amp['data']:
        fig.add_trace(trace, row=1, col=2)

        fig.update_xaxes(type='log', title_text='f (Hz)', row=1, col=2)
        fig.update_yaxes(title_text=r'[B(GP)/B(NP)] - 1', title_standoff=1, row=1, col=2)
        fig.update_layout(showlegend=True)
        fig.update_yaxes(range=[-1.05, 0.95], row=1, col=2)
        fig.update_xaxes(range=[math.log10(19), math.log10(150)], tickmode='array', tickvals=[math.log10(20), math.log10(40), math.log10(60), math.log10(80), math.log10(100), math.log10(140)], ticktext=['20', '40', '60', '80', '100', '140'], row=1, col=2)
        
    # Add traces and update axes for fig_phase
    for trace in fig_phase['data']:
        fig.add_trace(trace, row=1, col=3)

        fig.update_xaxes(type='log', title_text='f (Hz)', row=1, col=3)
        fig.update_yaxes(title_text=r'Φ(GP) - Φ(NP)', title_standoff=1, row=1, col=3)
        fig.update_layout(showlegend=False)
        fig.update_yaxes(range=[-3.25, 0.75], row=1, col=3)
        fig.update_xaxes(range=[math.log10(19), math.log10(150)], tickmode='array', tickvals=[math.log10(20), math.log10(40), math.log10(60), math.log10(80), math.log10(100), math.log10(140)], ticktext=['20', '40', '60', '80', '100', '140'], row=1, col=3)
        
    # Update layout
    fig.update_layout(
        title=None,
        margin=dict(t=45, b = 100),
        legend=dict(title='Legend', itemsizing='constant'),
        width=1350,
        height=450,
        xaxis=dict(tickmode='array', tickvals=[-np.pi, -7*np.pi/8, -3*np.pi/4, -5*np.pi/8, -np.pi/2, -3*np.pi/8, -np.pi/4, 0, np.pi/4, 3*np.pi/8, np.pi/2, 5*np.pi/8, 3*np.pi/4, 7*np.pi/8, np.pi], ticktext=['-π', '-7π/8', '-3π/4', '-5π/8', '-π/2', '-3π/8', '-π/4', '0', 'π/4', '3π/8', 'π/2', '5π/8', '3π/4', '7π/8', 'π']),
        xaxis2=dict(tickmode='array', tickvals=[20, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400], ticktext=['20', '40', '60', '80', '100', '150', '200', '250', '300', '350', '400']),
        xaxis3=dict(tickmode='array', tickvals=[20, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400], ticktext=['20', '40', '60', '80', '100', '150', '200', '250', '300', '350', '400']),
    )

    # Add annotation for the label
    fig.add_annotation(
        text=label,  # Set your caption text here
        xref="paper", yref="paper",
        x=0.5, y=-0.325,  # Position it further down to avoid cutoff
        showarrow=False,
        font=dict(
            family="Times New Roman, sans-serif",
            size=15,
            color="gray"
        ),
        align="center",
        opacity=0.95  # Semi-transparent caption
    )
    
    return fig

app = Dash(__name__, assets_folder='assets')

server = app.server

app.title = "Regular Precession"

app.layout = html.Div([
    
    html.Link(rel='icon', href='/assets/bbh_spin.jpg?v=1'),
    
    html.H3("Regular Precession", style={'margin-bottom': '1px'}),

    dcc.Graph(id='combined-plot', style={'margin-bottom': '10px'}),

    # Watermark Section
    html.Div(
        style={
            'position': 'absolute',
            'right': '60px',  # Adjust the right position to fit watermark
            'font-family': 'Times New Roman, sans-serif',
            'font-size': '15px',
            'color': 'gray',
            'padding': '1px',
            'opacity': 0.4,
            'text-align': 'center',
        },
        children="© Tamanjyot Singh @ UTD"  # Watermark text
    ),
    
    html.Div(
        [
            html.A(
                html.Img(
                    src='https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png',
                    style={
                        'width': '18px',
                        'height': '18px',
                    }
                ),
                href='https://github.com/singhtaman',
                target='_blank',
                style={
                    'font-size': '12px',
                    'color': 'gray',
                    'text-decoration': 'none',
                }
            )
        ],
        style={
            'position': 'absolute',
            'right': '35px',
            'background-color': 'rgba(255, 255, 255, 0.8)',
        }
    ),

    html.Div(
        [
            html.A(
                html.Img(
                    src='https://static.vecteezy.com/system/resources/thumbnails/003/731/316/small/web-icon-line-on-white-background-image-for-web-presentation-logo-icon-symbol-free-vector.jpg',
                    style={
                        'width': '18px',
                        'height': '18px',
                    }
                ),
                href='https://singhtaman.github.io/',
                target='_blank',
                style={
                    'font-size': '12px',
                    'color': 'gray',
                    'text-decoration': 'none',
                }
            )
        ],
        style={
            'position': 'absolute',
            'right': '18px',
            'background-color': 'rgba(255, 255, 255, 0.8)',
        }
    ),
     
    html.Div(style={'height': '20px'}),
    
    html.Div([
        html.Label(id='theta_tilde-label', children=f'Precession amplitude, θ̃ : {5.0}', style={'font-size': '15px', 'opacity': '0.9'}),
        dcc.Slider(id='theta_tilde-slider',
                   min=0, max=10, step=0.1, value=5,
                   marks={i: str(i) for i in range(11)},
                   tooltip={"placement": "bottom", "always_visible": False}),
    ], style={'width': '50%', 'display': 'inline-block'}),

    html.Div([
        html.Label(id='omega_tilde-label', children=f'Precession frequency, Ω̃ : {2.0}', style={'font-size': '15px', 'opacity': '0.9'}),
        dcc.Slider(id='omega_tilde-slider',
                   min=0, max=4.0, step=0.1, value=2.0,
                   marks={i: str(i) for i in [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]},
                   tooltip={"placement": "bottom", "always_visible": False}),
    ], style={'width': '25%', 'display': 'inline-block'}),

    html.Div([
        html.Label(id='gamma_P-label', children=r'Precessional phase, γ_P : {0:.2f}'.format(0.0), style={'font-size': '15px', 'opacity': '0.9'}),
        dcc.Slider(id='gamma_P-slider',
                   min=0, max=2*np.pi, step=np.pi/16, value=0.0,
                   marks={
                       0: '0',
                       np.pi/2: 'π/2',
                       np.pi: 'π',
                       3*np.pi/2: '3π/2',
                       2*np.pi: '2π'
                   },
                   tooltip={"placement": "bottom", "always_visible": False}),
    ], style={'width': '25%', 'display': 'inline-block'}),
    
    html.Div([
        html.Label(id='chirp_mass-label', children=f'Chirp mass, M_c : {10}', style={'font-size': '15px', 'opacity': '0.9'}),
        dcc.Slider(id='chirp_mass-slider',
                    min=5, max=50, step=1, value=10,
                    marks={i: str(i) for i in range(5, 51, 5)},
                    tooltip={"placement": "bottom", "always_visible": False}),
    ], style={'width': '50%', 'display': 'inline-block'}),
    
    html.Div([
        html.Label(id='mass_ratio-label', children=f'Mass ratio, q : {0.6}', style={'font-size': '15px', 'opacity': '0.9'}),
        dcc.Slider(id='mass_ratio-slider',
                   min=0.05, max=1, step=0.05, value=0.6,
                   marks={i: str(i) for i in [0.2, 0.4, 0.6, 0.8, 1]},
                   tooltip={"placement": "bottom", "always_visible": False}),
    ], style={'width': '25%', 'display': 'inline-block'}),
    
    html.Div([
        html.Label(id='redshift-label', children=f'Redshift, z : {0.3}', style={'font-size': '15px', 'opacity': '0.9'}),
        dcc.Slider(id='redshift-slider',
                     min=0.05, max=1, step=0.05, value=0.3,
                     marks={i: str(i) for i in [0.2, 0.4, 0.6, 0.8, 1]},
                     tooltip={"placement": "bottom", "always_visible": False}),
    ], style={'width': '25%', 'display': 'inline-block'}),
    
    html.Div([
        html.Label(id='theta_S-label', children=f'N polar angle, θ_S : {np.pi/4:.2f}', style={'font-size': '15px', 'opacity': '0.9'}),
        dcc.Slider(id='theta_S-slider',
                   min=0, max=np.pi, step=np.pi/16, value=np.pi/4,
                   marks={
                       0: '0',
                       np.pi/4: 'π/4',
                       np.pi/2: 'π/2',
                       3*np.pi/4: '3π/4',
                       np.pi: 'π'
                   },
                   tooltip={"placement": "bottom", "always_visible": False}),
    ], style={'width': '25%', 'display': 'inline-block'}),
    
    html.Div([
        html.Label(id='phi_S-label', children=f'N azimuthal angle, φ_S : {0}', style={'font-size': '15px', 'opacity': '0.9'}),
        dcc.Slider(id='phi_S-slider',
                     min=0, max=2*np.pi, step=np.pi/16, value=0,
                     marks={
                          0: '0',
                          np.pi/2: 'π/2',
                          np.pi: 'π',
                          3*np.pi/2: '3π/2',
                          2*np.pi: '2π'
                     },
                     tooltip={"placement": "bottom", "always_visible": False}),
    ], style={'width': '25%', 'display': 'inline-block'}),
    
    html.Div([
        html.Label(id='theta_J-label', children=f'J polar angle, θ_J : {np.pi/2:.2f}', style={'font-size': '15px', 'opacity': '0.9'}),
        dcc.Slider(id='theta_J-slider',
                     min=0, max=np.pi, step=np.pi/16, value=np.pi/2,
                     marks={
                          0: '0',
                          np.pi/4: 'π/4',
                          np.pi/2: 'π/2',
                          3*np.pi/4: '3π/4',
                          np.pi: 'π'
                     },
                     tooltip={"placement": "bottom", "always_visible": False}),
    ], style={'width': '25%', 'display': 'inline-block'}),
    
    html.Div([
        html.Label(id='phi_J-label', children=f'J azimuthal angle, φ_J : {np.pi/2:.2f}', style={'font-size': '15px', 'opacity': '0.9'}),
        dcc.Slider(id='phi_J-slider',
                    min=0, max=2*np.pi, step=np.pi/16, value=np.pi/2,
                    marks={
                        0: '0',
                        np.pi/2: 'π/2',
                        np.pi: 'π',
                        3*np.pi/2: '3π/2',
                        2*np.pi: '2π'
                    },
                    tooltip={"placement": "bottom", "always_visible": False}),
    ], style={'width': '25%', 'display': 'inline-block'}),
    
    html.Div(style={'height': '20px'}),
    
    html.Div([
        html.Label(id='cos_i_JN-label', children=f'Cosine of angle between J and N, cos(ι_JN) = {0:.2f}', style={'font-size': '15px', 'opacity': '0.9'})]),
    
    html.Div([
        html.Label(id='total_mass=label', children=f'Total mass of the binary is, M = {23.88:.2f} solar masses (redshifted, M_z = {31.04:.2f} solar masses), and the cut off frequency for this signal is, f_cut = {141.62:.2f} Hz', style={'font-size': '15px', 'opacity': '0.9'})])
    
])


@app.callback(
    Output('combined-plot', 'figure'),
    Input('chirp_mass-slider', 'value'),
    Input('mass_ratio-slider', 'value'),
    Input('redshift-slider', 'value'),
    Input('theta_S-slider', 'value'),
    Input('phi_S-slider', 'value'),
    Input('theta_J-slider', 'value'),
    Input('phi_J-slider', 'value'),
    Input('theta_tilde-slider', 'value'),
    Input('omega_tilde-slider', 'value'),
    Input('gamma_P-slider', 'value'),
)

def update_plot(chirp_mass, mass_ratio, redshift, theta_S, phi_S, theta_J, phi_J, theta_tilde, omega_tilde, gamma_P):
    return plot_combined(chirp_mass, mass_ratio, redshift, theta_S, phi_S, theta_J, phi_J, theta_tilde, omega_tilde, gamma_P)

@app.callback(
    Output('chirp_mass-label', 'children'),
    Input('chirp_mass-slider', 'value'))
def update_chirp_mass(value):
    return f'Chirp mass, M_c : {value:.1f}'

@app.callback(
    Output('mass_ratio-label', 'children'),
    Input('mass_ratio-slider', 'value'))
def update_mass_ratio(value):
    return f'Mass ratio, q : {value:.2f}'

@app.callback(
    Output('redshift-label', 'children'),
    Input('redshift-slider', 'value'))
def update_redshift(value):
    return f'Redshift, z : {value:.1f}'

@app.callback(
    Output('theta_S-label', 'children'),
    Input('theta_S-slider', 'value'))
def update_theta_S(value):
    return f'N polar angle, θ_S : {value:.2f}'

@app.callback(
    Output('phi_S-label', 'children'),
    Input('phi_S-slider', 'value'))
def update_phi_S(value):
    return f'N azimuthal angle, φ_S : {value:.2f}'

@app.callback(
    Output('theta_J-label', 'children'),
    Input('theta_J-slider', 'value'))
def update_theta_J(value):
    return f'J polar angle, θ_J : {value:.2f}'

@app.callback(
    Output('phi_J-label', 'children'),
    Input('phi_J-slider', 'value'))
def update_phi_J(value):
    return f'J azimuthal angle, φ_J : {value:.2f}'

@app.callback(
    Output('theta_tilde-label', 'children'),
    Input('theta_tilde-slider', 'value'))
def update_theta_tilde(value):
    return f'Precession amplitude, θ̃ : {value:.1f}'

@app.callback(
    Output('omega_tilde-label', 'children'),
    Input('omega_tilde-slider', 'value'))
def update_omega_tilde(value):
    return f'Precession frequency, Ω̃ : {value:.1f}'

@app.callback(
    Output('gamma_P-label', 'children'),
    Input('gamma_P-slider', 'value'))
def update_gamma_P(value):
    return f'Precessional phase, γ_P : {value:.2f}'

@app.callback(
    Output('cos_i_JN-label', 'children'),
    Input('theta_J-slider', 'value'),
    Input('theta_S-slider', 'value'),
    Input('phi_J-slider', 'value'),
    Input('phi_S-slider', 'value')
)
def update_cos_i_JN(theta_J, theta_S, phi_J, phi_S):
    cos_i_JN = (np.sin(theta_J) * np.sin(theta_S) * np.cos(phi_J - phi_S)) + (np.cos(theta_J) * np.cos(theta_S))
    return f'Cosine of angle between J and N, cos(ι_JN) : {cos_i_JN:.2f}'

@app.callback(
    Output('total_mass=label', 'children'),
    Input('chirp_mass-slider', 'value'),
    Input('mass_ratio-slider', 'value'),
    Input('redshift-slider', 'value'))
def update_total_mass(chirp_mass, mass_ratio, redshift):
    total_mass = chirp_mass/((mass_ratio/(1+mass_ratio)**2)**(3/5))
    total_mass_redshifted = total_mass*(1+redshift)
    f_cut = 1/(6**(3/2)*np.pi*total_mass_redshifted*solar_mass)
    return f'Total mass of the binary is, M = {total_mass:.2f} solar masses (redshifted, M_z = {total_mass_redshifted:.2f} solar masses), and the cut off frequency for this signal is, f_cut = {f_cut:.2f} Hz'

import os
if __name__ == '__main__':
    app.run_server(debug=True)  # defaults to http://127.0.0.1:8050
    # Get the port from the environment variable
    #port = int(os.environ.get('PORT', 5000))
    #app.run_server(host='0.0.0.0', port=port)
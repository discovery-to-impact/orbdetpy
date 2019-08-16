import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import os
import sys
from examples import astrotools
from examples import plotResults
import time
import json

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Get filenames of all files in data folder
dir_path    = os.path.dirname(os.path.realpath(__file__))
dataPath    = 'examples/data/'
fnDir       = os.listdir(dir_path+'/'+dataPath)
fnDirSimCfg = [fn for fn in fnDir if 'sim_cfg' in fn]
fnDirOdCfg  = [fn for fn in fnDir if 'od_cfg' in fn]
fnDirObs    = [fn for fn in fnDir if 'obs_data' in fn]
fnDirOd     = [fn for fn in fnDir if 'od_out' in fn]

# Begin UI design
app.layout = html.Div([
    
    # Top banner
    html.Div([
        html.H1("ORBDETPY", 
            style = {
                'color' : 'White',
                'margin-left' : '20px'
            }
        ),
    ],style = {'backgroundColor' : 'rgb(51, 153, 255)'}),

    # Everything below banner (create a margin on the left of screen)
    html.Div([

        # About, About info
        html.Div([
            # About
            html.H5(['ABOUT'], style = {'margin-bottom' : '5px'}),
            # About info
            html.P("""Welcome to the orbdetpy user interface!
            You can run various routines from here that access the files in the data folder.
            To begin, select the input files required for each of the processes below. Each 
            process has a short description and the input files it needs to run
            and the output files it creates or alters. After selecting the input files and 
            processes you want to execute, press submit. After the process is complete, press 
            CLEAR before running again. If at any point callback errors occur, refresh the page
            and try again. If the callback errors persist, make sure to see if the input files
            follow the correct format outlined in docs/filefmt.rst""")
        ], style = {'margin-bottom' : '1px', 'border' : '0px solid rgb(51, 153, 255)', 'padding' : '10px'}),
        
        # Process Selection, input dropdown, range slider, start end time, checklist, input bias, submit text, buttons
        html.Div([
            # Process Selection
                html.H5(['PROCESS SELECTION'], style = {'margin-bottom' : '10px'}),
            # Input dropdown
            dcc.Dropdown(
                id = 'input-sim_cfg',
                options = [{'label' : i,'value' : i} for i in fnDirSimCfg],
                value = '',
                style = {'margin-bottom' : "10px"},
                multi = False,
                placeholder="Select Simulation Configure File"
            ),
            dcc.Dropdown(
                id = 'input-od_cfg',
                options = [{'label' : i,'value' : i} for i in fnDirOdCfg],
                value = '',
                style = {'margin-bottom' : "10px"},
                multi = False,
                placeholder="Select OD Configure File"
            ),
            dcc.Dropdown(
                id = 'input-obs_data',
                options = [{'label' : i,'value' : i} for i in fnDirObs],
                value = '',
                style = {'margin-bottom' : "10px"},
                multi = False,
                placeholder="Select Observation File"
            ),
            dcc.Dropdown(
                id = 'input-od_out',
                options = [{'label' : i,'value' : i} for i in fnDirOd],
                value = '',
                style = {'margin-bottom' : "10px"},
                multi = False,
                placeholder="Select OD Output File"
            ),
            html.Div(id = 'start-display', style = {'margin-bottom' : '5px'}),
            html.Div(id = 'end-display', style = {'margin-bottom' : '10px'}),
            # Checklist
            dcc.Checklist(
                id = 'runChecks',
                options = [
                    {'label' : 'Run simulation (in-sim_cfg out-obs_data)', 'value': 'runSim'},                    
                    {'label' : 'Transform all data TEME to J2000 (in-obs_data,od_cfg out-obs_data,od_cfg)', 'value': 'transformTEME2J2000'},
                    {'label' : 'Add time bias, enter manually below or estimate (in-obs_data,od_cfg out-obs_data,od_cfg)', 'value': 'addTimeBias'},
                    {'label' : 'Estimate and correct observations for stellar aberration (in-obs_data out-obs_data)', 'value': 'correctAberr'},
                    {'label' : 'Estimate for time bias (in-obs_data,od_cfg out-obs_data,od_cfg)', 'value': 'estTimeBias'},
                    {'label' : 'Create plot of residuals between raw obs and reference obs (in-obs_data,od_cfg out-plots)', 'value': 'plotEstSensorError'},
                    {'label' : 'Run orbit determination (in-od_cfg,obs_data out-od_out)', 'value': 'runOD'},
                    {'label' : 'Create plots for orbit determination analysis (in-od_cfg,obs_data,od_out out-plots)', 'value': 'runPlotResults'},
                    {'label' : 'Run smoother (must run with above)', 'value': 'runSmoother'}
                ],
                value=[],
                labelStyle={'display': 'block'},
                style = {
                    'margin-bottom' : "10px"
                },
                inputStyle={"margin-right": "10px"},
            ),
            # Manually add time bias input
            dcc.Input(
                id = 'add-tbias',
                placeholder='Manually include time bias: Defualt = 0',
                type='number',
                value='',
                style = {
                    'margin-bottom' : "20px",
                    'margin-left' : "5px"
                },
            ),
            # Create buttons to run and reset processes
            html.Div(id = 'progress'),
            html.Div([
                html.Button('Submit', id='submit-button'),
                html.Button('Clear', id='clear-button', style = {'margin-left':'10px'})
            ], style = {'display' : 'inline-block'}),
        ], style = {'display' : 'inline-block', 'margin-right' : '30px', 'border' : '0px solid rgb(51, 153, 255)', 'padding' : '10px'}),

        html.Div([
            html.H5(['DATA CONFIGURATION'], style = {'margin-bottom' : '10px'}),
            html.P("""To change file contents, select a file from the dropdown 
            and edit contents. Content is automatically saved when altered"""),
            dcc.Dropdown(
                id = 'configure-dropdown-list',
                options = [{'label' : i,'value' : i} for i in fnDir],
                style = {
                    'margin-top' : '10px'
                },
            ),
            dcc.Textarea(
                id = 'view-file',
                style = {
                    'margin-top' : '20px',
                    'height' : '490px',
                    'width' : '100%',
                },
            ),
        ], style = {"margin-top" : "25px", 'display' : 'inline-block', 'border' : '0px solid rgb(51, 153, 255)', 'padding' : '10px'})

    ], style = {"margin-left" : "25px", "margin-right" : "25px"}), # Everything below banner

    html.Div(id='nothing', style={'display':'none'})

], style = {'margin' : '0', 'backgroundColor' : 'rgb(250,250,250'}) # entire app layout

# Begin callbacks

# Get start and end time of given file
@app.callback(
    [Output(component_id='start-display', component_property='children'),
     Output(component_id='end-display', component_property='children')],
    [Input(component_id='input-od_cfg', component_property='value')])
def update_output(od_cfg):
    
    if not od_cfg:
        od_cfg = ''
        startTime = 'Select sim_cfg or od_cfg file'
        endTime   = 'Select sim_cfg or od_cfg file'
    else:
        with open(dataPath+od_cfg, "r") as fp:
            inp = json.load(fp)
        startTime = inp['Propagation']['Start']
        endTime   = inp['Propagation']['End']

    return 'Start time = ' + startTime, 'End time = ' + endTime

# Load selected file into text area
@app.callback(
    Output(component_id='view-file', component_property='value'),
    [Input(component_id='configure-dropdown-list', component_property='value')]
)
def update_output_div(fn):

    if not fn:
     file_contents = ''
    else:
        file_contents = open(dataPath+fn, "r")
        file_contents = file_contents.read()

    return file_contents

# Save any changes made to file contents
@app.callback(
    Output(component_id='nothing', component_property = 'children'),
    [Input(component_id='configure-dropdown-list', component_property='value'),
     Input(component_id='view-file', component_property='value')]
)
def update_output_div(fn, fileData):
    if not fn:
     file_contents = ''
    else:
        file_contents = open(dataPath+fn, "w")
        file_contents.write(fileData)
    return

# Update dropdown lists for newly created files
@app.callback(
    [Output(component_id='input-sim_cfg', component_property='options'),
     Output(component_id='input-od_cfg', component_property='options'),
     Output(component_id='input-obs_data', component_property='options'),
     Output(component_id='input-od_out', component_property='options'),
     Output(component_id='configure-dropdown-list', component_property='options')],
    [Input(component_id='submit-button', component_property='n_clicks')]
)
def update_output_div(n_clicks):

    dir_path    = os.path.dirname(os.path.realpath(__file__))
    fnDir       = os.listdir(dir_path+'/'+dataPath)
    fnDirSimCfg = [fn for fn in fnDir if 'sim_cfg' in fn]
    fnDirOdCfg  = [fn for fn in fnDir if 'od_cfg' in fn]
    fnDirObs    = [fn for fn in fnDir if 'obs_data' in fn]
    fnDirOd     = [fn for fn in fnDir if 'od_out' in fn]
    options1 = [{'label' : i,'value' : i} for i in fnDirSimCfg]
    options2 = [{'label' : i,'value' : i} for i in fnDirOdCfg]
    options3 = [{'label' : i,'value' : i} for i in fnDirObs]
    options4 = [{'label' : i,'value' : i} for i in fnDirOd]
    options5 = [{'label' : i,'value' : i} for i in fnDir]

    return options1, options2, options3, options4, options5

# Execute chosen processes from checklist
@app.callback(
    Output(component_id='progress', component_property='children'),
    [Input(component_id='submit-button', component_property='n_clicks'),
     Input(component_id='input-sim_cfg', component_property='value'),
     Input(component_id='input-od_cfg', component_property='value'),
     Input(component_id='input-obs_data', component_property='value'),
     Input(component_id='input-od_out', component_property='value'),
     Input(component_id='add-tbias', component_property='value'),
     Input(component_id='runChecks', component_property='value')]
)
def update_output_div(n_clicks,sim_cfg,obs_data,od_cfg,od_out,tbias,runChecks):

    dataPath = ' examples/data/'

    if not sim_cfg:
        sim_cfg = ''
    if not obs_data:
        obs_data = ''
    if not od_cfg:
        od_cfg = ''
    if not od_out:
        od_out = ''

    if 'sim_cfg' in sim_cfg:
        idx = sim_cfg.find("sim_cfg.json")
        mT  = sim_cfg[:idx]
    if 'obs_data' in obs_data:
        idx = obs_data.find("obs_data.json")
        mT  = obs_data[:idx]
    if 'od_cfg' in od_cfg:
        idx = od_cfg.find("od_cfg.json")
        mT  = od_cfg[:idx]
    if 'od_out' in od_out:
        idx = od_out.find("od_out.json")
        mT  = od_out[:idx]
 
    if tbias == '':
        tbias = 0
    
    # EXECUTE DESIRED PROCESSES
    runGo = 'Press Submit When Ready: '
    if n_clicks == 1:
        runGo = 'Complete! Press CLEAR before next run'
        # Simulate Measurements
        if 'runSim' in runChecks:
            sim_cfg  = dataPath+mT+'sim_cfg.json'
            obs_data = dataPath+mT+'obs_data.json'
            os.system('python examples/testsim.py' + sim_cfg + obs_data)
        # Transform TEME data into J2000
        if 'transformTEME2J2000' in runChecks:
            obs_data = dataPath+mT+'obs_data.json'
            od_cfg   = dataPath+mT+'od_cfg.json'
            astrotools.transformDataTEME2J2000(obs_data[1:], od_cfg[1:])
        # Apply provided time bias 
        if 'addTimeBias' in runChecks:
            obs_data = dataPath+mT+'obs_data.json'
            od_cfg   = dataPath+mT+'od_cfg.json'
            astrotools.addTimeBias(obs_data[1:], od_cfg[1:], tbias)
        # Estimate and correct given obs for stellar aberration
        if 'correctAberr' in runChecks:
            obs_data = dataPath+mT+'obs_data.json'
            astrotools.aberrationCorrection(obs_data[1:], 1)
        # Estimate time bias given obs and reference data
        if 'estTimeBias' in runChecks:
            obs_data = dataPath+mT+'obs_data.json'
            od_cfg   = dataPath+mT+'od_cfg.json'
            a0 = astrotools.estTimeBiasRADEC(obs_data[1:], od_cfg[1:])
            if 'addTimeBias' in runChecks:
                astrotools.addTimeBias(obs_data[1:], od_cfg[1:], a0)
            print(a0)
        # Plot residuals of given obs and reference obs
        if 'plotEstSensorError' in runChecks:
            sim_cfg  = dataPath+mT+'sim_cfg.json'
            od_cfg   = dataPath+mT+'od_cfg.json'
            obs_data = dataPath+mT+'obs_data.json'
            astrotools.plotEstSensorError(od_cfg[1:], obs_data[1:])
        # Plot Simulation Results
        if 'runPlotSim' in runChecks:
            obs_data = dataPath+mT+'obs_data.json'
            os.system('python examples/plotsim.py' + sim_cfg + obs_data)
        # Run Orbit Determination
        if 'runOD' in runChecks:
            od_cfg   = dataPath+mT+'od_cfg.json'
            obs_data = dataPath+mT+'obs_data.json'
            od_out   = dataPath+mT+'od_out.json'
            os.system('python examples/testodet.py' + od_cfg + obs_data + od_out)
        # Plot Filtered Results
        if 'runPlotOD' in runChecks:
            os.system('python examples/plotodet.py ' + od_cfg + obs_data + od_out)
        # Plot All OD Results
        if 'runPlotResults' in runChecks:
            od_cfg   = dataPath+mT+'od_cfg.json'
            obs_data = dataPath+mT+'obs_data.json'
            od_out   = dataPath+mT+'od_out.json'
            if 'runSmoother' in runChecks:
                runS = 1
            else:
                runS = 0
            plotResults.plot(obs_data[1:], od_cfg[1:], od_out[1:], runS, 1)

    return runGo

# Reset checklist, resest number of clicks, clear run info
@app.callback(
    [Output(component_id='runChecks', component_property = 'value'),
     Output(component_id='submit-button', component_property = 'n_clicks'),
     Output(component_id='add-tbias', component_property = 'value')],
    [Input(component_id='clear-button', component_property='n_clicks')]
)
def update_output_div(clear):
    return [], 0, ''

if __name__ == '__main__':
    app.run_server(port = 3005, debug=True)

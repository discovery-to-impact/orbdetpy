import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import os
import astrotools
import plotResults
import time
import json

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

dir_path = os.path.dirname(os.path.realpath(__file__))
fnDir = os.listdir(dir_path+'/data')

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
                id = 'input-dropdown-list',
                options = [{'label' : i,'value' : i} for i in fnDir],
                value = [''],
                style = {'margin-bottom' : "10px", 'width' :'95%'},
                multi = True,
                placeholder="Select input files"
            ),
            # Range slider
            html.Div([
                dcc.RangeSlider(
                    id = 'date-rangeslider',
                    min = 0,
                    max = 100,
                    value = [0,100]
                ),
            ], style = {'margin-left' : '10px','margin-bottom' : '10px', 'width' : '95%'}),
            html.Div(id = 'start-display', style = {'margin-bottom' : '5px'}),
            html.Div(id = 'end-display', style = {'margin-bottom' : '10px'}),
            # Checklist
            dcc.Checklist(
                id = 'runChecks',
                options = [
                    {'label' : 'Run simulation (in-sim_cfg out-obs_data)', 'value': 'runSim'},
                    {'label' : 'Parse astro data and reference data (in-astro,ref out-obs_data,od_cfg)', 'value': 'parseAstroAndRef'},
                    {'label' : 'Transform all data TEME to J2000 (in-obs_data,od_cfg out-obs_data,od_cfg)', 'value': 'transformTEME2J2000'},
                    {'label' : 'Add time bias, enter manually below or estimate (in-obs_data,od_cfg out-obs_data,od_cfg)', 'value': 'addTimeBias'},
                    {'label' : 'Estimate and correct observations for stellar aberration (in-obs_data out-obs_data)', 'value': 'correctAberr'},
                    {'label' : 'Estimate for time bias (in-obs_data,od_cfg out-obs_data,od_cfg)', 'value': 'estTimeBias'},
                    {'label' : 'Create plot of residuals between raw obs and reference obs (in-obs_data,od_cfg out-plots)', 'value': 'plotEstSensorError'},
                    {'label' : 'Run orbit determination (in-od_cfg,obs_data out-od_out)', 'value': 'runOD'},
                    {'label' : 'Run smoother (run with below)', 'value': 'runSmoother'},
                    {'label' : 'Create plots for orbit determination analysis (in-od_cfg,obs_data,od_out out-plots)', 'value': 'runPlotResults'}
                ],
                value=[],
                labelStyle={'display': 'block'},
                style = {
                    'margin-bottom' : "10px",
                    'width' : '95%',
                },
                inputStyle={"margin-right": "10px", 'width' : '20px'},
            ),
            # Manually add time bias input
            dcc.Input(
                id = 'add-tbias',
                placeholder='Manually include time bias: Defualt = 0',
                type='number',
                value='',
                style = {
                    'margin-bottom' : "20px",
                    'margin-left' : "5px",
                    'width' : '400px'
                },
            ),
            # Create buttons to run and reset processes
            html.Div(id = 'progress'),
            html.Div([
                html.Button('Submit', id='submit-button'),
                html.Button('Clear', id='clear-button', style = {'margin-left':'10px'})
            ], style = {'display' : 'inline-block'}),
        ], style = {'display' : 'inline-block', 'margin-right' : '30px', 'width' : '40%', 'border' : '0px solid rgb(51, 153, 255)', 'padding' : '10px'}),

        html.Div([
            html.H5(['DATA CONFIGURATION'], style = {'margin-bottom' : '10px'}),
            html.P("""To change file contents, select a file from the dropdown 
            and edit contents. Content is automatically saved when altered""", style = {'width' : '95%'}),
            dcc.Dropdown(
                id = 'configure-dropdown-list',
                options = [{'label' : i,'value' : i} for i in fnDir],
                style = {
                    'width': '95%',
                    'margin-top' : '10px',
                },
            ),
            dcc.Textarea(
                id = 'view-file',
                style = {
                    'height' : '350px',
                    'width': '95%',
                    'margin-top' : '20px',
                },
            ),
        ], style = {"margin-top" : "25px", 'display' : 'inline-block', 'width' : '40%', 'border' : '0px solid rgb(51, 153, 255)', 'padding' : '10px'})

    ], style = {"margin-left" : "25px", "margin-right" : "25px"}), # Everything below banner

    html.Div(id='nothing', style={'display':'none'})

], style = {'margin' : '0', 'backgroundColor' : 'rgb(250,250,250'}) # entire app layout

# Get start and end time of given file
@app.callback(
    [Output(component_id='start-display', component_property='children'),
    Output(component_id='end-display', component_property='children')],
    [Input(component_id='input-dropdown-list', component_property='value')])
def update_output(inputList):

    startTime = 'Select sim_cfg or od_cfg file'
    endTime = 'Select sim_cfg or od_cfg file'

    # FILENAMES
    for fn in inputList:
        if 'sim_cfg' in fn:
            sim_cfg = 'data/'+fn 
            with open(sim_cfg, "r") as fp:
                inp = json.load(fp)
            startTime = inp['Propagation']['Start']
            endTime = inp['Propagation']['End']
        elif 'od_cfg' in fn:
            od_cfg = 'data/'+fn
            with open(od_cfg, "r") as fp:
                inp = json.load(fp)
            startTime = inp['Propagation']['Start']
            endTime = inp['Propagation']['End']

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
        file_contents = open('data/'+fn, "r")
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
        file_contents = open('data/'+fn, "w")
        file_contents.write(fileData)
    return

# Update dropdown lists for newly created files
@app.callback(
    [Output(component_id='configure-dropdown-list', component_property = 'options'),
    Output(component_id='input-dropdown-list', component_property = 'options')],
    [Input(component_id='submit-button', component_property='n_clicks')]
)
def update_output_div(n_clicks):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    fnDir = os.listdir(dir_path+'/data')
    options = [{'label' : i,'value' : i} for i in fnDir]
    return options, options

# Execute chosen processes from checklist
@app.callback(
    Output(component_id='progress', component_property='children'),
    [Input(component_id='submit-button', component_property='n_clicks'),
     Input(component_id='input-dropdown-list', component_property='value'),
     Input(component_id='add-tbias', component_property='value'),
     Input(component_id='runChecks', component_property='value')]
)
def update_output_div(n_clicks,inputList,tbias,runChecks):

    sim_cfg = ''
    obs_data = ''
    od_cfg = ''
    od_out = ''

    # FILENAMES
    for fn in inputList:
        if 'sim_cfg' in fn:
            sim_cfg = ' data/'+fn 
            mT = fn[:-13]
        elif 'obs_data' in fn:
            obs_data = ' data/'+fn
            mT = fn[:-14]
        elif 'od_cfg' in fn:
            od_cfg = ' data/'+fn
            mT = fn[:-12]
        elif 'od_out' in fn:
            od_out = ' data/'+fn
            mT = fn[:-12]
        elif 'RefsatData' in fn:
            ref_data = 'data/'+fn
        else:
            astro_data = 'data/'+fn 
            sensorID = fn[:4]
            day = fn[-12:-10]
            mT = sensorID + '_' + day

    if tbias == '':
        tbias = 0
    
    # EXECUTE DESIRED PROCESSES
    runGo = 'Press Submit When Ready: '
    if n_clicks == 1:
        runGo = 'Complete! Press CLEAR before next run'
        # Simulate Measurements
        if 'runSim' in runChecks:
            sim_cfg = ' data/'+mT+'_sim_cfg.json'
            obs_data = ' data/'+mT+'_obs_data.json'
            os.system('python testsim.py' + sim_cfg + obs_data)
        # Map astro format mat files to od_cfg file and obs_data file, apply time bias
        if 'parseAstroAndRef' in runChecks:
            obs_data = ' data/'+mT+'_obs_data.json'
            od_cfg = ' data/'+mT+'_od_cfg.json'
            astrotools.parseAstroAndRefData(
                astro_data, ref_data, obs_data[1:], od_cfg[1:], 0)
        if 'transformTEME2J2000' in runChecks:
            obs_data = ' data/'+mT+'_obs_data.json'
            od_cfg = ' data/'+mT+'_od_cfg.json'
            astrotools.transformDataTEME2J2000(obs_data[1:], od_cfg[1:])
        # Apply provided time bias 
        if 'addTimeBias' in runChecks:
            obs_data = ' data/'+mT+'_obs_data.json'
            od_cfg = ' data/'+mT+'_od_cfg.json'
            astrotools.addTimeBias(obs_data[1:], od_cfg[1:], tbias)
        # Estimate and correct given obs for stellar aberration
        if 'correctAberr' in runChecks:
            obs_data = ' data/'+mT+'_obs_data.json'
            astrotools.aberrationCorrection(obs_data[1:], 1)
        # Estimate time bias given obs and reference data
        if 'estTimeBias' in runChecks:
            obs_data = ' data/'+mT+'_obs_data.json'
            od_cfg = ' data/'+mT+'_od_cfg.json'
            a0 = astrotools.estTimeBiasRADEC(obs_data[1:], od_cfg[1:])
            if 'addTimeBias' in runChecks:
                astrotools.addTimeBias(obs_data[1:], od_cfg[1:], a0)
            print(a0)
        # Plot residuals of given obs and reference obs
        if 'plotEstSensorError' in runChecks:
            od_cfg = ' data/'+mT+'_od_cfg.json'
            obs_data = ' data/'+mT+'_obs_data.json'
            astrotools.plotEstSensorError(od_cfg[1:], obs_data[1:])
        # Plot Simulation Results
        if 'runPlotSim' in runChecks:
            os.system('python plotsim.py' + sim_cfg + obs_data)
        # Run Orbit Determination
        if 'runOD' in runChecks:
            od_cfg = ' data/'+mT+'_od_cfg.json'
            obs_data = ' data/'+mT+'_obs_data.json'
            od_out = ' data/'+mT+'_od_out.json'
            os.system('python testodet.py' + od_cfg + obs_data + od_out)
        # Plot Filtered Results
        if 'runPlotOD' in runChecks:
            os.system('python plotodet.py' + od_cfg + obs_data + od_out)
        # Plot All OD Results
        if 'runPlotResults' in runChecks:
            od_cfg = ' data/'+mT+'_od_cfg.json'
            obs_data = ' data/'+mT+'_obs_data.json'
            od_out = ' data/'+mT+'_od_out.json'
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
    app.run_server(debug=True)

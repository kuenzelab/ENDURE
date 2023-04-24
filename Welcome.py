import os
import sys
from functools import partial

import streamlit as st
from PIL import Image
from yaml import SafeLoader
import yaml

from st_pages import show_pages_from_config
from utility import add_logo, load_text

import warnings
warnings.filterwarnings("ignore")

# Constants
LOCAL_PATH = 'lib/rosetta_linux/source/bin/residue_energy_breakdown.static.linuxgccrelease'
RS_LOCAL_PATH = 'lib/rosetta_linux/source/bin/rosetta_scripts.static.linuxgccrelease'


# --- Initialization and Session State Functions ---

def clear_session() -> None:
    """
    Clear the Streamlit Session State
    :return: None
    """
    for key in st.session_state.keys():
        del st.session_state[key]


def ensure_state():
    sections = [
        "Home",
        "File Upload",
        "Interaction Analysis",
        "Residue Depth",
        "Energy Heatmap",
    ]
    sys.path.append(os.path.dirname(__file__)) #
    if 'root_dir' not in st.session_state.keys():# Save to the Session State the root directory of the project, remove the name of the file
        st.session_state['root_dir'] = os.path.dirname(__file__.replace('Welcome.py', ''))
        st.session_state['root_dir'] = '/app/ENDURE'
    for section in sections:
        if section not in st.session_state:
            st.session_state[section] = {}


# --- Rosetta Detection Functions ---

def check_local_rosetta() -> None:
    """
    Check if rosetta is included as part of the webserver
    :return:
    """
    exists = os.path.exists(LOCAL_PATH)
    rs_exists = os.path.exists(RS_LOCAL_PATH)
    if STATE:
        if 'rosetta_installed' not in STATE.keys():
            STATE['rosetta_installed'] = False
            STATE['rosetta_local'] = False
            st.session_state['Home']['rosetta_installed'] = False
            st.session_state['Home']['rosetta_local'] = False
            
        STATE['rosetta_local'] = exists
        
        if exists:
            STATE['rosetta_path'] = LOCAL_PATH
            STATE['rosetta_scripts_path'] = RS_LOCAL_PATH
            STATE['rosetta_installed'] = True
            st.session_state['rosetta_path'] = LOCAL_PATH
            st.session_state['rosetta_scripts_path'] = RS_LOCAL_PATH
            st.session_state['Home']['rosetta_installed'] = True
    else:
        if 'rosetta_installed' not in st.session_state.keys():
            st.session_state['rosetta_installed'] = False
            st.session_state['rosetta_local'] = False
        st.session_state['rosetta_local'] = exists
        if exists:
            st.session_state['rosetta_path'] = LOCAL_PATH
            st.session_state['rosetta_installed'] = True
        if rs_exists:
            st.session_state['rosetta_scripts_path'] = RS_LOCAL_PATH
            st.session_state['rosetta_installed'] = True

def check_user_rosetta(path: str) -> bool:
    """
    Validate the user-provided rosetta path
    :param path:
        The user-provided rosetta path
    :return:
    """
    valid_path = os.path.exists(path)
    STATE['rosetta_installed'] = \
        valid_path and 'residue_energy_breakdown' in path
    return valid_path and 'residue_energy_breakdown' in path


def path_input(container) -> None:
    """
    Callback function to dynamically update the status widget without having
    to wait for a page refresh.
    :param container:
        The container to write the status symbol to
    :return:
    """
    STATE['rosetta_path'] = st.session_state['rosetta_path']
    if check_user_rosetta(STATE['rosetta_path']):
        container.success('Successfully Found the Provided Executable')
    else:
        container.error('Unable to find provided filepath')


def detect_rosetta() -> None:
    """
    Ensure that the application knows where to find the Rosetta executable
    :return:
    """
    if STATE:
        rosetta_local = STATE.get('rosetta_local', False)
    else:
        rosetta_local = st.session_state.get('rosetta_local', False)

    if rosetta_local:
        st.success('Local Rosetta Installation Detected')
    else:
        status = st.container()
        if st.session_state.get('rosetta_installed', False):
            status.success('Successfully Found the Provided Executable')
        else:
            status.warning('Please Enter the Executable Path')
        st.text_input(
            label='',
            value=st.session_state.get('rosetta_path', 'main/source/bin/residue_energy_breakdown.static.linuxgccrelease'),
            key='rosetta_path',
            on_change=partial(path_input, status)
        )

# --- Status Display Functions ---

def file_status(
    name: str,
    error: str,
    success: str,
    warning: str
) -> None:
    """
    Creates a file status icon in the sidebar
    :param name:
        Name of the file stored in streamlit session state
    :param error:
        Message to display if file is not yet generated
    :param success:
        Message to display if file is up-to-date
    :param warning:
        Message to display if file is outdated
    :return: None
    """
    if name not in st.session_state['File Upload'].keys():
        st.error(error)
    elif st.session_state['File Upload'][name]:
        st.success(success)
    else:
        st.warning(warning)

# ---- Markdown Display Functions ---

def load_welcome_text() -> str:
    """
    Load the welcome text from the welcome.md file
    :return:
        The welcome text
    """ 
    
    with open("text/welcome/welcome.md", "r") as file:
        return file.read()

# --- Main Application Functions ---

def home() -> None:
    """
    Creates the Homepage Screen
    :return: None
    """
    left, center, right = st.columns([1, 2, 1])
    check_local_rosetta()
    with center:
        st.markdown(load_welcome_text())



# --- Main Application Execution ---

STATUS = {
    'cleaned': dict(
        error='PDB Files not Cleaned',
        success='PDB Files Cleaned',
        warning='PDB Files Changed, should re-clean'
    ),
    'mut_calc': dict(
        error='Mutations Not Calculated',
        success='Mutations Calculated',
        warning='PDB Files Changed, should re-calculate'
    ),
    'depth': dict(
        error='Residue Depth Not Calculated',
        success='Reside Depth Calculated',
        warning='PDB Files Changed, should re-calculate'
    ),
    'breakdown': dict(
        error='Energy Breakdown Not Calculated',
        success='Energy Breakdown Calculated',
        warning='PDB Files Changed, should re-calculate'
    )
}


st.set_page_config(
    page_title="Hello",
    page_icon="👋",
    layout="wide",
    initial_sidebar_state="expanded",
)

ensure_state()

add_logo("images/draft_logo_200.png")

with st.sidebar:
    st.title("ENDURE")
    STATE = st.session_state['Home']
    check_local_rosetta()
    detect_rosetta()
    for key, value in STATUS.items():
        file_status(name=key, **value)

show_pages_from_config(".streamlit/pages_sections.toml")
home()  
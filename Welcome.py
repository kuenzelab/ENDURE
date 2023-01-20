from pages.testing import c_Residue_Depth
import streamlit as st
import streamlit_authenticator as stauth
# Import yaml
import yaml
# Import SafeLoader
from yaml import SafeLoader
from functools import partial
import streamlit as st
import sys
import os
from PIL import Image
from pages import (
    a_File_Upload,   
    b_Interaction_Analysis, 
    d_Energy_Heatmap, 
    e_Structure_View,
    c_Residue_Depth_plotly#,Mutations
)
from utility import load_text
sys.path.append(os.path.dirname(__file__))

STATE: dict
LOCAL_PATH = 'lib/rosetta_linux/source/bin/residue_energy_breakdown' \
             '.static.linuxgccrelease'
RS_LOCAL_PATH = 'lib/rosetta_linux/source/bin/rosetta_scripts' \
             '.static.linuxgccrelease'

def clear_session() -> None:
    """
    Clear the Streamlit Session State
    :return: None
    """
    for key in st.session_state.keys():
        del st.session_state[key]


@st.cache
def load_logo() -> Image:
    """
    Load the Lab Logo to Display in the Sidebar
    :return: Image Object
    """
    logo = Image.open('images/lab_logo.png')
    return logo


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


def ensure_state() -> None:
    """
    Create Session State Dictionaries for Each Page
    :return:
    """
    for i in PAGES.keys():
        if i not in st.session_state.keys():
            st.session_state[i] = {}


def check_local_rosetta() -> None:
    """
    Check if rosetta is included as part of the webserver
    :return:
    """
    exists = os.path.exists(LOCAL_PATH)
    # If STATE is not empty, then we are on the homepage
    if STATE:
        if 'rosetta_installed' not in STATE.keys():
            STATE['rosetta_installed'] = False
            STATE['rosetta_local'] = False
            st.session_state['Home']['rosetta_installed'] = False
            st.session_state['Home']['rosetta_local'] = False
        STATE['rosetta_local'] = exists
        if exists:
            STATE['rosetta_path'] = LOCAL_PATH
            STATE['rosetta_installed'] = True
            st.session_state['rosetta_path'] = LOCAL_PATH
            st.session_state['Home']['rosetta_installed'] = True

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
        if STATE['rosetta_local']:
            st.success('Local Rosetta Installation Detected')
        else:
            status = st.container()
            if STATE['rosetta_installed']:
                status.success('Successfully Found the Provided Executable')
            else:
                status.warning('Please Enter the Executable Path')
            st.text_input(
                label='',
                value=STATE['rosetta_path'],
                key='rosetta_path',
                on_change=partial(path_input, status)
            )


def sidebar_title() -> None:
    st.markdown(
        body="<h1 style='text-align: center;'>Kuenze Lab</h1>",
        unsafe_allow_html=True
    )
    logo = load_logo()
    st.image(logo, use_column_width=True)


def home() -> None:
    """
    Creates the Homepage Screen
    :return: None
    """
    left, center, right = st.columns([1, 2, 1])
    check_local_rosetta()
    with center:
        st.title('Energetic Analysis Tools')
        st.header('Introduction')
        st.write(load_text('home', 'introduction'))
        st.subheader('Rosetta Requirement')
        st.write(load_text('home', 'rosetta'))
        if 'rosetta_path' not in STATE.keys():
            STATE['rosetta_path'] = 'main/source/bin/residue_energy_' \
                                    'breakdown.static.linuxgccrelease'
        detect_rosetta()


PAGES = {
    'Home': home,
    'File Upload': a_File_Upload.main,
    'Interaction Analysis': b_Interaction_Analysis.main,
    'Residue Depth': c_Residue_Depth_plotly.main,
    'Energy Heatmap': d_Energy_Heatmap.main,
    'Structure View': e_Structure_View.main,
    #'Mutations': Mutations.main,
}

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





# Initialization
# If user or username not in session state then set to None
if "name" not in st.session_state:
    st.session_state["name"] = None
if "username" not in st.session_state:
    st.session_state["username"] = None
ensure_state()
st.set_page_config(
    page_title="Hello",
    page_icon="👋",
    layout="wide",
)

# Create a new authenticator instance
with open('./config.yaml') as file:
    config = yaml.load(file, Loader=SafeLoader)
authenticator = stauth.Authenticate(
    config['credentials'],
    config['cookie']['name'],
    config['cookie']['key'],
    config['cookie']['expiry_days'],
    config['preauthorized']
    )

    
st.write("# Welcome to Protein Design Energy Analysis Tool 👋")



# Within the sidebar, we'll display a login form. If the user is not logged in,
# we'll show a login form. If the user is logged in, we'll show a logout button.
if st.session_state['name'] is None:

    name, authentication_status, username = authenticator.login('Login', 'sidebar')
    # Save authentication status to session state
    #st.session_state['authentication_status'] = authentication_status
    # Write the user's name to the session state
    #st.session_state['user'] = user
    #st.session_state['username'] = username
    if authentication_status:
        st.write("Welcome *%s*" % name)
    # your application
    elif authentication_status == False:
        st.error("Username/password is incorrect")
    elif authentication_status == None:
        st.warning("Please enter your username and password")


    
# Print in the sidebar the user's name with a smiley
if st.session_state["authentication_status"]:
    st.sidebar.write(f'Welcome **{st.session_state["name"]}**! 😃')
    authenticator.logout('Logout', 'sidebar')

    
# If the session state authentication status is false
elif st.session_state["authentication_status"] == None:
    try:
        if authenticator.register_user('Dont have a user yet?', preauthorization=False):
            st.success('User registered successfully')
    except Exception as e:
        st.error(e)


st.markdown(
    """

    **👈 Start exploring your protein desing in the pages!**

    This app contains a set of interactive tools to help determine key interactions
    within a protein structure that affect its stability and/or activity.
    This is designed to compare a single designed structure against its wild-type,
    not an array of designed structures. Furthermore, the analysis protocols do
    not support insertions and deletions. The variant structure should have the
    same length as the wild-type in order for the analysis to be accurate.
"""
)


with st.sidebar:
    STATE = st.session_state['Home']
    check_local_rosetta()
    if st.session_state['Home']:
        if 'rosetta_path' not in st.session_state['Home'].keys():
            st.session_state['Home']['rosetta_path'] = 'main/source/bin/residue_energy_' \
                                    'breakdown.static.linuxgccrelease'
    detect_rosetta()
    for key, value in STATUS.items():
        file_status(name=key, **value)




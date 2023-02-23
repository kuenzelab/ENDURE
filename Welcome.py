import streamlit as st
# Import yaml
import yaml
# Import SafeLoader
from yaml import SafeLoader
from functools import partial
import streamlit as st
import sys
import os
from PIL import Image
#from pages import (
#    a_File_Upload,   
#    b_Interaction_Analysis,
#    c_Residue_Depth_plotly, 
#    d_Energy_Heatmap,
#    #,Mutations
#)
from utility import load_text, add_logo
sys.path.append(os.path.dirname(__file__)) #
if 'root_dir' not in st.session_state.keys():# Save to the Session State the root directory of the project, remove the name of the file
    st.session_state['root_dir'] = os.path.dirname(__file__.replace('Welcome.py', ''))
    print("holiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii")
    print(__file__)
    print(st.session_state['root_dir'])
    # Add '/home/iwe30/Github/' to the session state
    st.session_state['root_dir'] = '/app/ENDURE'
    
st.set_page_config(
    page_title="Hello",
    page_icon="👋",
    layout="wide",
    initial_sidebar_state="expanded",
)

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




add_logo("images/draft_logo_200.png")


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



from st_pages import show_pages_from_config
show_pages_from_config(".streamlit/pages_sections.toml")



# Initialization
# If user or username not in session state then set to None
#if "name" not in st.session_state:
#    st.session_state["name"] = None
#if "username" not in st.session_state:
#    st.session_state["username"] = None
#ensure_state()


# Create a new authenticator instance
#with open('./config.yaml') as file:
#    config = yaml.load(file, Loader=SafeLoader)
#authenticator = stauth.Authenticate(
#    config['credentials'],
#    config['cookie']['name'],
#    config['cookie']['key'],
#    config['cookie']['expiry_days'],
#    config['preauthorized']
#    )
#
    



# Within the sidebar, we'll display a login form. If the user is not logged in,
# we'll show a login form. If the user is logged in, we'll show a logout button.
#if st.session_state['name'] is None:
#
#    name, authentication_status, username = authenticator.login('Login', 'sidebar')
#    # Save authentication status to session state
#    #st.session_state['authentication_status'] = authentication_status
#    # Write the user's name to the session state
#    #st.session_state['user'] = user
#    #st.session_state['username'] = username
#    if authentication_status:
#        st.write("Welcome *%s*" % name)
#    # your application
#    elif authentication_status == False:
#        st.error("Username/password is incorrect")
#    elif authentication_status == None:
#        st.warning("Please enter your username and password")
#
#
#    
## Print in the sidebar the user's name with a smiley
#if st.session_state["authentication_status"]:
#    st.sidebar.write(f'Welcome **{st.session_state["name"]}**! 😃')
#    authenticator.logout('Logout', 'sidebar')
#
#    
## If the session state authentication status is false
#elif st.session_state["authentication_status"] == None:
#    try:
#        if authenticator.register_user('Dont have a user yet?', preauthorization=False):
#            st.success('User registered successfully')
#    except Exception as e:
#        st.error(e)
#
#

st.write("# Welcome to ENDURE! 💻🔬")

st.markdown(
"""
**ENDURE is a user-friendly web application designed to help you analyze the energetic contributions of single and multi-variant protein designs. 
With ENDURE, you can quickly and easily evaluate the structural and thermodynamic changes that occur when mutations are introduced into your protein designs.**

## File Upload 📂

The first step in using ENDURE is to upload your protein models in PDB format. On the file upload page, simply select the PDB files you want to analyze and 
hit the "Upload" button. ENDURE will take care of the rest, running the necessary preprocessing steps to prepare your files for analysis.

## Interaction Analysis 🔍

This section provides an overview of the interactions between the residues in the uploaded protein structure. The user can select different types of interactions, 
such as salt bridges, sulfide bonds, and hydrogen bonds, and view their changes between a variant and wild-type structure.

## Residue Depth 📈

The residue depth page provides a visual representation of the change in residue depth that occurs when mutations are introduced. Residue depth is calculated 
using the Biopython library and is defined as the average distance (in angstroms) of the atoms in a residue from the solvent accessible surface. By analyzing 
the changes in residue depth, you can gain insights into how mutations affect the structural stability of your designs.

## Energy Heatmap 🔥

The energy heatmap page provides a visual representation of the changes in interaction energy that occur when mutations are introduced. 
The heatmap allows you to easily identify which residues are contributing the most to the changes in interaction energy, providing valuable 
insights into how to optimize your designs for stability and functionality.

## Example user workflow

In order to utilize the full capabilities of the ENDURE tool, it is important to understand the sequential nature of the tool and follow a specific workflow. 
A typical user would start by uploading their protein model in PDB format on the File Upload page. 

Once the model has been uploaded, the preprocessing steps including cleaning the PDB files, determining mutations, and calculating residue depth must be performed. 

Once these steps have been completed, the user can then proceed to the Interaction Analysis page, where the energetic contributions of single and multi-variant protein designs can be analyzed. 

The Residue Depth page provides a visual representation of the average distance of atoms in a residue from the solvent accessible surface, while the Energy Heatmap page displays the energy contribution of each residue in the protein model.

It is recommended for new users to follow this pipeline in order to effectively utilize the ENDURE tool:

    1. Upload protein model in PDB format on the File Upload page.
    2. Perform preprocessing steps including cleaning PDB files, determining mutations, and calculating residue depth.
    3. Analyze the energetic contributions of single and multi-variant protein designs on the Interaction Analysis page.
    4. View the Residue Depth and Energy Heatmap pages for additional insights into the protein model.
"""
)

# Add a streamlit logo to the sidebar
st.markdown("""
## Foreword 📖
We thank Paola Engelberger-Aliaga for the draft of the tool's logo! 🎨

With ENDURE, you have all the tools you need to make informed decisions about your protein designs. So why wait? Get started today and start exploring the exciting world of protein design! 🚀

We are constantly working to improve and expand the capabilities of ENDURE. We welcome suggestions for new functionalities and any feedback you may have about the tool. If you encounter any issues or would like to make a suggestion, please raise an issue in our GitHub repository. We are always looking for ways to make ENDURE the most helpful and efficient tool for analyzing protein designs. 

felipeengelberger@gmail.com

Thank you for using ENDURE! 💻🔬


"""
)





with st.sidebar:
    st.title("ENDURE")
    # Load the logo image function

    # initalize session state
    if 'Home' not in st.session_state:
        st.session_state['Home'] = {}
    if 'File Upload' not in st.session_state:
        st.session_state['File Upload'] = {}
    if 'Interaction Analysis' not in st.session_state:
        st.session_state['Interaction Analysis'] = {}
    if 'Residue Depth' not in st.session_state:
        st.session_state['Residue Depth'] = {}
    if 'Energy Heatmap' not in st.session_state:
        st.session_state['Energy Heatmap'] = {}
    
    # Hide all the python UserWarning

    STATE = st.session_state['Home']
    check_local_rosetta()
    if st.session_state['Home']:
        if 'rosetta_path' not in st.session_state['Home'].keys():
            st.session_state['Home']['rosetta_path'] = 'main/source/bin/residue_energy_' \
                                    'breakdown.static.linuxgccrelease'
    detect_rosetta()
    for key, value in STATUS.items():
        file_status(name=key, **value)




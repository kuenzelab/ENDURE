import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import streamlit as st
import streamlit_authenticator as stauth
# Import yaml
import yaml
# Import SafeLoader
from yaml import SafeLoader
from typing import Dict
from typing import Set
from streamlit_plotly_events import plotly_events
from st_aggrid import (
    AgGrid, DataReturnMode, GridUpdateMode, GridOptionsBuilder
)
from lib.visualization import WebViewer
from stmol import makeobj, showmol, render_pdb_resi, add_model, add_hover
import py3Dmol



import streamlit as st
from lib.visualization import WebViewer
from utility import load_text, add_logo

add_logo("images/draft_logo_200.png")

## TODO Move this to a stmol_mod module


def add_hover_res(obj,backgroundColor='white',fontColor='black'):
    """
    Adds a hover function to the Py3DMOL object to show the atom name when the mouse hovers over an atom.
    Example:
        obj = render_pdb()
        add_hover(obj)
        showmol(obj)
    Parameters
    ----------
    obj: Py3DMOL object
        Already existing Py3DMOL object, which can be created using the makeobj function.
    backgroundColor: String, default 'white'
        Is the background color of the hover text
    fontColor: String, default 'black'
        Is the color of the text
    Returns
    -------
    None.
    """

    js_script = """function(atom,viewer) {
                   if(!atom.label) {
                    atom.label = viewer.addLabel('Res # '+atom.resi,{position: atom, backgroundColor:"%s" , fontColor:"%s"});
                }
              }"""%(backgroundColor,fontColor)
    obj.setHoverable({},True,js_script,
               """function(atom,viewer) {
                   if(atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                   }
                }"""
               )


def add_model_sidechain(obj,molformat='pdb',model_style='sphere'):
    """
    Adds a model to an existing Py3DMOL object.
    
    Parameters
    ----------
    obj: Py3DMOL object
        Already existing Py3DMOL object.
    xyz: String
        Is the model to be added.
    molformat: String, default 'mol'    
        Is the format of the added model
    model_style: String, default 'stick'
        Is the style of the added model. Can be 'stick', 'sphere', 'cross', 'surface', 'ribbon', 'cartoon'
    
    Returns
    -------
    None.
    """
    #obj.addModel(xyz,molformat)
    obj.setStyle({'model':-1},{model_style:{}})
    obj.setStyle({'model':0},{'cartoon': {'color':'white'}},{model_style:'sphere'})
    obj.setStyle({'model':1},{'cartoon': {'color':'gray'}},{model_style:{}})
    obj.addStyle({'within':{'distance': 7, 'sel':{'resi':209}}},{'stick':{'colorscheme':'grayCarbon'}})
    obj.addStyle(
                    {
                        'model':0,'resi':209},
                    {
                        'stick':
                            {'colorscheme':'purpleCarbon'}
                    }
                )
    obj.addStyle(
                    {
                        'model':1,
                        'resi':209
                    },
                    {
                        'stick':
                        {
                        'colorscheme':'greenCarbon'
                        }
                    }
                )
    obj.zoomTo(
                {
                    'model':0,
                    'resi':209
                }
            )
    obj.setBackgroundColor('white')
    #Set a visualization style for chain B
    obj.setStyle({'chain':'B'},{'cartoon': {'color':'purple'}})
    #Add a visualization style for residue 158 in chain B
    obj.addStyle({'chain':'B','resi':158},{'stick':{'colorscheme':'grayCarbon'}})

def color_cartoon(obj, file_name: str, resi: int, color: str) -> None:
        """
        Set a particular location of the cartoon to a specific color
        :param file_name:
            A portion of the file name, such as "wild" or "variant", if the
            full name is "pdb_wild_clean" or "pdb_variant_clean"
        :param resi:
            The residue location to be colored
        :param color:
            The desired color
        :return:
            None
        """
        obj.__set_style(
            {
                'model': 0,
                'resi': resi
            },
            {
                'cartoon': {
                    'color': color
                }
            }
        )

STATE: dict


def check_files() -> bool:
    """
    Check that the required files exist in session state
    :return:
    """
    constraints = [
        'energy_wild' in st.session_state['File Upload'].keys(),
        'energy_variant' in st.session_state['File Upload'].keys(),
        'mutations' in st.session_state['File Upload'].keys()
    ]
    return all(constraints)

def parse_resi(user_input: str) -> None:
    """
    Read the residue selection string
    :param user_input:
    :return:
    """
    try:
        STATE['resi'] = [int(x) for x in user_input.split(',')]
    except ValueError:
        st.error('Invalid Residue String')



def color_select() -> None:
    """
    Create color pickers for changing colors
    :return:
    """
    # Background Color
    if 'background' not in STATE.keys():
        STATE['background'] = '#E2DFDF'
    color = st.color_picker('Background', STATE['background'])
    STATE['background'] = color

    # Cartoon Color
    if 'cartoon' not in STATE.keys():
        STATE['cartoon'] = '#858282'
    color = st.color_picker('Cartoon', STATE['cartoon'])
    STATE['cartoon'] = color

    # Wild Color
    if 'wild_color' not in STATE.keys():
        STATE['wild_color'] = '#04EEF3'
    color = st.color_picker('Wild Color', STATE['wild_color'])
    STATE['wild_color'] = color

    # Variant Color
    if 'variant_color' not in STATE.keys():
        STATE['variant_color'] = '#0FA81B'
    color = st.color_picker('Variant Color', STATE['variant_color'])
    STATE['variant_color'] = color




def changes() -> Dict[int, int]:
    """
    Determine if a residue is a mutation with better energy, mutation with
    worse energy, or not a mutation.
    :return:
    """
    wild: pd.DataFrame = st.session_state['File Upload']['energy_wild']

    variant: pd.DataFrame = st.session_state['File Upload']['energy_variant']


    mutations: pd.DataFrame = st.session_state['File Upload']['mutations']


    resi = wild['resi1'].values.tolist() + wild['resi2'].values.tolist()
    mut_idx = list(mutations.index)
    results = {}
    for i in range(1, max(resi) + 1):
        if i not in mut_idx:
            results[i] = 0
            continue
        query_wild = wild[
            (wild['resi1'] == i) |
            (wild['resi2'] == i)
        ]
        query_variant = variant[
            (variant['resi1'] == i) |
            (variant['resi2'] == i)
        ]
        net = query_variant['total'].sum() - query_wild['total'].sum()
        results[i] = 1 if net > 0 else 2
    return results


def create_source_wild(
    depth: Dict[int, float],
    net: Dict[int, float],
    mut_map: Dict[int, int]
) -> pd.DataFrame:
    """
    Create the Pandas DataFrame for the wild-type
    :param depth:
    :param net:
    :param mut_map:
    :return:
    """
    c1 = (9, 92, 224, 1)
    c2 = (224, 138, 9, 1)
    mut_list = [x != 0 for x in mut_map.values()]
    source = dict(
            x=list(depth.values()),
            y=list(net.values()),
            position=list(depth.keys()),
            color=[c2 if x else c1 for x in mut_list],
            label=['Mutated' if x else 'Conserved' for x in mut_list]
        )
    
    # Write as a CSV file from the dictionary
    # First save a dataframe
    df = pd.DataFrame(source)
    # Then save the dataframe as a CSV file
    df.to_csv('data/wild_dict_toDF.csv', index=False)



    return source


def create_source_variant(
    depth: Dict[int, float],
    net: Dict[int, float],
    mut_map: Dict[int, int]
) -> pd.DataFrame:
    """
    Create the Pandas DataFrame for the variant

    :param depth:
    :param net:
    :param mut_map:
    :return:
    """
    c1 = (9, 92, 224, 1)
    c2 = (224, 138, 9, 1)
    mut_list = [x != 0 for x in mut_map.values()]
    source = dict(
            x=list(depth.values()),
            y=list(net.values()),
            position=list(depth.keys()),
            color=[c2 if x else c1 for x in mut_list],
            label=['Mutated' if x else 'Conserved' for x in mut_list]
        )
    
    # Write as a CSV file from the dictionary
    # First save a dataframe
    df = pd.DataFrame(source)
    # Then save the dataframe as a CSV file
    df.to_csv('data/variant_dict_toDF.csv', index=False)




    return source


# A function lile create_plot but only for the data structures
# that are needed for the plot_master function
def create_data(file_name: str) -> pd.DataFrame:
    """
    A function lile create_plot but only for the data structures
    that are needed for the plot_master function
    :param file_name:
    :return:
    """
    # Fetch Data from Streamlit Session State
    inter: pd.DataFrame =\
        st.session_state['File Upload'][f'energy_{file_name}']


    # Create
    depth: Dict[int: float] =\
        st.session_state['File Upload'][f'depth_{file_name}']
    mut_map = changes()
    resi = inter['resi1'].values.tolist() + inter['resi2'].values.tolist()
    assert max(resi) == max(depth.keys())
    assert min(resi) == min(depth.keys()) == 1

    # Calculate Net Energy Changes
    entries = {}
    net = {}
    for i in depth.keys():
        query = inter[(inter['resi1'] == i) | (inter['resi2'] == i)]
        net[i] = query['total'].sum()
        entries[i] = query[['resi1', 'resi2', 'total']].values.tolist()
    # Convert the entries dictionary to a list of dataframes
    # Where for each key, the value is a dataframe
    entries_list = {}
    for key, value in entries.items():
        entries_list[key] = pd.DataFrame(value, columns=['Residue 1', 'Residue 2', 'Total energy[REU]']).sort_values(by=['Total energy[REU]'], ascending=True)
        # Round the column Residue 1 and Residue 2 to integers
        entries_list[key]['Residue 1'] = entries_list[key]['Residue 1'].astype(int)
        entries_list[key]['Residue 2'] = entries_list[key]['Residue 2'].astype(int)

    
    # Create Data Source
    if file_name == 'wild':
        source = create_source_wild(depth, net, mut_map)
        # Convert the dictionary to a dataframe
        source_df = pd.DataFrame(source)

    else:
        source = create_source_variant(depth, net, mut_map)
        # Convert the dictionary to a dataframe
        source_df = pd.DataFrame(source)

    return source_df, entries_list


@st.experimental_singleton
def load_data_variant() -> pd.DataFrame:
    # Load data\residue_depth_wild.csv
    # TODO update the source to the non static version 
    df = pd.read_csv("data/residue_depth_variant.csv")
    # Round the column x and y to two decimals
    df["x"] = df["x"].round(2)
    df["y"] = df["y"].round(2)
    # Remove the elements that are Conserved or Worse Energy in the column label
    df = df[df["label"] != "Conserved"]
    # Create a new df where the rows where label is Worse Energy are removed
    #df = df[df["label"] != "Worse Energy"]

    
    return df

@st.experimental_singleton
def load_data_wild() -> pd.DataFrame:
    # Load data\residue_depth_wild.csv
    # TODO update the source to the non static version
    df = pd.read_csv("data/residue_depth_wild.csv")
    # Round the column x and y to two decimals
    df["x"] = df["x"].round(2)
    df["y"] = df["y"].round(2)
    # Remove the elements that are Conserved in the column label
    df = df[df["label"] != "Conserved"]
    return df

def initialize_state():
    """Initializes all filters and counter in Streamlit Session State
    """
    for q in ["depth_wild", "depth_variant"]:
        if f"{q}_query" not in st.session_state:
            st.session_state[f"{q}_query"] = set()

    if "counter" not in st.session_state:
        st.session_state.counter = 0


def reset_state_callback():
    """Resets all filters and increments counter in Streamlit Session State
    """
    st.session_state.counter = 1 + st.session_state.counter

    for q in ["depth_variant", "depth_wild"]:
        st.session_state[f"{q}_query"] = set()


def query_data(df_wild: pd.DataFrame,
                df_variant: pd.DataFrame) -> pd.DataFrame:
    """Apply filters in Streamlit Session State
    to filter the input DataFrame
    """
    df_variant["depth_variant"] = (
        (100 * df_variant["x"]).astype(int).astype(str)
        + "-"
        + (100 * df_variant["y"]).astype(int).astype(str)
    )
    df_wild["depth_wild"] = (
        (100 * df_wild["x"]).astype(int).astype(str)
        + "-"
        + (100 * df_wild["y"]).astype(int).astype(str)
    )
    df_wild["selected"] = True
    df_variant["selected"] = True

    for q in ["depth_variant", "depth_wild"]:
        if st.session_state[f"{q}_query"]:
            if q == "depth_variant":
                df_wild.loc[~df_variant[q].isin(st.session_state[f"{q}_query"]), "selected"] = False
                df_variant.loc[~df_variant[q].isin(st.session_state[f"{q}_query"]), "selected"] = False
            if q == "depth_wild":
                df_wild.loc[~df_wild[q].isin(st.session_state[f"{q}_query"]), "selected"] = False
                df_variant.loc[~df_wild[q].isin(st.session_state[f"{q}_query"]), "selected"] = False



    return df_wild, df_variant


def build_depth_variant_figure(df: pd.DataFrame,
                            height : int,
                            width : int) -> go.Figure:
    # Create a streamlit slider from 1 to 20 with a step of 1
    slider_variant = st.slider("", 5, 20, 8, key="slider_variant")
    fig = px.scatter(
        df,
        "x",
        "y",
        color="selected",
        symbol="label",
        # 
        #labels={
        color_discrete_sequence=["rgba(99, 110, 250, 0.2)", "rgba(99, 110, 250, 1)"],
        category_orders={"selected": [False, True]},
        hover_data=[
            "position",
            "label",
            "x",
            "y",
            "selected"
        ],

    )
    fig.update_layout(paper_bgcolor="#FFFFFF", plot_bgcolor="#FFFFFF")
    fig.update_xaxes(gridwidth=0.1, gridcolor="#EDEDED")
    fig.update_yaxes(gridwidth=0.1, gridcolor="#EDEDED")
    # Make the dots bigger
    fig.update_traces(marker_size=slider_variant)
    # Add the x axis title and the y axis title depth in Å
    fig.update_xaxes(title_text="Depth in Å")
    # Net Interaction Energy in REU
    fig.update_yaxes(title_text="Net Interaction Energy [REU]")
    # Make the figure longer
    fig.update_layout(height=height, width=width)
    return fig


def build_depth_wild_figure(df: pd.DataFrame,
                            height : int,
                            width : int) -> go.Figure:
    # Create a streamlit slider from 1 to 20 with a step of 1
    slider_wild = st.slider("Select the size of the dots in the plot", 5, 20, 8, key="slider_wild")
    fig = px.scatter(
        df,
        "x",
        "y",
        color="selected",
        color_discrete_sequence=["rgba(99, 110, 250, 0.2)", "rgba(99, 110, 250, 1)"],
        category_orders={"selected": [False, True]},
        symbol="label",
        hover_data=[
            "position",
            "label",
            "x",
            "y",
            "selected"
        ],

    )
    fig.update_layout(paper_bgcolor="#FFFFFF", plot_bgcolor="#FFFFFF")
    fig.update_xaxes(gridwidth=0.1, gridcolor="#EDEDED")
    fig.update_yaxes(gridwidth=0.1, gridcolor="#EDEDED")
    # Make the dots bigger
    fig.update_traces(marker_size=slider_wild)
    # Add the x axis title and the y axis title depth in Å
    fig.update_xaxes(title_text="Depth in Å")
    # Net Interaction Energy in REU
    fig.update_yaxes(title_text="Net Interaction Energy [REU]")
    # Make the figure longer
    fig.update_layout(height=height, width=width)
    return fig




def render_preview_ui(df_variant: pd.DataFrame,
                        df_wild: pd.DataFrame):
    """Renders an expander with content of DataFrame and Streamlit Session State
    """
    with st.expander("Preview"):
        l, r = st.columns(2)
        l.dataframe(
            df_variant,
        )
        l.dataframe(
            df_wild,
        )
        r.json(
            {
                k: v
                for k, v in st.session_state.to_dict().items()
                if f'_{st.session_state["counter"]}' not in k
            }
        )


def render_plotly_ui(transformed_df_wild: pd.DataFrame,
                        transformed_df_variant: pd.DataFrame,
                            height : int,
                            width : int) -> Dict:
    """Renders all Plotly figures.
    Returns a Dict of filter to set of row identifiers to keep, built from the
    click/select events from Plotly figures.
    The return will be then stored into Streamlit Session State next.
    """
    c1, c2 = st.columns(2)

    with c1:
        # Write a subheader
        st.subheader("Depth Wild")
        depth_wild_selected = plotly_events(
            build_depth_wild_figure(transformed_df_wild, height, width),
            select_event=True,
            key=f"depth_wild_{st.session_state.counter}",
        )


    with c2:
        # Write a subheader
        st.subheader("Depth Variant")
        
        
        depth_variant_selected = plotly_events(
                build_depth_variant_figure(transformed_df_variant, height, width),
                select_event=True,
                key=f"depth_variant_{st.session_state.counter}",
            )



    current_query = {}
    current_query["depth_variant_query"] = {
        f"{int(100*el['x'])}-{int(100*el['y'])}" for el in depth_variant_selected
    }

    current_query["depth_wild_query"] = {
        f"{int(100*el['x'])}-{int(100*el['y'])}" for el in depth_wild_selected
    }


    return current_query


def update_state(current_query: Dict[str, Set]):
    """Stores input dict of filters into Streamlit Session State.
    If one of the input filters is different from previous value in Session State, 
    rerun Streamlit to activate the filtering and plot updating with the new info in State.
    """
    rerun = False
    for q in ["depth_wild","depth_variant"]:
        if current_query[f"{q}_query"] - st.session_state[f"{q}_query"]: # if there is a difference
            st.session_state[f"{q}_query"] = current_query[f"{q}_query"]
            rerun = True

    if rerun:
        st.experimental_rerun()




def main():
    # TODO Add a radio filter to the data by Mutated or Conserved and Better and Worse Energy for the varian
    
     # Create  a expandable text in the sidebar that explains the purpose of the app 
    with st.sidebar:
        st.expander("About", expanded=False).markdown("""The residue depth page provides a visual representation of the change in residue depth that occurs when mutations are introduced. Residue depth is calculated using the Biopython library and is defined as the average distance (in angstroms) of the atoms in a residue from the solvent accessible surface. By analyzing the changes in residue depth, you can gain insights into how mutations affect the structural stability of your designs.
    """)
    # Create a markdown text in the sidebar that indicate the molecule viewer style options
    st.sidebar.markdown("""
    ### Molecule Viewer Selection Options
    """)
    hl_resi_list = st.sidebar.multiselect(label="Extra residues to highlight in the structure viewer",
                                            options=list(range(1,1000)),
                                            # Write a detailed help message, step by step
                                            help="Select the residues to highlight in the structure viewer (e.g. 1, 2, 3), the residues will be highlighted in the structure viewer in addition to the residues selected in the table")
    st.sidebar.markdown("""
    ### Molecule Viewer Style Options
    """)
    hl_model = st.sidebar.text_input(label="Highlight Chain",
                                        value="A",
                                        # Write a detailed help message, step by step
                                        help="Enter the chain ID of the variant to highlight residues in the structure viewer (e.g. A), for the moment, only one chain can be highlighted at a time")

    label_resi = st.sidebar.checkbox(label="Label Residues", 
                                    value=True,
                                    # Write a detailed help message, step by step
                                    help="Do you want to label the residues in the structure viewer?, if yes, select the checkbox")
    
    surf_transp = st.sidebar.slider("Surface Transparency", 
                                        min_value=0.0, 
                                        max_value=1.0, 
                                        value=0.0,
                                        # Write a detailed help message, step by step
                                        help="Set the transparency of the surface in the structure viewer (0.0 is transparent, 1.0 is opaque)")

    # Top n energy interactions to show
    top_n = st.sidebar.slider("Top n interacting residue pairs", 
                                min_value=1, 
                                max_value=10, 
                                value=3,
                                # Write a detailed help message, step by step
                                help="Set the number of top interacting residue pairs to show in the structure viewer, these are the residue pairs with the lowest energy change meaning they are the most contributing to the energy change.")

    # Create a checkbox called zone selection and set it to False
    zone_selection = st.sidebar.checkbox(label="Zone Selection", 
                                        value=False,
                                        # Write a detailed help message, step by step
                                        help="Do you want to select a zone to view, if yes, select the checkbox and use the slider to select the zone in angstroms (1-10)")

    # If zone selection is True, show the zone slider
    if zone_selection:
        zone = st.sidebar.slider("Zone", min_value=1, max_value=10, value=7)
    hl_color = st.sidebar.text_input(label="Highlight Color",
                                    value="yellow",
                                    # Write a detailed help message, step by step
                                    help="Enter the color of the residue stick representation in the structure viewer (e.g. yellow) ")

    bb_color = st.sidebar.text_input(label="Backbone Color",
                                    value="lightgrey",
                                    # Write a detailed help message, step by step
                                    help="Enter the color of the backbone riboon representation in the structure viewer (e.g. lightgrey)")
    lig_color = st.sidebar.text_input(label="Ligand Color",
                                        value="white",
                                        # Write a detailed help message, step by step
                                        help="Enter the color of the ligand representation in the structure viewer (e.g. white), if no ligand is present, this will be ignored")
    # Create a sidebar slider to select the width of the structure viewer
    width = st.sidebar.slider("Width", min_value=300, max_value=1000, value=590)

    # Create a sidebar slider to select the height of the structure viewer
    height = st.sidebar.slider("Height", min_value=300, max_value=1000, value=590)

    cartoon_radius = 0.2
    stick_radius = 0.2


    df_variant = load_data_variant()
    # Create a list with the positions that have Worse Energy
    worse_energy_list = df_variant[df_variant.label == "Worse Energy"].position.to_list()
    
    
    df_wild = load_data_wild()
    # Drop the rows that have Worse Energy in the df_wild
    df_wild = df_wild[~df_wild.position.isin(worse_energy_list)]
    # Drop the rows that have Worse Energy in the df_variant
    df_variant = df_variant[~df_variant.position.isin(worse_energy_list)]

    transformed_df_wild, transformed_df_variant = query_data(df_wild,df_variant)
    st.button("Reset filters", on_click=reset_state_callback)
    # Run the create_data function to create the data
    source_wild, entries_wild = create_data("wild")
    source_variant, entries_variant = create_data("variant")
    st.title("Energy Breakdown by Residue Depth")

    #render_preview_ui(transformed_df_wild, transformed_df_variant)

    current_query = render_plotly_ui(transformed_df_wild,transformed_df_variant,height,width)
    update_state(current_query)

    # Write a subheader
    st.subheader("Selected Residue Pairwise Energy Breakdown")
    # Write a text
    st.write("Click on a residue in the scatter plot to see an energy breakdown of all the pairwise interactions of that residue")
    c1, c2 = st.columns(2)
    wild_structure = st.session_state['File Upload'][f'pdb_wild_clean'].getvalue()
    variant_structure = st.session_state['File Upload'][f'pdb_variant_clean'].getvalue()
    
    if len(transformed_df_wild[transformed_df_wild["selected"]]["position"]) > 1 and len(transformed_df_variant[transformed_df_variant["selected"]]["position"]) > 1:
        st.write("Please select only one residue in the scatter plot")
    else:
                # With c1
                with c1:
                    position_wild = int(transformed_df_wild[transformed_df_wild["selected"]]["position"].astype(int))
                    # Create a sum of the total energy for the wild type selected residue save it as a variable
                    total_energy_wild = entries_wild[position_wild]["Total energy[REU]"].sum()
                    wild_df = entries_wild[position_wild]
                    wild_df['Residue pair'] = wild_df.apply(lambda x: x['Residue 2'] if x['Residue 1'] == position_wild else x['Residue 1'], axis=1)
                    # Reset the index
                    wild_df = wild_df.reset_index(drop=True)
                    
                    # Convert to int
                    wild_df['Residue pair'] = wild_df['Residue pair'].astype(int)
                    # Drop the Residue 1 and Residue 2 columns
                    wild_df = wild_df.drop(['Residue 1', 'Residue 2'], axis=1)
                    # Reorder the columns
                    wild_df = wild_df[['Residue pair', 'Total energy[REU]']]
                    st.metric(f"Residue {position_wild} Sum Total Energy [REU]", f"{total_energy_wild:.2f}")
                    st.dataframe(wild_df.style.background_gradient(axis=0, subset='Total energy[REU]', cmap='coolwarm_r'), use_container_width=True)
                    # Create a list of the top n energy interactions Residue pairs
                    # Head returns the top n rows, extract the Residue pair column and convert to a list
                    top_n_wild = wild_df.head(top_n)["Residue pair"].tolist()

                    view_wt = makeobj(wild_structure,molformat='pdb',style='stick',background='white')

                    view_wt.setStyle({"cartoon": {"style": "oval","color": bb_color,"thickness": cartoon_radius}})

                    view_wt.addSurface(py3Dmol.VDW, {"opacity": surf_transp},
                                                        {"hetflag": False})

                    view_wt.addStyle({"elem": "C", "hetflag": True},
                                    {"stick": {"color": lig_color, "radius": stick_radius}})

                    view_wt.addStyle({"hetflag": True},
                                        {"stick": {"radius": stick_radius}})
                    # Append the position_variant to the h1_resi_list
                    hl_resi_list.append(position_wild)
                    # Add the top n energy interactions to the h1_resi_list
                    hl_resi_list.append(top_n_wild)
                    for hl_resi in hl_resi_list:
                        view_wt.addStyle({"chain": hl_model, "resi": hl_resi, "elem": "C"},
                                        {"stick": {"color": hl_color, "radius": stick_radius}})
                        

                        view_wt.addStyle({"chain": hl_model, "resi": hl_resi},
                                            {"stick": {"radius": stick_radius}})

                    if label_resi:
                        for hl_resi in hl_resi_list:
                            view_wt.addResLabels({"chain": hl_model,"resi": hl_resi},
                            {"backgroundColor": "lightgray","fontColor": "black","backgroundOpacity": 0.5})
                    view_wt.zoomTo(
                        {
                            'model':0,
                            'resi': position_wild
                        }
                    )
                    # If zone selection is True, add a style to the view_wt
                    if zone_selection:
                        view_wt.addStyle({'within':{'distance': zone, 'sel':{'resi':position_wild}}},{'stick':{'colorscheme':'grayCarbon'}})
                    # Add hover
                    add_hover_res(view_wt)
                    showmol(view_wt,height=height+20, width=width)                

                with c2:
                    position_variant = int(transformed_df_variant[transformed_df_variant["selected"]]["position"].astype(int))
                    # Create a sum of the total energy for the variant selected residue save it as a variable
                    total_energy_variant = entries_variant[position_variant]["Total energy[REU]"].sum()
                    variant_df = entries_variant[position_variant]
                    variant_df['Residue pair'] = variant_df.apply(lambda x: x['Residue 2'] if x['Residue 1'] == position_variant else x['Residue 1'], axis=1)
                    # Reset the index
                    variant_df = variant_df.reset_index(drop=True)
                    # Convert to int
                    variant_df['Residue pair'] = variant_df['Residue pair'].astype(int)
                    # Drop the Residue 1 and Residue 2 columns
                    variant_df = variant_df.drop(['Residue 1', 'Residue 2'], axis=1)
                    # Reorder the columns
                    variant_df = variant_df[['Residue pair', 'Total energy[REU]']]

                    st.metric(f"Residue {position_variant} Sum Total Energy [REU]", f"{total_energy_variant:.2f}")
                    st.dataframe(variant_df.style.background_gradient(axis=0, 
                                                                                            subset='Total energy[REU]', 
                                                                                            cmap='coolwarm_r'), 
                                                                                            use_container_width=True)
                    # Create a list of the top n energy interactions Residue pairs
                    # Head returns the top n rows, extract the Residue pair column and convert to a list
                    top_n_variant = variant_df.head(top_n)["Residue pair"].tolist()

                    view_variant = makeobj(variant_structure,molformat='pdb',style='stick',background='white')
                    view_variant.setStyle({"cartoon": {"style": "oval","color": bb_color,"thickness": cartoon_radius}})
                    view_variant.addSurface(py3Dmol.VDW, {"opacity": surf_transp},
                                                        {"hetflag": False})
                    view_variant.addStyle({"elem": "C", "hetflag": True},
                                    {"stick": {"color": lig_color, "radius": stick_radius}})
                    view_variant.addStyle({"hetflag": True},
                                        {"stick": {"radius": stick_radius}})
                    # Append the position_variant to the h1_resi_list
                    hl_resi_list.append(position_variant)
                    # Add the top n energy interactions to the h1_resi_list
                    hl_resi_list.append(top_n_variant)
                    for hl_resi in hl_resi_list:
                        view_variant.addStyle({"chain": hl_model, "resi": hl_resi, "elem": "C"},
                                        {"stick": {"color": hl_color, "radius": stick_radius}})
                        view_variant.addStyle({"chain": hl_model, "resi": hl_resi},
                                            {"stick": {"radius": stick_radius}})

                        
                    if label_resi:
                        for hl_resi in hl_resi_list:
                            view_variant.addResLabels({"chain": hl_model,"resi": hl_resi},
                            {"backgroundColor": "lightgray","fontColor": "black","backgroundOpacity": 0.5})
                    view_variant.zoomTo(
                        {
                            'model':0,
                            'resi': position_variant
                        }
                    )

                    add_hover_res(view_variant)
                    showmol(view_variant,height=height+20, width=width)


                



if __name__ == "__main__":

    # Initialize the session state
    if 'Structure View' not in st.session_state:
        st.session_state['Structure View'] = {}
    # if the pre-requisites are met, run the app

    STATE = st.session_state['Structure View']
    initialize_state()
    #If all the pre-requisites are not calculated, then return
    if not check_files():
        st.error('Error: Not all Pre-requisites are calculated')
    else:
        main()
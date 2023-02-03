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

## If user or username not in session state then set to None
if "name" not in st.session_state:
    st.session_state["name"] = None
if "username" not in st.session_state:
    st.session_state["username"] = None


import streamlit as st
from utility import load_text
from lib.visualization import WebViewer

STATE: dict


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


def create_viewer() -> None:
    """
    Create the PDB structure viewer
    :return:
    """
    from stmol import showmol
    import py3Dmol


    viewer = py3Dmol.view()
    for i in ['wild', 'variant']:
        viewer.add_model(i)
        viewer.show_cartoon(i, STATE['cartoon'])
        viewer.show_sc(
            i, STATE['resi'], STATE[f'{i}_color'], STATE['cartoon']
        )
    viewer.set_background(STATE['background'])
    viewer.center('wild', STATE['resi'])
    viewer.show()
    showmol(viewer, height = 500,width=800)


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


def tools() -> None:
    """
    Create toolbar for selecting residues
    :return:
    """
    # Residue Selection
    if 'resi' not in STATE.keys():
        STATE['resi'] = [1, 2]
    resi = st.text_input(
        label='Show Residue Side Chains:',
        value=','.join([str(x) for x in STATE['resi']])
    )
    parse_resi(resi)


def check_files() -> bool:
    """
    Check that necessary files exist in session state
    :return:
    """
    constraints = [
        'pdb_wild_clean' in st.session_state['File Upload'].keys(),
        'pdb_variant_clean' in st.session_state['File Upload'].keys()
    ]
    return all(constraints)



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
    # Remove the elements that are Conserved in the column label
    df = df[df["label"] != "Conserved"]
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


def build_depth_variant_figure(df: pd.DataFrame) -> go.Figure:
    # Create a streamlit slider from 1 to 20 with a step of 1
    slider_variant = st.slider("", 5, 20, 9, key="slider_variant")
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
        height=800,
        width=1100,
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
    fig.update_layout(height=500, width=600)
    return fig


def build_depth_wild_figure(df: pd.DataFrame) -> go.Figure:
    # Create a streamlit slider from 1 to 20 with a step of 1
    slider_wild = st.slider("Select the size of the dots in the plot", 5, 20, 9, key="slider_wild")
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
        height=800,
        width=1100,
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
    fig.update_layout(height=500, width=600)
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


def render_plotly_ui(transformed_df_variant: pd.DataFrame,
                        transformed_df_wild: pd.DataFrame) -> Dict:
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
            build_depth_wild_figure(transformed_df_variant),
            select_event=True,
            key=f"depth_wild_{st.session_state.counter}",
        )


    with c2:
        # Write a subheader
        st.subheader("Depth Variant")
        depth_variant_selected = plotly_events(
                build_depth_variant_figure(transformed_df_wild),
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
    # TODO Add a radio filter to the data by Mutated or Conserved and Better and Worse Energy for the variant
    """
    st.subheader("Wild Type Filter")
        # Create a selectbox with a label
        wild_type_filter = st.radio(
            # Use Conserved and Mutated as label
            "Wild Type",
            # List of options
            ("Mutated", "Conserved"),
            key = "wild_type_filter"
        )
        df_wild = load_data_wild(wild_type_filter)
        st.subheader("Variant Filter")
        # Create a selectbox with a label
        variant_filter = st.radio(
            # Use Conserved, Better Energy and Worse Energy as label
            "Variant",
            # List of options
            ( "Better Energy", "Conserved", "Worse Energy"),
            key = "variant_filter"
        )
        df_variant = load_data_variant(variant_filter)
    """
    df_variant = load_data_variant()
    df_wild = load_data_wild()
    transformed_df_wild, transformed_df_variant = query_data(df_wild,df_variant)
    st.button("Reset filters", on_click=reset_state_callback)
    # Run the create_data function to create the data
    source_wild, entries_wild = create_data("wild")
    source_variant, entries_variant = create_data("variant")
    st.title("Energy Breakdown by Residue Depth")
    #render_preview_ui(transformed_df_wild, transformed_df_variant)

    current_query = render_plotly_ui(transformed_df_wild,transformed_df_variant)
    update_state(current_query)

    
    c1, c2 = st.columns(2)
    if len(transformed_df_wild[transformed_df_wild["selected"]]["position"]) > 1 and len(transformed_df_variant[transformed_df_variant["selected"]]["position"]) > 1:
        st.write("Please select only one residue in the scatter plot")
    else:
                # With c1
                with c1:
                    position_wild = int(transformed_df_wild[transformed_df_wild["selected"]]["position"].astype(int))
                    position_wild
                    st.dataframe(entries_wild[position_wild].style.background_gradient(axis=0, subset='Total energy[REU]', cmap='coolwarm_r'), use_container_width=True)


                with c2:
                    position_variant = int(transformed_df_variant[transformed_df_variant["selected"]]["position"].astype(int))
                    position_variant
                    st.dataframe(entries_variant[position_wild].style.background_gradient(axis=0, subset='Total energy[REU]', cmap='coolwarm_r'), use_container_width=True)



# if the user is None then ask them to login
if st.session_state["name"] is None:
    # Please login
    st.write("Please go to the homepage to login")
else:
    if __name__ == "__main__":
        STATE = st.session_state['Structure View']

        initialize_state()
        main()
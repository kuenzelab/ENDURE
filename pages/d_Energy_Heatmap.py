import streamlit as st
import pandas as pd
import numpy as np
import colorcet as cc
from typing import Dict

from lib.visualization import WebViewer
import json
from typing import List
from Bio.PDB.PDBParser import PDBParser
import plotly.express as px
import plotly.graph_objs as go
from streamlit_plotly_events import plotly_events
# Import yaml
import yaml
# Import SafeLoader
from yaml import SafeLoader
from typing import Dict
from typing import Set
from stmol import makeobj, showmol, render_pdb_resi, add_model, add_hover
import py3Dmol

from utility import load_text, add_logo

add_logo("images/draft_logo_200.png")
STATE: dict

ROWS = [
    'fa_atr', 'fa_rep', 'fa_sol', 'fa_intra_rep', 'fa_intra_sol_xover4',
    'lk_ball_wtd', 'fa_elec', 'pro_close', 'hbond_sr_bb', 'hbond_lr_bb',
    'hbond_bb_sc', 'hbond_sc', 'dslf_fa13', 'omega', 'fa_dun', 'p_aa_pp',
    'yhh_planarity', 'ref', 'rama_prepro', 'total'
]

with open('lib/aa_map.json', 'r') as my_file:
    aa_map = json.load(my_file)

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

def show_cartoon_test(obj, color: str) -> None:
        """
        Show the entire structure as a cartoon
        :param file_name:
            A portion of the file name, such as "wild" or "variant", if the
            full name is "pdb_wild_clean" or "pdb_variant_clean"
        :param color:
            The color of the cartoon
        :return:
            None
        """
        obj.setStyle(
            {
                'model': 0
            },
            {
                'cartoon': {
                    'color': color
                }
            }
        )

def color_cartoon_scale(obj, resi: int, color: str) -> None:
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
        obj.setStyle(
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



@st.cache
def read_js(version: int) -> str:
    """
    Read JS Code from file
    :param version:
        The serial number of the code file
    :return:
        The contents of the file
    """
    with open(f'js/energy_heatmap_{version}.js', 'r') as file:
        return file.read()


def check_files() -> bool:
    """
    Check to see if necessary files are in session state
    :return:
    """
    constraints = [
        'energy_wild' in st.session_state['File Upload'].keys(),
        'energy_variant' in st.session_state['File Upload'].keys(),
        'mutations' in st.session_state['File Upload'].keys()
    ]
    return all(constraints)


def cool_stuff(
    data: pd.DataFrame,
    column_1: str,
    column_2: str,
    value: str
) -> pd.DataFrame:
    """
    Experiment with pd.melt
    :param data:
    :param column_1:
    :param column_2:
    :param value:
    :return:
    """
    pivot = pd.pivot_table(data, value, column_1, column_2, np.sum)
    pivot.fillna(value=0, inplace=True)
    pivot.reset_index(level=0, inplace=True)
    return pivot.melt(id_vars=column_1, var_name=column_2, value_name=value)


def fill_holes() -> Dict[str, pd.DataFrame]:
    """
    Ensure that the wild-type and variant interaction energy dataframes
    are of the same length. Important for JS Code when syncing selections
    :return:
    """
    wild = st.session_state['File Upload']['energy_wild']
    variant = st.session_state['File Upload']['energy_variant']
    cw = pd.DataFrame(wild[ROWS + ['resi1', 'resi2']], copy=True)
    cv = pd.DataFrame(variant[ROWS + ['resi1', 'resi2']], copy=True)
    pairs_w = cw[['resi1', 'resi2']].values.tolist()
    pairs_v = cv[['resi1', 'resi2']].values.tolist()

    new_values = []
    for resi1, resi2 in [x for x in pairs_v if x not in pairs_w]:
        row = {'resi1': resi1, 'resi2': resi2}
        row.update({x: 0 for x in ROWS})
        new_values.append(row)
    cw = pd.concat([cw, pd.DataFrame(new_values)])

    new_values = []
    for resi1, resi2 in [x for x in pairs_w if x not in pairs_v]:
        row = {'resi1': resi1, 'resi2': resi2}
        row.update({x: 0 for x in ROWS})
        new_values.append(row)
    cv = pd.concat([cv, pd.DataFrame(new_values)])

    cw.sort_values(by=['resi1', 'resi2'], inplace=True)
    cv.sort_values(by=['resi1', 'resi2'], inplace=True)
    cw.reset_index(inplace=True, drop=True)
    cv.reset_index(inplace=True, drop=True)
    return {'wild': cw, 'variant': cv}


def create_heatmap(
    file_name: str,
    data: pd.DataFrame,
    extrema: int = 5
) -> dict:
    """
    Create the Bokeh Components
    :param file_name:
        The partial file name
    :param data:
        The residue energy breakdown dataframe
    :param extrema:
        The extrema to use for the color bar, which is centered at 0
    :return:
    """
    # Setup Bokeh Plot
    reset = ResetTool()
    wheel_zoom = WheelZoomTool()
    pan_tool = PanTool()
    tap_tool = TapTool()
    poly = BoxSelectTool()
    save = SaveTool()
    tool_tips = [
        ('Resi1', '@x'),
        ('Resi2', '@y'),
        ('Total Energy', '@total{0.000}')
    ]
    plot = figure(
        tools=[reset, wheel_zoom, pan_tool, tap_tool, poly, save],
        tooltips=tool_tips
    )
    plot.title = f'Interaction Energy Pairs for {file_name.capitalize()}'
    plot.xaxis.axis_label = 'Position 1'
    plot.yaxis.axis_label = 'Position 2'
    plot.title.align = 'center'
    plot.title.text_font_size = '25px'

    # Create Data Source
    source_data = {
        'x': data['resi1'].values.tolist(),
        'y': data['resi2'].values.tolist()
    }
    source_data.update({x: data[x].values.tolist() for x in ROWS})

    source = ColumnDataSource(data=source_data)

    # Create Heatmap
    mapper = LinearColorMapper(
        palette=cc.b_linear_bmy_10_95_c78,
        low=-extrema,
        high=extrema
    )
    plot.rect(
        source=source,
        width=1,
        height=1,
        fill_color=transform('total', mapper),
        line_color=None
    )

    # Final Plot Configuration
    color_bar = ColorBar(color_mapper=mapper)
    plot.add_layout(color_bar, 'right')
    plot.toolbar.active_scroll = wheel_zoom
    return {
        'plot': plot,
        'source': source
    }


def plot_side() -> None:
    """
    Creates side by side heatmaps with linked axes for wild-type and variant
    :return:
    """
    # Create Heatmaps
    df = fill_holes()
    wild = create_heatmap('wild', df['wild'])
    variant = create_heatmap('variant', df['variant'])
    wild['plot'].width = 575
    variant['plot'].width = 575

    # Link Pan and Scroll
    wild['plot'].x_range = variant['plot'].x_range
    wild['plot'].y_range = variant['plot'].y_range

    # Bokeh Table
    source_table = ColumnDataSource(
        data=dict(
            energy=ROWS,
            wild=[0] * len(ROWS),
            variant=[0] * len(ROWS)
        )
    )
    table = DataTable(
        source=source_table,
        columns=[
            TableColumn(field='energy', title='Energy Term'),
            TableColumn(field='wild', title='Wild-Type'),
            TableColumn(field='variant', title='Variant'),
        ],
        index_position=None,
        width=250,
        height=535
    )

    # JS Code Linking Selections
    wild['source'].selected.js_on_change(
        'indices',
        CustomJS(
            args=dict(
                source=wild['source'],
                other=variant['source']
            ),
            code=read_js(1)
        )
    )
    variant['source'].selected.js_on_change(
        'indices',
        CustomJS(
            args=dict(
                source=variant['source'],
                other=wild['source']
            ),
            code=read_js(1)
        )
    )

    # JS Code Linking Selection to Table
    wild['source'].selected.js_on_change(
        'indices',
        CustomJS(
            args=dict(
                wild=wild['source'],
                variant=variant['source'],
                rows=ROWS,
                table=source_table
            ),
            code=read_js(2)
        )
    )

    # Show Bokeh Chart
    st.bokeh_chart(
        gridplot([
            [wild['plot'], variant['plot'], table],
        ])
    )

# A function to create the the dataset for the difference heatmap
def difference_dataset() -> list:
    """
    Create the dataset for the difference heatmap
    :return:
    """
    # Create Data and Heatmaps
    df = fill_holes()
    data = df['variant'] - df['wild']
    data['resi1'] = df['wild']['resi1']
    data['resi2'] = df['variant']['resi2']

    source_data = {
        'x': data['resi1'].values.tolist(),
        'y': data['resi2'].values.tolist()
    }
    source_data.update({x: data[x].values.tolist() for x in ROWS})
    # Create a dataframe from the source data
    dataframe = pd.DataFrame(source_data)

    minima = data['total'].min()
    maxima = data['total'].max()
    sequence_length = len(fasta('pdb_wild_clean'))
    # Currently source_data is a dictionary with keys 'x', 'y', 'total', 'vdw', 'elec', 'hbond', 'desolv', 'hba', 'hbd', 'ss', 'polar', 'apolar'
    # moreover the lists only contain the non zero values of the matrix, so we need to create a matrix with the right dimensions and fill it with the values
    # the matrix should be of size sequence_length x sequence_length
    # the values should be in the right position, so we need to iterate over the 'x' and 'y' and fill with the corresponding value in the 'total' list
    # we also need to fill the diagonal with 0
    matrix_list = []
    for i in range(sequence_length):
        matrix_list.append([0] * sequence_length)
    for i in range(len(source_data['x'])):
        matrix_list[source_data['x'][i] - 1][source_data['y'][i] - 1] = source_data['total'][i]
    # Replace the 0 with None
    for i in range(sequence_length):
        for j in range(sequence_length):
            if matrix_list[i][j] == 0:
                matrix_list[i][j] = None


    return matrix_list, minima, maxima, dataframe



def select_mode() -> str:
    """
    User selection of heatmap display mode
    :return:
    """
    if 'mode' not in STATE.keys():
        STATE['mode'] = 0
    radio = st.radio(
        label='Select Mode',
        options=['Difference', 'Side-by-Side'],
        horizontal=True
    )
    return radio


def resi_energy_map(
    wild: pd.DataFrame,
    variant: pd.DataFrame,
    colormap: list[str],
    min_value: float,
    max_value: float
) -> Dict[int, str]:
    """
    Create a colormap for the residues of the 3D structure based on their
    change in interaction energy from wild-type to variant
    :param wild:
    :param variant:
    :param colormap:
    :param min_value:
    :param max_value:
    :return:
    """

    resi_max = max(
        wild['resi1'].values.tolist() + wild['resi2'].values.tolist()
    ) # resi_max is the maximum residue number in the wild type

    energy: dict[int, float] = {}
    for i in range(1, resi_max + 1): # iterate over all residues
        energy_wild = wild[
            (wild['resi1'] == i) | (wild['resi2'] == i)
        ]['total'].sum() # sum the interaction energy of all interactions involving residue i in the wild type
        energy_variant = variant[
            (variant['resi1'] == i) | (variant['resi2'] == i)
        ]['total'].sum() # sum the interaction energy of all interactions involving residue i in the variant
        energy[i] = energy_variant - energy_wild # calculate the difference in interaction energy
    results = {}
    for key, value in energy.items():
        if value < min_value:
            results[key] = colormap[0]
        elif value > max_value:
            results[key] = colormap[-1]
        else:
            index = (value - min_value) / (max_value - min_value)
            results[key] = colormap[round(index * len(colormap))]
    return results


def fasta(file_name: str) -> List[str]:
    """
    Generate the FASTA sequence of a PDB File
    :param file_name:
        The name of the PDB file as stored in streamlit session state
    :return:
        The FASTA sequence as a list of strings
    """
    parser = PDBParser()
    parser.QUIET = True
    pdb_file = st.session_state['File Upload'][file_name]
    structure = parser.get_structure(0, pdb_file)
    pdb_file.seek(0)
    return [aa_map[x.resname] for x in structure.get_residues()]

def view_difference() -> None:
    """
    Create a 3D WebViewer to identify where energy changes occur in the
    conformation
    :return:
    """
    st.header('3D Structure Heatmap')
    viewer = WebViewer()
    viewer.add_model('wild')
    viewer.show_cartoon('wild', '#858282')
    cartoon_color = resi_energy_map(
        st.session_state['File Upload']['energy_wild'],
        st.session_state['File Upload']['energy_variant'],
        cc.b_linear_bmy_10_95_c78,
        -4, 4
    )
    for resi, color in cartoon_color.items():
        viewer.color_cartoon(viewer, resi, color)
    viewer.set_background('#E2DFDF')
    viewer.show()

### Interactive plotly functions ###

def initialize_state():
    """Initializes all filters and counter in Streamlit Session State
    """
    for q in ["difference_heatmap"]:
        if f"{q}_query" not in st.session_state:
            st.session_state[f"{q}_query"] = set()

    if "counter" not in st.session_state:
        st.session_state.counter = 0


def reset_state_callback():
    """Resets all filters and increments counter in Streamlit Session State
    """
    st.session_state.counter = 1 + st.session_state.counter

    for q in ["difference_heatmap"]:
        st.session_state[f"{q}_query"] = set()


def query_data(df: pd.DataFrame) -> pd.DataFrame:
    """Apply filters in Streamlit Session State
    to filter the input DataFrame
    """


    for q in ["difference_heatmap"]:
        if st.session_state[f"{q}_query"]:
            df.loc[~df[q].isin(st.session_state[f"{q}_query"]), "selected"] = False

    return df


def build_difference_heatmap_figure(df: pd.DataFrame) -> go.Figure:
    # Print st.session_state keys create a list
    #print (list(st.session_state['Energy Heatmap']))
    #st.write(len(fasta('pdb_wild_clean')))

    

    import plotly.express as px

    # Create a list from
    fig = px.imshow(df, x=df.columns, y=df.index,
                    labels=dict(x="Position 1", y="Position 2", color="[REU]"),
                    # Change the color scale invert the color scale
                    color_continuous_scale=px.colors.sequential.Bluered,
                )
    fig.update_xaxes(side="top")
    # Make the color legend vertical 90 degrees
    fig.update_layout(coloraxis_colorbar=dict(
        title="[REU]",
        thicknessmode="pixels", thickness=10,
        lenmode="pixels", len=300,
        yanchor="top", y=1,
        xanchor="left", x=1.05
    ))

    
    # Make the x-axis legend bigger
    fig.update_xaxes(title_font_size=15)
    # Make the y-axis legend bigger
    fig.update_yaxes(title_font_size=15)
    # Make the x axis labels bigger
    fig.update_xaxes(tickfont_size=10)
    fig.update_yaxes(tickfont_size=10)
    # Set the ticks every 10 and rotate the labels
    fig.update_xaxes(tickangle=-45, tickmode='linear', tick0=0, dtick=20)
    fig.update_yaxes(tickangle=-45, tickmode='linear', tick0=0, dtick=20)
    # Make the color legend bigger
    fig.update_layout(coloraxis_colorbar=dict(title_font_size=8))
    # make the height of the figure bigger
    fig.update_layout(width=450, height=450)
    if st.session_state["difference_heatmap_query"]:
        fig.add_annotation(
            x=st.session_state["difference_heatmap_query"][0],
            y=st.session_state["difference_heatmap_query"][1],
            xref="x",
            yref="y",
            text=f"Selected Pair",
            showarrow=True,
            font=dict(
                family="Courier New, monospace",
                size=16,
                color="#ffffff"
                ),
            align="center",
            arrowhead=2,
            arrowsize=1,
            arrowwidth=2,
            arrowcolor="#636363",
            ax=20,
            ay=-30,
            bordercolor="#c7c7c7",
            borderwidth=2,
            borderpad=4,
            bgcolor="#ff7f0e",
            opacity=0.8
            )
    return fig


def render_plotly_ui(transformed_df: pd.DataFrame) -> Dict:
    """Renders all Plotly figures.
    Returns a Dict of filter to set of row identifiers to keep, built from the
    click/select events from Plotly figures.
    The return will be then stored into Streamlit Session State next.
    """
    
    difference_heatmap_figure = build_difference_heatmap_figure(transformed_df)


    difference_heatmap_selected = plotly_events(
            difference_heatmap_figure,
            select_event=True,
            key=f"difference_heatmap_{st.session_state.counter}",
        )


    current_query = {}
    if difference_heatmap_selected:
        current_query["difference_heatmap_query"] = difference_heatmap_selected[0]['x'], difference_heatmap_selected[0]['y']
    


    return current_query


def update_state(current_query: Dict[str, Set]):
    """Stores input dict of filters into Streamlit Session State.
    If one of the input filters is different from previous value in Session State, 
    rerun Streamlit to activate the filtering and plot updating with the new info in State.
    """
    rerun = False
    for q in ["difference_heatmap"]:
        # If the current query is different from the previous one, update the state
        # If the current query contains elements
        if current_query.get(f"{q}_query"):
            if current_query[f"{q}_query"] != st.session_state[f"{q}_query"]:
                st.session_state[f"{q}_query"] = current_query[f"{q}_query"]
                rerun = True

    if rerun:
        st.experimental_rerun()


def main():
    """
    Create the Energy Heatmap Main Page
    :return:
    """
    global STATE
    STATE = st.session_state['Energy Heatmap']
    if check_files():
        
        # Create 3 columns
        col1, col2, col3 = st.columns(3)
        with col1:
            current_query = render_plotly_ui(df)
            update_state(current_query)
            st.button("Reset filters", on_click=reset_state_callback)
        with col2:
            if st.session_state["difference_heatmap_query"]:
                st.write("")
                st.write("")
                st.write("")
                st.write("")
                view_variant = makeobj(variant_structure,molformat='pdb',style='stick',background='white')
                view_variant.setStyle({"cartoon": {"style": "oval","color": bb_color,"thickness": cartoon_radius}})
                view_variant.addSurface(py3Dmol.VDW, {"opacity": surf_transp, "color": bb_color},{"hetflag": False})
                cartoon_color = resi_energy_map(
                    st.session_state['File Upload']['energy_wild'],
                    st.session_state['File Upload']['energy_variant'],
                    cc.b_linear_bmy_10_95_c78,
                    -4, 4
                )
                for resi, color in cartoon_color.items():
                    color_cartoon_scale(view_variant, resi, color)
                view_variant.addStyle({"elem": "C", "hetflag": True},
                                            {"stick": {"color": lig_color, "radius": stick_radius}})
                view_variant.addStyle({"hetflag": True},
                                                {"stick": {"radius": stick_radius}})

                if st.session_state["difference_heatmap_query"]:
                    hl_resi_list.append(st.session_state["difference_heatmap_query"][0])
                    hl_resi_list.append(st.session_state["difference_heatmap_query"][1])

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
                                    'resi': hl_resi_list
                                }
                            )
                add_hover_res(view_variant)
                showmol(view_variant,width=width, height=height)
            else:
                st.warning("Please select a residue pair to view the pair in 3D")
        with col3:

            # If st.session_state["difference_heatmap_query"]
            if st.session_state["difference_heatmap_query"]:
                st.write("")
                st.write("")
                st.write("")
                st.write("")
                # Filter the table based on st.session_state["difference_heatmap_query"][0] on column x
                filtered_table = table[table['x'] == st.session_state["difference_heatmap_query"][1]]
                # Filter the table based on st.session_state["difference_heatmap_query"][1] on column y
                filtered_table = filtered_table[filtered_table['y'] == st.session_state["difference_heatmap_query"][0]]
                # Transpose the table
                filtered_table = filtered_table.T
                # Drop first two rows
                filtered_table = filtered_table.drop(filtered_table.index[0:2])
                # Rename the first column to 'Energy Difference'
                filtered_table = filtered_table.rename(columns={filtered_table.columns[0]: 'Energy Difference'})
                # Create a new column with the index
                filtered_table['Rosetta Score Term'] = filtered_table.index
                # Make the index from 0 to len(filtered_table)
                filtered_table = filtered_table.reset_index(drop=True)
                # Rename the index to 'Rosetta Score Term'
                filtered_table.index.name = 'Rosetta Score Term'
                # Reorder the columns
                filtered_table = filtered_table[['Rosetta Score Term', 'Energy Difference']]
                # Convert the 'Energy Difference' column to float
                filtered_table['Energy Difference'] = filtered_table['Energy Difference'].astype(float)
                # Sort the table based on the 'Energy Difference' column
                filtered_table = filtered_table.sort_values(by=['Energy Difference'], ascending=True)
                # Drop Energy Difference rows with value 0
                filtered_table = filtered_table[filtered_table['Energy Difference'] != 0]
                #st.dataframe(filtered_table,use_container_width = True)
                # Create a plotly px bar chart
                fig = px.bar(filtered_table, x='Energy Difference', y='Rosetta Score Term', orientation='h')
                # Color each bar based with a color in a Pastel1 discrete color scale
                fig.update_traces(marker_color=px.colors.qualitative.Pastel1)
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("Select a residue pair in the heatmap to see the per score term energy difference.")
    else:
        st.error('Not all Pre-Requisites are Calculated')


if __name__ == '__main__':
    # If all the pre-requisites are not calculated, then return
    if not check_files():
        st.error('Error: Not all Pre-requisites are calculated')
    else:
        st.title("Energy Heatmap")
        wild_structure = st.session_state['File Upload'][f'pdb_wild_clean'].getvalue()
        variant_structure = st.session_state['File Upload'][f'pdb_variant_clean'].getvalue()
        # Create  a expandable text in the sidebar that explains the purpose of the app 
        with st.sidebar:
            st.expander("About", expanded=False).markdown("""The energy heatmap page provides a visual representation of the changes in interaction energy that occur when mutations are introduced. The heatmap allows you to easily identify which residues are contributing the most to the changes in interaction energy, providing valuable insights into how to optimize your designs for stability and functionality.
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
        width = st.sidebar.slider("Width", min_value=300, max_value=1000, value=400)

        # Create a sidebar slider to select the height of the structure viewer
        height = st.sidebar.slider("Height", min_value=300, max_value=1000, value=400)
        cartoon_radius = 0.2
        stick_radius = 0.2
        
        initialize_state()

        STATE = st.session_state['Energy Heatmap']


        # Print st.session_state keys create a list
        #print (list(st.session_state['Energy Heatmap']))
        #st.write(len(fasta('pdb_wild_clean')))
        matrix_list, minima, maxima, table = difference_dataset()

        # Create a st.slider for with a range from minima to maxima
        threshold = st.slider('Threshold', minima, maxima, 0.0)
        # Based on the threshold convert to None if an is above the threshold, if the value is already None skip it!!
        matrix_list = [[None if x is None or x > threshold else x for x in row] for row in matrix_list]
        # Create a dataframe from the matrix_list
        df = pd.DataFrame(matrix_list)
        wild_residues = fasta('pdb_wild_clean')
        variant_residues = fasta('pdb_variant_clean')
        # Create a list from 1 to len(fasta('pdb_wild_clean')
        wild_residues = list(range(1, len(fasta('pdb_wild_clean')) + 1))
        variant_residues = list(range(1, len(fasta('pdb_variant_clean')) + 1))
        # Add the variant residues as the column names
        df.columns = variant_residues
        # Add the wild residues as the index
        df.index = wild_residues
        df.fillna(0)
        main()
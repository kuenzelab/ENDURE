import pandas as pd
import streamlit as st
from typing import Dict
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import CustomJS, DataTable, TableColumn
from bokeh.models.tools import (
    WheelZoomTool, ResetTool, PanTool, TapTool, SaveTool
)
from bokeh.layouts import gridplot
from bokeh.models.widgets import Slider, CheckboxGroup, TextInput

import plotly.express as px
import plotly.graph_objs as go

import streamlit_authenticator as stauth
# Import yaml
import yaml
# Import SafeLoader
from yaml import SafeLoader
from typing import Dict
from typing import Set
from streamlit_plotly_events import plotly_events


STATE: dict


def check_files() -> bool:
    """
    Check that necessary data exists in session state
    :return:
    """
    constraints = [
        'depth_wild' in st.session_state['File Upload'].keys(),
        'energy_wild' in st.session_state['File Upload'].keys(),
        'depth_variant' in st.session_state['File Upload'].keys(),
        'energy_variant' in st.session_state['File Upload'].keys()
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


@st.cache
def read_js(version: int) -> str:
    """
    Read JS Code from file
    :param version:
        The serial number of the file to be read
    :return:
    """
    with open(f'js/residue_depth_{version}.js', 'r') as file:
        return file.read()


def create_source_wild(
    depth: Dict[int, float],
    net: Dict[int, float],
    mut_map: Dict[int, int]
) -> ColumnDataSource:
    """
    Create the Bokeh ColumnDataSource for the wild-type
    :param depth:
    :param net:
    :param mut_map:
    :return:
    """
    c1 = (9, 92, 224, 1)
    c2 = (224, 138, 9, 1)
    mut_list = [x != 0 for x in mut_map.values()]
    source = ColumnDataSource(
        data=dict(
            x=list(depth.values()),
            y=list(net.values()),
            position=list(depth.keys()),
            color=[c2 if x else c1 for x in mut_list],
            label=['Mutated' if x else 'Conserved' for x in mut_list]
        )
    )
    # Write as a CSV file
    df = pd.DataFrame(source.data)
    df.to_csv('data/residue_depth_wild.csv', index=False)


    return source


def create_source_variant(
    depth: Dict[int, float],
    net: Dict[int, float],
    mut_map: Dict[int, int]
) -> ColumnDataSource:
    """
    Create the Bokeh ColumnDataSource for the Variant
    :param depth:
    :param net:
    :param mut_map:
    :return:
    """
    color_map = [(9, 92, 224, 1), (230, 9, 9, 1), (69, 214, 95, 1)]
    label_map = ['Conserved', 'Worse Energy', 'Better Energy']
    source = ColumnDataSource(
        data=dict(
            x=list(depth.values()),
            y=list(net.values()),
            position=list(depth.keys()),
            color=[color_map[x] for x in mut_map.values()],
            label=[label_map[x] for x in mut_map.values()]
        )
    )
    # Write as a CSV file
    df = pd.DataFrame(source.data)
    df.to_csv('data/residue_depth_variant.csv', index=False)
    return source


def create_plot(file_name: str) -> dict:
    """
    Create Bokeh Scatter plot components and tables to be assembled
    in the plot_master function
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

    # Setup Bokeh Plot
    reset = ResetTool()
    wheel_zoom = WheelZoomTool()
    pan_tool = PanTool()
    tap_tool = TapTool()
    save = SaveTool()
    tool_tips = [
        ('Position', '@position'),
        ('Depth', '@x'),
        ('Net Energy', '@y')
    ]
    plot = figure(
        tools=[reset, wheel_zoom, pan_tool, tap_tool, save],
        tooltips=tool_tips
    )
    plot.title = f'{file_name.capitalize()}: Net Energy vs Depth'
    plot.xaxis.axis_label = 'Residue Depth from Surface'
    plot.yaxis.axis_label = 'Net Interaction Energy'

    # Create Data Source
    if file_name == 'wild':
        source = create_source_wild(depth, net, mut_map)

    else:
        source = create_source_variant(depth, net, mut_map)
    circles = plot.circle(
        source=source,
        legend_group='label',
        color='color',
        size=6
    )

    # Setup Bokeh Table
    source_table = ColumnDataSource(
        data=dict(resi1=[], resi2=[], total=[])
    )
    table = DataTable(
        source=source_table,
        columns=[
            TableColumn(field='resi1', title='Position 1'),
            TableColumn(field='resi2', title='Position 2'),
            TableColumn(field='total', title='Total Score')
        ]
    )

    # Final Plot Configuration
    plot.legend.location = 'top_right'
    plot.toolbar.active_scroll = wheel_zoom
    plot.title.align = 'center'
    plot.title.text_font_size = '25px'
    return {
        'plot': plot,
        'source': source,
        'table': table,
        'source_table': source_table,
        'entries': entries,
        'circles': circles
    }

# A function lile create_plot but only for the data structures
# that are needed for the plot_master function
def create_plot_data(file_name: str) -> dict:
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

    # Create Data Source
    if file_name == 'wild':
        source = create_source_wild(depth, net, mut_map)

    else:
        source = create_source_variant(depth, net, mut_map)

    return {
        'source': source,
        'entries': entries
    }

def plot_master() -> None:
    """
    Create and show the Bokeh Figure
    :return:
    """
    # Create Scatter Plots
    wild = create_plot('wild')
    variant = create_plot('variant')
    wild['plot'].x_range = variant['plot'].x_range
    wild['plot'].y_range = variant['plot'].y_range

    # Slider to control scatter point sizes
    slider = Slider(
        title='Circle Size',
        start=1,
        end=20,
        value=6,
        step=0.5
    )
    slider.js_link('value', wild['circles'].glyph, 'size')
    slider.js_link('value', variant['circles'].glyph, 'size')

    # Checkbox to toggle legend groups
    check_wild = CheckboxGroup(
        labels=['Conserved', 'Mutated'],
        active=[0, 1]
    )
    check_wild.js_on_click(
        CustomJS(
            args=dict(
                source=wild['source'],
                circles=wild['circles']
            ),
            code=read_js(2)
        )
    )
    check_variant = CheckboxGroup(
        labels=['Conserved', 'Better Energy', 'Worse Energy'],
        active=[0, 1, 2]
    )
    check_variant.js_on_click(
        CustomJS(
            args=dict(
                source=variant['source'],
                circles=variant['circles']
            ),
            code=read_js(3)
        )
    )

    var_maps = dict(
        wild=dict(
            source=wild['source'],
            s_e=wild['entries'],
            sc=wild['source_table'],
            oc=variant['source_table'],
            o_e=variant['entries'],
            o_s=variant['source']
        ),
        variant=dict(
            source=variant['source'],
            s_e=variant['entries'],
            sc=variant['source_table'],
            oc=wild['source_table'],
            o_e=wild['entries'],
            o_s=wild['source']
        )
    )
    # Create a df from the dict
    dicts=dict(
            s_e=wild['entries'],
            o_e=variant['entries'],
            
        )
    # Create df from wild dict
    dfj = pd.DataFrame.from_dict(dicts, orient='index')
    st.write(dfj)
    # JS Code Linking Graph Selection to Datatables
    for item in var_maps.values():
        item['source'].selected.js_on_change(
            'indices',
            CustomJS(
                args={**item},
                code=read_js(1)
            )
        )

    # Text Box to Select Residue by Entering number
    text = TextInput(title='Enter Residue ID:')
    text.js_on_change(
        'value',
        CustomJS(
            args=dict(source=wild['source']),
            code="""
            source.selected.indices = [parseInt(this.value) - 1];
            """
        )
    )

    # JS Code Linking Table Selection
    wild['source_table'].selected.js_on_change(
        'indices',
        CustomJS(
            args=dict(
                source=wild['source_table'],
                other=variant['source_table']
            ),
            code=read_js(4)
        )
    )
    variant['source_table'].selected.js_on_change(
        'indices',
        CustomJS(
            args=dict(
                source=variant['source_table'],
                other=wild['source_table']
            ),
            code=read_js(4)
        )
    )

    # Organize all components in Bokeh Grid
    st.bokeh_chart(
        gridplot([
            [slider, text],
            [check_wild, check_variant],
            [wild['plot'], variant['plot']],
            [wild['table'], variant['table']]
        ]),
        use_container_width=False
    )


# Create plotly figure
def create_plotly(df : pd.DataFrame) -> go.Figure:
    

    # Create Plotly Figure
    fig = px.scatter(
        df,
        "x",
        "y",
        color="label",
        color_discrete_sequence=["rgba(99, 110, 250, 0.2)", "rgba(99, 110, 250, 1)"],
        category_orders={"selected": ["Conserved", "Mutated"]},
        hover_data=[
            "position",
            "x",
            "y",
        ],
        height=800,
    )
    fig.update_layout(paper_bgcolor="#FFFFFF", plot_bgcolor="#FFFFFF")
    fig.update_xaxes(gridwidth=0.1, gridcolor="#EDEDED")
    fig.update_yaxes(gridwidth=0.1, gridcolor="#EDEDED")
    # Make the figure longer
    #fig.update_layout(height=2000)
    return fig


def initialize_state():
    """Initializes all filters and counter in Streamlit Session State
    """
    for q in ["bill_to_tip"]:
        if f"{q}_query" not in st.session_state:
            st.session_state[f"{q}_query"] = set()

    if "counter" not in st.session_state:
        st.session_state.counter = 0


# A function that returns the wild and variant source.data dataframes
def generate_dataframe(
    file_name : str
) -> pd.DataFrame:
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

    # Create Data Source
    if file_name == 'wild':
        source = create_source_wild(depth, net, mut_map)

    else:
        source = create_source_variant(depth, net, mut_map)
    return source.data

def render_plotly_ui() -> Dict:
    """Renders all Plotly figures.
    Returns a Dict of filter to set of row identifiers to keep, built from the
    click/select events from Plotly figures.
    The return will be then stored into Streamlit Session State next.
    """
    
    # Generate wild_df
    wild_df = generate_dataframe('wild')
    # Generate variant_df
    variant_df = generate_dataframe('variant')


    
    c1, c2 = st.columns(2)
    residue_depth_plotly_wild = create_plotly(wild_df)
    residue_depth_plotly_variant = create_plotly(variant_df)

    with c1:
        residue_depth_plotly_wild_selected = plotly_events(
            residue_depth_plotly_wild,
            select_event=True,
            key=f"residue_depth_variant_{st.session_state.counter}",
        )
    with c2:
        residue_depth_plotly_variant_selected = plotly_events(
            residue_depth_plotly_variant,
            select_event=True,
            key=f"residue_depth_variant_{st.session_state.counter}",
        )


def main():
    """
    Create the Residue Depth Main Page
    :return:
    """
    global STATE
    STATE = st.session_state['Residue Depth']
    st.title('Residue Depth in the Protein')
    if check_files():
        plot_master()
    else:
        st.error('Not all Pre-Requisites are Calculated')

if __name__ == '__main__':
    STATE = st.session_state['Residue Depth']
    # Print st.session_state keys create a list
    print (list(st.session_state['Residue Depth']))
    main()

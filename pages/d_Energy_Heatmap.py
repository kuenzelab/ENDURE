import streamlit as st
import pandas as pd
import numpy as np
import colorcet as cc
from typing import Dict
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import LinearColorMapper, ColorBar, DataTable, TableColumn
from bokeh.transform import transform
from bokeh.models.tools import (
    WheelZoomTool, ResetTool, PanTool, TapTool, BoxSelectTool, SaveTool
)
from bokeh.models.callbacks import CustomJS
from bokeh.layouts import gridplot, column
from bokeh.models.widgets import Slider, Button
from lib.visualization import WebViewer

STATE: dict

ROWS = [
    'fa_atr', 'fa_rep', 'fa_sol', 'fa_intra_rep', 'fa_intra_sol_xover4',
    'lk_ball_wtd', 'fa_elec', 'pro_close', 'hbond_sr_bb', 'hbond_lr_bb',
    'hbond_bb_sc', 'hbond_sc', 'dslf_fa13', 'omega', 'fa_dun', 'p_aa_pp',
    'yhh_planarity', 'ref', 'rama_prepro', 'total'
]


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
        'energy_variant' in st.session_state['File Upload'].keys()
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


def plot_difference() -> None:
    """
    Create a heatmap showing the difference in interaction energy from
    wild-type to variant
    :return:
    """
    # Create Data and Heatmaps
    df = fill_holes()
    data = df['variant'] - df['wild']
    data['resi1'] = df['wild']['resi1']
    data['resi2'] = df['variant']['resi2']
    st.write(data)
    diff = create_heatmap('Difference', data, extrema=2)

    # Bokeh Table
    source_table = ColumnDataSource(
        data=dict(
            energy=ROWS,
            diff=[0] * len(ROWS)
        )
    )
    # Create a dictionary of the data
    source_data_dict = dict(
            energy=ROWS,
            diff=[0] * len(ROWS)
        )
    st.write(diff['source'].data)
    # Create a dataframe from the diff['source'].data
    df = pd.DataFrame.from_dict(diff['source'].data)
    st.write(df)
    table = DataTable(
        source=source_table,
        columns=[
            TableColumn(field='energy', title='Energy Term'),
            TableColumn(field='diff', title='Difference')
        ],
        index_position=None,
        width=200,
        height=535
    )

    # JS Code Linking Selection to Table
    diff['source'].selected.js_on_change(
        'indices',
        CustomJS(
            args=dict(
                diff=diff['source'],
                rows=ROWS,
                table=source_table
            ),
            code=read_js(3)
        )
    )

    # JS Code to Jump to Specific Area of Map
    slider_x = Slider(
        title='Resi1',
        start=1,
        end=258,
        value=50,
        step=1
    )
    slider_y = Slider(
        title='Resi2',
        start=1,
        end=258,
        value=50,
        step=1
    )
    submit = Button(
        label='Go To Interaction Area'
    )
    submit.js_on_click(
        CustomJS(
            args=dict(
                slider_x=slider_x,
                slider_y=slider_y,
                x=diff['plot'].x_range,
                y=diff['plot'].y_range
            ),
            code="""
            x.start = slider_x.value - 15;
            x.end = slider_x.value + 15;
            y.start = slider_y.value - 15;
            y.end = slider_y.value + 15;
            """
        )
    )

    st.bokeh_chart(
        gridplot([
            [diff['plot'], table, column([slider_x, slider_y, submit])]
        ])
    )


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
    )
    energy: dict[int, float] = {}
    for i in range(1, resi_max + 1):
        energy_wild = wild[
            (wild['resi1'] == i) | (wild['resi2'] == i)
        ]['total'].sum()
        energy_variant = variant[
            (variant['resi1'] == i) | (variant['resi2'] == i)
        ]['total'].sum()
        energy[i] = energy_variant - energy_wild
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
        viewer.color_cartoon('wild', resi, color)
    viewer.set_background('#E2DFDF')
    viewer.show()


def main():
    """
    Create the Energy Heatmap Main Page
    :return:
    """
    global STATE
    STATE = st.session_state['Energy Heatmap']
    st.title('Energy Heatmap')
    if check_files():
        if select_mode() == 'Side-by-Side':
            plot_side()
        else:
            plot_difference()
            view_difference()
    else:
        st.error('Not all Pre-Requisites are Calculated')


if __name__ == '__main__':
    STATE = st.session_state['Energy Heatmap']
    # Print st.session_state keys create a list
    print (list(st.session_state['Energy Heatmap']))
    main()
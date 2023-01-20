import json
from functools import partial

import pandas as pd
import streamlit as st
from utility import load_text
from lib.energy_breakdown import energy_calc
from st_aggrid import (
    AgGrid, DataReturnMode, GridUpdateMode, GridOptionsBuilder
)

STATE: dict

KEY = 1
categories = {
    'salt_changes': 'Salt Bridges',
    'sulfide_changes': 'Sulfide Bonds',
    'hbonds_sc_sc': 'Hydrogen Bonds: Side-Chain to Side-Chain',
    'hbonds_bb_sc': 'Hydrogen Bonds: Side-Chain to Backbone',
    'hbonds_bb_bb_sr': 'Hydrogen Bonds: Backbone to Backbone Short Range',
    'hbonds_bb_bb_lr': 'Hydrogen Bonds: Backbone to Backbone Long Range',
    'all_changes': 'All Interactions'
}


def change_types() -> None:
    """
    Display information explaining the 6 categories of changes
    :return:
    """
    data = pd.DataFrame([
        ['A', True, True, False],
        ['B', True, True, True],
        ['C', True, False, False],
        ['D', True, False, True],
        ['E', False, True, False],
        ['F', False, True, True]
    ])
    data.columns = [
        'Type of Change',
        'Positions Interacting in Wild-Type',
        'Positions Interacting in Variant',
        'Involving Mutated Residues'
    ]
    data.set_index('Type of Change', inplace=True)
    st.dataframe(data,use_container_width=True)
    with open(f'text/interaction_changes/changes.json', 'r') as file:
        descriptions = json.load(file)
    data_2 = pd.DataFrame(descriptions['data'])
    data_2.columns = ['Change Type', 'Description']
    data_2.set_index('Change Type', inplace=True)
    data_2.columns = ['Description']
    #st.dataframe(data_2,use_container_width=True)


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


# TODO: Put this in a separate thread to prevent GUI Freeze
def start_calculations(progress_bar) -> None:
    """
    Execute the interaction analysis
    :param progress_bar:
        Progress bar to push updates to
    :return:
    """
    STATE['results'] = energy_calc(
        variant=st.session_state['File Upload']['energy_variant'],
        wild_type=st.session_state['File Upload']['energy_wild'],
        mutations=st.session_state['File Upload']['mutations'],
        bar=progress_bar
    )
    print(STATE['results'])
    STATE['complete'] = True


def build_options(data: pd.DataFrame) -> dict:
    """
    Configure options for the AgGrid Widget
    :param data:
    :return:
    """
    gb = GridOptionsBuilder.from_dataframe(data)
    gb.configure_default_column(resizable=False)
    gb.configure_pagination(
        paginationAutoPageSize=False,
        paginationPageSize=30
    )
    gb.configure_column(
        field='Resi2',
        menuTabs=['generalMenuTab', 'filterMenuTab']
    )
    options = gb.build()
    return options


def display_changes(label: str) -> None:
    """
    Display the results in an AgGrid Widget
    :param label:
    :return:
    """
    data = []
    for key, df in STATE['results'][label].items():
        if len(df) > 0:
            df: pd.DataFrame
            new = df[['resi1', 'resi2', 'total']].values.tolist()
            new = [x + [key.upper()] for x in new]
            data.extend(new)
    data = pd.DataFrame(data)
    data.columns = ['Resi1', 'Resi2', 'Total', 'Change Type']
    data['Total'] = data['Total'].round(3)
    data.insert(0, 'Change Type', data.pop('Change Type'))
    options = build_options(data)
    data = data.sort_values(by=['Total'], ascending=True)
    # Convert the Resi1 and Resi2 columns to integers
    data['Resi1'] = data['Resi1'].astype(int)
    data['Resi2'] = data['Resi2'].astype(int)
    st.dataframe(data.style.background_gradient(axis=0, subset='Total', cmap='coolwarm_r'), use_container_width=False)
    #TODO integrate AgGrid for interactive selection and view of residues with st-mol
    #grid_response = AgGrid(
    #    dataframe=data,
    #    data_return_mode=DataReturnMode.AS_INPUT,
    #    update_mode=GridUpdateMode.NO_UPDATE,
    #    gridOptions=options,
    #    try_to_convert_back_to_original_types=False,
    #    theme='streamlit'
    #)


def total_energy_changes(all_changes) -> float:
    print('here')
    sum = 0
    interactionList = [
                'a',
                'b',
                'c',
                'd',
                'e',
                'f'
            ]
    for interaction in interactionList:
            # Add the sum to sum
        sum += all_changes[f'{interaction}']['total'].sum()
    return sum

def significant_changes(all_changes) -> float:
    print('here')
    sum = 0
    interactionList = [
                'a',
                'b',
                'c',
                'd',
                'e',
                'f'
            ]
    for interaction in interactionList:
            # Add the sum to sum
        sum += all_changes[f'{interaction}'][(all_changes[f'{interaction}']['total'] >= 1) | (all_changes[f'{interaction}']['total'] <= -1 )]['total'].sum()
    return sum


def main():
    """
    Create the Interaction Analysis Main Page
    :return:
    """
    global STATE
    STATE = st.session_state['Interaction Analysis']
    
    if 'complete' not in STATE.keys():
        STATE['complete'] = False
    if 'check' not in STATE.keys():
        STATE['check'] = {
            x: {
                y: False for y in ['a', 'b', 'c', 'd', 'e', 'f']
            } for x in categories.keys()
        }

    st.title('Changes in Pairwise Interactions')
    st.header('Execute the Analysis')
    if not check_files():
        st.error('Error: Not all Pre-requisites are calculated')
        return
    st.success('Checking Pre-requisites: Success!')
    container = st.container()
    analysis_bar = st.progress(100 if STATE['complete'] else 0)
    container.button(
        label='Start Calculations',
        on_click=partial(
            start_calculations,
            progress_bar=analysis_bar
        ),
        disabled=STATE['complete']
    )
    if not STATE['complete']:
        return
    st.success('Successfully Executed Analysis and Stored the Results')
    st.header('Introduction')
    st.write(load_text('interaction_changes', 'introduction'))
    change_types()

    st.success('Successfully Executed Analysis and Stored the Results')
    tsum = total_energy_changes(STATE['results']['all_changes'])
    ssum = significant_changes(STATE['results']['all_changes'])
    # round to 1
    tsum = round(tsum, 1)
    ssum = round(ssum, 1)
    print(ssum)
    # Create a st.metric
    st.metric(label="Total Change in Energy", value=f"{tsum} REU")
    st.metric(label="Significant Change in Energy", value=f"{ssum} REU")
    st.header('Results')
    with st.expander('Summary', expanded=True):

        summary_df = pd.DataFrame(STATE['results']['summary'])
        # Convert to int
        summary_df = summary_df.astype(int).T
        with open(f'text/interaction_changes/changes.json', 'r') as file:
            descriptions = json.load(file)
        data_2 = pd.DataFrame(descriptions['data'])
        data_2.columns = ['Change Type', 'Description']
        data_2.set_index('Change Type', inplace=True)
        data_2.columns = ['Description']
        # Create a merged dataframe
        summary_df_m = data_2.merge(summary_df, left_index=True, right_index=True)
        # Add a new line in the column Description every 20 characters
        import textwrap
        summary_df_m['Description'] = summary_df_m['Description'].apply(lambda x: textwrap.fill(x, 20))
        st.dataframe(summary_df_m[['Description','Changes', 'Sum']].style.background_gradient(axis=1, vmax=100,cmap='coolwarm_r'),use_container_width=True)
        # st.dataframe the summary_df excluding the columns ['Changes', 'Sum']
        st.dataframe(summary_df.drop(['Changes', 'Sum'], axis=1).style.background_gradient(axis=1, cmap='coolwarm_r'),use_container_width=True)

        # 
    #Save the number of items in categories.items()
    num_items = 8
    num_items_per_col = int(num_items / 2)
    # Create two columns
    col1, col2 = st.columns(2)
    # With col1 add the first 4 items
    with col1:
        for key, value in list(categories.items())[:num_items_per_col]:
            with st.expander(value):
                display_changes(key)
    # With col2 add the last 4 items
    with col2:
        for key, value in list(categories.items())[num_items_per_col:]:
            with st.expander(value):
                display_changes(key)

if __name__ == '__main__':
    STATE = st.session_state['Interaction Analysis']
    # Print st.session_state keys create a list
    print (list(st.session_state['Interaction Analysis']))
    main()
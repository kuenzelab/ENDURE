import json
from functools import partial

import pandas as pd
import streamlit as st
from lib.energy_breakdown import energy_calc
from st_aggrid import (
    AgGrid, DataReturnMode, GridUpdateMode, GridOptionsBuilder
)
from lib.visualization import WebViewer
from stmol import makeobj, showmol, render_pdb_resi, add_model, add_hover
import py3Dmol


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


from utility import load_text, add_logo

add_logo("images/draft_logo_200.png")
STATE: dict

KEY = 1
categories_full = {
    'salt_changes': 'Salt Bridges',
    'sulfide_changes': 'Sulfide Bonds',
    'hbonds_sc_sc': 'Hydrogen Bonds: Side-Chain to Side-Chain',
    'hbonds_bb_sc': 'Hydrogen Bonds: Side-Chain to Backbone',
    'hbonds_bb_bb_sr': 'Hydrogen Bonds: Backbone to Backbone Short Range',
    'hbonds_bb_bb_lr': 'Hydrogen Bonds: Backbone to Backbone Long Range',
    'all_changes': 'All Interactions'
}

# Create a list of the values of the categories
category_values = list(categories_full.values())


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
    # Add the selection column
    gb.configure_selection(
        selection_mode = 'single', 
        pre_selected_rows = [0],
        use_checkbox=True,

    )
    options = gb.build()
    return options

def show_grid(data: pd.DataFrame):
    """
    Show the grid
    :return:
    """
    options = build_options(data)
    grid_response = AgGrid(
        dataframe=data,
        data_return_mode=DataReturnMode.AS_INPUT,
        update_mode=GridUpdateMode.SELECTION_CHANGED,
        gridOptions=options,
        try_to_convert_back_to_original_types=False,
        theme='streamlit'
    )
    selected = grid_response['selected_rows']
    #
    # Extract Resi1 and Resi2 from the selected rows and convert to a list where resi1 is the first element and resi2 is the second
    selected = [[x['Resi1'], x['Resi2']] for x in selected]

    return selected


def display_changes(label: str,
                    variant_structure: str ,
                    wild_structure: str,
                    top_n : int, 
                    hl_resi_list
                    , hl_model, label_resi, surf_transp, hl_color, bb_color, lig_color, width, height, cartoon_radius, stick_radius, change_type
                    

                    ) -> None:
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
    # Drop the rows that are not in the selected list change_type
    data = data[data['Change Type'].isin(change_type)]
    # Convert the Resi1 and Resi2 columns to integers
    data['Resi1'] = data['Resi1'].astype(int)
    data['Resi2'] = data['Resi2'].astype(int)
    # In the Change Type column, if the row contains 'A' append (WT+Variant) to the end of the string
    data['Change Type'] = data['Change Type'].apply(lambda x: x + ' (WT+Variant)' if 'A' in x else x)
    # In the Change Type column, if the row contains 'B' append (WT+Variant) to the end of the string
    data['Change Type'] = data['Change Type'].apply(lambda x: x + ' (WT+Variant)' if 'B' in x else x)
    # In the Change Type column, if the row contains 'C' append WT-Only to the end of the string
    data['Change Type'] = data['Change Type'].apply(lambda x: x + ' (WT-Only)' if 'C' in x else x)
    # In the Change Type column, if the row contains 'D' append WT-Only to the end of the string
    data['Change Type'] = data['Change Type'].apply(lambda x: x + ' (WT-Only)' if 'D' in x else x)
    # In the Change Type column, if the row contains 'E' append Variant-Only to the end of the string
    data['Change Type'] = data['Change Type'].apply(lambda x: x + ' (Variant-Only)' if 'E' in x else x)
    # In the Change Type column, if the row contains 'F' append Variant-Only to the end of the string
    data['Change Type'] = data['Change Type'].apply(lambda x: x + ' (Variant-Only)' if 'F' in x else x)
    
    # Create two columns 
    col1, col2 = st.columns(2)
    with col1:
        #st.dataframe(data.style.background_gradient(axis=0, subset='Total', cmap='coolwarm_r'), use_container_width=False)
        #TODO integrate AgGrid for interactive selection and view of residues with st-mol
        # If the data is empty, display a warning message
        if data.empty:
            st.warning(f'No {label} interactions found')
            selected = None
        else:
            selected = show_grid(data)
            # If selected is empty, cre


    with col2:
        if data.empty:
            st.warning(f'No {label} interactions found')
            selected = None
        else:

            view_variant = makeobj(variant_structure,molformat='pdb',style='stick',background='white')
            view_variant.setStyle({"cartoon": {"style": "oval","color": bb_color,"thickness": cartoon_radius}})
            view_variant.addSurface(py3Dmol.VDW, {"opacity": surf_transp, "color": bb_color},{"hetflag": False})
            view_variant.addStyle({"elem": "C", "hetflag": True},
                                        {"stick": {"color": lig_color, "radius": stick_radius}})
            view_variant.addStyle({"hetflag": True},
                                            {"stick": {"radius": stick_radius}})
            # Append the last two selected residues to the list of residues to highlight
            # If selected is not None, then append the last two selected residues to the list of residues to highlight
            if selected:

                for resi in selected:

                    hl_resi_list.append(resi[0])
                    hl_resi_list.append(resi[1])

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

    
    wild_structure = st.session_state['File Upload'][f'pdb_wild_clean'].getvalue()
    variant_structure = st.session_state['File Upload'][f'pdb_variant_clean'].getvalue()

    with st.sidebar:
        st.expander("About", expanded=False).markdown("""Once your files have been uploaded, you can dive into the interaction analysis page to get a detailed look at the changes in pairwise interactions that occur when mutations are introduced. ENDURE uses the Rosetta energy breakdown protocol to compare the output of the wild-type and variant designs, giving you a clear picture of how the interactions have changed.Once your files have been uploaded, you can dive into the interaction analysis page to get a detailed look at the changes in pairwise interactions that occur when mutations are introduced. ENDURE uses the Rosetta energy breakdown protocol to compare the output of the wild-type and variant designs, giving you a clear picture of how the interactions have changed.
    """)
        # Create a markdown text in the sidebar that indicate the molecule viewer style options
    st.sidebar.markdown("""
    ### Dataset filtering Options
    """)
    # Create multiselect in the sidebar called Change Type ['A', 'B', 'C', 'D', 'E', 'F'] set 'F' as default
    change_type = st.sidebar.multiselect(label="Change Type", 
                                        options=['A', 'B', 'C', 'D', 'E', 'F'], 
                                        default=['E','F'],
                                        # Write a detailed help message, step by step
                                        help="Select the type of interaction to show in the table and structure viewer (e.g. E and F which are the interactions present in the variant but not in the wild type, for more information about the interaction types, please refer to the table in the Introduction section)")
    





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
    width = st.sidebar.slider("Width", min_value=300, max_value=1000, value=1000)

    # Create a sidebar slider to select the height of the structure viewer
    height = st.sidebar.slider("Height", min_value=300, max_value=1000, value=1000)
    cartoon_radius = 0.2
    stick_radius = 0.2
    # Create a dictionary to store the selection and style options
    # the idea is to sync the selection and style options to a variable in the session state
    # so that the selection and style options are preserved when the user refreshes the page
    # or when the user comes back to the page

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

    categories_select = st.selectbox(
        'Interaction Categories',
        category_values,
        key=KEY,
        help='Select the interaction categories to display',


    )
    # Map each selected category to its key
    selected_categories = [k for k, v in categories_full.items() if v in categories_select]

    # Create a dictionary of the selected categories
    categories = {k: categories_full[k] for k in selected_categories}


    for key, value in list(categories.items()):
            with st.expander(value, expanded=True):
                
                display_changes(key, variant_structure, wild_structure, top_n, hl_resi_list, hl_model, label_resi, surf_transp, hl_color, bb_color, lig_color, width, height, cartoon_radius, stick_radius,change_type)
if __name__ == '__main__':


    STATE = st.session_state['Interaction Analysis']
    #If all the pre-requisites are not calculated, then return
    if not check_files():
        st.error('Error: Not all Pre-requisites are calculated')
    else:
        main()
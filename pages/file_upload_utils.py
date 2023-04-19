import json
import os
import re
import requests
import streamlit as st
import inspect
import pandas as pd
from io import StringIO, BytesIO
from typing import List, Callable
from functools import partial
from threading import Thread
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import ResidueDepth
from lib.energy_breakdown import run as eb_run, convert_outfile as eb_convert_outfile
from lib.fast_relax import run as fr_run

with open('lib/aa_map.json', 'r') as my_file:
    aa_map = json.load(my_file)


from streamlit.runtime.scriptrunner.script_run_context import (
    add_script_run_ctx
)

def add_logo(image_path: str) -> None:
    st.sidebar.image(image_path, use_column_width=True)


def file_status(name: str, error: str, success: str, warning: str) -> None:
    if name not in st.session_state['File Upload']:
        st.error(error)
    elif st.session_state['File Upload'][name]:
        st.success(success)
    else:
        st.warning(warning)


def new_files() -> None:
    if 'cleaned' in STATE.keys():
        STATE['cleaned'] = False


def file_uploader(key: str, value: dict) -> None:
    STATE[key] = st.file_uploader(**value, on_change=new_files)


def renumber_pdb(file_name: str) -> None:
    pdb_file = STATE[file_name]
    temp_file = StringIO()
    raw = pdb_file.read()
    text = raw.decode('utf-8').split('\n')
    rows = [x for x in text if x[0:4] == 'ATOM']
    offset = 1 - find_start(rows[0])
    temp_file.write('\n'.join(rows))
    temp_file.seek(0)
    data = pd.read_fwf(temp_file, header=None, infer_nrows=2000)
    data.columns = ['ATOM', 'atom_number', 'PDB_atom', 'resi_3', 'resi_1', 'resi_number', 'A', 'B', 'C', 'D', 'E', 'atom']
    data['resi_number'] = data['resi_number'] + offset
    data['PDB_atom'] = data['PDB_atom'].apply(adjust_pdb)
    rows = []
    spacing = [4, 7, 5, 4, 2, 4, 12, 8, 8, 6, 6, 12]
    for row in data.iterrows():
        row_txt = ''
        counter = 0
        for column in row[1].values.tolist():
            row_txt += f"{column:>{spacing[counter]}}"
            counter += 1
        rows.append(row_txt)
    result = StringIO()
    result.write('\n'.join(rows))
    result.seek(0)
    pdb_file.seek(0)
    STATE[f"{file_name}_clean"] = result


def fasta(file_name: str) -> List[str]:
    parser = PDBParser()
    parser.QUIET = True
    pdb_file = STATE[file_name]
    structure = parser.get_structure(0, pdb_file)
    pdb_file.seek(0)
    return [aa_map[x.resname] for x in structure.get_residues()]


def mutations() -> pd.DataFrame:
    variant = fasta('pdb_variant_clean')
    wild = fasta('pdb_wild_clean')
    counter = 1
    results = []
    for (v, w) in zip(variant, wild):
        if v != w:
            results.append([counter, v, w])
        counter += 1
    results = pd.DataFrame(results)
    results.columns = ['Position', 'Mutated', 'Wild']
    results.sort_values(by='Position', inplace=True)
    results.set_index(keys='Position', inplace=True)
    return results


def adjust_pdb(x: str) -> str:
    initial = len(x)
    if initial == 1:
        x += '  '
    if initial == 2:
        x += ' '
    if x[0].isnumeric() and initial == 3:
        x += ' '
    return x


def find_start(data: str) -> int:
    start = data.find(re.findall('\\w\\w\\w \\w', data)[0]) + 5
    number = ''
    final = False
    while len(number) == 0 or not final:
        if data[start].isnumeric():
            number += data[start]
        elif len(number) > 0:
            final = True
        start += 1
    return int(number)


def clean_pdb() -> None:
    for i in ['pdb_wild', 'pdb_variant']:
        if STATE[i] is not None:
            renumber_pdb(i)
    STATE['cleaned'] = True
    if 'mut_calc' in STATE.keys():
        STATE['mut_calc'] = False
    if 'depth' in STATE.keys():
        STATE['depth'] = False
    if 'breakdown' in STATE.keys():
        STATE['breakdown'] = False


def find_mutations() -> None:
    clean = ['pdb_wild_clean', 'pdb_variant_clean']
    if any([x not in STATE.keys() for x in clean]):
        return
    data = mutations()
    STATE['mutations'] = data
    STATE['mut_calc'] = True


def re_upload(key: str) -> None:
    STATE[key] = None


def calculate_depth(file_name: str) -> None:
    parser = PDBParser()
    parser.QUIET = True
    pdb_file = STATE[f"pdb_{file_name}_clean"]
    structure = parser.get_structure(0, pdb_file)
    pdb_file.seek(0)
    rd = ResidueDepth(model=structure[0], msms_exec='lib/msms_linux/msms.x86_64Linux2.2.6.1')
    results = {x[1][1]: y[0] for (x, y) in rd.property_dict.items()}
    STATE[f"depth_{file_name}"] = results
    STATE['depth'] = True


def calculate_energy(file_type: str) -> None:
    check_local_rosetta()
    pdb_file: StringIO = STATE[f"pdb_{file_type}_clean"]
    with open(f"lib/storage/{file_type}.pdb", 'w') as file:
        file.write(pdb_file.read())
    pdb_file.seek(0)
    if STATE['rosetta_local']:
        this_dir = st.session_state['root_dir']
        this_dir += '/ENDURE'
        root = '/'.join(this_dir.split('/')[:-1])
        executable = f"{root}/{STATE['rosetta_path']}"
    else:
        executable = st.session_state['Home']['rosetta_path']
    eb_run(file_name=f"lib/storage/{file_type}.pdb", save_path=f"lib/storage/energy_{file_type}.out",
           log_path=f"lib/storage/log_{file_type}.txt", executable=executable)
    eb_convert_outfile(file_name=f"lib/storage/energy_{file_type}.out", save_path=f"lib/storage/energy_{file_type}.csv")
    energy = pd.read_csv(f"lib/storage/energy_{file_type}.csv")
    energy.drop(energy[energy['resi2'] == '--'].index, inplace=True)
    energy['resi2'] = energy['resi2'].astype(int)
    STATE[f"energy_{file_type}"] = energy
    STATE['breakdown'] = True
    os.remove(f"lib/storage/energy_{file_type}.csv")
    os.remove(f"lib/storage/energy_{file_type}.out")


def calculate_relax() -> None:
    check_local_rosetta()
    for i in ['wild', 'variant']:
        if f"pdb_{i}_clean" in STATE.keys():
            pdb_file: StringIO = STATE[f"pdb_{i}_clean"]
            with open(f"lib/storage/{i}.pdb", 'w') as file:
                file.write(pdb_file.read())
            pdb_file.seek(0)
    if STATE['rosetta_local']:
        this_dir = st.session_state['root_dir']
        this_dir += '/ENDURE'
        root = '/'.join(this_dir.split('/')[:-1])
        executable = f"{root}/{STATE['rosetta_path']}"
    else:
        executable = st.session_state['Home']['rosetta_path']
    fr_run(file_name='lib/storage/wild.pdb', save_score_path='lib/storage/score_wild_relaxed.out',
           save_pdb_path='lib/storage/pdb_wild_relaxed.pdb', log_path='lib/storage/log_wild_relaxed.txt',
           executable=executable, nstruct=1)
    fr_run(file_name='lib/storage/variant.pdb', save_score_path='lib/storage/score_variant.out',
           save_pdb_path='lib/storage/pdb_variant.pdb', log_path='lib/storage/log_variant.txt',
           executable=executable, nstruct=1)
    STATE['relax'] = True


def check_rosetta() -> bool:
    if STATE['rosetta_local']:
        return True
    return STATE['rosetta_installed']


def check_local_rosetta() -> None:
    LOCAL_PATH = 'lib/rosetta_linux/source/bin/residue_energy_breakdown.static.linuxgccrelease'
    RS_LOCAL_PATH = 'lib/rosetta_linux/source/bin/rosetta_scripts.static.linuxgccrelease'
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
            STATE['rosetta_installed'] = True
            STATE['rosetta_scripts_path'] = RS_LOCAL_PATH
            st.session_state['rosetta_path'] = LOCAL_PATH
            st.session_state['rosetta_scripts_path'] = RS_LOCAL_PATH
            st.session_state['Home']['rosetta_installed'] = True


def show_action(header: str, text_file_name: str, button_label: str, callback: Callable) -> None:
    st.caption(f"**{header}**")
    with st.expander('Explanation'):
        st.write(load_text('file_upload', text_file_name))
    if text_file_name == 'energy_files' and not check_rosetta():
        st.error('No Rosetta Executable Available!')
        return
    status = st.container()
    status.write('')
    if len(inspect.signature(callback).parameters.keys()):
        st.button(label=button_label, on_click=partial(callback, status))
    else:
        st.button(label=button_label, on_click=callback)


def use_example(file_name: str) -> None:
    file = BytesIO()
    with open(f"lib/example_{file_name[4:]}.pdb", 'rb') as stream:
        file.write(stream.read())
    file.seek(0)
    file.name = f"example_{file_name[4:]}.pdb"
    STATE[file_name] = file

def main() -> None:
    """
    Creates the File Upload Page
    :return:
    """
    global STATE
    STATE = st.session_state['File Upload']
    # Create a reload function
    def reload() -> None:
        """
        Reloads the page
        :return:
        """
        st.experimental_rerun()
    # Add a reload button
    st.sidebar.button(
        label='Reload',
        on_click=reload,
        # add info about the reload button reload the page to check if the calculations are done
        # by default, the page will not reload
        help='Reload the page to check if the calculations are done, by default, the page will not reload'
    )
    # Add warning about the reload button
    st.sidebar.warning(
        # add a arrow up emoji
        '⬆️⬆️⬆️ \n Reload the page to check if the calculations are done(green), by default, the page will not reload'
    )
    #left, center, right = st.columns([1, 2, 1])
    # Create a container to hold the file uploaders
    file_uploaders = st.container()
    check_local_rosetta()
    with file_uploaders:

        st.header('Upload PDB Files or paste your sequences')

        file_uploader_widgets()
        st.header('Run pre-processing actions')

    

def second() -> None:
    check_local_rosetta()
    # If pdb_wild and pdb_variant are in the state, then we can
    # display the actions
    # Create a container to hold the actions
    col1, col2, col3, col4, col5 = st.columns(5)
    actions_container = st.container()
    if 'pdb_wild' in STATE.keys() and 'pdb_variant' in STATE.keys():
        with actions_container:
            counter = 0
            for col in [col1, col2, col3, col4, col5]:
                # with col:
                with col:
                    show_action(list(actions.keys())[counter],**list(actions.values())[counter])
                    counter += 1
# Create a third section to hold a general overview after the clean pdb files and Determining mutations is done
def third() -> None:
    """
    Creates the third section of the page, which displays a general overview of the session state mutations
    :return:
    """
    # Print this as a data frame
    #  st.session_state['File Upload']['mutations']
    # Create three columns to hold the data
    df = pd.DataFrame(st.session_state['File Upload']['mutations'])
    col1, col2, col3= st.columns([1, 1, 1])
    # Create a dictionary with aminoacid one letter codes as keys and the three letter codes as values and the type of amino acid
    # as the value
    amino_acids = {
        'A': {'three_letter': 'ALA', 'type': 'Aliphatic', 'type_emoji': '🟢'},
        'R': {'three_letter': 'ARG', 'type': 'Basic', 'type_emoji': '🔴'},
        'N': {'three_letter': 'ASN', 'type': 'Polar', 'type_emoji': '🟡'},
        'D': {'three_letter': 'ASP', 'type': 'Acidic', 'type_emoji': '🟤'},
        'C': {'three_letter': 'CYS', 'type': 'Sulfur', 'type_emoji': '🟣'},
        'E': {'three_letter': 'GLU', 'type': 'Acidic', 'type_emoji': '🟤'},
        'Q': {'three_letter': 'GLN', 'type': 'Polar', 'type_emoji': '🟡'},
        'G': {'three_letter': 'GLY', 'type': 'Aliphatic', 'type_emoji': '🟢'},
        'H': {'three_letter': 'HIS', 'type': 'Basic', 'type_emoji': '🔴'},
        'I': {'three_letter': 'ILE', 'type': 'Aliphatic', 'type_emoji': '🟢'},
        'L': {'three_letter': 'LEU', 'type': 'Aliphatic', 'type_emoji': '🟢'},
        'K': {'three_letter': 'LYS', 'type': 'Basic', 'type_emoji': '🔴'},
        'M': {'three_letter': 'MET', 'type': 'Sulfur', 'type_emoji': '🟣'},
        'F': {'three_letter': 'PHE', 'type': 'Aromatic', 'type_emoji': '🟠'},
        'P': {'three_letter': 'PRO', 'type': 'Aliphatic', 'type_emoji': '🟢'},
        'S': {'three_letter': 'SER', 'type': 'Polar', 'type_emoji': '🟡'},
        'T': {'three_letter': 'THR', 'type': 'Polar', 'type_emoji': '🟡'},
        'W': {'three_letter': 'TRP', 'type': 'Aromatic', 'type_emoji': '🟠'},
        'Y': {'three_letter': 'TYR', 'type': 'Aromatic', 'type_emoji': '🟠'},
        'V': {'three_letter': 'VAL', 'type': 'Aliphatic', 'type_emoji': '🟢'},
    }

    # In the first column, display the file status
    with col1:
        st.header('Overview')
        # Write in human language a in depth summary of the mutations
        # Create a list of strings from df that follows the format of the following
        # index, "Mutated" --> "Wild"
        list_of_strings = []
        list_of_strings_indexes = []
        # Count how many mutations are polar, aliphatic, basic, acidic, aromatic or sulfur in the Mutated column of df
        # Use a list comprehension to to get how many of the rows in the Mutated column of df are of a certain type
        # Use the dictionary amino_acids to get the type of amino acid
        counts_polar = len([row for index, row in df.iterrows() if amino_acids[row['Mutated']]['type'] == 'Polar'])
        counts_aliphatic = len([row for index, row in df.iterrows() if amino_acids[row['Mutated']]['type'] == 'Aliphatic'])
        counts_basic = len([row for index, row in df.iterrows() if amino_acids[row['Mutated']]['type'] == 'Basic'])
        counts_acidic = len([row for index, row in df.iterrows() if amino_acids[row['Mutated']]['type'] == 'Acidic'])
        counts_aromatic = len([row for index, row in df.iterrows() if amino_acids[row['Mutated']]['type'] == 'Aromatic'])
        counts_sulfur = len([row for index, row in df.iterrows() if amino_acids[row['Mutated']]['type'] == 'Sulfur'])


            

        for index, row in df.iterrows():
            # Use the three letter code of the amino acid to display the type of mutation
            # use the dictionary amino_acids to get the three letter code
            list_of_strings.append(f"{index},{row['Wild']} ➡️ {row['Mutated']}")
        st.write(f"""
                        With respect to the wild type PDB file, we found
                        **{len(st.session_state['File Upload']['mutations'])}** mutations
                """)
        st.write(f"""
                        Of which:
                        - {counts_polar} are {amino_acids['N']['type_emoji']} 
                        - {counts_aliphatic} are {amino_acids['G']['type_emoji']}
                        - {counts_basic} are {amino_acids['R']['type_emoji']}
                        - {counts_acidic} are {amino_acids['D']['type_emoji']}
                        - {counts_aromatic} are {amino_acids['F']['type_emoji']}
                        - {counts_sulfur} are {amino_acids['C']['type_emoji']}
                    """)
        
        

    with col2:
        st.header('In detail')
        # Create a expandable section to display the mutations in detail, by default it should be collapsed
        with st.expander('All mutations', expanded=True):
            counter = 0
            for string in list_of_strings:
                # use the type emoji to display the type of mutation
                st.write(f" - Position {df.index[counter]} : {string.split(',')[1]}  mutation")
                counter += 1

    # In the second column, display the mutations
    with col3:
        st.header('Legend')
        # Write a legend for the aminoacid types
        # Use the st.write function to write the legend
        st.write("")
        st.write("")
        st.write("")
        st.write("")
        st.write(f"""
                        - Polar {amino_acids['N']['type_emoji']}
                        - Aliphatic {amino_acids['A']['type_emoji']}
                        - Basic {amino_acids['R']['type_emoji']}        
                        - Acidic {amino_acids['D']['type_emoji']}
                        - Aromatic {amino_acids['F']['type_emoji']}
                        - Sulfur {amino_acids['C']['type_emoji']}
        """)


def esmfold_api_request_wild() -> tuple:
    api_url = 'https://api.esmatlas.com/foldSequence/v1/pdb/'
    seq_name = 'pdb_wild'
    if st.session_state.esmfold_prediction_wild is not None:
        prediction = requests.post(api_url, data=st.session_state.esmfold_prediction_wild)
        if prediction.status_code == 200:
            structure = prediction.content
            file = BytesIO()
            file.write(structure)
            file.seek(0)
            file.name = f"esm_{seq_name[4:]}.pdb"
            STATE[seq_name] = file
            message = 'Prediction for is now complete and ready for download'
        else:
            message = 'Sorry try again - the ESMFold API might be having problems'
            structure = None
        return message, structure


def esmfold_api_request_mutant() -> tuple:
    api_url = 'https://api.esmatlas.com/foldSequence/v1/pdb/'
    seq_name = 'pdb_variant'
    if st.session_state.esmfold_prediction_mutant is not None:
        prediction = requests.post(api_url, data=st.session_state.esmfold_prediction_mutant)
        if prediction.status_code == 200:
            structure = prediction.content
            file = BytesIO()
            file.write(structure)
            file.seek(0)
            file.name = f"esm_{seq_name[4:]}.pdb"
            STATE[seq_name] = file
            message = 'Prediction for is now complete and ready for download'
        else:
            message = 'Sorry try again - the ESMFold API might be having problems'
            structure = None
        return message, structure
    
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
        
def run(
    file_name: str,
    save_score_path: str,
    save_pdb_path: str,
    log_path: str,
    executable: str,
    nstruct: int
) -> None:
    pass

def file_uploader(key: str, value: dict) -> None:
    """
    Create the file uploader widget to accept uploaded PDB files
    :param key:
        The file name that will be stored in session state
    :param value:
        A dictionary of arguments for the streamlit widget
    :return:
    """
    STATE[key] = st.file_uploader(**value, on_change=new_files)

def file_uploader_widgets() -> None:
    """
    Create the File Uploaders and the associated functionality, including
    an option to re-upload a file and use an example file.
    :return:
    """
    global KEY
    # Create two columns to hold the file uploaders
    left, right = st.columns([1, 1])
    for key, value in files.items():
        # If key is pdb_wild 
        if key == 'pdb_wild':
            with left:
                if key not in STATE.keys() or STATE[key] is None:
                    file_uploader(key, value)
                    st.button(
                        label='Use Example File',
                        key=KEY,
                        on_click=partial(use_example, file_name=key)
                    )
                    KEY += 1
                    st.text_input(
                        value="",
                        max_chars=300,
                        label='Paste wildtype sequence for ESMFold Prediction',
                        key='esmfold_prediction_wild',
                        on_change=esmfold_api_request_wild
                    )
                    KEY += 1
                else:
                    st.success(
                        f'{key} is uploaded --- {STATE[key].name}'
                    )
                    st.button(
                        label='Re-upload?',
                        key=KEY,
                        on_click=partial(re_upload, key=key)
                    )
                    KEY += 1
        if key == 'pdb_variant':
            with right:
                if key not in STATE.keys() or STATE[key] is None:
                    file_uploader(key, value)
                    st.button(
                        
                        label='Use Example File',
                        key=KEY,
                        on_click=partial(use_example, file_name=key)
                    )
                    KEY += 1
                    st.text_input(
                        value="",
                        max_chars=300,
                        label='Paste mutant sequence for ESMFold Prediction',
                        key='esmfold_prediction_mutant',
                        on_change=esmfold_api_request_mutant
                    )
                    KEY += 1
                else:
                    st.success(
                        f'{key} is uploaded --- {STATE[key].name}'
                    )
                    st.button(
                        label='Re-upload?',
                        key=KEY,
                        on_click=partial(re_upload, key=key)
                    )
                    KEY += 1

# Find relaxed
def find_relaxed(container) -> None:
    """
    Call the relax_pdb function in a separate thread, and monitor
    this thread using add_script_run_ctx for the wild-type and variant
    :return:
    """
    for i in ['wild', 'variant']:
        if f'pdb_{i}_relaxed' in STATE.keys():
            task = Thread(target=partial(calculate_relax))
            add_script_run_ctx(task)
            task.start()
            container.warning(
                f'Calculations for {i} initiated in separate thread'
            )


def find_depth(container) -> None:
    """
    Call the calculate_depth function in a separate thread, and monitor
    this thread using add_script_run_ctx
    :return:
    """
    for i in ['wild', 'variant']:
        if f'pdb_{i}_clean' in STATE.keys():
            task = Thread(target=partial(calculate_depth, file_name=i))
            add_script_run_ctx(task)
            task.start()
            container.warning(
                f'Calculations for {i} initiated in separate thread'
            )


def find_energy(container) -> None:
    """
    Call the calculate_energy function in a separate thread, and monitor
    this thread using add_script_run_ctx
    :return:
    """
    for i in ['wild', 'variant']:
        if f'pdb_{i}_clean' in STATE.keys():
            task = Thread(target=partial(calculate_energy, i))
            add_script_run_ctx(task)
            task.start()
            container.warning(
                f'Calculations for {i} initiated in separate thread'
            )

def load_text(folder_name: str, file_name: str) -> str:
    """
    Load text from file
    :param folder_name:
    :param file_name:
    :return:
    """
    with open(f'text/{folder_name}/{file_name}.txt', 'r') as file:
        data = file.read()
    return data.strip('\n')


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

files = {
    'pdb_wild': {
        'label': 'PDB Structure: Wild-Type',
        'type': ['pdb'],
    },
    'pdb_variant': {
        'label': 'PDB Structure: Variant',
        'type': ['pdb'],
    }
}


actions = {
    '1) Cleaning PDB Files': dict(
        text_file_name='pdb_files',
        button_label='Clean PDB Files',
        callback=clean_pdb
    ),
    '2) Relaxing PDB Files': dict(
        text_file_name='pdb_files',
        button_label='Relax PDB Files',
        callback=find_relaxed
    ),
    '3) Determining Mutations': dict(
        text_file_name='mutations',
        button_label='Find Mutations',
        callback=find_mutations
    ),
    '4) Residue Depth': dict(
        text_file_name='residue_depth',
        button_label='Calculate Depth',
        callback=find_depth
    ),
    '5) Rosetta Energy Breakdown': dict(
        text_file_name='energy_files',
        button_label='Calculate Energy',
        callback=find_energy
    )
}


KEY = 1



import sys
import os
import pandas as pd
from typing import List, Dict, Tuple, Hashable
from streamlit.elements.progress import ProgressMixin
sys.path.append(os.path.dirname(__file__))
from rosetta import rosetta_simple


def run(
    file_name: str,
    save_score_path: str,
    save_pdb_path: str,
    log_path: str,
    executable: str,
    nstruct: int
) -> None:
    """
    Run the Rosetta Fast Relax Protocol in a Subprocess
    :param file_name:
        The name of the input file
    :param save_score_path:
        The path to save the output score file
    :param save_pdb_path:
        The path to save the output pdb file
    :param log_path:
        The path to save the log file
    :param executable:
        The path to the Rosetta executable
    :return: None
    """
    print(f'Started Job for {file_name} with nstruct: {nstruct}')
    executable = 'rosetta_scripts.static.linuxgccrelease'

    options = [
        f'-s {file_name}',
        f'-native {file_name}',
        f'-parser:protocol lib/XML/fast_relax.xml',
        f'-out:pdb',
        f'-out:path:pdb {save_pdb_path}',
        f'-nstruct {nstruct}',
        f'-out:file:scorefile {save_score_path}'
    ]
    # Add the options Unconstrained and Constrained
    #options.append('-relax:constrain_relax_to_start_coords')

    log = rosetta_simple(executable, options)
    with open(log_path, 'w') as file:
        file.write(log)

def convert_outfile(
    file_name: str,
    save_path: str,
) -> None:
    """
    Convert the output of Rosetta Energy Breakdown to a CSV File
    :param file_name:
        Location of the output file
    :param save_path:
        Location to save the CSV file
    :return: None
    """
    spacer = len(file_name) - 17
    with open(file_name, 'r') as file:
        data = file.read()
    lines = data.split('\n')
    lines[0] = lines[0][0:6] + ' '*spacer + lines[0][6:]

    with open(f'{file_name[:-3]}txt', 'w') as file:
        file.write('\n'.join(lines))

    data = pd.read_fwf(f'{file_name[:-3]}txt')
    data.drop(columns=['description', 'SCORE:', 'pose_id'], inplace=True)
    data.to_csv(save_path, index=False)
    os.remove(f'{file_name[:-3]}txt')
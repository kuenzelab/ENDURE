import json, os, re, inspect, streamlit as st, pandas as pd
from io import StringIO, BytesIO
from typing import List, Callable
from functools import partial
from threading import Thread
from pages.file_upload_utils import (load_text, add_logo, file_status, new_files, file_uploader, renumber_pdb, fasta, mutations, adjust_pdb, find_start, clean_pdb, find_mutations,
                     re_upload, calculate_depth, calculate_energy, run, calculate_relax, check_rosetta, check_local_rosetta, show_action, use_example, esmfold_api_request_wild, esmfold_api_request_mutant,
                     file_uploader_widgets, main, second, third, STATUS)
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import ResidueDepth
from lib.energy_breakdown import run as eb_run, convert_outfile as eb_convert_outfile
from lib.fast_relax import run as fr_run

STATE: dict
with open('lib/aa_map.json', 'r') as my_file:
    aa_map = json.load(my_file)

if 'esmfold_prediction_wild' not in st.session_state:
    st.session_state['esmfold_prediction_wild'] = ''
if 'esmfold_prediction_mutant' not in st.session_state:
    st.session_state['esmfold_prediction_mutant'] = ''

files = {'pdb_wild': {'label': 'PDB Structure: Wild-Type', 'type': ['pdb']},
         'pdb_variant': {'label': 'PDB Structure: Variant', 'type': ['pdb']}}
KEY = 1

if __name__ == '__main__':
    add_logo('images/draft_logo_200.png')

    with st.sidebar:
        for (key, value) in STATUS.items():
            file_status(name=key, **value)
    main()
    second()
    if 'mutations' not in st.session_state['File Upload']:
        st.warning('Please upload a PDB file to see the mutations and run the clean and calculate mutations functions')
    else:
        third()
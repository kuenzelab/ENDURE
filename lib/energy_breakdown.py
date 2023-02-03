import sys
import os
import pandas as pd
from typing import List, Dict, Tuple, Hashable
from streamlit.elements.progress import ProgressMixin
sys.path.append(os.path.dirname(__file__))
from rosetta import rosetta_simple
# Print the path variable
print(sys.path)

def run(
    file_name: str,
    save_path: str,
    log_path: str,
    executable: str,
) -> None:
    """
    Run the Rosetta Energy Breakdown Protocol in a Subprocess
    :param file_name:
        The input PDB file on disk
    :param save_path:
        Location to output the result file
    :param log_path:
        Location to save the log file
    :param executable:
        Filepath of the Rosetta executable to run
    :return: None
    """
    options = [
        f'-in:file:s {file_name}',
        f'-out:file:silent {save_path}'
    ]
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


def energy_calc(
    variant: pd.DataFrame,
    wild_type: pd.DataFrame,
    mutations: pd.DataFrame,
    bar: ProgressMixin
) -> Dict[str, pd.DataFrame or Dict[str, pd.DataFrame]]:
    """
    Identify Important Changes in Interaction Energies
    :param variant:
        The residue energy breakdown for the variant structure
    :param wild_type:
        The residue energy breakdown for the wild-type structure
    :param mutations:
        The mutations between the wild-type and variant
    :param bar:
        The streamlit progressbar to push updates to
    :return:
        A dictionary containing assorted dataframes
    """
    bar.progress(5)
    results = interaction_analysis(
        variant, wild_type, mutations.index.tolist()
    )
    bar.progress(25)
    salt_changes = interaction_analysis(
        salt_bridges(variant),
        salt_bridges(wild_type),
        mutations.index.tolist()
    )
    bar.progress(35)
    sulfide_changes = interaction_analysis(
        sulfide_bonds(variant),
        sulfide_bonds(wild_type),
        mutations.index.tolist()
    )
    bar.progress(45)
    hbonds_v = hydrogen_bonds(variant)
    hbonds_w = hydrogen_bonds(wild_type)
    hbonds_sc_sc = interaction_analysis(
        hbonds_v['sc_sc'],
        hbonds_w['sc_sc'],
        mutations.index.tolist()
    )
    bar.progress(55)
    hbonds_bb_sc = interaction_analysis(
        hbonds_v['bb_sc'],
        hbonds_w['bb_sc'],
        mutations.index.tolist()
    )
    bar.progress(65)
    hbonds_bb_bb_sr = interaction_analysis(
        hbonds_v['bb_bb_sr'],
        hbonds_w['bb_bb_sr'],
        mutations.index.tolist()
    )
    bar.progress(75)
    hbonds_bb_bb_lr = interaction_analysis(
        hbonds_v['bb_bb_lr'],
        hbonds_w['bb_bb_lr'],
        mutations.index.tolist()
    )
    bar.progress(85)
    data = [
        {x.upper(): len(y) for x, y in results.items()},
        {x.upper(): round(y['total'].sum(), 4) for x, y in results.items()},
        {x.upper(): len(y[y['total'] >= 1]) for x, y in results.items()},
        {x.upper(): len(y[y['total'] <= -1]) for x, y in results.items()},
        {x.upper(): len(y) for x, y in salt_changes.items()},
        {x.upper(): len(y) for x, y in sulfide_changes.items()},
        {x.upper(): len(y) for x, y in hbonds_sc_sc.items()},
        {x.upper(): len(y) for x, y in hbonds_bb_sc.items()},
        {x.upper(): len(y) for x, y in hbonds_bb_bb_sr.items()},
        {x.upper(): len(y) for x, y in hbonds_bb_bb_lr.items()}
    ]
    bar.progress(95)
    data = pd.DataFrame(data)
    # Round to 1 decimal place
    data = data.round(1)
    data.index = [
        'Changes',
        'Sum',
        'Energy > 1',
        'Energy < -1',
        'Salt Bridge',
        'Sulfide Bonds',
        'SC-SC HBonds',
        'BB-SC HBonds',
        'BB-BB-SR HBonds',
        'BB-BB-LR HBonds'
    ]
    bar.progress(99)
    return {
        'summary': data,
        'all_changes': results,
        'salt_changes': salt_changes,
        'sulfide_changes': sulfide_changes,
        'hbonds_sc_sc': hbonds_sc_sc,
        'hbonds_bb_sc': hbonds_bb_sc,
        'hbonds_bb_bb_sr': hbonds_bb_bb_sr,
        'hbonds_bb_bb_lr': hbonds_bb_bb_lr
    }


def salt_bridges(interactions: pd.DataFrame) -> pd.DataFrame:
    """
    Find Salt Bridges
    :param interactions:
        The residue energy breakdown dataframe
    :return:
        A dataframe of the salt bridge interactions
    """
    query = interactions[
        (interactions['hbond_sc'] < 0) &
        (
            (
                (interactions['restype1'].isin(['ASP', 'GLU'])) &
                (interactions['restype2'].isin(['ARG', 'HIS', 'LYS']))
            ) |
            (
                (interactions['restype1'].isin(['ARG', 'HIS', 'LYS'])) &
                (interactions['restype2'].isin(['ASP', 'GLU']))
            )
        )
    ]
    return query


def sulfide_bonds(interactions: pd.DataFrame) -> pd.DataFrame:
    """
    Find Disulfide Bonds
    :param interactions:
        The residue energy breakdown dataframe
    :return:
        A dataframe of the disulfide interactions
    """
    query = interactions[
        (
            (interactions['restype1'] == 'CYS') &
            (interactions['restype2'] == 'CYS')
        ) |
        (interactions['dslf_fa13'] < 0)
    ]
    return query


def hydrogen_bonds(interactions: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Find Hydrogen Bonds
    :param interactions:
        The residue energy breakdown dataframe
    :return:
        A dictionary containing a dataframe of interactions
         for each type of hydrogen bond
    """
    sc_sc = interactions[interactions['hbond_sc'] < 0]
    bb_sc = interactions[interactions['hbond_bb_sc'] < 0]
    bb_bb_sr = interactions[interactions['hbond_sr_bb'] < 0]
    bb_bb_lr = interactions[interactions['hbond_lr_bb'] < 0]
    return {
        'sc_sc': sc_sc,
        'bb_sc': bb_sc,
        'bb_bb_sr': bb_bb_sr,
        'bb_bb_lr': bb_bb_lr
    }


def interaction_analysis(
    variant: pd.DataFrame,
    wild_type: pd.DataFrame,
    mutated: List[int]
) -> Dict[str, pd.DataFrame]:
    """
    Classify interactions into 6 categories
    :param variant:
        The residue energy breakdown for the variant
    :param wild_type:
        The residue energy breakdown for the wild-type
    :param mutated:
        The mutated residues between the variant and wild-type
    :return:
        A dictionary containing a dataframe of interactions for
        each category
    """
    position = ['resi1', 'resi2']

    position_v = variant[position].values.tolist()
    position_w = wild_type[position].values.tolist()

    change_a = []
    change_b = []
    change_c = []
    change_d = []
    change_e = []
    change_f = []

    for row in variant.iterrows():
        if (
            row[1][position].values.tolist() in position_w and
            (
                row[1]['resi1'] not in mutated and
                row[1]['resi2'] not in mutated
            )
        ):
            change_a.append(subtract_rows(row, wild_type))

        if (
            row[1][position].values.tolist() in position_w and
            (
                row[1]['resi1'] in mutated or
                row[1]['resi2'] in mutated
            )
        ):
            change_b.append(subtract_rows(row, wild_type))

        if (
            row[1][position].values.tolist() not in position_w and
            (
                row[1]['resi1'] not in mutated and
                row[1]['resi2'] not in mutated
            )
        ):
            change_e.append(row[1].values)

        if (
            row[1][position].values.tolist() not in position_w and
            (
                row[1]['resi1'] in mutated or
                row[1]['resi2'] in mutated
            )
        ):
            change_f.append(row[1].values)

    for row in wild_type.iterrows():
        if (
            row[1][position].values.tolist() not in position_v and
            (
                row[1]['resi1'] not in mutated and
                row[1]['resi2'] not in mutated
            )
        ):
            change_c.append(row[1].values)

        if (
            row[1][position].values.tolist() not in position_v and
            (
                row[1]['resi1'] in mutated or
                row[1]['resi2'] in mutated
            )
        ):
            change_d.append(row[1].values)

    change_a = pd.DataFrame(change_a)
    change_b = pd.DataFrame(change_b)
    change_c = pd.DataFrame(change_c)
    change_d = pd.DataFrame(change_d)
    change_e = pd.DataFrame(change_e)
    change_f = pd.DataFrame(change_f)
    dfs = [change_a, change_b, change_c, change_d, change_e, change_f]
    for i in dfs:
        if len(i) > 0:
            i.columns = variant.columns

    return {
        'a': change_a,
        'b': change_b,
        'c': change_c,
        'd': change_d,
        'e': change_e,
        'f': change_f
    }


def subtract_rows(
    row: Tuple[Hashable, pd.Series],
    df: pd.DataFrame
) -> list:
    """
    Subtract the numerical values of corresponding interactions
    :param row:
        The row whose values are being subtracted from
    :param df:
        The dataframe that contains the matching row whose values will
        be used for subtraction
    :return:
        The adjusted values of the row, to be added into a new dataframe
    """
    query = df[
        (df['resi1'] == row[1]['resi1']) &
        (df['resi2'] == row[1]['resi2'])
    ].values
    data = list(row[1][0:6].values)
    data.extend(list(row[1][6:] - query[0, 6:]))
    return data


def buried_hbonds(
    resi_depth: dict[int, float],
    threshold: float,
    hbonds_sc_sc: Dict[str, pd.DataFrame],
    hbonds_bb_sc: Dict[str, pd.DataFrame],
    hbonds_bb_bb_lr: Dict[str, pd.DataFrame],
    hbonds_bb_bb_sr: Dict[str, pd.DataFrame],
    unsatisfied: bool = True,
    **kwargs
) -> None:
    """
    An incomplete attempt to identify buried and unsatisfied
     hydrogen bonds using residue depth
    :param resi_depth:
    :param threshold:
    :param hbonds_sc_sc:
    :param hbonds_bb_sc:
    :param hbonds_bb_bb_lr:
    :param hbonds_bb_bb_sr:
    :param unsatisfied:
    :param kwargs:
    :return:
    """
    frames = {
        f'sc_sc_{"c" if unsatisfied else "e"}':
            hbonds_sc_sc['c' if unsatisfied else 'e'],
        f'sc_sc_{"d" if unsatisfied else "f"}':
            hbonds_sc_sc['d' if unsatisfied else 'f'],
        f'bb_sc_{"c" if unsatisfied else "e"}':
            hbonds_bb_sc['c' if unsatisfied else 'e'],
        f'bb_sc_{"d" if unsatisfied else "f"}':
            hbonds_bb_sc['d' if unsatisfied else 'f'],
        f'bb_bb_sr_{"c" if unsatisfied else "e"}':
            hbonds_bb_bb_sr['c' if unsatisfied else 'e'],
        f'bb_bb_sr_{"d" if unsatisfied else "f"}':
            hbonds_bb_bb_sr['d' if unsatisfied else 'f'],
        f'bb_bb_lr_{"c" if unsatisfied else "e"}':
            hbonds_bb_bb_lr['c' if unsatisfied else 'e'],
        f'bb_bb_lr_{"d" if unsatisfied else "f"}':
            hbonds_bb_bb_lr['d' if unsatisfied else 'f']
    }

    for key, df in frames.items():
        for row in df.iterrows():
            data_1 = resi_depth[row[1]['resi1']]
            data_2 = resi_depth[row[1]['resi2']]
            if data_1 > threshold or data_2 > threshold:
                print(f"{key:<10}", end=' ')
                print(row[1][['resi1', 'resi2', 'total']].values, end=' ')
                print(max([data_1, data_2]))

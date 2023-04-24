import sys
from multiprocessing.pool import Pool
from functools import partial
sys.path.insert(0, f'/media/data/zakaryjd/Python_Projects/Leipzig')
from linux_handles import rosetta_simple


def run(filename: str, nstruct: int) -> None:
    print(f'Started Job for {filename} with nstruct: {nstruct}')
    executable = 'rosetta_scripts.static.linuxgccrelease'

    options = [
        f'-s data/inputs/{filename}.pdb',
        f'-native data/inputs/{filename}.pdb',
        f'-parser:protocol XML/fast_relax.xml',
        f'-out:pdb',
        f'-out:path:pdb data/results',
        f'-nstruct {nstruct}',
        f'-out:file:scorefile data/scores/{filename}.txt'
    ]
    log = rosetta_simple(executable, options)
    with open(f'logs/relax_{filename}.txt', 'w') as file:
        file.write(log)


if __name__ == '__main__':
    """
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Fast relax')
    parser.add_argument('--nstruct', type=int, required=True)
    parser.add_argument('--name', type=str, required=True)
    args = parser.parse_args()
    """
    # t5 d1 = To protMPNN desing at temp 0.1
    # t5 d2 = To protMPNN desing at temp 0.2

    struct_names = [f't6d{i}' for i in range(1, 3)]
    worker = partial(run, nstruct=5)
    pool = Pool(processes=10)
    pool.map(worker, struct_names)

from subprocess import Popen, PIPE


def simple_run(command: list) -> str:
    """
    A wrapper for running linux programs via python
    :param command:
        The command line to be run
    :return:
        The console output
    """
    process = Popen(command, shell=False, stdout=PIPE)
    output = list()
    while True:
        output.append(process.stdout.readline().strip())
        if process.poll() is not None:
            break
    return '\n'.join([x.decode() for x in output])


def rosetta_simple(executable: str, args: list = None) -> str:
    """
    A wrapper for running Rosetta via python
    :param executable:
        The path to the rosetta executable that will be run
    :param args:
        Any command line flags
    :return:
        The console output
    """
    command = [executable]
    if args:
        command += args
    log = simple_run(command)
    return log

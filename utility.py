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

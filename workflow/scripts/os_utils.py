import os


def check_path(file_path):
    """
    Check to see if path exists

    :param file_path: Path to check
    :return: True if path exists, false if not
    """
    if not os.path.exists(file_path):
        return False
    else:
        return True


def remove_file(file_path):
    """
    Delete file if it exists
    :param file_path: Path to file to delete
    """
    if check_path(file_path):
        os.remove(file_path)

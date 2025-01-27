import os

def check_create_dir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
    return
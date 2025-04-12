import os

def file_exists(file):
    if not os.path.isfile(file):
        print(f"Error: File '{file}' does not exist.")
        return False
    return True

def can_read_file(file):
    if not os.access(file, os.R_OK):
        print(f"Error: You do not have read permission for '{file}'.")
        return False
    return True

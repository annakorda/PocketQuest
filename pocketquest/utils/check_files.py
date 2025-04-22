import os
import argparse

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

def positive_float(value):
    fvalue = float(value)
    if fvalue <= 0:
        raise argparse.ArgumentTypeError(f"Distance must be a positive float, got {value}")
    return fvalue

def positive_int(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"Value must be a positive integer, got {value}")
    return ivalue

def probability_float(value):
    fvalue = float(value)
    if not (0 <= fvalue <= 1):
        raise argparse.ArgumentTypeError(f"Probability must be between 0 and 1, got {value}")
    return fvalue
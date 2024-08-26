# pyvibtools/filetypes.py

# Define default file types and check formats for I/O stuff

import os
from ase.io import read
from ase.io.formats import UnknownFileTypeError

class FileFormatChecker:
    def __init__(self):
        self.format_checks = {}

    def register_format(self, format_name, check_functions):
        """
        Register a new file format with one or more check functions.

        Parameters:
        - format_name: The name of the format (e.g., "vibspectrum").
        - check_functions: A single function or a list of functions that take a filename
                           and return True if the file matches the format.
        """
        if not isinstance(check_functions, list):
            check_functions = [check_functions]
        self.format_checks[format_name] = check_functions

    def check_format(self, filename):
        """
        Check the format of a given file by running all registered checks.

        Parameters:
        - filename: The path to the file.

        Returns:
        - The name of the format if a match is found, otherwise None.
        """
        # Check if the file exists
        if not os.path.isfile(filename):
            return None
        # Check format
        for format_name, check_functions in self.format_checks.items():
            if any(check_function(filename) for check_function in check_functions):
                return format_name
        return None

    @staticmethod
    def check_extension(filename, extension):
        """
        Check if a file has a specific extension.
        """
        return filename.lower().endswith(extension.lower())

    @staticmethod
    def check_contains_keyword(filename, keyword):
        """
        Check if a file contains a specific keyword in its contents.
        """
        with open(filename, 'r') as file:
            for line in file:
                if keyword in line:
                    return True
        return False

    @staticmethod
    def check_exact_name(filename, exact_name):
        """
        Check if a file has an exact name (without the path).
        """
        return os.path.basename(filename) == exact_name

    @staticmethod
    def check_is_plain_text_numbers(filename):
        """
        Check if the entire file contains only plain text numbers.
        """
        try:
            with open(filename, 'r') as file:
                for line in file:
                    if not all(FileFormatChecker.is_number(part) for part in line.split()):
                        return False
            return True
        except Exception:
            return False

    @staticmethod
    def is_number(s):
        """
        Helper function to check if a string is a number.
        """
        try:
            float(s)
            return True
        except ValueError:
            return False


def register_default_formats(checker):
    # Register a check for Turbomole vibspectrum files with multiple alternative conditions
    checker.register_format("TM_vibspectrum", [
        lambda filename: (
            FileFormatChecker.check_exact_name(filename, "vibspectrum")
        ),
        lambda filename: (
            FileFormatChecker.check_contains_keyword(filename, "$vibrational spectrum")
        )
    ])

    # Register a check for Turbomole hessian files
    checker.register_format("TM_hessian", [
        lambda filename: FileFormatChecker.check_contains_keyword(filename, "$hessian")
    ])

    # Register a check for plain text files containing only numbers
    checker.register_format("plain", [
        lambda filename: FileFormatChecker.check_is_plain_text_numbers(filename)
    ])
 



def check_ASE_readable(filename):
    try:
        structure = read(filename)
        return True
    except (UnknownFileTypeError, IOError, ValueError) as e:
        return False

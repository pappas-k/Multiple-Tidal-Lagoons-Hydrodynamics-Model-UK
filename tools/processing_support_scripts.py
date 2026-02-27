from thetis import *
import pickle


def output_field_h5(output_directory, field, name):
    """
    Simply outputs a field for further processing
    :param output_directory: output directory
    :param field: Firedrake function field
    :param name: Name of Function to be used in the output file
    :return:
    """
    checkpoint_file = checkpointing.DumbCheckpoint(output_directory + "/" + name)
    checkpoint_file.store(field)
    checkpoint_file.close()
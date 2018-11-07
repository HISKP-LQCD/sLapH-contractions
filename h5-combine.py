#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import h5py


def main():
    options = _parse_args()

    # Open the output file.
    with h5py.File(options.output, 'w') as file_output:
        first_iteration = True

        # Iterate through all inputs.
        for ip_id, input_path in enumerate(options.input):
            print(input_path)
            # Open the input file.
            with h5py.File(input_path, 'r') as file_input:
                # On further files, assert that the available datasets match exactly.
                # Anything else is an error.
                if not first_iteration:
                    keys_input = list(file_input.keys())
                    keys_output = list(file_output.keys())
                    assert keys_input == keys_output, \
                            "Contained datasets must be exactly the same for all files"

                # Retrieve a list of all available data sets.
                for ds_id, datasetname in enumerate(file_input.keys()):
                    # If we are on the first file, create the needed datasets as groups
                    # in the output.
                    if first_iteration:
                        file_output.create_group(datasetname)


                    # Copy the dataset over into the according group.
                    file_input.copy(datasetname, file_output['/{}'.format(datasetname)])
                    file_output.move('/{}/{}'.format(datasetname, datasetname),
                                     '/{}/{:04d}'.format(datasetname, ip_id))

            file_output.flush()

            first_iteration = False


def _parse_args():
    '''
    Parses the command line arguments.

    :return: Namespace with arguments.
    :rtype: Namespace
    '''
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input', nargs='+')
    parser.add_argument('output')

    options = parser.parse_args()

    return options


if __name__ == '__main__':
    main()

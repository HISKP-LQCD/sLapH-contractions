#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright © 2019 Martin Ueding <mu@martin-ueding.de>

import argparse
import itertools
import json
import pprint
import re


def diagram_name(correlator_name, pattern=re.compile(r'^([^_]+)')):
    m = pattern.search(correlator_name)
    if m:
        return m.group(1)
    else:
        return None


def flatten(forest):
    '''
    Taken from a comment to https://stackoverflow.com/a/952952/653152.
    '''
    return [leaf for tree in forest for leaf in tree]


def extract_single(path):
    with open(path) as f:
        data = json.load(f)

    correlators = extract_recursive(data)

    while len(correlators) > 0 and isinstance(correlators[0], list):
        correlators = flatten(correlators)

    return correlators


def extract_recursive(data):
    if isinstance(data, list):
        result = list(map(extract_recursive, data))
    elif isinstance(data, dict):
        if 'datasetname' in data:
            result = data['datasetname']
        else:
            result = list(map(extract_recursive, data.values()))

    return result


def main():
    options = _parse_args()

    print('Processing prescriptions …')
    correlators = set()
    for prescription in options.prescription:
        print(' ', prescription)
        result = extract_single(prescription)
        correlators |= set(result)

    print()
    print('Found {} unique correlators.'.format(len(correlators)))

    grouped = {}
    for group, values in itertools.groupby(sorted(correlators), diagram_name):
        grouped[group] = list(values)

    with open(options.o, 'w') as f:
        json.dump(grouped, f, indent=2, sort_keys=True)


def _parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('prescription', nargs='+')
    parser.add_argument('-o', required=True)
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    main()

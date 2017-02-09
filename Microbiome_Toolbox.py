#!/usr/bin/python3
"""
This is the start point of a series (developing) of microbiome analysis tools.
For Vinatzer Lab only.
"""

# IMPORT
import argparse
import sys
import os

# FUNCTIONS
def get_parsed_args():
    parser = argparse.ArgumentParser(
        description="Microbiome analysis tool box"
    )
    parser.add_argument("-i", dest="input", help="Directory of the QIIME workspace")
    args = parser.parse_args()
    return args

# MAIN
if __name__ == "__main__":
    tool = sys.argv[1]
    tools = {"origin_track":"python origin_track.py {0}"}
    argv = sys.argv[2:]
    args = get_parsed_args()
    input_dir = args.input
    os.system(tools[tool])
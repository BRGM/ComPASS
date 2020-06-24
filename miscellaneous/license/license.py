# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 15:34:29 2017

@author: lopez
"""

import os
import sys
import glob
from optparse import OptionParser

parser = OptionParser()
parser.add_option(
    "-a",
    "--add",
    action="store_true",
    dest="add_all",
    default=False,
    help="add licence disclaimer to all registered files",
)
parser.add_option(
    "-c",
    "--check",
    action="store_true",
    dest="check_all",
    default=False,
    help="check licence disclaimer of all registered files",
)
options, args = parser.parse_args()


license_file = "license-disclaimer.txt"

with open(license_file) as f:
    license_lines = f.readlines()


def get_extension(filename):
    root, ext = os.path.splitext(filename)
    return ext.lower()


def header_size(filename):
    ext = get_extension(filename)
    preserved_header = 0
    with open(filename) as f:
        loop = True
        while loop:
            loop = False
            line = f.readline()
            if ext == ".py":
                if line.startswith("#!"):
                    preserved_header += 1
                    loop = True
                if line.startswith("# -*-"):
                    preserved_header += 1
                    loop = True
    return preserved_header


def license_as_comment(filename):
    ext = get_extension(filename)
    comment_symbol = {
        ".c": "//",
        ".h": "//",
        ".cc": "//",
        ".cpp": "//",
        ".f90": "!",
        ".py": "#",
    }[ext] + " "
    return [comment_symbol + l.strip() for l in license_lines]


def check_license(filename):
    comment = license_as_comment(filename)
    with open(filename) as f:
        preserved_header = header_size(filename)
        for _ in range(preserved_header):
            f.readline()
        for l in comment:
            if not f.readline().strip() == l.strip():
                return False
    return True


def add_license(filename):
    if check_license(filename):
        return False
    preserved_header = header_size(filename)
    if preserved_header > 0:
        print(preserved_header, filename)
    with open(filename) as f:
        lines = f.readlines()
        header = lines[:preserved_header]
        body = lines[preserved_header:]
    comment = license_as_comment(filename)
    with open(filename, "w") as f:
        for line in header + comment:
            print(line.strip(), file=f)
        print(file=f)
        for line in body:
            print(line.rstrip(), file=f)
    return True


check_dirs = ["../../src", "../../python"]


def perform_action_on_files(action):
    for dirname in check_dirs:
        for filename in glob.iglob(dirname + "/**/*", recursive=True):
            root, ext = os.path.splitext(filename)
            ext = ext.lower()
            if ext in [".cc", ".cpp", ".h", ".py", ".f90"]:
                action(filename)


if options.add_all:
    print("Adding lincense disclaimer...")

    def add_license_action(filename):
        if add_license(filename):
            print("processed:", filename)

    perform_action_on_files(add_license_action)

if options.check_all:
    print("Checking lincense disclaimer...")

    def check_license_action(filename):
        if not check_license(filename):
            print(filename, "is missing license dislaimer!")
            sys.exit(-1)

    perform_action_on_files(check_license_action)

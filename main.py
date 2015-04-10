#!/usr/bin/env python

import sys, os, re
import argparse
import imp
topology = imp.load_source("topo", "./src/topology.py")

def main():
    topology.main()

if __name__ == "__main__":
    main()

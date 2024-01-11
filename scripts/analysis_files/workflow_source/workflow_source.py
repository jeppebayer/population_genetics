#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def test():
    print(f'{glob.glob(f"{os.path.dirname(os.path.realpath(__file__))}/scripts/popoolation2*")[0]}')
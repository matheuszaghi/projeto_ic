# -*- coding: utf-8 -*-
"""

Created on 09/10/17

@author: Carlos Eduardo Barbosa

Program to handle workspaces in different computers in the collaboration

"""
from __future__ import division, print_function

import os
import getpass

username = getpass.getuser()

if username == "kadu":
    home_dir = "/home/kadu/Dropbox/matheusic"
    data_dir = os.path.join(home_dir, "data")


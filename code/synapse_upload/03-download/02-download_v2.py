#   The PEC dataset has been updated since the original download. This script
#   pulls all the new files from the staging area for this dataset, since it
#   will take a bit of time for the staged files to be moved to their final
#   location.

import sys
import synapseclient
import synapseutils
import yaml
import pyhere
import pandas as pd
from pathlib import Path

credential_path = '/users/neagles/synapse_credentials.yaml'
download_dir = pyhere.here('raw-data', 'psychENCODE', 'version2')
synapse_id = 'syn30106435'

Path(download_dir).mkdir(parents=True, exist_ok=True)

with open(credential_path, 'r') as f:
    cred_dict = yaml.safe_load(f)

#   Log in to synapse
syn = synapseclient.Synapse()
syn.login(authToken=cred_dict['token'])

files = synapseutils.syncFromSynapse(syn, synapse_id, path = download_dir)


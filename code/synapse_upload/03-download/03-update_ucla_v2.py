#   Some h5ad files for the UCLA dataset were modified/fixed in a separate
#   location on synapse than the other datasets. This script pulls the fixed
#   files for the UCLA dataset and places them with the other PEC files

import sys
import synapseclient
import synapseutils
import yaml
import pyhere
import pandas as pd
from pathlib import Path

credential_path = '/users/neagles/synapse_credentials.yaml'
download_dir = pyhere.here('raw-data', 'psychENCODE', 'version2', 'UCLA-ASD')
synapse_id = 'syn30106539'

Path(download_dir).mkdir(parents=True, exist_ok=True)

with open(credential_path, 'r') as f:
    cred_dict = yaml.safe_load(f)

#   Log in to synapse
syn = synapseclient.Synapse()
syn.login(authToken=cred_dict['token'])

files = synapseutils.syncFromSynapse(syn, synapse_id, path = download_dir)

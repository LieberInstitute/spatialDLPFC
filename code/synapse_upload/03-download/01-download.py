import sys
import synapseclient
import synapseutils
import yaml
import pyhere
import pandas as pd

credential_path = '/users/neagles/synapse_credentials.yaml'
table_path = pyhere.here('raw-data', 'psychENCODE', 'version1', 'exported_table.csv')
download_dir = pyhere.here('raw-data', 'psychENCODE', 'version1')

with open(credential_path, 'r') as f:
    cred_dict = yaml.safe_load(f)

#   Log in to synapse
syn = synapseclient.Synapse()
syn.login(authToken=cred_dict['token'])

#   Read in table of files for the study and download
exported_table = pd.read_csv(table_path)
for x in exported_table.id:
    syn.get(x, downloadLocation = download_dir)

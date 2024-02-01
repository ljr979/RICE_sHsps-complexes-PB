
import os
import os, re
import pandas as pd
import numpy as np
import json
import networkx
import pandas as pd
import tarfile
import gzip, shutil, os, re
import time
from loguru import logger

logger.info(f'Import OK')

resource_folder = 'resources/bioinformatics_databases/'


def gz_unzipper(filename, input_path=resource_folder, output_path=resource_folder):
    with gzip.open(f'{input_path}{filename}.gz', 'rb') as f_in:
        with open(f'{output_path}{filename}', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def tar_file_to_folder(input_path, output_path):
    tar = tarfile.open(f'{input_path}', 'r')
    tar.extractall(f'{output_path}')
    tar.close()



    uniprot_db = create_uniprot_db(f'{resource_folder}{tax_id}.tab.gz', genes)
    uniprot_map = create_uniprot_map(f'{resource_folder}{tax_id}_idmapping.tab.gz', genes)
    merged = pd.merge(uniprot_db, uniprot_map, left_on='Entry', right_on='UniProtKB-AC', how='outer')
    merged['name'] = merged['Entry name'].str.split('_').str[0]
    if reviewed:
        merged = merged[merged['Status'] == 'reviewed']
    return merged
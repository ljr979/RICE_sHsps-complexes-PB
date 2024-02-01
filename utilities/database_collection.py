import gzip, shutil, os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import OrderedDict
import json
from contextlib import closing
import urllib.request as request
import time
from utilities.database_map_and_filter import gz_unzipper
from loguru import logger


logger.info('Import OK')


def download_resources(filename, url, resource_folder):
    """
    Worker function to download and save file from URL.
    
    inputs
    ======
    filename: (str) name of output file (including extension)
    url: (str) complete location of file to be downloaded
    output_path: (str) relative or complete path to directory where folder will be saved.

    returns:
    ======
    None

    """
    if not os.path.exists(resource_folder):
        os.makedirs(resource_folder)

    try:
        with closing(request.urlopen(url)) as r:
            with open(f'{resource_folder}{filename}', 'wb') as f:
                shutil.copyfileobj(r, f)
        logger.info(f'Downloaded {filename}')
    except:
        logger.info(f'Downloaded failed for {filename}.')


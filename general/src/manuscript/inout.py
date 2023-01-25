import os
import socket

import pandas as pd

from pathlib import Path


def gustav_path():
    """
    Get path of gustav
    """
    
    p = os.path.join(
        _settings_from_file()['gustav_path'],
        'src'
    )
    return p


def get_internal_path(extension=None):
    '''
    Returns subfolder within input folder of pstm.
    Input:
        extension   str, optional, subfolder
    Output:
        outpath     str, folder within internal part of pstm
    '''

    if extension is not None:
        extension = _adjust_to_current_file_separator(extension)
    
    settings = _settings_from_file()
    
    
    base_folder = settings['internal_path']
    
    if not os.path.exists(base_folder):
        raise AssertionError(
            'Could not find input folder {}. Please ensure that specified.'.format(
                base_folder))

    if extension is not None:
        extension = str.replace(extension, '\\', os.path.sep)
        extension = str.replace(extension, '/', os.path.sep)

        outpath = os.path.join(base_folder, extension)
    else:
        outpath = base_folder


    return outpath


def ensure_presence_of_directory(directory_path=None, ):
    '''
    Ensure that the directory of exists. Creates dictionary with cognate
    name in case that directory does not exist. Anticipates that files have
    extensions separated by '.' symbol (can not create directory with . in
    its name); If file does not have an extension separated by '.' a folder
    will with its filname will be created, a behavior that can be avoided
    by calling os.path.dirname prior this function.
    Input:
        directory_path      str; Name of a directory or the full path of a file
    '''
    if directory_path is None:
        raise ValueError('No input specfied for ensure_presence_of_directory')

    directory_path_n, ext = os.path.split(directory_path)

    if '.' in ext:
        directory_path = directory_path_n

    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

        
        
def _settings_from_file():
    """
    Loads settings from settings file.
    
    Output:
    settings       dict
    
    
    """

    home = str(Path.home())


    path_to_settings = os.path.join(
        str(Path.home()), 'Documents', 'data_paths')

    if not os.path.exists(path_to_settings):
        raise EnvironmentError(rf"""
            Could not find directory reserved for settings:
            {path_to_settings}
        """)

    path_to_settings = os.path.join(
        path_to_settings,
        'promising_genes.csv'
    )

    if not os.path.exists(path_to_settings):
        raise EnvironmentError(rf"""
            Could not find promising_genes.csv file:
            {path_to_settings}
            This file needs to be UTF-8 formatted
            csv file with two columns: key, value
            Also see readme of repository for
            further guidance.
        """)
        
        

    settings = pd.read_csv(
        path_to_settings
    )
    
    
    if not all(settings.columns[:2] == ['key', 'value']):
        raise EnvironmentError(rf"""
            promising_genes.csv must start with the two
            columns: key and value.
            
            {path_to_settings}
        """)
        
        
    settings = settings.drop_duplicates()
    
    if any(settings['key'].duplicated()):
        raise EnvironmentError(rf"""
            At least one key within pstm.csv is
            duplicated and therefore ambiguous
            
            {path_to_settings}
        """)
    
    
    settings = settings.set_index(
        'key', 
        verify_integrity=True
    )['value'].to_dict()
    
    return settings

def _adjust_to_current_file_separator(x):
    '''
    Replaces backslashes and forward slashes
    by file separtor used on current operating
    system.
    '''
    x = x.replace('\\', os.path.sep)
    x = x.replace('/', os.path.sep)

    return x
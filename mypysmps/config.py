# -*- coding: utf-8 -*-
#################
import os
import imp
import traceback
import warnings
#################

"""
mypysmps.config
============

PySMPS configuration.
    load_config
    get_metadata
    get_field_name
    get_fillvalue

Created on Fri Jul 10 09:54 2020

@author: flovan / fvanden

Revision history:   10.07.2020 - Created
                    
"""
## -------------------------------------------------------------------------- ##


# the path to the default configuration file
_dirname = os.path.dirname(__file__)
_DEFAULT_CONFIG_FILE = os.path.join(_dirname, 'default_config.py')

def load_config(filename = None):
    """
    Load a PySMPS configuration from a config file.
    
    The default values for a number of PySMPS parameters and metadata is
    controlled by a single Python configuration file. An example of this 
    file can be found in the PySMPS source directory named 
    default_config.py. 
    
    The recommended method for changing these defaults is for users to
    copy this file, rename it and point either this scipt to the new
    configuration file or execute load_config with the new configuration
    file as an input parameter.
    
    NOTE: this method and script is adapted from .default_config.py
    available at https://arm-doe.github.io/pyart/

    Parameters
    ----------
    filename : str
        Filename of configuration file.  If None the default configuration
        file is loaded from the Py-ART source code directory.

    """
    if filename is None:
        filename = _DEFAULT_CONFIG_FILE
        
    # these are private since they should not be accessed by users or other
    # modules, use the get_ functions.
    global _DEFAULT_FIELD_NAMES
    global _DEFAULT_METADATA
    global _DEFAULT_FIELD_COLORMAP
    global _DEFAULT_FIELD_LIMITS
    global _FILL_VALUE
    global _DEFAULT_VARIABLES
    global _FIELD_MAPPING
    global _DEFAULT_FIGURE_SETTINGS
    global _INSTRUMENT_SETTINGS
    global _INSTRUMENT_HEADERS
    global _CONVERSIONS
    
    #global _FILE_SPECIFIC_METADATA
    #global _FIELD_MAPPINGS
    
    cfile = imp.load_source('metadata_config', filename)
    _DEFAULT_METADATA = cfile.DEFAULT_METADATA
    #_FILE_SPECIFIC_METADATA = cfile.FILE_SPECIFIC_METADATA
    #_FIELD_MAPPINGS = cfile.FIELD_MAPPINGS
    _FILL_VALUE = cfile.FILL_VALUE
    _DEFAULT_FIELD_NAMES = cfile.DEFAULT_FIELD_NAMES
    _DEFAULT_FIELD_COLORMAP = cfile.DEFAULT_FIELD_COLORMAP
    _DEFAULT_FIELD_LIMITS = cfile.DEFAULT_FIELD_LIMITS
    _DEFAULT_VARIABLES = cfile.DEFAULT_VARIABLES
    _FIELD_MAPPING = cfile.FIELD_MAPPING
    _DEFAULT_FIGURE_SETTINGS = cfile.DEFAULT_FIGURE_SETTINGS
    _INSTRUMENT_SETTINGS = cfile.INSTRUMENT_SETTINGS
    _INSTRUMENT_HEADERS = cfile.INSTRUMENT_HEADERS
    _CONVERSIONS = cfile.CONVERSIONS
    return

# load the configuration from the enviromental parameter if it is set
# if the load fails issue a warning and load the default config.
_config_file = os.environ.get('PYART_CONFIG')
if _config_file is None:
    load_config(_DEFAULT_CONFIG_FILE)
else:
    try:
        load_config(_config_file)
    except:
        msg = ("\nLoading configuration from PYART_CONFIG enviromental "
               "variable failed:"
               "\n--- START IGNORED TRACEBACK --- \n" +
               traceback.format_exc() +
               "\n --- END IGNORED TRACEBACK ---"
               "\nLoading default configuration")
        warnings.warn(msg)
        load_config(_DEFAULT_CONFIG_FILE)
        
def get_metadata(p):
    """
    Return a dictionary of metadata for a given parameter, p.

    An empty dictionary will be returned in no metadata dictionary exists for
    parameter p.
    """
    if p in _DEFAULT_METADATA:
        return _DEFAULT_METADATA[p].copy()
    else:
        return {}
    
def get_fillvalue():
    """
    Return the current fill value.
    """
    return _FILL_VALUE


def get_field_name(field):
    """
    Return the field name from the configuration file for a given field.
    """
    return str(_DEFAULT_FIELD_NAMES[field])


def get_field_colormap(field):
    """
    Return the colormap name from the configuration file for a field name.
    """
    if field in _DEFAULT_FIELD_COLORMAP:
        return _DEFAULT_FIELD_COLORMAP[field]
    else:
        import matplotlib.cm
        return matplotlib.cm.get_cmap().name
    
def get_field_limits(field):
    """
    Return the data limits from the configuration file for a given field

    Parameters
    ----------
    field: str
        Field name.

    Returns
    -------
    vmin, vmax: 2-tuplet of float
        Minimun and Maximun theorical value for field, if field is not
        in the configuration file returns (None, None).
    """
    if field in _DEFAULT_FIELD_LIMITS:
        limits = _DEFAULT_FIELD_LIMITS[field]
        return limits
    else:
        return None, None
    
def get_figure_settings(figtype):
    """
    Return figure settings from the configuration file for a given figure type
    
    Parameters
    ----------
    figtype : str
        type of figure
        
    Returns
    -------
    figdict: dict
        dictionary with figure settings
    """
    if figtype in _DEFAULT_FIGURE_SETTINGS:
        figdict = _DEFAULT_FIGURE_SETTINGS[figtype]
        return figdict
    else:
        return None
    
def get_instrument_header(instrument):
    """
    Return instrument header from the configuration file for a given instrument
    
    Parameters
    ----------
    instrument : str
        instrument name
        
    Returns
    -------
    header : list
        list of instrument header names
    """
    if instrument in _INSTRUMENT_HEADERS:
        header = _INSTRUMENT_HEADERS[instrument]
        return header
    else:
        return None
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#################
import os
import gzip

from ..core.met import MET
from .csv_read import read_aim_csv, read_opc_csv, read_csv
from .txt_read import read_aim_txt, read_opc_txt, read_grimm_txt
from .nc_read import read_mpl
#################
"""
mypysmps.io.read
================

Automatic reading of files by detecting format:
    read
    determine_filetype
    
Created on Thu Jul 9 14:37 2020

@author: flovan / fvanden

Revision history:   09.07.2020 - Created
                    20.07.2020 - filetype added to allow for different
                    file organisations
                    30.09.2020 - metdata added
                    11.08.2022 - read_grimm_txt added



"""
## -------------------------------------------------------------------------- ##


def read(filename, fileorg = 'AIM', **kwargs):
    """
    Read a SMPS file and return a SMPS object
    
    Parameters
    ----------
    filename : str
        path and name of file to read
    
    fileorg : str
        refers to the organisation of the file
        
        
        
    Returns
    -------
    smps : smps
        mysmps.core.smps object
        
    """
    
    filetype = determine_filetype(filename)
    
    # Gzip, uncompress and see if we can determine the type
    if filetype == 'GZ':
        gzfile = gzip.open(filename, 'rb')
        try:
            smps = read(gzfile, **kwargs)
        except:
            raise ValueError(
                'Gzip file cannot be read compressed, '
                'uncompress and try again')
        finally:
            gzfile.close()
        return smps
    
    # CSV
    if filetype == 'CSV':
        if fileorg == 'AIM':
            return read_aim_csv(filename, fileorg = fileorg, **kwargs)
        elif fileorg == 'OPC':
            return read_opc_csv(filename, fileorg = fileorg, **kwargs)
        elif fileorg == 'MET':
            vardict = read_csv(filename, fileorg = fileorg, **kwargs)
            MET(**vardict)
            
        else:
            try:
                vardict = read_csv(filename, fileorg = fileorg, **kwargs)
                return vardict
            except:
                raise TypeError('Unknown or unsupported file organisation: ' + fileorg)
    
    # TXT
    if filetype in ['TXT','.gtm','.g6s']:
        if fileorg == 'AIM':
            return read_aim_txt(filename, fileorg = fileorg, **kwargs)
        if fileorg == 'OPC':
            return read_opc_txt(filename, fileorg = fileorg, **kwargs)
        elif fileorg == 'Grimm':
            return read_grimm_txt(filename, fileorg = fileorg, **kwargs)
        else:
            raise TypeError('Unknown or unsupported file organisation: ' + fileorg)
            
    if filetype == "NETCDF3" or filetype == "NETCDF4":
        if fileorg == 'MPL':
            return read_mpl(filename, fileorg = fileorg, **kwargs)

        
    raise TypeError('Unknown or unsupported file format: ' + filetype)
    
    
def determine_filetype(filename):
    """
    Return the filetype of a given file by examining the first few bytes.
    
    Adapted from pyart.io.auto_read.py script by : https://arm-doe.github.io/pyart/ 

    The following filetypes are detected:

    * 'csv'
    * 'txt'
    * 'excel'
    * 'NETCDF3'
    * 'NETCDF4'
    * 'HDF4'
    * 'gzip'

    Parameters
    ----------
    filename : str
        Name of file to examine.

    Returns
    -------
    filetype : str
        Type of file.
        
    """
    
    # read the first 12 bytes from the file
    try:
        f = open(filename, 'rb')
        begin = f.read(12)
        f.close()
    except TypeError:
        f = filename
        begin = f.read(12)
        f.seek(-12, 1)
        
    # CSV - no file signature as far as I know
    csv_signature = "csv"
    if filename[-3:] == csv_signature:
        return "CSV"
    
    # txt
    txt_signature = "txt"
    if filename[-3:] == txt_signature:
        return "TXT"
    
    # txt
    txt_signature = "TXT"
    if filename[-3:] == txt_signature:
        return "TXT"
        
    # xlsx
    xlsx_signature = b'PK\x03\x04\x14\x00\x08\x08\x08\x00ss'
    if begin == xlsx_signature:
        return "XLSX"
        
    # NetCDF3, read with read_cfradial
    if begin[:3] == b"CDF":
        return "NETCDF3"

    # NetCDF4, read with read_cfradial, contained in a HDF5 container
    # HDF5 format signature from HDF5 specification documentation
    hdf5_signature = b'\x89\x48\x44\x46\x0d\x0a\x1a\x0a'
    if begin[:8] == hdf5_signature:
        return "NETCDF4"
    
    # HDF4 file
    # HDF4 format signature from HDF4 specification documentation
    hdf4_signature = b'\x0e\x03\x13\x01'
    if begin[:4] == hdf4_signature:
        return "HDF4"
    
    # gzip filetype
    gzip_signature = b'\x1f\x8b'
    if begin[:2] == gzip_signature:
        return 'GZ'
    
    # zip filetype
    zip_signature = b'PK\x03\x04\x14\x00\x08\x00\x08\x00\x84y'
    if begin == zip_signature:
        return 'ZIP'
    
    # try os
    try:
        return os.path.splitext(filename)[1]
    except:
        # Cannot determine filetype
        return "UNKNOWN"
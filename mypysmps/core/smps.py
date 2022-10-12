# -*- coding: utf-8 -*-
#################
import warnings
import datetime as dt
import numpy as np
import math

from ..config import get_metadata
from ..util.timetransform import TimeTransform
tt = TimeTransform()
#################

"""
mysmps.core.smps
================

A general particle sizer instrument class
    ParticleSizer
    SMPS

Created on Thu Jul 9 14:37 2020

@author: flovan / fvanden

Revision history:   09.07.2020 - Created
                    10.07.2020 - Addittdions to ParticleSizer class
                    20.07.2020 - Added self.__dict__.update(kwargs)
                        for open ended kwargs addition
                    13.08.2020 - Added Empty
                    18.08.2020 - Added createTimeDate and findSample functions
                    13.11.2020 - Extended createTimeDate to also create a time
                        and date list from a datetime list


"""
## -------------------------------------------------------------------------- ##


class ParticleSizer(object):
    """
    A class storing particle sizer data

    Parameters - required
    ---------------------
    time : dict
        Time of each sample
    date: dict
        Date of each sample
    sample: dict
        Sample number
    data : dict
        Data fields of sample measurements
    diameter : dict
        Diameters for data fields

    Parameters - optional
    ---------------------
    temperature : dict
        Temperature of each sample
    pressure : dict
        Pressure for each sample
    relative_humidity : dict
        Relative humidity for each sample
    mean_free_path : dict
        Mean Free Path for each sample
    viscosity : dict
        Gas viscosity for each sample
    scan_time : dict
        Scan duration
    retrace_time : dict
        Retrace duration
    scan_resolution : dict
        Scan resolution
    scans_per_sample : dict
        Number of scans per sample
    aerosol_flow : dict
        Aerosol flow in instrument
    lower_size : dict
        Lower size limit in instrument
    upper_size : dict
        Upper size limit in instrument
    density : dict
        Flow density
    td+05 : dict
        Aerosol delay time
    tf : dict
        Time for aerosol to flow through the sample column of the classifier
    D50 : dict
        Cut point diameter of the impactor.
    median : dict
        Median value of sample
    mean : dict
        Mean value of sample
    geo_mean : dict
        Geometric mean value of sample ?
    mode : dict
        Mode of sample
    geo_std_dev : dict
        Geometric standard deviation of sample?
    total_concentration : dict
        Total concentration of samples
    title : dict
        Title of measurement
    user_name : dict
        User name of instrument operator
    sample_id : dict
        ID of sample
    instrument_id : dict
        ID of instrument
    lab_id : dict
        ID of laboratory
    leak_test_rate : dict
        Leak test and leakage rate
    instrument_errors : dict
        Errors reported by instrument
    comment : dict
        Sample comments

    metadata : dict
        Metadata describing the instrument and data
    pddata : pandas table
        data organised in pandas table

    #latitude : dict
        Latitude of the instrument
    #longitude : dict
        Longitude of the instrument
    #altitude : dict
        Altitude of the instrument


    """
    ## ------------------------------------------------------------------ ##
    ## Constructors/Destructors                                           ##
    ## ------------------------------------------------------------------ ##
    def __init__(self, time, sample, data, diameter, metadata, header, **kwargs):

        if time['standard_name'] == 'time':
            self.time = time
        elif time['standard_name'] == 'datetime':
            self.datetime = time

        #self.date = date
        self.sample = sample
        self.data = data
        self.diameter = diameter
        self.metadata = metadata
        #self.pddata = pddata
        self.header = header

        self.__dict__.update(kwargs)

        # run this automatically when a ParticleSizer instance is created:
        try:
            self.createTimeDate()
        except AttributeError:
            pass

    def __del__(self):
        pass


    ## ------------------------------------------------------------------ ##
    ## Methods                                                            ##
    ## ------------------------------------------------------------------ ##

    # private:
    def checkAttr(self, attr):
        """
        Checks if an attribute is in the instance

        Parameters
        ----------
        attr : string
            attribute for which to check if it is present in the instance

        Returns
        -------
        bool :
            True if present, False if not present

        """
        if hasattr(self, attr):
            return True
        else:
            return False

    def createTimeDate(self, output = False, **kwargs):
        """
        Creates a dict which combines date and time information (if these are
        separate). Creates a time and date dict out of a datetime dict

        Parameters
        ----------
        output : bool
            if True, the list is returned, otherwise it is stored
            in the class instance

        outformat : str ot list of str
            (optional) date time format for output list(s)
            i.e. for datetime: '%d.%m.%Y %H:%M:%S
            for date and time: ['%d.%m.%Y','%H:%M:%S']
            default the units from the default config file are used

        Returns
        -------
        list : list of str
            date and time combined in a single time format
        """

        if self.checkAttr('time'):

            tdlist = []
            timeformat = self.date['units'] + ' ' + self.time['units']
            self.datetime = get_metadata('datetime')
            outformat = kwargs.get("outformat",self.datetime['units'])

            for i in range(len(self.date['data'])):
                nd = self.date['data'][i] + ' ' + self.time['data'][i]
                tdlist.append(dt.datetime.strftime(dt.datetime.strptime(nd,timeformat), outformat))


            self.datetime['units'] = outformat
            self.datetime['data'] = tdlist

            if output:
                return tdlist

        elif self.checkAttr('datetime'):

            tlist = []
            dlist = []
            timeformat = self.datetime['units']
            self.time = get_metadata('time')
            self.date = get_metadata('date')
            outformatD = kwargs.get("outformat",[self.date['units'], self.time['units']])[0]
            outformatT = kwargs.get("outformat",[self.date['units'], self.time['units']])[1]

            for i in range(len(self.datetime['data'])):
                nd = dt.datetime.strptime(self.datetime['data'][i],timeformat)
                tlist.append(dt.datetime.strftime(dt.datetime(2020,1,1,nd.hour,nd.minute,nd.second), outformatT))
                dlist.append(dt.datetime.strftime(dt.datetime(nd.year, nd.month, nd.day), outformatD))

            self.time['units'] = outformatT
            self.time['data'] = tlist
            self.date['units'] = outformatD
            self.date['data'] = dlist


    # public:

    def findSample(self,sample, **kwargs):
        """
        Finds the sample index or time

        Parameters
        ----------
        sample : int or str
            either an integer (sample number) or a string (date)

        dtformat : str
            date time format for output

        Returns
        -------
        loc : str or int
            depending on the type of input: the date of the sample
            or the index of the date (or closest date)
        """
        dtformat = kwargs.get("dtformat",self.datetime['units'])
        if isinstance(sample, int):
            # given sample is an integer,
            # the date and time for the sample
            # will be returned
            inputformat = self.date['units'] + ' ' + self.time['units']
            self.createTimeDate(outputformat = [self.date['units'], self.time['units']])
            return dt.datetime.strftime(dt.datetime.strptime(self.date['data'][sample] + ' ' + self.time['data'][sample],inputformat), dtformat)
        elif isinstance(sample, str):
            # given sample is a string,
            # the indice for the given time
            # will be returned

            self.createTimeDate(outformat = dtformat)
            idx, date = tt.findNearestDate(self.datetime['data'], sample)
            return idx
        else:
            print( ("%s is not a relevant input format"%(type)) )
            pass


    def add_field(self, fieldname, data, method = None, metadata = None):
        """
        Adds a new field to ParticleSizer object

        Parameters
        ----------
        fieldname : str
            name of the field to add

        data : list or numpy.ma.core.MaskedArray
            a warning is emitted if the data does not have
            the same dimensions as the data present in the
            ParticleSizer object

        method : str
            in case of filtered or adapted data, this allows
            to describe the technique applied if the user
            does not wish to add a whole metadata dictionary

        metadata : dict
            a dictionary with metadata for the added field,
            if None, a dictionary from the default config file
            is loaded

        """
        # check if field already exists in object

        if fieldname in self.data['variables']:
            warnings.warn( ("field name %s already exists, field is overwritten...")%(fieldname) )

        # get metadata for field

        if metadata is None:
            metadata = get_metadata(fieldname)
            if metadata == {}:
                metadata = {'units': '-', 'standard_name': fieldname, 'axis': '-', 'valid_min': None, 'valid_max': None, 'comment': None}

        if method is not None:
            metadata['method'] = method

        # check if data has the right dimensions

        objectshape = self.data[self.data['variables'][0]]['data'].shape
        datashape = data.shape

        if datashape == objectshape:
            pass
        elif (datashape[1], datashape[0]) == objectshape:
            # try to match shape
            warnings.warn("field data is translated to match existing field shape")
            data = data.T
        elif datashape[0] == objectshape[1]:
            warnings.warn("only one dimension matches existing field shape: these dimensions are matched")
            data = data.T
        elif datashape[1] == objectshape[0]:
            warnings.warn("only one dimension matches existing field shape: these dimensions are matched")
            data = data.T
        else:
            warnings.warn("field data does not match existing field shape")

        metadata['data'] = data

        # add new field to object data

        self.data[fieldname] = metadata
        self.data['variables'].append(fieldname)

    def create_bins(self,):
        """
        Creates bins for histogram plot

        Parameters
        ----------

        See Also
        --------

        """

        if self.instrument_type == 'SMPS' or self.instrument_type == 'SMPS_concatenated':
            # create empty matrix
            bins = np.empty((np.asarray(self.diameter['data']).shape[0], 3))

            # obtain number of channels per decade
            cpd = float(self.metadata['Channels/Decade'])

            # fill bins with NaNs so that bins fail if method fails
            bins.fill(np.NaN)

            # fill bins with diameters
            bins[:, 1] = np.asarray(self.diameter['data'])
            bins[0, 0] = self.diameter['valid_min']
            bins[-1, -1] = self.diameter['valid_max']


            for i in range(bins.shape[0] - 1):
                bins[i, 2] = round(math.pow(10, np.log10(bins[i, 0]) + 1./cpd), 4)
                bins[i+1, 0] = bins[i, 2]


            # take only bin widths
            dwidths = bins[:,2] - bins[:,0]

            # keep only a list of bins
            binmatrix = bins
            bins = bins[:,1]

        elif self.instrument_type ==  'OPC_concatenated':
            bins = np.append(self.diameter['data'], np.asarray([40.0])) # very dirty fix.. TODO: something about this..
            dwidths = np.diff(bins)
            bins = bins[:-1]
            binmatrix = None

        else: #if self.instrument_type == 'OPC':
            bins = self.diameter['data']
            dwidths = np.diff(bins)
            bins = bins[:-1]
            binmatrix = None

        return bins, dwidths, binmatrix



class SMPS(ParticleSizer):
    """
    A class storing SMPS particle sizer data.

    Parameters
    ----------
    inherited from ParticleSizer class

    sheath_flow : dict
        Sheath Flow in instrument
    bypass_flow : dict
        Bypass flow in instrument
    low_voltage : dict
        Low voltage value in instrument
    high_voltage : dict
        High voltage value in instrument
    neutralizer_status : dict
        Status of the neutralizer

    See Also
    --------
    mysmps.core.read.ParticleSizer

    """
    ## ------------------------------------------------------------------ ##
    ## Constructors/Destructors                                           ##
    ## ------------------------------------------------------------------ ##
    def __init__(self,time, sample, data, diameter, metadata, header, **kwargs):
        ParticleSizer.__init__(self, time, sample, data, diameter, metadata, header, **kwargs)
        # other stuff if necessary
        self.instrument_type = 'SMPS'

    def __del__(self):
        pass

    ## ------------------------------------------------------------------ ##
    ## Methods                                                            ##
    ## ------------------------------------------------------------------ ##

    # private:

class Empty(ParticleSizer):
    """
    An empty class for storing particle sizer data.

    Parameters
    ----------
    None : class can be initiated empty

    See Also
    --------
    mysmps.core.read.ParticleSizer

    """
    ## ------------------------------------------------------------------ ##
    ## Constructors/Destructors                                           ##
    ## ------------------------------------------------------------------ ##

    def __init__(self, **kwargs):
        time = get_metadata('time')
        time['data'] = []

        date = get_metadata('date')
        date['data'] = []

        sample = get_metadata('sample')
        sample['data'] = []

        data = {}
        data['variables'] = []

        diameter = get_metadata('diameter')
        diameter['data'] = []

        metadata = {}

        pddata = []

        header = []
        ParticleSizer.__init__(self,  time, sample, data, diameter, metadata, header, **kwargs)
        # other stuff if necessary

    def __del__(self):
        pass

    ## ------------------------------------------------------------------ ##
    ## Methods                                                            ##
    ## ------------------------------------------------------------------ ##

    # private:

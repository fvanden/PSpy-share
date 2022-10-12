# -*- coding: utf-8 -*-
#################
import numpy as np
import datetime as dt
import copy

from .basic import conversion
from .timetransform import TimeTransform
tt = TimeTransform()
#################

"""
mysmps.util.stack_ps
mysmps.util.convert_units
mysmps.util.OPCtimetransform
mysmps.util.OPCtimetransformOld
mysmps.util.check_diameters
==================

Functions for working with particle sizer instances
    

Created on Tue Sep 22 13:56 2020

@author: flovan / fvanden

Revision history:   22.09.2020 - Created
                               - stack_ps added
                    23.09.2020 - stacking finished, added fill_time to stack_ps 
                                debugging of stack_ps (sample problem)
                
"""

def stack_ps(ps1, ps2, keep_unique = False, fill_time = False, message = True):
    """
    Combine two particle sizer instances into one.
    
    Parameters
    ----------
    ps1 : ParticleSizer
        mypysmps.core.smps.ParticleSizer object
        
    ps2 : ParticleSizer
        mypysmps.core.smps.ParticleSizer object
        
    keep_unique : bool
        if set to True, attributes which are present in
        one instance but not in the other are kept, if
        False, only attributes common to both instances
        are preserved
        
    fill_time : bool
        if set to True, a time gap between the two 
        particle sizer instances is filled with NaNs
        
    message : bool
        if set to True, helpful warning messages will
        be printed
        
        
    Returns
    -------
    new_ps : ParticleSizer
        mypysmps.core.smps.ParticleSizer object with
        data from both input objects
    """
    # create deepcopies to avoid changing original instances
    
    ps1 = copy.deepcopy(ps1)
    ps2 = copy.deepcopy(ps2)
    
    # create datetime information in PS instances
    
    try:
        _ = getattr(ps1, "datetime")
    except AttributeError:
        ps1.createTimeDate()
        
    try: 
        _ = getattr(ps2, "datetime")
    except AttributeError:
        ps2.createTimeDate()
        
    # check time resolutions
    res1 = (dt.datetime.strptime(ps1.datetime['data'][1], ps1.datetime['units']) - dt.datetime.strptime(ps1.datetime['data'][0], ps1.datetime['units'])).seconds
    res2 = (dt.datetime.strptime(ps2.datetime['data'][1], ps2.datetime['units']) - dt.datetime.strptime(ps2.datetime['data'][0], ps2.datetime['units'])).seconds
    
    if abs(res1-res2) > 60:
        if message:
            print( ("warning: resolutions differ %d seconds")%(abs(res1-res2)) )
         
    # check if ps1 is "older" than ps2
    
    keeponlyoneinstance = False
    reversed_order = False
    cut = None
    
    if dt.datetime.strptime(ps1.datetime['data'][-1], ps1.datetime['units']) < dt.datetime.strptime(ps2.datetime['data'][0], ps2.datetime['units']):
        # ps2 starts after ps1 ends
        timediff = (dt.datetime.strptime(ps2.datetime['data'][0], ps2.datetime['units']) - dt.datetime.strptime(ps1.datetime['data'][-1], ps1.datetime['units'])).total_seconds()
    elif dt.datetime.strptime(ps2.datetime['data'][-1], ps2.datetime['units']) < dt.datetime.strptime(ps1.datetime['data'][0], ps1.datetime['units']):
        # ps1 starts after ps2 ends (user has inadvertently switched the order of the instances)
        reversed_order = True
        timediff = (dt.datetime.strptime(ps1.datetime['data'][0], ps1.datetime['units']) - dt.datetime.strptime(ps2.datetime['data'][-1], ps2.datetime['units'])).total_seconds()
    elif dt.datetime.strptime(ps2.datetime['data'][-1], ps2.datetime['units']) < dt.datetime.strptime(ps1.datetime['data'][-1], ps1.datetime['units']):
        # ps1 ends after ps2 ends
        # The particle sizer instances have overlapping data, but the 
        # first instance covers the whole time period of the second instance.
        # the second instance is ignored
        keeponlyoneinstance = True
        timediff = (dt.datetime.strptime(ps2.datetime['data'][0], ps2.datetime['units']) - dt.datetime.strptime(ps1.datetime['data'][-1], ps1.datetime['units'])).total_seconds()
        if message:
            print("second instance overlaps with first, only first is preserved")
    else:
        # yikes! The particle sizer instances have overlapping data
        # it is assumed that ps2 data replaces ps1 data starting 
        # from the overlapping time
        cut, cutdate = tt.findNearestDate(ps1.datetime['data'], ps2.datetime['data'][0])  
        timediff = (dt.datetime.strptime(ps2.datetime['data'][0], ps2.datetime['units']) - dt.datetime.strptime(ps1.datetime['data'][-1], ps1.datetime['units'])).total_seconds() 
        fill_time = False
        
    #print(timediff, 1.5*res1)
    # check if filling is required
    if fill_time is True:
        # check time difference
        if reversed_order:
        # ps1 starts after ps2 ends
            if timediff > 1.5*res2:
            # the time gap between two instances has to be
            #  larger than twice the normal resolution
                numdates = int(np.ceil(timediff/res2))
                base = dt.datetime.strptime(ps1.datetime['data'][0], ps1.datetime['units'])
                date_list = [base - dt.timedelta(seconds=res2*x) for x in range(numdates)]
                date_list = list(reversed(date_list[1:]))# because numdates starts at 0, first date on date_list is the same as the startdate from the second instance
                datetimelist = [dt.datetime.strftime(dl, ps2.datetime['units']) for dl in date_list]
                ps2.datetime['data'] = np.append(ps2.datetime['data'], datetimelist)
                timelist =  [dt.datetime.strftime(dl, ps2.time['units']) for dl in date_list]
                ps2.time['data'] = np.append(ps2.time['data'], timelist)
                datelist = [dt.datetime.strftime(dl, ps2.date['units']) for dl in date_list]
                ps2.date['data'] = np.append(ps2.date['data'], datelist)
            else:
                fill_time = False
        else:
            if timediff > 1.5*res1:
            # the time gap between two instances has to be
            #  larger than twice the normal resolution
                numdates = int(np.ceil(timediff/res1))
                base = dt.datetime.strptime(ps2.datetime['data'][0], ps2.datetime['units'])
                date_list = [base - dt.timedelta(seconds=res1*x) for x in range(numdates)]
                date_list = list(reversed(date_list[1:])) # because numdates starts at 0, first date on date_list is the same as the startdate from the second instance
                datetimelist = [dt.datetime.strftime(dl, ps1.datetime['units']) for dl in date_list]
                ps1.datetime['data'] = np.append(ps1.datetime['data'], datetimelist)
                timelist =  [dt.datetime.strftime(dl, ps1.time['units']) for dl in date_list]
                ps1.time['data'] = np.append(ps1.time['data'], timelist)
                datelist = [dt.datetime.strftime(dl, ps1.date['units']) for dl in date_list]
                ps1.date['data'] = np.append(ps1.date['data'], datelist)
            else:
                fill_time = False
        
    if message:
        print("reversed order:", reversed_order)
    # check which attributes are similar in both instances
    if reversed_order:
        # ps1 starts after ps2 ends
        new_ps = copy.deepcopy(ps2)
        for attribute in ps1.__dict__.keys():
            if attribute in ps2.__dict__.keys():
                afield = getattr(new_ps, attribute)
                if attribute == 'diameter':
                    st11, st12, st21, st22, diamlist = check_diameters(ps1.diameter['data'], ps2.diameter['data'], ps1.instrument_type)
                              
                    for var in new_ps.data['variables']:
                        if fill_time is True:
                            add = np.ma.zeros((ps2.data[var]['data'].shape[0],len(date_list))) 
                            add[:] = np.nan
                            newdata = np.append(ps2.data[var]['data'],add,axis=1)
                            ps2.data[var]['data'] = newdata
                            
                        sh1 = ps1.data[var]['data'].shape
                        sh2 = ps2.data[var]['data'].shape
                        if cut:
                            sh2 = (sh2[0],cut)
                        newfields = (len(diamlist) ,sh1[1] + sh2[1])
                        new_field = np.ma.zeros(newfields)
                        new_field[:] = np.ma.masked
                        
                        new_field[st21:st22, 0:ps2.data[var]['data'][:,:cut].shape[1]] = ps2.data[var]['data'][:,:cut]
                        new_field[st11:st12, ps2.data[var]['data'][:,:cut].shape[1]:] = ps1.data[var]['data']
                        
                        new_ps.data[var]['data'] = new_field
                        
                    afield['data'] = diamlist
                    
                elif attribute == 'data':
                    # data has been appended with diameters
                    pass
                else:
                    try:
                        field_ps2 = getattr(ps2, attribute)
                        field_ps1 = getattr(ps1, attribute)
                    except TypeError:
                        if attribute == 'header':
                            pass
                        else:
                            if message:
                                print( ("Could not append %s attribute")%(attribute) )
                    try:
                        data_ps2 = field_ps2['data']
                        data_ps1 = field_ps1['data']
                        if attribute in ['date', 'datetime', 'time']: # these have already been extended with the correct data
                            afield['data'] = np.append(data_ps2[:cut], data_ps1)
                        elif fill_time:
                            add = np.ma.zeros(len(date_list))
                            add[:] = np.nan
                            afield['data'] = np.append(np.append(data_ps2[:cut],add), data_ps1)
                        else:
                            afield['data'] = np.append(data_ps2[:cut], data_ps1)
                    except:
                        if message:
                            print( ("Could not append %s attribute")%(attribute) )
                   
            else:
                if keep_unique:
                    newattribute = getattr(ps1,attribute)
                    newattribute['time'] = ps1['datetime']['data']
                    setattr(new_ps, attribute, newattribute)
                else:
                    pass
        if keep_unique is False:
            # get rid of attributes which were in ps2 but not in ps1
            for attribute in ps2.__dict__.keys():
                if attribute in ps1.__dict__.keys():
                    pass
                else:
                    delattr(new_ps, attribute)
            
            
    else:
        # ps2 starts after ps1 ends
        new_ps = copy.deepcopy(ps1)
        if keeponlyoneinstance:
            pass
        else:
            for attribute in ps2.__dict__.keys():
                if attribute in ps1.__dict__.keys():
                    afield = getattr(new_ps, attribute)
                    if attribute == 'diameter':
                        st11, st12, st21, st22, diamlist = check_diameters(ps1.diameter['data'], ps2.diameter['data'], ps1.instrument_type)
                        for var in new_ps.data['variables']:
                            if fill_time is True:
                                add = np.ma.zeros((ps1.data[var]['data'].shape[0],len(date_list))) 
                                add[:] = np.nan
                                newdata = np.append(ps1.data[var]['data'],add,axis=1)
                                ps1.data[var]['data'] = newdata
                                
                            sh1 = ps1.data[var]['data'].shape
                            sh2 = ps2.data[var]['data'].shape
                            if cut:
                                sh1 = (sh1[0],cut)
                            newfields = (len(diamlist) ,sh1[1] + sh2[1])
                            new_field = np.ma.zeros(newfields)
                            new_field[:] = np.ma.masked
                    
                            new_field[st11:st12, 0:ps1.data[var]['data'][:,:cut].shape[1]] = ps1.data[var]['data'][:,:cut]
                            new_field[st21:st22, ps1.data[var]['data'][:,:cut].shape[1]:] = ps2.data[var]['data']
                            
                            new_ps.data[var]['data'] = new_field
                            
                        afield['data'] = diamlist
                        
                    elif attribute == 'data':
                        # data has been appended with diameters
                        pass
                    else:
                        try:
                            field_ps2 = getattr(ps2, attribute)
                            field_ps1 = getattr(ps1, attribute)
                        except TypeError:
                            if attribute == 'header':
                                pass
                            else:
                                if message:
                                    print( ("Could not append %s attribute")%(attribute) )
                        try:
                            data_ps2 = field_ps2['data']
                            data_ps1 = field_ps1['data']
                            if attribute in ['date', 'datetime', 'time']: # these have already been extended with the correct data
                                afield['data'] = np.append(data_ps1[:cut], data_ps2)
                            elif fill_time:
                                add = np.ma.zeros(len(date_list))
                                add[:] = np.nan
                                afield['data'] = np.append(np.append(data_ps1[:cut],add), data_ps2)
                            else:
                                afield['data'] = np.append(data_ps1[:cut], data_ps2)
                        except:
                            if message:
                                print( ("Could not append %s attribute")%(attribute) )
                        
                else:
                    if keep_unique:
                        newattribute = getattr(ps2,attribute)
                        newattribute['time'] = ps2['datetime']['data']
                        setattr(new_ps, attribute,newattribute)
                    else:
                        pass
            if keep_unique is False:
                # get rid of attributes which were in ps2 but not in ps1
                for attribute in ps1.__dict__.keys():
                    if attribute in ps2.__dict__.keys():
                        pass
                    else:
                        delattr(new_ps, attribute)

    if keeponlyoneinstance:
        pass
    else:               
        new_ps.sample['data'] = np.arange(1.0, len(new_ps.datetime['data'])+1)
        new_ps.instrument_type = ps1.instrument_type.split('_')[0] + '_concatenated'
    
    if message:
        print('filltime: ', fill_time)
    
    return new_ps

def convert_units(data, typec,fromm, to):
    """
    Converts data from one unit to another
    
    Parameters
    ----------
    data : list
        list of data values
        
    typec : string
        type of conversion (i.e. temporal or units)
        
    fromm : string
        units to convert from
        
    to : string
        units to convert to
        
    Returns
    -------
    cdata : list
        converted data
    """
    cdata = []
    
    if typec == 'OPC_temporal':
        cdata = OPCtimetransform(data, to)
    elif typec == 'temporal':
        if isinstance(data[0], float):
            data = [str(i) for i in data]
        cdata = tt.convertTime(data, fromm, to)
    elif typec == 'units':
        for measure in data:
            cdata.append(conversion(measure, fromm, to))
    else:
        raise TypeError('Unknown or unsupported conversion tyoe: ' + typec)
        cdata = data
            
    return cdata

def OPCtimetransform(data, to):
    """
    OPC times go from minutes 1..60 and seconds 1..60
    Python can't handle times with minutes 60+ or 
    floating point values with leading zeros. Here
    values are rounded (i.e. 146001 to 150001.)
    
    Parameters
    ----------
    data : list of str
        list of data values to be converted
        
    to : str
        time format (i.e. '%H:%M:%S')
        
    Returns
    -------
    outtimes : list of str
        list of time values with transformed formatting
    """
    
    remove_times = []
    outtimes = []
    times = {'ms':[],'SS':[],'MM':[],'HH':[]}

    for i in range(0, len(data)):
        times['HH'] = 0
        times['MM'] = 0
        times['SS'] = 0
        times['ms'] = 0

        item = data[i]
        
        try:
            if len(item.split('.')[1]) < 2:
                item += '0'
        except IndexError:
            item += '.00'
        if len(item) < 9:
            item = item.zfill(9)
        if int(item[:2]) > 23:
            item = '0' + item
            
        # remove items with extra zero (2319010.00 to 231910)
        if len(item) > 9:
            olditem = item
            newitem = item[:4] + item[5:]
            print( ('Repairing strange value %s into %s')%(olditem, newitem) )
            item = newitem
        else:
            pass
        try:
            md = dt.datetime.strptime(item, "%H%M%S.%f")
            
        # round off items which exceed 59 minutes or 59 seconds 
        # (i.e. 146001 to 150001.)
        except ValueError:
            
            try:
                times['HH'] = int(item[0:2])
                times['MM'] = int(item[2:4])
                times['SS'] = int(item[4:6])
                times['ms'] = int(item[7:9])
            except ValueError:
                print(i, item)

            if times['SS'] > 59:
                times['MM'] += 1
                times['SS'] = 0
            if times['MM'] > 59:
                times['HH'] += 1
                times['MM'] = 0
            # discard items which exceed 23 hours
            if times['HH'] > 23:
                times['HH'] = 23
                print( ('resetting value %s')%(item) )
            

            md = dt.datetime(1900,1,1,times['HH'], times['MM'], times['SS']) 

                
        outtimes.append( dt.datetime.strftime(md, to) )

    return outtimes

def OPCtimetransformOld(data, to):
    """
    OPC times go from minutes 1..60 and seconds 1..60
    Python can't handle times with minutes 60+ or 
    floating point values with leading zeros. Here
    one minute and one second is subtracted from each
    value.
    
    Parameters
    ----------
    data : list of str
        list of data values to be converted
        
    to : str
        time format (i.e. '%H:%M:%S')
        
    Returns
    -------
    outtimes : list of str
        list of time values with transformed formatting
    """
    outtimes = []
    
    times = {
        'ms':[],
        'SS':[],
        'MM':[],
        'HH':[]
    }
    for i in range(0, len(data)):
        item = data[i]
        try:       
            times['HH'].append(int(item[0:2]))
            times['MM'].append(int(item[2:4]))
            times['SS'].append(int(item[4:6]))
            times['ms'].append(int(item[7:9]))
        except ValueError:
            # strange value 2319010.00 in 201129 file...
            olditem = item
            newitem = item[:4] + item[4+1:]
            print( ('Repairing strange value %s into %s')%(olditem, newitem) )
            try:
                times['HH'].append(int(newitem[0:2]))
                times['MM'].append(int(newitem[2:4]))
                times['SS'].append(int(newitem[4:6]))
                times['ms'].append(int(newitem[7:9]))
            except ValueError:
                print(newitem)

    # OPC times go up to 60 minutes. This is corrected by moving one minute
    times['MM'] = [max(0,x-1) for x in times['MM']]
    times['SS'] = [max(0,x-1) for x in times['SS']]

    for i in range(0, len(data)):
        md = dt.datetime(1900,1,1,times['HH'][i], times['MM'][i], times['SS'][i]) 
        outtimes.append( dt.datetime.strftime(md, to))

    return outtimes

def check_diameters(diameter1, diameter2, instrument_type):
    """
    """
    
    """
    nlist1 = []
    nlist2 = []
    olist = []
    
    for i in range(0, len(diameter1)):
        if diameter1[i] in diameter2:
            olist.append(diameter1[i])
        else:
            nlist2.append(diameter1[i])

    for i in range(0, len(diameter2)):
        if diameter2[i] in diameter1:
            olist.append(diameter2[i])
        else:
            nlist1.append(diameter2[i])
            
    """
    diamlist = np.append(diameter1, diameter2)
    diamlist.sort()
    diamlist = np.unique(diamlist)
    
    if instrument_type == 'OPC' or instrument_type == 'OPC_concatenated':
        diamlist = diamlist[:-1]
    
    st11 = abs(diamlist - diameter1[0]).argmin()
    st12 = abs(diamlist - diameter1[-1]).argmin()+1
    st21 = abs(diamlist - diameter2[0]).argmin()
    st22 = abs(diamlist - diameter2[-1]).argmin()+1
    
    if st12 >= len(diamlist):
        st12 = None
    if st22 >= len(diamlist):
        st22 = None
            
    return st11, st12, st21, st22, diamlist


def renameFiles(mypath):
    """
    renames SMPS files and orders them in seperate
    folders for seperate variables

    Note: all files will be copied, renamed and 
    archived in the folder structure:
        mypath
            converted
                raw_counts
                dwdlogdDp
                etc
    
    Parameters
    ----------
    mypath : string
        path to smps files
        
        
    Returns
    -------
    None
    """
    # read in file paths in folder
    files =  sorted( filter( os.path.isfile, glob.glob(mypath + '*.TXT') ) )
    if files == []:
        files =  sorted( filter( os.path.isfile, glob.glob(mypath + '*.csv') ) )
        
    # organise files according to units
    filedict = {}

    filepath = []
    tstart = []
    tend = []
    units = []

    for file in files:
        print(file)
        SMPS = read(file, fileorg = 'AIM')
        filepath.append(file)
        tstart.append(SMPS.datetime['data'][0])
        tend.append(SMPS.datetime['data'][1])
        units.append(SMPS.metadata['Units'])

    for un in np.unique(units):
        filedict[un] = {'start':[],'end':[],'filepath':[]}
    for f in range(0, len(filepath)):
        filedict[units[f]]['start'].append(tstart[f])
        filedict[units[f]]['end'].append(tend[f])
        filedict[units[f]]['filepath'].append(filepath[f])

    # create output folder if necessary
    outpath = os.path.join(mypath,'converted')
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    for un in list(filedict.keys()):
        subpath = os.path.join(outpath,un.replace('/', ''))
        if not os.path.exists(subpath):
            os.makedirs(subpath)
    
    # copy and rename files
    for un in list(filedict.keys()):
        subpath = os.path.join(outpath,un.replace('/', ''))
        for file in filedict[un]['filepath']:
            SMPS = read(file, fileorg = 'AIM')
            ext = os.path.splitext(file)[1]
            newname = SMPS.datetime['data'][0].translate(str.maketrans('', '', string.punctuation)).replace(' ','_')
            newname = newname + ext
            shutil.copy(file, os.path.join(subpath, str(newname)))


def concatenateFilesInFolder(mypath, fileorg):
    """
    Uses stack_ps to concatenate all files in a
    given folder

    Parameters
    ----------
    mypath : string
        path to folder containing particle sizer files

    fileorg : string
        i.e. AIM or OPC

    Returns
    -------
    None
    
    """
    # read in file paths in folder
    files =  sorted( filter( os.path.isfile, glob.glob(mypath + '*.TXT') ) )
    if files == []:
        files =  sorted( filter( os.path.isfile, glob.glob(mypath + '*.csv') ) )
    
    SMPS = read(files[0], fileorg = fileorg)
       
    for file in files[1:]:
        smps2 = read(file, fileorg = fileorg)
        SMPS = stack_ps(SMPS, smps2, message = False, fill_time = True)
        
    return SMPS
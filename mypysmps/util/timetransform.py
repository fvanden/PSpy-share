# -*- coding: utf-8 -*-
#################
from datetime import datetime
import datetime as dt
import re
import sys
#################
"""
Module containing class TimeTransform

Created on Fri Jun 17 17:02:46 2016

@author: fvanden

Content:
    Class TimeTransform
        convertTime
        commonTimeFormat
        findNearestDate
        findSameDate
        findSameDate_Mult
        findCommonDate
        withinDate
        groupbyDate
        seperateDateTime

Revision history:   17.06.2016 - Created
                    07.07.2016 - findSameDate and findSameDate_Mult added,
                    output_format added to all findDate functions so user can
                        decide what type of timeformat is output by these
                    11.07.2016 - findCommonDate added
                    26.10.2017 - groupByDate added

"""

## -------------------------------------------------------------------------- ##
class TimeTransform():
    """
    A class containing transformations related to date and time stamps

    Parameters
    ----------
    timearray : list, array with date or time stamps

    Attributes
    ----------
    timearray : list
        Array with date or time stamps

    """

    ## ------------------------------------------------------------------ ##
    ## Constructors/Destructors                                           ##
    ## ------------------------------------------------------------------ ##

    def __init__(self):
        pass

    def __del__(self):
        pass


    ## ------------------------------------------------------------------ ##
    ## Methods                                                            ##
    ## ------------------------------------------------------------------ ##

    # public:

    def convertTime(self,time_list,inputtimeform, outputtimeform):
        """
        converts one time format into another

        Parameters
        ----------
        time_list : list, of strings, time data units

        inputtimeform : str, i.e. '%d%m%Y'

        outputtimeform : str, i.e. '%d-%m-%Y'

        Returns
        -------
        time_list_conv : list, of strings, of converted times

        See Also
        --------
        datetime.datetime
        """

        time_list_re = [datetime.strptime(x, inputtimeform) for x in time_list]
        output_time = [datetime.strftime(x, outputtimeform) for x in time_list_re]



        return output_time


    def commonTimeFormat(self,time_1, time_2, time_1_Format = None, time_2_Format = None, outputFormat = None):
        """
        Transforms different date formats into one single format

        Parameters
        ----------
        time_1 : list, of strings, time units can be seperated by any
            character in ' -.:' and must contain at least yyyymmddHHMM

        time_2 : same as time_1

        Returns
        -------
        time_1 : list, of dates either in datetime format or in outputFormat

        time_2 : list, of dates either in datetime format or in outputFormat
        """
        if time_1_Format:
            # check if a format is given
            pass
        else:
            # if no given format, try to estimate format
            time_1 = [re.sub(r"[ -.:]", r"", x) for x in time_1]
            if len(time_1[0]) == 8:
                time_1_Format = '%Y%m%d'
            elif len(time_1[0]) == 12:
                time_1_Format = '%Y%m%d%H%M'
            elif len(time_1[0]) == 14:
                time_1_Format = '%Y%m%d%H%M%S'
            else:
                print('cant recognise time format, retry with time_1_Format string')
                return

        # the same thing for list2

        if time_2_Format:
            # check if a format is given
            pass
        else:
            # if no given format, try to estimate format
            time_2 = [re.sub(r"[ -.:]", r"", x) for x in time_2]
            if len(time_2[0]) == 8:
                time_2_Format = '%Y%m%d'
            elif len(time_2[0]) == 12:
                time_2_Format = '%Y%m%d%H%M'
            elif len(time_2[0]) == 14:
                time_2_Format = '%Y%m%d%H%M%S'
            else:
                print('cant recognise time format, retry with time_2_Format string')
                return

        # transform date lists into datetime formats

        time_1_d = [datetime.strptime(x, time_1_Format) for x in time_1]
        time_2_d = [datetime.strptime(x, time_2_Format) for x in time_2]

        # if an outputformat is given, transform into this format, else keep
        # datetime format

        if outputFormat:
            output_time1 = [datetime.strftime(x, outputFormat) for x in time_1_d]
            output_time2 = [datetime.strftime(x, outputFormat) for x in time_2_d]
        else:
            output_time1 = time_1_d
            output_time2 = time_2_d

        return output_time1, output_time2


    def findNearestDate(self, time_list, base_date, output_format = None):
        """
        Finds the closest date within a list of dates. Note: if dates in list
        are equally close, the earliest date from the list is output.

        Parameters
        ----------
        time_list : list, of strings, time units can be seperated by any
            character in ' -.:' and must contain at least yyyymmddHHMM

        base_date : string, time units can be seperated by any
            character in ' -.:' and must contain at least yyyymmddHHMM

        output_format : str, of desired output format (i.e. '%Y.%m.%d %H:%M:%S')
            if None, datetime class is output

        Returns
        -------
        idx : integer, the indice of the date in the timelist which is closest
            to the given base date

        date : datetime class or given output format, the date in the timelist
            which is closest to the given base date

        """
        basedate = [base_date]

        time_d, base_d = self.commonTimeFormat(time_list, basedate)

        base_d = base_d[0]

        # calculate absolute time difference between time_list and base_date
        delta =  [abs(x - base_d) for x in time_d]

        # finds the indice of the date in time_list which is closest to the
        # base_date
        idx = delta.index(min(delta))
        date = time_d[idx]

        if output_format is None:
            pass
        else:
            date = datetime.strftime(date, output_format)


        return idx, date

    def findSameDate(self,time_list, base_date, output_format = None):
        """
        Finds the date similar to the given base date within a list of dates.
        Only takes the last similar date in the list. If multiple dates from
        the list are expected to be the same, use findSameDate_Mult.

        Parameters
        ----------
        time_list : list, of strings, time units can be seperated by any
            character in ' -.:' and must contain at least yyyymmddHHMM

        base_date : string, time units can be seperated by any
            character in ' -.:' and must contain at least yyyymmddHHMM

        output_format : str, of desired output format (i.e. '%Y.%m.%d %H:%M:%S')
            if None, datetime class is output

        Returns
        -------
        idx : integer, the indice of the date in the timelist which is similar
            to the given base date

        date : datetime class or given output format, the date in the timelist
            which is similar to the given base date

        """
        basedate = [base_date]

        time_d, base_d = self.commonTimeFormat(time_list, basedate)

        base_d = base_d[0]

        idx = None
        dates = None

        #
        for i in range(0,len(time_d)):
            date = time_d[i]
            if date == base_d:
                idx = i
                dates = date
            else:
                pass

        if output_format is None:
            pass
        else:
            dates = datetime.strftime(dates, output_format)


        return idx, dates


    def findSameDate_Mult(self,time_list, base_date, output_format = None):
        """
        Finds the date similar to the given base date within a list of dates.
        Allows for multiple same dates in the list. Else use findSameDate.

        Parameters
        ----------
        time_list : list, of strings, time units can be seperated by any
            character in ' -.:' and must contain at least yyyymmddHHMM

        base_date : string, time units can be seperated by any
            character in ' -.:' and must contain at least yyyymmddHHMM

        output_format : str, of desired output format (i.e. '%Y.%m.%d %H:%M:%S')
            if None, datetime class is output

        Returns
        -------
        idx : list, the indices of the dates in the timelist which are similar
            to the given base date

        date : list of datetime class or given output format, the dates in the
            timelist which are similar to the given base date

        """
        basedate = [base_date]

        time_d, base_d = self.commonTimeFormat(time_list, basedate)

        base_d = base_d[0]

        idx = []
        dates = []

        #
        for i in range(0,len(time_d)):
            date = time_d[i]
            if date == base_d:
                idx.append(i)
                if output_format is None:
                    pass
                else:
                    date = datetime.strftime(date, output_format)

                dates.append(date)
            else:
                pass


        return idx, dates


    def findCommonDate(self, time_list, base_date, time_diff, output_format = None):
        """
        Finds the date which is within a certain range of the the base date
        within a list of dates. If multiple dates, it takes the last date in
        the list

        Parameters
        ----------
        time_list : list, of strings, time units can be seperated by any
            character in ' -.:' and must contain at least yyyymmddHHMM

        base_date : string, time units can be seperated by any
            character in ' -.:' and must contain at least yyyymmddHHMM

        time_diff : int, accepted time difference (+ and -) in seconds

        output_format : str, of desired output format (i.e. '%Y.%m.%d %H:%M:%S')
            if None, datetime class is output

        Returns
        -------
        idx : int, the indices of the dates in the timelist which are similar
            to the given base date

        date : list of datetime class or given output format, the dates in the
            timelist which are similar to the given base date
        """

        basedate = [base_date]

        time_d, base_d = self.commonTimeFormat(time_list, basedate)

        base_d = base_d[0]

        idx = []
        dates = []

        #
        for i in range(0,len(time_d)):
            date = time_d[i]
            if date >= base_d:
                t_delta = date - base_d
            else:
                t_delta = base_d - date

            if t_delta <= dt.timedelta(0,time_diff):
                idx.append(i)
                if output_format is None:
                    pass
                else:
                    date = datetime.strftime(date, output_format)

                dates.append(date)
            else:
                pass

        return idx, dates

    def withinDate(self, date, timeperiod):
        """
        Finds whether a date is within a given time window

        Parameters
        ----------
        date : str, time units can be seperated by any character in ' -.:'
            and must contain at least yyyymmddHHMM

        timeperiod : 2-tuple of strings, time window. units can be seperated
            by any character in ' -.:' and must contain at least yyyymmddHHMM

        Returns
        -------
        bool : returns True if the data is within the given time window
            (inclusive), else False

        """

        date, _, = self.commonTimeFormat([date], [date])

        timeperiod_start, timeperiod_end = self.commonTimeFormat([timeperiod[0]], [timeperiod[1]])

        if timeperiod_start <= date <= timeperiod_end:
            return True
        else:
            return False

    def groupbyDate(self, data, dates, res, return_timelist = False):
        """
        Groups data by date increment

        Parameters
        ----------
        data : list, list with data

        dates : list of strings, time units can be seperated by any character
            in ' -.:' and must contain at least yyyymmdd

        res : int, time resolution of increment in minutes

        return_timelist : bool, if set to True, the grouped timelist is also
            returned

        Returns
        -------
        grouped_data : list of lists, the data list grouped according to the
            time increment

        """
        # check if data and dates are of the same length

        if len(data) == len(dates):
            pass
        else:
            print("data and dates must be of the same length")

        # sort the data chronologically

        data = [x for _,x in sorted(zip(dates,data))]
        dates = sorted(dates)

        # try to estimate the time format

        if len(dates[0]) == 8:
            time_fmt = '%Y%m%d'
        elif len(dates[0]) == 12:
            time_fmt = '%Y%m%d%H%M'
        elif len(dates[0]) == 14:
            time_fmt = '%Y%m%d%H%M%S'
        else:
            print('cant recognise time format')

        # transform dates into datetime

        tsfm_dates = []

        for i in range(0, len(dates)):
            tsfm_dates.append(datetime.strptime(dates[i], time_fmt))

        grouped_data = []
        grouped_time = []

        i = 0
        while i  < len(tsfm_dates):
            enddate = tsfm_dates[i] + dt.timedelta(minutes = res)
            nn = self.findNearestDate(dates,dt.datetime.strftime(enddate, time_fmt))[0]
            grouped_data.append(data[i:nn+1])
            grouped_time.append(dates[i:nn+1])
            i = nn + 1
            if nn == len(dates)-1:
                break

        if return_timelist is True:
            return grouped_data, grouped_time
        else:
            return grouped_data


    def seperateDateTime(self, datetime, inputformat, outputtimeformat, outputdatformat):
        """
        Separates a DateTime list into a date list and a time list

        Parameters
        ----------

        datetime : list of strings
            input list with datetimes

        inputformat : str
            the datetime format of the input list (i.e. '%Y.%m.%d %H:%M:%S')

        outputtimeformat : str
            the time format of the output str (i.e. '%H:%M:%S')

        outputdateformat : str
            the date format of the output str (i.e. '%Y.%m.%d')

        Returns
        -------

        times : list of strings
            list with times

        dates : list of strings
            list with dates

        """

        return times, dates

    def decimalToDate(self, number, year, outformat = "%Y-%m-%d %H:%M:%S"):
        """
        Transforms decimal value to Date

        Parameters
        ----------

        number : float
            float part of decimal date

        year : int
            year part of decimal date (in four numbers)

        outformat : string
            out format of datetime

        Returns
        -------

        outdate : string
            datetime string in outformat

        """

        # get date
        date = dt.datetime(year, 1, 1) + dt.timedelta(int(number) - 1)

        # get hours
        dec_doy = float(str(number-int(number))[1:])
        hours = dec_doy * 24
        outhours = int(hours)

        # get minutes
        dec_minutes = float(str(hours-int(hours))[1:])
        minutes = dec_minutes * 60
        outminutes = int(minutes)

        # get seconds
        dec_seconds = float(str(minutes-int(minutes))[1:])
        seconds = dec_seconds * 60
        outseconds = round(seconds)

        # correct outseconds

        if outseconds == 60:
            outseconds = 0

        #print(date.year,date.month,date.day, outhours, outminutes,outseconds-1)
        adate = dt.datetime(date.year,date.month,date.day, outhours, outminutes,outseconds) # seconds seem to be from 1 to 60
        outdate = dt.datetime.strftime(adate, outformat)

        return outdate

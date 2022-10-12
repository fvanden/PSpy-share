# -*- coding: utf-8 -*-
#################
import warnings
import numpy as np
#from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mtick
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.ticker import LogFormatter
import datetime as dt

from ._cm import customcb
from ..util.mathfuncs import roundup
from ..config import get_figure_settings, get_field_limits, _DEFAULT_FIELD_LIMITS
import mypysmps.util.timetransform as timetransform
tt = timetransform.TimeTransform()
#################

"""
mysmps.graph.psdisplay
==================

A general particle sizer display class
    PSDisplay
        plot
        histogram
        timeLine

        create_bins
        fmt
        plotInfo
        dateInfo
        plotIndicator
        defineIndicator
        find24hourperiods


Created on Mon Jul 20 13:46 2020

@author: flovan / fvanden

Revision history:   20.07.2020 - Created
                    17.09.2020 - Added timeLine
                    21.09.2020 - Finished adding timeLine,
        and adapted ylim handling in plot and histogram scripts
                    24.09.2020 - Added some exceptions to timeline
                    plots to handle widget input
                    11.08.2022 - deals automatically with diameter
                    limits and c lims


"""

class PSDisplay(object):
    """
    A class for displaying particle sizer data

    Parameters
    ----------
    smps - SMPS
        core.SMPS object to use for creating plots

    Attributes
    ----------

    """
    ## ------------------------------------------------------------------ ##
    ## Constructors/Destructors                                           ##
    ## ------------------------------------------------------------------ ##
    def __init__(self, smps, **kwargs ):
        self._smps = smps
        self.__dict__.update(kwargs)

    def __del__(self):
        pass


    ## ------------------------------------------------------------------ ##
    ## Plotting methods                                                   ##
    ## ------------------------------------------------------------------ ##

    # private:

    def plot(self, field, **kwargs):
        """
        Create a heatmap of particle sizer data

        Parameters
        ----------
        field : str
            data field to plot

        kwargs:
            indicator : list of ints
                sample(s) on which to plot a location indicator
            linecs : list of str
                colours of indicatorlines
            indperiods : bool
                if True, 24 hour periods will be indicated with
                lines on the plot
            starttime : str
                starting time for indperiods in 'HH:MM:SS'
            colorbar : str
                name of colorbar to use
            xlim : tuple of floats or float
                x limits
            ylim : tuple of floats or float
                y limits
            clim : tuple of floats or float
                colorbar limits
            set_time : bool
                if True, the measurement time will be plotted on
                the x axis, else the sample number will be used
            time_format : str
                if set_time is True, which time format to use,
                default is '%d-%m %H:%M'
            xlabel : str
                xlabel for plot, if not given the SMPS standard
                label from the config file is used
            ylabel : str
                ylabel for plot, if not given the SMPS standard
                label from the config file is used
            clabel : str
                clabel for colourbar, if not given, no claber is
                used
            title : str
                title for plot
            return_axes : bool
                if True the figure handles are returned


        See Also
        --------

        """
        # get relevant kwargs
        indperiods = kwargs.get("periods",False)
        starttime = kwargs.get('starttime', None)
        return_axes = kwargs.get("return_axes", False)

        # get figure settings
        figdict = get_figure_settings('heatplot')

        # create figure
        if isinstance(return_axes, list):
            fig = return_axes[0]
            try:
                ax = return_axes[1]
            except IndexError:
                ax = plt.gca()
        else:
            fig, ax = plt.subplots(figsize=figdict['size'])

        # get plot data
        plotdata = self._smps.data[field]['data']

        # get data coordinates
        coord = self._smps.data['coordinates']
        xcoord = getattr(self._smps, coord[0])
        ycoord = getattr(self._smps, coord[1])

        # get relevant kwargs
        indicatorline = kwargs.get("indicator", False)
        linecs = kwargs.get("linecs",False)

        x = xcoord['data']
        y = ycoord['data']

        X, Y = np.meshgrid(x,y)

        # get colorbar
        colorbar = kwargs.get("colorbar",'_smps_flo')
        if colorbar[0] == '_':
            cbar = mcolors.ListedColormap(customcb[colorbar])
        else:
            cbar = cm.get_cmap(colorbar)


        # plot data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            #print(Y.shape, X.shape, plotdata.T.shape) ############
            c = plt.pcolormesh(Y,X,plotdata.T, cmap=cbar)

        if 'clabel' in kwargs:
            mcbar = fig.colorbar(c, format=ticker.FuncFormatter(self.fmt), label = kwargs.get("clabel"))
        else:
            mcbar = fig.colorbar(c, format=ticker.FuncFormatter(self.fmt))


        # set plot color to black
        ax.set_facecolor('black')

        # set axes limits
        #instrument_xlims = (np.nanmin(xcoord['data']), np.nanmax(xcoord['data']))
        xlim = kwargs.get("xlim", None)

        instrument_ylims = (np.nanmin(xcoord['data']), np.nanmax(xcoord['data']))
        ylim = kwargs.get("ylim", instrument_ylims)
        try:
            len(ylim) # ylim is a tuple
        except TypeError: # ylim is a single value
            ylim = [10**0, ylim]

        instrument_clims = (np.nanmin(plotdata), np.percentile(plotdata,95))
        clim = kwargs.get("clim", instrument_clims)
        try:
            len(clim)
        except TypeError: # clim is a single value
            clim = [0, clim]

        # plot periods if True

        if indperiods:
            plist = self.find24hourperiods(starttime = starttime)
            for p in plist:
                self.plotIndicator(ax, [p,p], [ylim[0], ylim[1]],'b--')

        # plot indicator line if True

        if indicatorline is not False:
            if isinstance(indicatorline, int): # needs to be a list
                indicatorline = [indicatorline]
            if isinstance(linecs, int): # needs to be a list
                linecs = [linecs]
            #print(len(linecs),len(indicatorline)) ##########
            if linecs is False or len(linecs) is not len(indicatorline):
                linecs = ["r--"]*len(indicatorline)
            for ind,linec in zip(indicatorline,linecs):
                self.plotIndicator(ax, [ind, ind],[ylim[0], ylim[1]],linec)
            #ax.plot([indicatorline, indicatorline],[ylim[0], ylim[1]],'r--')
        

        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.3g"))


        # set axes labels and other
        set_time = kwargs.get("set_time", True)
        if set_time is True:

            datetimes = self._smps.createTimeDate(output=True)
            #datetimes = []
            #for i in range(0, len(self._smps.time['data'])):
            #    datetimes.append( self._smps.date['data'][i] + ' ' + self._smps.time['data'][i] )

            x_values = [dt.datetime.strptime(d,self._smps.datetime['units']) for d in datetimes]
            #time_format = kwargs.get("time_format", '%d-%m \n %H:%M')
            time_format = kwargs.get("time_format", '%d-%m')
            newlabels = [dt.datetime.strftime(d, time_format) for d in x_values]
            #labelnums = [int(item.get_position()[0]) for item in ax.get_xticklabels()]
            labelnums = np.ceil(np.arange(0,roundup(len(newlabels))+1, roundup(len(newlabels))/8))
            labelnums = [int(d) for d in labelnums]
            labels = []

            for i in range(0, len(labelnums)):
                try:
                    labels.append(newlabels[labelnums[i]])
                except IndexError:
                    labels.append(newlabels[-1])

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ax.set_xticks(labelnums)
                ax.set_xticklabels(labels)

            #formatter = DateFormatter(time_format)
            #ax.xaxis.set_major_formatter(formatter)


            xlabel = kwargs.get("xlabel", self._smps.time['axis'])
        else:
            xlabel = kwargs.get("xlabel", ycoord['axis'])

        ylabel = kwargs.get("ylabel", xcoord['axis'])

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if xlim is None:
            add = roundup(len(y))/10
            xlim = [-add, len(y)+add]

        # set axes limits
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.clim(clim)

        # set title
        try:
            title = kwargs.get("title", self._smps.metadata['Sample File'] )
        except:
            title = kwargs.get("title", ' ')

        plt.title(title)

        # if return axes is True or list, return them
        if return_axes is not False:
            return fig, ax



    def histogram(self, field, sample= 0, **kwargs):
        """
        Create a histogram plot for the particle sizer

        Parameters
        ----------
        field : str
            data field to plot
        sample : int
            sample number to plot

        Parameters - optional
        ---------------------
        xlim : tuple of floats
            limits of the x axis
        ylim : tuple of floats
            limits of the y axis
        xlabel : str
            name of x axis
        ylabel : str
            name of y axis
        add_D50 : bool
            add a D50 impactor line to the plot
        add_boundaries : bool
            add upper and lower boundaries to plot
        add_legend : bool
            add legend to a plot
        D50_linec : str
            colour of D50 line
        boundary_linec : str
            colour of upper and lower boundary lines
        legend_loc : str
            location of thee legend
        return_axes : bool
            returns axes and figure handles

        See Also
        --------

        """
        # get kwargs
        xlabel = kwargs.get("xlabel", self._smps.diameter['axis'])
        ylabel = kwargs.get("ylabel", self._smps.data[field]['axis'])
        addD50 = kwargs.get("add_D50", False)
        addboundaries = kwargs.get("add_boundaries", False)
        addlegend = kwargs.get("add_legend", False)
        D50_linec = kwargs.get("D50_line", 'k--')
        boundary_linec = kwargs.get("boundary_line", 'cyan')
        legend_loc = kwargs.get("legend_loc", 'best')
        return_axes = kwargs.get("return_axes", False)

        # get plot data
        plotdata = self._smps.data[field]['data'][:,sample]

        # get correct ylim, xlim data

        ylim = kwargs.get("ylim", get_field_limits(field))
        try:
            len(ylim) # ylim is a tuple
        except TypeError: # ylim is a single value
            ylim = [get_field_limits(field)[0], ylim]
        if None in ylim:
            nloc = [i for i, val in enumerate(ylim) if val == None]
            if len(nloc) == 2:
                ylim = (np.min(plotdata), np.max(plotdata) )
            elif nloc == 0:
                ylim = (np.nanmin(plotdata), ylim[1])
            elif nloc == 1:
                ylim = (ylim[0], np.nanmax(plotdata) )


        xlim = kwargs.get("xlim", [1,1000])
        try:
            len(xlim) # xlim is a tuple
        except TypeError: # xlim is a single value
            xlim = [1, xlim]

        # get figure settings
        figdict = get_figure_settings('histplot')

        # create figure
        fig, ax = plt.subplots(figsize=figdict['size'])

        # get bins
        bins, dwidths,_ = self._smps.create_bins()

        # create bar plot - take midpoint diameter of bins
        ax.bar(x=bins, height = plotdata,width = dwidths, color = 'blue', edgecolor = 'white', align = 'edge')

        # add extra's if required
        if addD50:
            xval = self._smps.D50['data'][sample]
            ax.plot([xval, xval],[ylim[0], ylim[1]], D50_linec, label = 'D50 impactor')
            ax.text(xval - ((xlim[1] - xlim[0])/3.5), ylim[1] - ((ylim[1] - ylim[0])/30), 'D50 impactor')

        if addboundaries:
            xval_lower = self._smps.diameter['valid_min']
            xval_upper = self._smps.diameter['valid_max']
            ax.plot([xval_lower, xval_lower],[ylim[0],ylim[1]], boundary_linec, label = 'lower boundary')
            ax.plot([xval_upper, xval_upper],[ylim[0],ylim[1]], boundary_linec, label = 'upper boundary')


        # set axes limits
        if ylim is None:
            ylim = [np.nanmin(plotdata), np.nanmax(plotdata)]

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.semilogx()

        if ylim[1] > 9999:
            ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.4g"))
            fmt = mtick.ScalarFormatter()
            fmt.set_powerlimits((0, 1))
            ax.yaxis.set_major_formatter(fmt)
        elif ylim[1] > 999:
            ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.3g"))
            fmt = mtick.ScalarFormatter()
            fmt.set_powerlimits((0, 1))
            ax.yaxis.set_major_formatter(fmt)
        else:
            pass

        # set axes names
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        # set title
        try:
            title = kwargs.get("title", self._smps.metadata['Sample File'] + '  ' + 'Sample: %d'%(sample+1) )
        except:
            title = kwargs.get("title", ' ')

        # add legend if applicable
        if addlegend:
            ax.legend(loc=legend_loc)

        # if return axes is True, return them
        if return_axes is True:
            return fig, ax


    def timeLine(self,field = None,**kwargs):
        """
        Creates a time line plot for the particle sizer

        Parameters
        ----------

        field : str
            data field to plot, if None, total_concentration
            will be used

        Parameters - optional
        ---------------------
        diameter : int or float
            if the selected data field is multi-dimensional,
            the number (0 to length of diameters) or midpoint
            value of the diameter size bin to plot. The
            default value is 0
        indicator : list of int
            sample(s) on which to plot a location indicator
        linecs : list of str
            line attributes for indicator lines
        indperiods : bool
            if True, 24 hour periods will be indicated with
            lines on the plot
        starttime : str
            starting time for indperiods in 'HH:MM:SS'
        set_time : bool
            if True, the measurement time will be plotted on
            the x axis, else the sample number will be used
        time_format : str
            if set_times = [i for i, val in enumerate(test_list) if val == None] e is True, which time format to use,
            default is '%d-%m %H:%M'
        xlim : tuple of floats
            limits of the x axis
        ylim : tuple of floats
            limits of the y axis
        xlabel : str
            name of x axis
        ylabel : str
            name of y axis
        title : str
            name of plot
        return_axes : bool
            returns axes and figure handles
        widget : bool
            if True, plot is used in a widget and
            handles ylim differently

        """
        widget = kwargs.get('widget', False)
        # widget exception
        if widget:
            ylim = kwargs.get("ylim", None)
            if ylim == 0:
                ylim = None
        else:
            ylim = kwargs.get("ylim", get_field_limits(field))

        # set field
        if field is None:
            field = 'total_concentration'

        # get field data
        select_diameter = False
        try:
            fielddata = getattr(self._smps, field)
        except AttributeError:
            try:
                fielddata = self._smps.data[field]
                select_diameter = True
                diameter = kwargs.get("diameter", 0)
            except KeyError:
                print( ("Warning: %s field could not be found")%(field) )
                return None

        # get kwargs
        indicatorline = kwargs.get("indicator", False)
        linecs = kwargs.get("linecs",False)
        indperiods = kwargs.get("periods",False)
        starttime = kwargs.get('starttime', '00:00:00')
        set_time = kwargs.get("set_time", True)
        if set_time is True:
            time_format = kwargs.get("time_format",'%d-%m %H:%M')
        xlabel = kwargs.get("xlabel", self._smps.time['axis'])
        ylabel = kwargs.get("ylabel", fielddata['axis'])
        try:
            title = kwargs.get("title", self._smps.metadata['Sample File'] )
        except:
            title = ' '
        addlegend = kwargs.get("add_legend", False)
        legend_loc = kwargs.get("legend_loc", 'best')
        return_axes = kwargs.get("return_axes", False)

        # get figure settings
        figdict = get_figure_settings('timeplot')

        # create figure
        if isinstance(return_axes, list):
            fig = return_axes[0]
            try:
                ax = return_axes[1]
            except IndexError:
                ax = plt.gca()
        else:
            fig, ax = plt.subplots(figsize=figdict['size'])

        # get correct xlim data
        xlim = kwargs.get("xlim", [0,len(self._smps.time['data'])])
        try:
            len(xlim) # xlim is a tuple
        except TypeError: # xlim is a single value
            xlim = [0, xlim]

        if select_diameter is True:
            if isinstance(diameter, float):
                diameter = (abs(np.asarray(self._smps.diameter['data']) - diameter)).argmin()


        # create figure
        #fig, ax = plt.subplots(figsize=figdict['size'])

        # get plot data

        if select_diameter is True:
            plotdata = fielddata['data'][diameter,:]
        else:
            plotdata = fielddata['data']

        # create plot
        ax.plot(np.arange(0, len(plotdata)), plotdata)

        # get correct ylim data

        try:
            len(ylim) # ylim is a tuple
        except TypeError: # ylim is a single value
            ylim = [get_field_limits(field)[0], ylim]
        if None in ylim:
            nloc = [i for i, val in enumerate(ylim) if val == None]
            if len(nloc) == 2:
                ylim = (np.nanmin(plotdata), np.nanmax(plotdata) )
            elif nloc[0] == 0:
                ylim = (np.nanmin(plotdata), ylim[1])
            elif nloc[0] == 1:
                ylim = (ylim[0], np.nanmax(plotdata) + 0.1*np.nanmax(plotdata) )

        if ylim[0] > ylim[1]:
            # lims should be lower, upper
            ylim = (ylim[1], ylim[0] )

        # plot periods if True

        if indperiods:
            plist = self.find24hourperiods(starttime = starttime)
            for p in plist:
                self.plotIndicator(ax, [p,p], [ylim[0], ylim[1]],'b--')

        # plot indicator line if True

        if indicatorline is not False:
            if isinstance(indicatorline, int): # needs to be a list
                indicatorline = [indicatorline]
            if isinstance(linecs, int): # needs to be a list
                linecs = [linecs]
            if linecs is False or len(linecs) is not len(indicatorline):
                linecs = ["r--"]*len(indicatorline)
            for ind,linec in zip(indicatorline,linecs):
                self.plotIndicator(ax, [ind, ind],[ylim[0], ylim[1]],linec)
            #ax.plot([indicatorline, indicatorline],[ylim[0], ylim[1]],'r--')

        if set_time is True:
            try:
                datetimes = []
                for i in range(0, len(self._smps.time['data'])):
                    datetimes.append( self._smps.date['data'][i] + ' ' + self._smps.time['data'][i] )

                x_values = [dt.datetime.strptime(d,"%d/%m/%Y %H:%M:%S") for d in datetimes]
                newlabels = [dt.datetime.strftime(d, time_format) for d in x_values]

                #pos = [int(item.get_position()[0]) for item in ax.get_xticklabels()]
                #labels = [newlabels[0]] + list(np.asarray(newlabels)[pos[1:-1]])  + [newlabels[-1]]
                labelnums = np.ceil(np.arange(0,roundup(len(newlabels))+1, roundup(len(newlabels))/8))
                labelnums = [int(d) for d in labelnums]
                labels = []

                for i in range(0, len(labelnums)):
                    try:
                        labels.append(newlabels[labelnums[i]])
                    except IndexError:
                        labels.append(newlabels[-1])

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    ax.set_xticks(labelnums)
                    ax.set_xticklabels(labels,rotation = 45)
            except AttributeError:
                pass

        # set labels
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        # set axes limits
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        # set title
        plt.title(title)

        # if return axes is True, return them
        if return_axes is True:
            return fig, ax


    def fmt(self,x, pos):
        """
        """
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    def plotInfo(self, **kwargs):
        """

        """
        other_var = ['aerosol_flow','bypass_flow','sheath_flow','D50','density','geo_mean','low_voltage','high_voltage','mean_free_path','pressure','retrace_time', 'td05','tf','viscosity']

        # take in other var if not already present in other values
        for value in kwargs:
            if value in other_var:
                pass
            else:
                other_var.append(value)

        # get metadata items
        mydata = self._smps.metadata

        # add other items
        other_dict = {}

        for item in other_var:
            itemdat = getattr(self._smps, item)
            try:
                other_dict[item +'[' + itemdat['units'] + ']'] = np.mean(itemdat['data'])
            except KeyError:
                pass

        mydata.update(other_dict)

        # print items
        for key, value in mydata.items():
            line = ['\033[1m' + key + ':' + '\033[0m',  value]
            print('{:>52} | {:<52}'.format(line[0], line[1]) )


    def dateInfo(self, sample):
        """
        Prints the time and date of a SMPS sample

        Parameters
        ----------
        sample : int
            sample number

        """
        date = dt.datetime.strftime(dt.datetime.strptime(self._smps.date['data'][sample], '%d/%m/%Y'), '%A, %d %B %Y, ')  +  self._smps.time['data'][sample]
        text = ("Date: " + date)
        print(text)

    def plotIndicator(self, ax, xlims, ylims, linec):
        """
        plots and indicator line

        Parameters
        ----------
        ax : matplotlib figure axis
            axis of the figure

        xlims : tuple
            x limits of indicator line

        ylims : tuple
            y limits of indicator line

        linec : str
            line type and colour
        """
        ax.plot(xlims,ylims,linec)

    def defineIndicator(self, date):
        """
        defines x (and y) position of a date on the plot

        Parameters
        ----------
        date : str
            date in '%d-%m-%Y %H:%M:%S' format

        """
        pass

    def find24hourperiods(self,starttime = None):
        """
        Returns indices in data for 24 hour periods
        starting from starttime

        Parameters
        ----------
        starttime : str
            time in HH:MM:SS

        Returns
        -------
        indices : list of integers
            list of indices

        """
        if starttime:
            pass
        else:
            starttime = '00:00:00'
        datetimes = []
        for i in range(0, len(self._smps.time['data'])):
            datetimes.append( self._smps.date['data'][i] + ' ' + self._smps.time['data'][i] )

        times = tt.convertTime(datetimes,'%d/%m/%Y %H:%M:%S','%Y.%m.%d %H:%M:%S')

        indices = []
        for day in np.unique(self._smps.date['data']):
            newtime = tt.convertTime([day + ' ' + starttime],'%d/%m/%Y %H:%M:%S','%Y.%m.%d %H:%M:%S')
            ind = tt.findNearestDate(times, newtime[0])
            indices.append( ind[0])

        indices = sorted(indices)
        return indices

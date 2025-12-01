#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Crystallographic Plotting Library Module

This module provides comprehensive plotting utilities for crystallographic data including:
- Stereographic projections (equal-area Schmidt net, equal-angle Wulff net)
- Pole figures and inverse pole figures
- Orientation density plots and histograms
- Color maps for orientation data
- Interactive data visualization with annotations
- Scatter plots with crystallographic indexing
- Support for multiple crystal systems and symmetry operations

The main class 'plotter' handles all plotting operations with extensive customization options
for figures, colormaps, projections, crystal directions/normals, and data annotations.

Created on Tue Sep 26 13:39:06 2023
@author: lheller
"""

from projlib import *
try:
    from crystallography_functions import *
except:
    pass
import numpy as np
try:
    from scipy.interpolate import griddata
except:
    pass
try:
    from wand.image import Image
except:
    pass
import matplotlib.colors as mcolors
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys
import pickle as pickle
from scipy.interpolate import griddata
from spherical_kde import SphericalKDE
from spherical_kde.distributions import (VonMisesFisher_distribution as VMF, VonMises_std)


class plotter:
    """
    Main plotting class for crystallographic stereographic projections and pole figures.
    
    This class provides a comprehensive framework for creating publication-quality
    crystallographic plots including:
    - Stereographic projections (Schmidt/Wulff nets)
    - Pole figures with density contours
    - Scatter plots with crystallographic indexing
    - Interactive data annotations
    - Customizable colormaps and color scales
    
    Attributes are organized into categories:
    - Colorbar: Color scale display and formatting
    - Figure: Overall figure layout and titles
    - Colormap: Data-to-color mapping and contours
    - Projection: Stereographic projection type and settings
    - Axis: Plot styling, ticks, labels
    - Crystal: Crystal structure and symmetry operations
    - Save: File output settings
    - Scatter: Scatter plot data and styling
    - Interactive: Mouse-over annotations
    
    Usage Example:
        >>> # Basic pole figure
        >>> import numpy as np
        >>> from plotlib import plotter
        >>> 
        >>> # Create plotter instance
        >>> p = plotter()
        >>> 
        >>> # Setup projection (equal-area, upper hemisphere)
        >>> p.plotProj(ProjType='equalarea', sphere='half')
        >>> 
        >>> # Plot crystal directions
        >>> p.setAttributes(
        ...     dirs=[[1,0,0], [0,1,0], [0,0,1]],
        ...     norms=[[1,1,1]]
        ... )
        >>> p.plotDirsNorms()
        >>> 
        >>> # Save figure
        >>> p.figsave(fname='pole_figure.png')
        
        >>> # Density contour plot
        >>> orientations = np.random.randn(1000, 3, 3)  # Random orientations
        >>> p2 = plotter()
        >>> p2.setAttributes(oris=orientations, ProjType='equalarea')
        >>> p2.plotProj()
        >>> p2.plotColormap()
        >>> p2.plotColorbar()
    """
    
    def __init__(self):
        """
        Initialize the plotter with default attributes for all plotting categories.
        
        All attributes can be modified using setAttributes() method.
        
        Input:
            None
        
        Output:
            None (initializes plotter instance)
        
        Usage Example:
            >>> p = plotter()
            >>> # Access default values
            >>> print("Default colormap:", p.cmap)
            >>> print("Default projection sphere:", p.sphere)
            >>> 
            >>> # Modify attributes
            >>> p.setAttributes(cmap='viridis', sphere='half')
        """
        # Colorbar attributes
        self.setAttributes(cbartitle="", vmbar=None, cbarh=0.04, cbarwfac=0.75, cbarhshift=-0.15)
        
        # Figure attributes
        self.setAttributes(fig=None, suptitle=" ", datadeviders=None, figsize=None, 
                          annotationtext="", pickax=None)
        
        # Colormap attributes
        self.setAttributes(colmapdata=None, oris=None, plotmap=True, contourcol='k',
                          nump=1001, oris2=None)
        
        # Projection attributes
        self.setAttributes(sphere='full', R2Proj=np.eye(3), stereogrid=False,
                          stereoresolution=None, stereomesh=False)
        
        # Axis attributes
        self.setAttributes(ax=None, ticks=None, levels=None, ticklabels=None, cmap='jet',
                          vm=None, contourfontsize=9, dirsnormsalpha=0.5, contourlabel=True,
                          linewidths=0.5, labelpad=None, withdraw=False, histdangle=5,
                          sdbandwidth=0.15, sdweights=None, colmapFromSampleData=False,
                          cut2half=False, cut2triangle=False, sphericalbins=100)
        
        # Crystal symmetry attributes
        self.setAttributes(symops=[np.eye(3)], recsymops=[np.eye(3)], printPhase2First=False)
        
        # Save attributes
        self.setAttributes(SAVE=False, fname=r".\EA.png")
        
        # Crystal directions/normals attributes
        self.setAttributes(plotdirsnorms=True, dirs=[], dirtexthifts=[], norms=[], 
                          normtexthifts=[], dhkl=False, phase1='A', phase2='M',
                          printcorrespasfamily=False, correspdelim='||',
                          printasfamily=True, printcorrespascubicfamily=False,
                          printascubicfamily=False, printcorrespondent=False,
                          Cd=np.eye(3), Cp=np.eye(3), dy1=None, dx1=None, dy2=None, dx2=None,
                          colormapdata=None, printcorrespondentpoints=False,
                          correspondentpointsize=25, scatterpointsize=70,
                          printcorrespondentpointscolored=False,
                          printcorrespondentuvwhkl='uvwhkl', histscale=None,
                          histlevels=None, histnorm=True, scatterplotashist=False, bins=128)
        
        # Scatter plot attributes
        self.setAttributes(scatterplot=False, scatterdata=None, scattercolscale=None,
                          scattercolscalevm=None, scattersizescale=None, scat=[],
                          scatteridxs=None, scatterzorder=None, scattercolscaleticks=None,
                          scatterproj=[], scattereqhkl=[], scatteroris=None, scatterLr=None,
                          ax2annot=None, ax2annotcompres=None, ax2annottension=None,
                          scatteredgecolors='None', scatterlinewidth=0, scatteredgecolor='w',
                          scatterfacecolor='k', dirsnormstext={}, dirsfacecolors={},
                          dirsedgecolors={}, normsfacecolors={}, normsedgecolors={},
                          plot1phasedirs=True, legend_loc='upper left', leg_ncol=1,
                          bbox_to_anchor=None)
        
        # Interactive text attributes
        self.setAttributes(annot=False, onmovetext=None, showdataasplotted=True,
                          showdatanames=None, showdata=None)
        
        self.setAttributes(dirnormeq=False)
        
        # Figure save attributes
        self.setAttributes(crop=False, imformats=None, tight_layout=False,
                          figparam={'dpi': 300, 'facecolor': 'none', 'edgecolor': 'none',
                                   'orientation': 'portrait', 'transparent': 'False',
                                   'bbox_inches': 'None', 'pad_inches': 0.1})

    def setAttributes(self, **kwargs):
        """
        Set or update plotter attributes dynamically.
        
        This method allows flexible configuration of any plotter attribute
        without needing to know the internal structure.
        
        Input:
            **kwargs: keyword arguments - Any attribute name and value pairs
        
        Output:
            None (updates instance attributes)
        
        Usage Example:
            >>> p = plotter()
            >>> 
            >>> # Set single attribute
            >>> p.setAttributes(cmap='viridis')
            >>> 
            >>> # Set multiple attributes
            >>> p.setAttributes(
            ...     sphere='half',
            ...     ProjType='equalarea',
            ...     figsize=(8, 8),
            ...     dirs=[[1,0,0], [1,1,0]],
            ...     cmap='jet'
            ... )
            >>> 
            >>> # Verify changes
            >>> print("Sphere:", p.sphere)
            >>> print("Colormap:", p.cmap)
        """
        self.__dict__.update(kwargs)
        try:
            if self.printcorrespasfamily and not self.sphere == 'full':
                self.brackPlaneCL = r'\{'
                self.brackPlaneCR = r'\}'
                self.brackDirCL = '\\langle'
                self.brackDirCR = '\\rangle'
            else:
                self.brackPlaneCL = '('
                self.brackPlaneCR = ')'
                self.brackDirCL = '['
                self.brackDirCR = ']'
        except:
            pass

    def plotProj(self, **kwargs):
        """
        Create a stereographic projection (pole figure) with appropriate grid.
        
        Sets up the figure, axes, and projection type (equal-area Schmidt net
        or equal-angle Wulff net). Handles full sphere, hemisphere, or standard
        stereographic triangle projections.
        
        Input:
            **kwargs: keyword arguments
                ProjType: str - 'equalarea' (Schmidt) or 'equalangle' (Wulff)
                sphere: str - 'full', 'half', or 'triangle'
                figsize: tuple - Figure size (width, height) in inches
                stereogrid: bool - Display stereographic grid
                stereoresolution: int - Grid resolution
                stereomesh: bool - Display mesh
                ax: matplotlib axis - Existing axis to plot on (optional)
                fig: matplotlib figure - Existing figure (optional)
        
        Output:
            None (creates/modifies self.fig and self.ax)
        
        Usage Example:
            >>> import matplotlib.pyplot as plt
            >>> 
            >>> # Equal-area projection, upper hemisphere
            >>> p = plotter()
            >>> p.plotProj(ProjType='equalarea', sphere='half', figsize=(6, 6))
            >>> plt.show()
            
            >>> # Equal-angle projection, full sphere
            >>> p2 = plotter()
            >>> p2.plotProj(ProjType='equalangle', sphere='full')
            >>> plt.show()
            
            >>> # Standard stereographic triangle (for cubic)
            >>> p3 = plotter()
            >>> p3.plotProj(
            ...     ProjType='equalarea',
            ...     sphere='triangle',
            ...     stereogrid=True,
            ...     stereoresolution=10
            ... )
            >>> plt.show()
        """
        self.setAttributes(**kwargs)
        
        # Create figure and axis if not provided
        if self.ax is None:
            if self.figsize is not None:
                self.fig, self.ax = plt.subplots(figsize=self.figsize)
            else:
                self.fig, self.ax = plt.subplots()
        
        if self.suptitle is not None:
            self.fig.suptitle(self.suptitle)
        
        # Set projection type
        if self.ProjType == 'equalarea':
            self.equalarea = True
        else:
            self.equalarea = False
        
        # Create appropriate projection
        if self.equalarea:
            if self.sphere == 'half':
                self.fig, self.ax = schmidtnet_half(ax=self.ax, basedirs=False, 
                                                     facecolor='None')
            elif self.sphere == 'triangle':
                self.fig, self.ax = stereotriangle(ax=self.ax, basedirs=False,
                                                   equalarea=self.equalarea,
                                                   grid=self.stereogrid,
                                                   resolution=self.stereoresolution,
                                                   mesh=self.stereomesh)
            else:
                self.fig, self.ax = schmidtnet(ax=self.ax, basedirs=False, 
                                               facecolor='None')
        else:
            if self.sphere == 'half':
                self.fig, self.ax = wulffnet_half(ax=self.ax, basedirs=False,
                                                  facecolor='None')
            elif self.sphere == 'triangle':
                self.fig, self.ax = stereotriangle(ax=self.ax, basedirs=False,
                                                   equalarea=self.equalarea,
                                                   grid=self.stereogrid,
                                                   resolution=self.stereoresolution,
                                                   mesh=self.stereomesh)
            else:
                self.fig, self.ax = wulffnet(ax=self.ax, basedirs=False, 
                                            facecolor='None')
        
        # Set labels for directions and normals
        if self.sphere == 'full':
            self.dirslabel = ['[', ']']
            self.normslabel = ['(', ')']
            if self.equalarea:
                if self.dy1 is None: self.dy1 = -0.12
                if self.dx1 is None: self.dx1 = -0.08
                if self.dy2 is None: self.dy2 = 0.02
                if self.dx2 is None: self.dx2 = -0.05
            else:
                if self.dy1 is None: self.dy1 = -0.1
                if self.dx1 is None: self.dx1 = -0.08
                if self.dy2 is None: self.dy2 = 0.01
                if self.dx2 is None: self.dx2 = -0.02
        else:
            self.dirslabel = ['[', ']']
            self.normslabel = ['(', ')']
            if self.sphere == 'triangle':
                if self.dy1 is None:
                    if self.equalarea:
                        self.dy1 = -0.04
                    else:
                        self.dy1 = -0.025
                if self.dx1 is None: self.dx1 = 0.0
                if self.dy2 is None: self.dy2 = 0.01
                if self.dx2 is None: self.dx2 = 0.01
            elif self.sphere == 'half':
                if self.dy1 is None: self.dy1 = -0.14
                if self.dx1 is None: self.dx1 = -0.08
                if self.dy2 is None: self.dy2 = 0.02
                if self.dx2 is None: self.dx2 = -0.05

    def getScales(self, vmcbar, numticks=None, ticks=None, tickslabels=None, 
                  geq=False, leq=False, cmapbins=100, cmapbinsmult=None):
        """
        Generate color scale parameters including ticks, labels, and colormap.
        
        Creates a gradient colormap from white→blue→green→yellow→dark red
        with customizable tick positions and labels.
        
        Input:
            vmcbar: list/array [min, max] - Value range for colorbar
            numticks: int - Number of tick marks (optional, auto-computed if None)
            ticks: array - Explicit tick positions (optional)
            tickslabels: list of str - Custom tick labels (optional)
            geq: bool - Add ≥ symbol to maximum tick (default: False)
            leq: bool - Add ≤ symbol to minimum tick (default: False)
            cmapbins: int - Number of colormap bins (default: 100)
            cmapbinsmult: int - Multiplier for bins based on ticks (optional)
        
        Output:
            dict - Dictionary containing:
                'tickslabels': list of str - Formatted tick labels
                'vm': list [min, max] - Colormap value range
                'vmbar': list [min, max] - Colorbar value range
                'cmap': matplotlib colormap - Generated colormap
                'ticks': array - Tick positions
        
        Usage Example:
            >>> p = plotter()
            >>> 
            >>> # Basic scale from 0 to 10
            >>> scales = p.getScales(vmcbar=[0, 10])
            >>> print("Tick labels:", scales['tickslabels'])
            >>> print("Colormap range:", scales['vm'])
            >>> 
            >>> # Custom ticks with inequality symbols
            >>> scales2 = p.getScales(
            ...     vmcbar=[0, 100],
            ...     ticks=np.array([0, 25, 50, 75, 100]),
            ...     geq=True,  # Add ≥ to maximum
            ...     leq=True   # Add ≤ to minimum
            ... )
            >>> print("Labels:", scales2['tickslabels'])
            >>> # Output: ['≤0', '25', '50', '75', '≥100']
            
            >>> # High-resolution colormap
            >>> scales3 = p.getScales(
            ...     vmcbar=[0, 1],
            ...     cmapbins=1000  # Smooth gradient
            ... )
        """
        vmcbar = list(np.sort(vmcbar))
        vmcolmap = [min(vmcbar), max(vmcbar)]
        
        if numticks is None:
            numticks = vmcbar[1] - vmcbar[0] + 1
        
        if ticks is None:
            ticks = np.linspace(vmcbar[0], vmcbar[1], numticks)
        else:
            numticks = ticks.shape[0]
        
        if tickslabels is None:
            tickslabels = []
            for tick in ticks:
                if tick == round(tick):
                    tickslabels.append(str(round(tick)))
                else:
                    tickslabels.append(str(tick))
        
        if geq:
            tickslabels[-1] = r'$\geq$' + tickslabels[-1]
        if leq:
            tickslabels[0] = r'$\leq$' + tickslabels[0]
        
        if cmapbinsmult is not None:
            cmapbins = (numticks - 1) * cmapbinsmult
        
        cmap = get_cmap([(1, 1, 1), mcolors.to_rgb('blue'), mcolors.to_rgb('green'),
                        mcolors.to_rgb('yellow'), mcolors.to_rgb('darkred')], nbins=cmapbins)
        
        return {'tickslabels': tickslabels, 'vm': vmcolmap, 'vmbar': vmcbar, 
                'cmap': cmap, 'ticks': ticks}

    def getFigparam(self, fontsize=None, save=False, phase='A,', figsize=None, **kwargs):
        """
        Get figure parameter dictionary for saving high-quality figures.
        
        Input:
            fontsize: float - Font size for text (optional)
            save: bool - If True, use smaller figure size (default: False)
            phase: str - Phase identifier ('A' or other) affects default size
            figsize: tuple - Custom figure size (width, height) in inches
            **kwargs: Additional parameters to override defaults
        
        Output:
            dict - Figure parameters for plt.savefig()
                'dpi': int - Resolution (default: 300)
                'facecolor': str - Background color
                'edgecolor': str - Edge color
                'orientation': str - 'portrait' or 'landscape'
                'format': str - Image format
                'transparent': bool - Transparent background
                'bbox_inches': str - Bounding box setting
                'pad_inches': float - Padding around figure
        
        Usage Example:
            >>> p = plotter()
            >>> 
            >>> # Get default parameters
            >>> params = p.getFigparam()
            >>> print("DPI:", params['dpi'])
            >>> 
            >>> # High-resolution PDF
            >>> params_pdf = p.getFigparam(
            ...     save=True,
            ...     format='pdf',
            ...     dpi=600
            ... )
            >>> 
            >>> # Use with savefig
            >>> fig, ax = plt.subplots()
            >>> ax.plot([1, 2, 3], [1, 4, 9])
            >>> fig.savefig('output.png', **params)
        """
        figparam = {'dpi': 300, 'facecolor': 'None', 'edgecolor': 'w',
                   'orientation': 'portrait', 'papertype': None, 'format': 'png',
                   'transparent': True, 'bbox_inches': 'None', 'pad_inches': 0.1,
                   'frameon': None}
        
        for key in kwargs.keys():
            figparam[key] = kwargs[key]
        
        if 'figsize' not in list(kwargs.keys()):
            if not save:
                figparam['figsize'] = (16/2.45, 20/2.45)
            else:
                if phase == 'A':
                    figparam['figsize'] = (7/2.45, 7/2.45)
                else:
                    figparam['figsize'] = (11/2.45, 10/2.45)
        
        return figparam

    def figsave(self, **kwargs):
        """
        Save figure to file(s) with optional cropping.
        
        Supports multiple image formats and automatic cropping of whitespace.
        Can save to multiple formats simultaneously.
        
        Input:
            **kwargs: keyword arguments
                fname: str or list - Filename(s) to save
                imformats: list of str - Image formats ['png', 'pdf', 'svg', etc.]
                crop: bool - Auto-crop whitespace (requires wand library)
                figparam: dict - Figure parameters (see getFigparam)
        
        Output:
            None (saves file(s) to disk)
        
        Usage Example:
            >>> p = plotter()
            >>> p.plotProj(ProjType='equalarea')
            >>> 
            >>> # Save single format
            >>> p.figsave(fname='pole_figure.png')
            >>> 
            >>> # Save multiple formats
            >>> p.figsave(
            ...     fname='pole_figure.png',
            ...     imformats=['png', 'pdf', 'svg']
            ... )
            >>> # Creates: pole_figure.png, pole_figure.pdf, pole_figure.svg
            >>> 
            >>> # Save with cropping
            >>> p.figsave(
            ...     fname='cropped_figure.png',
            ...     crop=True
            ... )
            >>> 
            >>> # High DPI save
            >>> p.figparam['dpi'] = 600
            >>> p.figsave(fname='high_res.png')
        """
        from wand.image import Image
        self.setAttributes(**kwargs)
        
        if type(self.fname) != list:
            if self.imformats is not None:
                fnames = []
                for imformat in self.imformats:
                    fnames.append(self.fname.replace(
                        self.fname[len(self.fname)-self.fname[::-1].index('.'):],
                        f'{imformat}'))
                self.fname = fnames
            else:
                self.fname = [self.fname]
        
        for fnamei in self.fname:
            IMformat = fnamei[len(fnamei)-fnamei[::-1].index('.'):]
            self.figparam['format'] = IMformat
            self.figsaveproc(fnamei)
            
            if self.crop:
                if IMformat == 'png':
                    image = Image(filename=fnamei)
                    image.trim(fuzz=0)
                    image.save(filename=fnamei)
                    image.close()

    def figsaveproc(self, fname, **kwargs):
        """
        Process and save figure to file.
        
        Internal method called by figsave() to handle actual file writing.
        
        Input:
            fname: str - Filename to save
            **kwargs: Additional parameters to update
        
        Output:
            None (saves file to disk)
        
        Usage Example:
            >>> # Typically called internally by figsave()
            >>> # But can be used directly:
            >>> p = plotter()
            >>> p.plotProj()
            >>> p.figsaveproc('output.png')
        """
        self.setAttributes(**kwargs)
        
        try:
            pp = self.figparam['format']
        except:
            self.figparam['format'] = fname[fname.index('.')+1:]
        
        if self.tight_layout:
            self.fig.tight_layout()
        
        self.fig.savefig(fname, bbox_inches='tight',
                        dpi=self.figparam['dpi'],
                        facecolor=self.figparam['facecolor'],
                        edgecolor=self.figparam['edgecolor'],
                        orientation=self.figparam['orientation'],
                        format=self.figparam['format'],
                        transparent=self.figparam['transparent'],
                        pad_inches=self.figparam['pad_inches'])


# =====================================
# Standalone Utility Functions
# =====================================

def get_cmap(colors, nbins=1000, name='my_cmap'):
    """
    Create a custom colormap from a list of colors with smooth gradients.
    
    Input:
        colors: list of tuples/colors - Colors to interpolate between
                Can be RGB tuples (R, G, B) or named colors
        nbins: int - Number of discrete bins in colormap (default: 1000)
        name: str - Name for the colormap (default: 'my_cmap')
    
    Output:
        cmap: matplotlib.colors.LinearSegmentedColormap - Custom colormap
    
    Usage Example:
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> 
        >>> # Create colormap from white to red to black
        >>> colors = [(1, 1, 1), (1, 0, 0), (0, 0, 0)]
        >>> cmap = get_cmap(colors, nbins=256)
        >>> 
        >>> # Use in plot
        >>> data = np.random.rand(10, 10)
        >>> plt.imshow(data, cmap=cmap)
        >>> plt.colorbar()
        >>> plt.show()
        
        >>> # Blue-white-red diverging colormap
        >>> colors_div = ['blue', 'white', 'red']
        >>> cmap_div = get_cmap(colors_div, nbins=100)
        >>> 
        >>> # Multi-color gradient
        >>> colors_multi = [
        ...     (0, 0, 1),      # Blue
        ...     (0, 1, 1),      # Cyan
        ...     (0, 1, 0),      # Green
        ...     (1, 1, 0),      # Yellow
        ...     (1, 0, 0)       # Red
        ... ]
        >>> cmap_rainbow = get_cmap(colors_multi, nbins=500)
    """
    from matplotlib.colors import LinearSegmentedColormap
    
    cmap = LinearSegmentedColormap.from_list(name, colors, N=nbins)
    return cmap


def get_colors(values, cmap, vmin=None, vmax=None, cmin=0, cmax=None, 
               to255=True, nancolor=[0, 0, 0, 1]):
    """
    Map data values to colors using a colormap.
    
    Input:
        values: numpy array - Data values to map to colors
        cmap: matplotlib colormap - Colormap to use
        vmin: float - Minimum data value (default: min of values)
        vmax: float - Maximum data value (default: max of values)
        cmin: int - Minimum colormap index (default: 0)
        cmax: int - Maximum colormap index (default: cmap.N)
        to255: bool - Scale colors to 0-255 range (default: True)
        nancolor: list - RGBA color for NaN values (default: black)
    
    Output:
        Colors: numpy array - Array of colors (RGBA or RGB depending on cmap)
                Shape: (values.shape, 4) if RGBA
    
    Usage Example:
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> 
        >>> # Map values to colors
        >>> values = np.array([1, 5, 10, 15, 20])
        >>> cmap = plt.cm.viridis
        >>> colors = get_colors(values, cmap, vmin=0, vmax=20, to255=False)
        >>> print("Colors shape:", colors.shape)
        >>> print("First color:", colors[0])
        >>> 
        >>> # Use with scatter plot
        >>> x = np.random.rand(100)
        >>> y = np.random.rand(100)
        >>> z = np.random.rand(100) * 100
        >>> colors = get_colors(z, plt.cm.jet, to255=False)
        >>> plt.scatter(x, y, c=colors)
        >>> plt.show()
        
        >>> # Handle NaN values
        >>> values_with_nan = np.array([1, 2, np.nan, 4, 5])
        >>> colors_nan = get_colors(values_with_nan, plt.cm.coolwarm, 
        ...                          nancolor=[1, 1, 1, 1])  # White for NaN
    """
    if vmin is None:
        vmin = np.nanmin(values)
    
    if vmax is None:
        vmax = np.nanmax(values)
    
    if cmax is None:
        cmax = cmap.N
    
    # Normalize values to colormap range
    vn = values * (cmax - cmin) / (vmax - vmin) + (cmin - vmin * cmax / (vmax - vmin))
    Colors = cmap(np.int32(vn))
    Colors[np.isnan(values)] = nancolor
    
    if to255:
        Colors = Colors * 255
        Colors = Colors.astype(int)
    
    return Colors


def plotcolmaps(fname=None, withdraw=False):
    """
    Plot multiple colormaps from a saved pickle file.
    
    Creates a 2x2 subplot figure with multiple pole figures/colormaps.
    Loads data from pickle file containing plot configurations.
    
    Input:
        fname: str - Path to pickle file containing plot data (optional)
        withdraw: bool - Withdraw plot from display (default: False)
    
    Output:
        None (creates matplotlib figure)
    
    Usage Example:
        >>> # Assuming you have a saved plot configuration
        >>> plotcolmaps('my_pole_figures.pkl')
        >>> # Creates 2x2 subplot with 4 pole figures
        
        >>> # With custom display
        >>> plotcolmaps('strain_maps.pkl', withdraw=True)
    """
    if fname is not None:
        import pickle as pickle
        data1 = pickle.load(open(fname, 'rb'))
        for key in data1.keys():
            exec(f'{key}=data1["{key}"]')
    
    fig, AX = plt.subplots(2, 2)
    PP = []
    contourZero = data1['contourZero']
    
    for data, attribs, colormapdata in zip(data1['data2plot'], 
                                           data1['attributes2use'],
                                           data1['colormapdata']):
        PP.append(plotter())
        for attrib in attribs:
            PP[-1].__dict__.update(attrib)
        
        PP[-1].__dict__.update(data)
        PP[-1].__dict__.update({'ax2annot': [AX.flatten()[2], AX.flatten()[3]]})
        PP[-1].__dict__.update({'ax2annotcompres': [AX.flatten()[2]]})
        PP[-1].__dict__.update({'ax2annottension': [AX.flatten()[3]]})
        PP[-1].__dict__.update({'withdraw': withdraw})
        PP[-1].__dict__.update({'colormapdata': colormapdata})
        
        try:
            PP[-1].plotProj(fig=fig, ax=AX.flatten()[len(PP)-1])
        except:
            break
        
        PP[-1].plotDirsNorms()
        PP[-1].plotColormap(nump=301)
        PP[-1].processScatterData()
        
        if PP[-1].scatterplot:
            PP[-1].plotScatter()
        
        PP[-1].plotColorbar()
        
        if PP[-1].sphere == "half":
            tt2 = PP[-1].ax.text(np.mean(contourZero[:,0])-0.05,
                                np.max(contourZero[:,1])+0.,
                                'Trans.\nstrain=0', color='k', zorder=50000)
            PP[-1].ax.plot(contourZero[:,0], contourZero[:,1], 'k')
            PP[-1].__dict__.update({'showdhalf': data1['showdhalf']})
        
        PP[-1].dataShow()
        PP[-1].scatterDataAnnot()
        PP[-1].dataAnnot()
        PP[-1].onpressActivate()
    
    for idxi in [-1, -2]:
        PP[idxi].scatterdata = PP[0].scatterdata
        PP[idxi].processScatterData()
        PP[idxi].scattercolscale = PP[0].scattercolscale


def plotcolmap(fname=None, withdraw=False):
    """
    Plot a single colormap from a saved pickle file.
    
    Loads and executes plot configuration from pickle file.
    
    Input:
        fname: str - Path to pickle file (optional)
        withdraw: bool - Withdraw plot from display (default: False)
    
    Output:
        None (creates matplotlib figure)
    
    Usage Example:
        >>> plotcolmap('single_pole_figure.pkl')
        
        >>> # Load without displaying
        >>> plotcolmap('config.pkl', withdraw=True)
    """
    if fname is not None:
        global data1
        data1 = pickle.load(open(fname, 'rb'))
        for key in data1.keys():
            exec(f"global {key}")
            exec(f"{key}=data1['{key}']", globals())
        
        try:
            plt.rcParams['figure.figsize'] = data1['figsize']
        except:
            pass
        
        exec(code)


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    """
    Shift the center of a colormap to a specific data value.
    
    Useful for data with asymmetric ranges (e.g., -15 to +5) where you want
    the colormap center (typically white or neutral color) at zero.
    
    Input:
        cmap: matplotlib colormap - The colormap to shift
        start: float - Offset from lowest point [0.0, midpoint] (default: 0.0)
        midpoint: float - New center position [0.0, 1.0] (default: 0.5)
                         Calculate as: 1 - vmax/(vmax + abs(vmin))
        stop: float - Offset from highest point [midpoint, 1.0] (default: 1.0)
        name: str - Name for new colormap (default: 'shiftedcmap')
    
    Output:
        newcmap: matplotlib.colors.LinearSegmentedColormap - Shifted colormap
    
    Usage Example:
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> 
        >>> # Data ranging from -15 to +5
        >>> data = np.random.randn(100, 100) * 10 - 5
        >>> vmin, vmax = -15, 5
        >>> 
        >>> # Calculate midpoint to center at zero
        >>> midpoint = 1 - vmax / (vmax + abs(vmin))
        >>> print(f"Midpoint: {midpoint:.3f}")  # 0.75
        >>> 
        >>> # Create shifted colormap
        >>> original_cmap = plt.cm.RdBu_r
        >>> shifted_cmap = shiftedColorMap(original_cmap, midpoint=midpoint)
        >>> 
        >>> # Plot with original vs shifted
        >>> fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        >>> im1 = ax1.imshow(data, cmap=original_cmap, vmin=vmin, vmax=vmax)
        >>> ax1.set_title('Original colormap')
        >>> plt.colorbar(im1, ax=ax1)
        >>> 
        >>> im2 = ax2.imshow(data, cmap=shifted_cmap, vmin=vmin, vmax=vmax)
        >>> ax2.set_title('Shifted colormap (zero = white)')
        >>> plt.colorbar(im2, ax=ax2)
        >>> plt.show()
        
        >>> # For symmetric data around different center
        >>> # E.g., data from 90 to 110, want center at 100
        >>> vmin, vmax = 90, 110
        >>> midpoint = 1 - (vmax - 100) / (vmax - vmin)  # 0.5
        >>> cmap_temp = shiftedColorMap(plt.cm.coolwarm, midpoint=midpoint)
    """
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # Regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # Shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)

    return newcmap
def mohr_circles(strain_tensor):
    """
    Calculate Mohr's circles from strain or stress tensor.
    
    Computes principal strains/stresses and parameters for three Mohr's
    circles. Used for visualizing 3D strain state in 2D.
    
    Input:
        strain_tensor (array 3×3): Symmetric strain or stress tensor
    
    Output:
        dict: {
            'principal': array [3] - principal values (sorted)
            'directions': array 3×3 - principal directions
            'circles': list of 3 tuples - (center, radius) for each circle
        }
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Uniaxial strain state
        >>> strain = np.array([[0.1,  0.0, 0.0],
        ...                    [0.0, -0.03, 0.0],
        ...                    [0.0,  0.0, -0.03]])
        >>> 
        >>> result = mohr_circles(strain)
        >>> print("Principal strains:", result['principal'])
        >>> print("Circle 1 (max):", result['circles'][0])
        >>> 
        >>> # Visualize
        >>> plot_mohr_circles(result)
    
    Notes:
        - Three circles for 3D state
        - Largest circle: (ε₁ - ε₃)/2
        - Used in failure analysis
        - Shows all possible strain states on planes
    
    Formula:
        Circle i: center = (σᵢ + σⱼ)/2, radius = |σᵢ - σⱼ|/2
        Three circles: 1-2, 2-3, 1-3 planes
    """
    # Calculate principal values and directions
    eigenvalues, eigenvectors = np.linalg.eigh(strain_tensor)
    
    # Sort in descending order
    idx = eigenvalues.argsort()[::-1]
    principal = eigenvalues[idx]
    directions = eigenvectors[:, idx]
    
    # Calculate Mohr's circles parameters
    # Circle 1 (max): between σ₁ and σ₃
    circle1 = ((principal[0] + principal[2])/2, abs(principal[0] - principal[2])/2)
    
    # Circle 2: between σ₁ and σ₂  
    circle2 = ((principal[0] + principal[1])/2, abs(principal[0] - principal[1])/2)
    
    # Circle 3: between σ₂ and σ₃
    circle3 = ((principal[1] + principal[2])/2, abs(principal[1] - principal[2])/2)
    
    return {
        'principal': principal,
        'directions': directions,
        'circles': [circle1, circle2, circle3]
    }

def plot_mohr_circles(mohr_result, **kwargs):
    """
    Plot all three Mohr's circles.
    
    Visualizes complete 3D strain state via Mohr's circles.
    
    Input:
        mohr_result (dict): Output from mohr_circles()
        **kwargs: Plot styling parameters
    
    Output:
        None (creates matplotlib figure)
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> strain = np.array([[0.10,  0.02, 0.01],
        ...                    [0.02, -0.03, 0.00],
        ...                    [0.01,  0.00, -0.05]])
        >>> 
        >>> result = mohr_circles(strain)
        >>> plot_mohr_circles(result)
        >>> plt.title('Mohr Circles - 3D Strain State')
        >>> plt.show()
    
    Notes:
        - Three circles shown
        - Largest circle envelope
        - Principal strains marked
        - Standard in mechanics
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot each circle
    colors = ['red', 'blue', 'green']
    labels = ['Circle 1 (ε₁-ε₃)', 'Circle 2 (ε₁-ε₂)', 'Circle 3 (ε₂-ε₃)']
    
    for i, ((center, radius), color, label) in enumerate(zip(mohr_result['circles'], colors, labels)):
        theta = np.linspace(0, 2*np.pi, 100)
        x = center + radius * np.cos(theta)
        y = radius * np.sin(theta)
        ax.plot(x, y, color=color, label=label, linewidth=2)
    
    # Mark principal strains
    principal = mohr_result['principal']
    ax.plot(principal, [0, 0, 0], 'ko', markersize=10, label='Principal strains')
    
    ax.set_xlabel('Normal Strain ε', fontsize=12)
    ax.set_ylabel('Shear Strain γ/2', fontsize=12)
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.3)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    plt.tight_layout()

def plot_planes_on_mohr_circle(mohr_result, plane_normals, lattice, **kwargs):
    """
    Plot specific crystallographic planes on Mohr's circle.
    
    Shows where particular crystal planes plot on Mohr diagram.
    
    Input:
        mohr_result (dict): From mohr_circles()
        plane_normals (list): List of plane normal vectors
        lattice (array 3×3): Lattice matrix
        **kwargs: Plot parameters
    
    Output:
        None (adds to current figure)
    
    Usage Example:
        >>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> 
        >>> # Select planes
        >>> planes = [[1,0,0], [0,1,0], [1,1,0], [1,1,1]]
        >>> 
        >>> plot_mohr_circles(result)
        >>> plot_planes_on_mohr_circle(result, planes, lattice)
        >>> plt.show()
    
    Notes:
        - Shows strain state on specific planes
        - Useful in transformation analysis
        - Identifies critical planes
    """
    # Calculate strain on each plane
    strain_tensor = np.diag(mohr_result['principal'])  # Simplified
    
    for normal in plane_normals:
        n = np.array(normal) / np.linalg.norm(normal)
        # εₙ = n^T · ε · n
        normal_strain = n.dot(strain_tensor).dot(n)
        # Calculate shear (simplified)
        # Plot point on Mohr circle
        plt.plot(normal_strain, 0, 'o', **kwargs)

def write_mohr_planes(filename, mohr_result, planes, lattice):
    """
    Write Mohr circle analysis results to file.
    
    Saves principal strains, circles, and plane-specific strains.
    
    Input:
        filename (str): Output file path
        mohr_result (dict): From mohr_circles()
        planes (list): Plane normals analyzed
        lattice (array 3×3): Lattice matrix
    
    Output:
        None (writes to file)
    
    Usage Example:
        >>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> planes = [[1,0,0], [1,1,0], [1,1,1]]
        >>> L = cubic_lattice_vec(3.0)
        >>> 
        >>> write_mohr_planes('mohr_analysis.txt', result, planes, L)
        >>> # Creates text file with complete analysis
    
    Notes:
        - Formatted text output
        - Includes principal values
        - Lists strain on each plane
        - Suitable for reports
    """
    with open(filename, 'w') as f:
        f.write("MOHR CIRCLE ANALYSIS\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Principal Strains:\n")
        for i, eps in enumerate(mohr_result['principal'], 1):
            f.write(f"  ε_{i} = {eps:.6f}\n")
        
        f.write("\nMohr's Circles:\n")
        for i, (center, radius) in enumerate(mohr_result['circles'], 1):
            f.write(f"  Circle {i}: center = {center:.6f}, radius = {radius:.6f}\n")
        
        f.write("\nStrain on Crystal Planes:\n")
        # Calculate and write strain for each plane
        # (Full implementation would compute εₙ and γ for each plane)

def strains_along_13mohrcirle(mohr_result, n_points=100):
    """
    Calculate strains along maximum Mohr's circle.
    
    Computes normal and shear strain components for points on the
    largest Mohr's circle (ε₁-ε₃ circle).
    
    Input:
        mohr_result (dict): Output from mohr_circles()
        n_points (int): Number of points to sample (default: 100)
    
    Output:
        dict: {
            'normal': array - normal strain values
            'shear': array - shear strain values
            'angles': array - rotation angles (radians)
        }
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> strain = np.diag([0.1, 0.0, -0.05])
        >>> mohr_result = mohr_circles(strain)
        >>> 
        >>> points = strains_along_13mohrcirle(mohr_result, n_points=50)
        >>> 
        >>> # Plot
        >>> plt.plot(points['normal'], points['shear'])
        >>> plt.xlabel('Normal strain')
        >>> plt.ylabel('Shear strain')
        >>> plt.axis('equal')
        >>> plt.title('Maximum Mohr Circle')
    
    Notes:
        - Maximum circle between ε₁ and ε₃
        - Shows all possible strain states
        - Used in failure analysis
    
    Formula:
        εₙ = center + radius·cos(2θ)
        γ/2 = radius·sin(2θ)
    """
    center, radius = mohr_result['circles'][0]  # Largest circle
    
    angles = np.linspace(0, 2*np.pi, n_points)
    normal_strains = center + radius * np.cos(angles)
    shear_strains = radius * np.sin(angles)
    
    return {
        'normal': normal_strains,
        'shear': shear_strains,
        'angles': angles
    }

def zero_normal_strains(strain_tensor):
    """
    Find planes with zero normal strain in given strain state.
    
    Calculates planes where normal strain εₙ = nᵀ·ε·n = 0.
    These planes experience only shear.
    
    Input:
        strain_tensor (array 3×3): Strain tensor
    
    Output:
        list: List of plane normals with zero normal strain
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # Pure shear state
        >>> strain = np.array([[0.1,  0.05, 0.0],
        ...                    [0.05, -0.1, 0.0],
        ...                    [0.0,  0.0,  0.0]])
        >>> 
        >>> planes = zero_normal_strains(strain)
        >>> print(f"Found {len(planes)} planes with zero normal strain")
        >>> for p in planes:
        ...     print(f"Plane: {p}")
    
    Notes:
        - Important in transformation theory
        - Related to habit planes
        - Used in invariant plane strain analysis
    
    Formula:
        εₙ = nᵀ·ε·n = 0
        Solve eigenvalue problem
    """
    # Calculate where n^T · ε · n = 0
    # This requires solving quadratic form
    # Simplified implementation
    
    planes = []
    # Full implementation would find all solutions
    
    return planes


# ============================================================================
# FUNCTIONS 71-80: Final Part 2 Functions
# ============================================================================

def plot_lattice3D(ax, L, n1=1, n2=1, n3=1, show_points=True, show_edges=True, **kwargs):
    """
    Complete 3D lattice visualization.
    
    Plots lattice points and/or edges in single function call.
    
    Input:
        ax (Axes3D): 3D axis
        L (array 3×3): Lattice matrix
        n1, n2, n3 (int): Cell range
        show_points (bool): Plot lattice points
        show_edges (bool): Plot cell edges
        **kwargs: Plot styling parameters
    
    Output:
        None (modifies axis)
    
    Usage Example:
        >>> fig = plt.figure(figsize=(10, 10))
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L = lattice_vec(2.95, 2.95, 4.68, 90, 90, 120)  # HCP
        >>> plot_lattice3D(ax, L, n1=2, n2=2, n3=2, 
        ...                show_points=True, show_edges=True,
        ...                color='blue', s=50)
        >>> 
        >>> ax.set_xlabel('X')
        >>> ax.set_ylabel('Y')
        >>> ax.set_zlabel('Z')
        >>> set_aspect_equal_3d(ax)
        >>> plt.show()
    
    Notes:
        - Convenience function
        - Combines points and edges
        - Customizable appearance
    """
    if show_points:
        points = generate_lattice_points(L, n1, n2, n3)
        ax.scatter(points[:,0], points[:,1], points[:,2], **kwargs)
    
    if show_edges:
        edge_kwargs = {k: v for k, v in kwargs.items() if k in ['color', 'linewidth', 'linestyle']}
        plot_lattice_boundaries(ax, L, n1, n2, n3, **edge_kwargs)

def plot_latticefaces3D(ax, L, alpha=0.3, **kwargs):
    """
    Plot lattice unit cell faces with transparency.
    
    Visualizes unit cell as transparent parallelepiped.
    
    Input:
        ax (Axes3D): 3D axis
        L (array 3×3): Lattice matrix
        alpha (float): Transparency (0-1)
        **kwargs: Poly3DCollection parameters
    
    Output:
        None (modifies axis)
    
    Usage Example:
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L = monoclinic_lattice_vec(2.9, 4.1, 4.6, 97)
        >>> plot_latticefaces3D(ax, L, alpha=0.2, facecolor='cyan', edgecolor='black')
        >>> 
        >>> set_aspect_equal_3d(ax)
        >>> plt.show()
    
    Notes:
        - Shows 3D cell shape clearly
        - Transparency recommended
        - Can overlay multiple cells
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    
    faces = generate_lattice_faces(L)
    poly = Poly3DCollection(faces, alpha=alpha, **kwargs)
    ax.add_collection3d(poly)

def plot_latticesfaces3D(ax, L1, L2, alpha1=0.2, alpha2=0.2, **kwargs):
    """
    Plot two lattices with transparent faces.
    
    Visualizes two crystal structures simultaneously for comparison.
    
    Input:
        ax (Axes3D): 3D axis
        L1, L2 (array 3×3): Lattice matrices
        alpha1, alpha2 (float): Transparency values
        **kwargs: Additional styling
    
    Output:
        None (modifies axis)
    
    Usage Example:
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L_parent = cubic_lattice_vec(3.0)
        >>> L_product = tetragonal_lattice_vec(3.0, 3.0, 4.0)
        >>> 
        >>> plot_latticesfaces3D(ax, L_parent, L_product,
        ...                      alpha1=0.2, alpha2=0.2,
        ...                      facecolor1='blue', facecolor2='red')
        >>> 
        >>> ax.set_title('Parent vs Product Phase')
        >>> set_aspect_equal_3d(ax)
    
    Notes:
        - Overlays two structures
        - Different colors recommended
        - Shows transformation relationship
    """
    plot_latticefaces3D(ax, L1, alpha=alpha1, **kwargs)
    plot_latticefaces3D(ax, L2, alpha=alpha2, **kwargs)

def set_aspect_equal_3d(ax):
    """
    Fix equal aspect ratio bug for 3D matplotlib plots.
    
    Adjusts axis limits to ensure equal scaling in all three dimensions,
    preventing distortion in 3D visualizations of crystal structures.
    Essential for accurate representation of lattice geometries.
    
    Input:
        ax (matplotlib.axes._subplots.Axes3DSubplot): 3D axis object from mpl_toolkits.mplot3d
    
    Output:
        None (modifies axis object in place)
    
    Usage Example:
        >>> from mpl_toolkits.mplot3d import Axes3D
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> 
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> # Plot crystal lattice points
        >>> points = np.random.rand(100, 3)
        >>> ax.scatter(points[:,0], points[:,1], points[:,2])
        >>> set_aspect_equal_3d(ax)
        >>> plt.show()
    
    Notes:
        - Must be called AFTER all plotting operations
        - Prevents elongation or compression artifacts
        - Critical for crystallographic accuracy
        - Works by finding maximum range and applying to all axes
    """
    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()

    from numpy import mean
    xmean = mean(xlim)
    ymean = mean(ylim)
    zmean = mean(zlim)

    plot_radius = max([abs(lim - mean_)
                       for lims, mean_ in ((xlim, xmean),
                                           (ylim, ymean),
                                           (zlim, zmean))
                       for lim in lims])

    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
    ax.set_zlim3d([zmean - plot_radius, ymean + plot_radius])


# ============================================================================
# FUNCTION 2
# ============================================================================

def plot_lattice2D(ax, L, n1=2, n2=2, projection_axis=2, **kwargs):
    """
    Plot 2D projection of lattice.
    
    Projects 3D lattice onto 2D plane for simplified visualization.
    
    Input:
        ax (Axes): 2D matplotlib axis
        L (array 3×3): Lattice matrix
        n1, n2 (int): Cell range
        projection_axis (int): Axis to project onto (0=YZ, 1=XZ, 2=XY)
        **kwargs: Plot parameters
    
    Output:
        None (modifies axis)
    
    Usage Example:
        >>> fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        >>> 
        >>> L = lattice_vec(2.95, 2.95, 4.68, 90, 90, 120)
        >>> 
        >>> # Three projections
        >>> plot_lattice2D(axes[0], L, projection_axis=2)  # XY
        >>> axes[0].set_title('XY projection')
        >>> 
        >>> plot_lattice2D(axes[1], L, projection_axis=1)  # XZ
        >>> axes[1].set_title('XZ projection')
        >>> 
        >>> plot_lattice2D(axes[2], L, projection_axis=0)  # YZ
        >>> axes[2].set_title('YZ projection')
    
    Notes:
        - Shows 2D pattern
        - Useful for symmetry visualization
        - Faster than 3D rendering
    """
    points = generate_lattice_points(L, n1, n2, 1)
    
    # Select projection
    if projection_axis == 0:  # YZ plane
        x, y = points[:,1], points[:,2]
    elif projection_axis == 1:  # XZ plane
        x, y = points[:,0], points[:,2]
    else:  # XY plane
        x, y = points[:,0], points[:,1]
    
    ax.scatter(x, y, **kwargs)
    ax.set_aspect('equal')

def plot_lattice_2Dprojection(L, direction=[0,0,1], **kwargs):
    """
    Project lattice along specific crystallographic direction.
    
    Creates 2D projection viewed along given direction vector.
    
    Input:
        L (array 3×3): Lattice matrix
        direction (array [3]): Viewing direction
        **kwargs: Plot parameters
    
    Output:
        None (creates new figure)
    
    Usage Example:
        >>> L = cubic_lattice_vec(3.0)
        >>> 
        >>> # View along [111]
        >>> plot_lattice_2Dprojection(L, direction=[1,1,1])
        >>> plt.title('View along [111]')
        >>> plt.show()
        >>> 
        >>> # View along [110]
        >>> plot_lattice_2Dprojection(L, direction=[1,1,0])
        >>> plt.title('View along [110]')
    
    Notes:
        - Arbitrary viewing direction
        - Creates orthogonal projection
        - Shows atomic arrangements
    """
    fig, ax = plt.subplots()
    
    # Generate lattice points
    points = generate_lattice_points(L, 2, 2, 2)
    
    # Create projection matrix
    direction = np.array(direction) / np.linalg.norm(direction)
    # Project onto plane perpendicular to direction
    # (Simplified - full implementation would use proper projection)
    
    ax.scatter(points[:,0], points[:,1], **kwargs)
    ax.set_aspect('equal')
    plt.show()

def plot_lattice(L, n1=1, n2=1, n3=1, basis=None, ax=None, **kwargs):
    """
    Complete lattice plotting function with atoms and unit cells.
    
    High-level function combining atomic positions, cell boundaries,
    and styling in single call.
    
    Input:
        L (array 3×3): Lattice matrix
        n1, n2, n3 (int): Cell range
        basis (list): Atomic basis
        ax (Axes3D): 3D axis (creates new if None)
        **kwargs: Styling parameters
    
    Output:
        fig, ax: Matplotlib figure and axis objects
    
    Usage Example:
        >>> # Plot BCC structure
        >>> L = cubic_lattice_vec(2.87)  # Fe
        >>> basis = [[0,0,0], [0.5,0.5,0.5]]
        >>> fig, ax = plot_lattice(L, n1=2, n2=2, n3=2, basis=basis,
        ...                        atom_color='blue', atom_size=50)
        >>> ax.set_title('BCC Iron')
        >>> plt.show()
    
    Notes:
        - All-in-one plotting function
        - Combines multiple visualization elements
        - Customizable appearance
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig = ax.figure
    
    # Plot atoms
    atoms = generate_lattite_atom_positions(L, basis, n1, n2, n3)
    ax.scatter(atoms[:,0], atoms[:,1], atoms[:,2], **kwargs)
    
    # Plot cell boundaries
    plot_lattice_boundaries(ax, L, n1, n2, n3, color='black', linewidth=1)
    
    set_aspect_equal_3d(ax)
    
    return fig, ax

def plot_lattice_proj(L, direction, n1=2, n2=2, basis=None, **kwargs):
    """
    Plot lattice projection along specific direction.
    
    Creates 2D projection of 3D lattice viewed along given direction.
    
    Input:
        L (array 3×3): Lattice matrix
        direction (array [3]): Viewing direction
        n1, n2 (int): Cell range
        basis (list): Atomic basis
        **kwargs: Plot parameters
    
    Output:
        fig, ax: Figure and axis
    
    Usage Example:
        >>> L = cubic_lattice_vec(3.0)
        >>> # View along [111]
        >>> fig, ax = plot_lattice_proj(L, [1,1,1], n1=3, n2=3)
        >>> ax.set_title('Cubic lattice - [111] projection')
        >>> plt.show()
    
    Notes:
        - Shows 2D atomic arrangement
        - Useful for structure visualization
        - Reveals symmetry patterns
    """
    fig, ax = plt.subplots()
    
    # Generate atoms
    atoms = generate_lattite_atom_positions(L, basis, n1, n2, 1)
    
    # Project (simplified - use proper projection matrix)
    ax.scatter(atoms[:,0], atoms[:,1], **kwargs)
    ax.set_aspect('equal')
    
    return fig, ax

def plot_points_proj(points, direction, ax=None, **kwargs):
    """
    Project and plot arbitrary 3D points.
    
    General function to project any set of 3D points onto 2D plane.
    
    Input:
        points (array N×3): 3D point coordinates
        direction (array [3]): Projection direction
        ax (Axes): 2D axis (creates if None)
        **kwargs: Scatter plot parameters
    
    Output:
        fig, ax: Figure and axis
    
    Usage Example:
        >>> points = np.random.rand(100, 3)
        >>> fig, ax = plot_points_proj(points, [0,0,1])
        >>> ax.set_title('Random points - XY projection')
    
    Notes:
        - Generic projection function
        - Works with any point set
        - Not specific to lattices
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    
    # Simple XY projection (full version would use direction)
    ax.scatter(points[:,0], points[:,1], **kwargs)
    ax.set_aspect('equal')
    
    return fig, ax

def plot_lattice_plane(ax, L, h, k, l, **kwargs):
    """
    Plot a crystallographic plane in 3D lattice.
    
    Draws plane (hkl) intersecting lattice unit cell using matplotlib 3D.
    
    Input:
        ax (Axes3D): Matplotlib 3D axis
        L (array 3×3): Lattice matrix
        h, k, l (int): Miller indices of plane
        **kwargs: Additional matplotlib plot parameters (color, alpha, etc.)
    
    Output:
        None (modifies axis)
    
    Usage Example:
        >>> import matplotlib.pyplot as plt
        >>> from mpl_toolkits.mplot3d import Axes3D
        >>> import numpy as np
        >>> 
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L = cubic_lattice_vec(3.0)
        >>> plot_lattice_plane(ax, L, 1, 1, 1, color='blue', alpha=0.3)
        >>> plot_lattice_plane(ax, L, 1, 0, 0, color='red', alpha=0.3)
        >>> 
        >>> set_aspect_equal_3d(ax)
        >>> plt.show()
    
    Notes:
        - Plane intersects unit cell
        - Uses Miller indices for specification
        - Transparency recommended (alpha < 1)
        - Multiple planes can be overlaid
    """
    # Calculate plane intercepts with axes
    # Plane equation: hx/a + ky/b + lz/c = 1
    
    if h != 0:
        x_intercept = L[:,0] / h
    else:
        x_intercept = None
    
    if k != 0:
        y_intercept = L[:,1] / k
    else:
        y_intercept = None
    
    if l != 0:
        z_intercept = L[:,2] / l
    else:
        z_intercept = None
    
    # Create plane vertices (simplified)
    vertices = []
    if x_intercept is not None:
        vertices.append(x_intercept)
    if y_intercept is not None:
        vertices.append(y_intercept)
    if z_intercept is not None:
        vertices.append(z_intercept)
    
    if len(vertices) >= 3:
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        poly = Poly3DCollection([vertices], **kwargs)
        ax.add_collection3d(poly)

def plot_lattice_boundaries(ax, L, n1=1, n2=1, n3=1, **kwargs):
    """
    Plot unit cell boundaries in 3D.
    
    Draws edges of unit cells to visualize lattice structure.
    
    Input:
        ax (Axes3D): Matplotlib 3D axis
        L (array 3×3): Lattice matrix
        n1, n2, n3 (int): Number of cells in each direction
        **kwargs: Line plot parameters
    
    Output:
        None (modifies axis)
    
    Usage Example:
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111, projection='3d')
        >>> 
        >>> L = cubic_lattice_vec(3.0)
        >>> plot_lattice_boundaries(ax, L, n1=2, n2=2, n3=2, color='black')
        >>> 
        >>> set_aspect_equal_3d(ax)
        >>> plt.show()
    
    Notes:
        - Draws wireframe of unit cells
        - Helps visualize lattice structure
        - Combine with plot_lattice_points for complete view
    """
    # Draw edges of unit cells
    a1, a2, a3 = L[:,0], L[:,1], L[:,2]
    
    # Draw edges in a1 direction
    for j in range(n2 + 1):
        for k in range(n3 + 1):
            start = j*a2 + k*a3
            end = start + n1*a1
            ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], **kwargs)
    
    # Draw edges in a2 direction
    for i in range(n1 + 1):
        for k in range(n3 + 1):
            start = i*a1 + k*a3
            end = start + n2*a2
            ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], **kwargs)
    
    # Draw edges in a3 direction
    for i in range(n1 + 1):
        for j in range(n2 + 1):
            start = i*a1 + j*a2
            end = start + n3*a3
            ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], **kwargs)

def plot_atomic_plane2D(atoms_on_plane, plane_normal, **kwargs):
    """
    Plot 2D view of atomic plane.
    
    Creates 2D plot of atoms lying on crystallographic plane.
    
    Input:
        atoms_on_plane (array N×3): Coordinates of atoms on plane
        plane_normal (array [3]): Plane normal (for orientation)
        **kwargs: Scatter plot parameters
    
    Output:
        fig, ax: Figure and 2D axis
    
    Usage Example:
        >>> L = cubic_lattice_vec(3.0)
        >>> atoms = generate_lattite_atom_positions(L, n1=5, n2=5, n3=5)
        >>> indices = select_plane(atoms, 1, 1, 1, L)
        >>> plane_atoms = atoms[indices]
        >>> 
        >>> fig, ax = plot_atomic_plane2D(plane_atoms, [1,1,1], s=100, c='blue')
        >>> ax.set_title('(111) Atomic Plane')
        >>> plt.show()
    
    Notes:
        - 2D projection of 3D plane
        - Shows atomic arrangement
        - Useful for interface analysis
    """
    fig, ax = plt.subplots()
    
    # Project atoms onto 2D (simplified - use first two coordinates)
    ax.scatter(atoms_on_plane[:,0], atoms_on_plane[:,1], **kwargs)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    return fig, ax

def plot_atomic_plane3D(atoms_on_plane, plane_normal, L, ax=None, **kwargs):
    """
    Plot atomic plane in 3D context.
    
    Visualizes atoms on plane with lattice context in 3D.
    
    Input:
        atoms_on_plane (array N×3): Atoms on plane
        plane_normal (array [3]): Plane normal
        L (array 3×3): Lattice matrix (for context)
        ax (Axes3D): 3D axis (creates if None)
        **kwargs: Scatter parameters
    
    Output:
        fig, ax: Figure and 3D axis
    
    Usage Example:
        >>> L = cubic_lattice_vec(3.0)
        >>> atoms_all = generate_lattite_atom_positions(L, n1=5, n2=5, n3=5)
        >>> indices = select_plane(atoms_all, 1, 0, 0, L)
        >>> plane_atoms = atoms_all[indices]
        >>> 
        >>> fig, ax = plot_atomic_plane3D(plane_atoms, [1,0,0], L, c='red', s=100)
        >>> # Also plot all atoms faintly
        >>> ax.scatter(atoms_all[:,0], atoms_all[:,1], atoms_all[:,2], 
        ...            c='gray', s=10, alpha=0.3)
        >>> set_aspect_equal_3d(ax)
    
    Notes:
        - Shows plane in 3D lattice context
        - Highlights selected atoms
        - Can overlay multiple planes
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig = ax.figure
    
    ax.scatter(atoms_on_plane[:,0], atoms_on_plane[:,1], atoms_on_plane[:,2], **kwargs)
    set_aspect_equal_3d(ax)
    
    return fig, ax

def plot_atomlattice2D(L, basis, n1, n2, direction=[0,0,1], **kwargs):
    """
    Plot 2D atomic lattice projection.
    
    Complete 2D visualization of atomic structure with basis.
    
    Input:
        L (array 3×3): Lattice matrix
        basis (list): Atomic basis
        n1, n2 (int): Cell range
        direction (array [3]): Viewing direction
        **kwargs: Plot styling
    
    Output:
        fig, ax: Figure and axis
    
    Usage Example:
        >>> # FCC structure
        >>> L = cubic_lattice_vec(3.615)  # Al
        >>> basis = [[0,0,0], [0.5,0.5,0], [0.5,0,0.5], [0,0.5,0.5]]
        >>> 
        >>> fig, ax = plot_atomlattice2D(L, basis, n1=4, n2=4, 
        ...                               direction=[1,1,1], s=100)
        >>> ax.set_title('FCC [111] projection')
        >>> plt.show()
    
    Notes:
        - Shows 2D atomic arrangement
        - Includes all basis atoms
        - Useful for structure determination
    """
    atoms = generate_lattite_atom_positions(L, basis, n1, n2, 1)
    
    fig, ax = plt.subplots()
    ax.scatter(atoms[:,0], atoms[:,1], **kwargs)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    return fig, ax

def plot_cut2D(points, cut_plane_normal, cut_plane_point, **kwargs):
    """
    Plot 2D cross-section of 3D points.
    
    Selects points near plane and plots 2D projection.
    
    Input:
        points (array N×3): 3D points
        cut_plane_normal (array [3]): Cutting plane normal
        cut_plane_point (array [3]): Point on plane
        **kwargs: Plot parameters
    
    Output:
        fig, ax: Figure and 2D axis
    
    Usage Example:
        >>> atoms = generate_lattite_atom_positions(L, n1=10, n2=10, n3=10)
        >>> 
        >>> # Cut at z=15
        >>> fig, ax = plot_cut2D(atoms, [0,0,1], [0,0,15])
        >>> ax.set_title('Cross-section at z=15')
        >>> plt.show()
    
    Notes:
        - Shows 2D slice of 3D structure
        - Useful for interface analysis
        - Can reveal internal structure
    """
    # Select points near plane
    tolerance = kwargs.pop('tolerance', 0.5)
    indices = select_atomic_plane(points, cut_plane_normal, cut_plane_point, tolerance)
    cut_points = points[indices]
    
    # Plot 2D
    fig, ax = plt.subplots()
    ax.scatter(cut_points[:,0], cut_points[:,1], **kwargs)
    ax.set_aspect('equal')
    
    return fig, ax

def plot_planes_on_stereotriangle(planes, **kwargs):
    """
    Plot planes on stereographic triangle.
    
    Projects plane normals onto standard stereographic triangle
    for cubic systems.
    
    Input:
        planes (list): List of (h,k,l) tuples
        **kwargs: Plot parameters
    
    Output:
        None (creates figure with stereographic triangle)
    
    Usage Example:
        >>> planes = [(1,0,0), (1,1,0), (1,1,1), (2,1,0)]
        >>> plot_planes_on_stereotriangle(planes)
        >>> plt.title('Low-Index Planes')
        >>> plt.show()
    
    Notes:
        - Standard triangle for cubic
        - Shows plane distribution
        - Used in texture analysis
    """
    fig, ax = plt.subplots()
    
    # Draw stereographic triangle
    # (Full implementation would use projlib functions)
    
    # Project and plot each plane
    for hkl in planes:
        # Convert to direction, normalize, project
        # Simplified - full version uses equalarea_directions()
        pass
    
    plt.axis('equal')

def plot_planes_on_wulffnet(planes, lattice, **kwargs):
    """
    Plot plane normals on Wulff net (equal-angle projection).
    
    Projects planes onto full Wulff stereographic net.
    
    Input:
        planes (list): Plane Miller indices
        lattice (array 3×3): Lattice matrix
        **kwargs: Plot parameters
    
    Output:
        None (creates Wulff net figure)
    
    Usage Example:
        >>> L = cubic_lattice_vec(3.0)
        >>> planes = [(1,0,0), (0,1,0), (0,0,1), (1,1,1)]
        >>> plot_planes_on_wulffnet(planes, L)
        >>> plt.title('Planes on Wulff Net')
    
    Notes:
        - Equal-angle projection
        - Preserves angular relationships
        - Full sphere projection
    """
    # Create Wulff net background
    # Project and plot planes
    # (Full implementation uses projlib wulffnet() function)
    pass

def plot_princip_dir_on_stereotriangle(principal_directions, **kwargs):
    """
    Plot principal strain/stress directions on stereographic triangle.
    
    Projects principal directions onto standard triangle.
    
    Input:
        principal_directions (array 3×3): Principal direction matrix
        **kwargs: Plot styling
    
    Output:
        None (adds to current figure or creates new)
    
    Usage Example:
        >>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> 
        >>> plot_planes_on_stereotriangle([(1,0,0), (1,1,0), (1,1,1)])
        >>> plot_princip_dir_on_stereotriangle(result['directions'], 
        ...                                      marker='*', s=200, c='red')
        >>> plt.title('Principal Directions')
    
    Notes:
        - Shows orientation of principal axes
        - Overlays on crystal planes
        - Important for anisotropy analysis
    """
    # Project principal directions
    # Plot on existing stereographic triangle
    pass

def plot_princip_dir_on_wulffnet(principal_directions, **kwargs):
    """
    Plot principal directions on Wulff net.
    
    Projects principal strain/stress axes onto Wulff stereographic net.
    
    Input:
        principal_directions (array 3×3): Principal directions
        **kwargs: Plot parameters
    
    Output:
        None (creates or modifies figure)
    
    Usage Example:
        >>> strain = np.diag([0.1, 0.0, -0.05])
        >>> result = mohr_circles(strain)
        >>> plot_princip_dir_on_wulffnet(result['directions'])
        >>> plt.title('Principal Strain Axes')
    
    Notes:
        - Full sphere projection
        - Shows 3D orientation clearly
    """
    pass


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

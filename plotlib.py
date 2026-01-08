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
import numpy as np
try:
    from scipy.interpolate import griddata
except:
    pass
try:
    from wand.image import Image
except:
    pass
from crystlibs import  *
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


    # ============================================================
    # Additional Plotter Methods (Integrated from plotlib_ini5)
    # ============================================================

    def dataAnnot(self,**kwargs):
        """Annotate data points with crystallographic information."""

        self.setAttributes(**kwargs)
        self.annot=True
        if self.onmovetext is None:
            #print('o')
            self.onmovetext=self.fig.text(0.1,0.85,'',fontsize=12)
        self.fig.texts[1].set_visible(True)
        #print("okkkkkkkk")
        #self.fig.canvas.mpl_connect('motion_notify_event', self.onmove)
        #self.fig.canvas.mpl_connect('button_press_event', self.onmove)


    def dataShow(self,**kwargs):
        """Display data information in plot or console."""
        self.setAttributes(**kwargs)
        self.ax.format_coord = self.format_coord

    def format_annot(self,x, y,**kwargs):
        """Format annotation text for data points."""
        
        self.setAttributes(**kwargs)
        HKL=np.linalg.inv(self.LrPhase1)
        UVW=np.linalg.inv(self.LPhase1)
        if self.scatterdata is not None:
            if self.scatteridxs is None:
                scatteridxs=list(np.where(~np.isnan(self.colmapdata))[0])     
            else:
                scatteridxs=self.scatteridxs
        #on pointer motion vizualization defnition
        if True:#try:
            dp=np.sqrt((self.poris[0,:]-x)**2+(self.poris[1,:]-y)**2)
            idx=np.where(dp==min(dp))[0][0]
            basictext=""
            if min(dp)<0.05:
                #print(idx)
                tol=100
                hkl=HKL.dot(self.R2Proj.T.dot(self.oris[:,idx]))
                #npwhere
                hkl=np.round(hkl/min(hkl[np.abs(hkl)>1/tol])*tol)/tol
                dhkl=1/np.sqrt(hkl.T.dot(self.GrA.dot(hkl)))
                uvw=UVW.dot((self.R2Proj.T.dot(self.oris[:,idx])))
                uvw=np.round(uvw/min(uvw[np.abs(uvw)>1/tol])*tol)/tol
                #basictext=f"(h,k,l)$^{{{self.phase1}}}$={str(hkl).replace('[','(').replace(']',')')}$^{{{self.phase1}}}$, [u,v,w]$^{{{self.phase1}}}$={str(uvw)}$^{{{self.phase1}}}$\n"
                basictext=f"(h,k,l)$^{{{self.phase1}}}$=({hkl[0]},{hkl[1]},{hkl[2]})$^{{{self.phase1}}}$, [u,v,w]$^{{{self.phase1}}}$=[{uvw[0]},{uvw[1]},{uvw[2]}]$^{{{self.phase1}}}$\n"
                if self.dhkl:
                    basictext+=f"dhkl={dhkl}"
                if self.printcorrespondent:
                    HKLc=np.linalg.inv(self.LrPhase2)
                    UVWc=np.linalg.inv(self.LPhase2)

                    #self.R2Proj.dot(self.T[:,:,self.varsel])
                    #self.R2Proj.T.dot(self.oris[:,idx]
                    #hklc=vector2miller(HKLc.dot(self.T[:,:,self.varsel].T.dot(self.R2Proj.T.dot(self.oris[:,idx]))),decimals=2,MIN=True)
                    #uvwc=vector2miller(UVWc.dot(self.T[:,:,self.varsel].T.dot(self.R2Proj.T.dot(self.oris[:,idx]))),decimals=2,MIN=True)
                    hklc=HKLc.dot(self.T[:,:,self.varsel].T.dot(self.R2Proj.T.dot(self.oris[:,idx])))
                    uvwc=UVWc.dot(self.T[:,:,self.varsel].T.dot(self.R2Proj.T.dot(self.oris[:,idx])))
                    hklc=np.round(hklc/min(hklc[np.abs(hklc)>1/tol])*tol)/tol
                    uvwc=np.round(uvwc/min(uvwc[np.abs(uvwc)>1/tol])*tol)/tol
                    if self.printcorrespascubicfamily:
                        hklc=np.sort(np.abs(hklc))[::-1]
                        uvwc=np.sort(np.abs(uvwc))[::-1]
                        
                    basictext+=f"${self.brackPlaneCL}$h,k,l${self.brackPlaneCR}$$^{{{self.phase2}}}={self.brackPlaneCL}{hklc[0]},{hklc[1]},{hklc[2]}{self.brackPlaneCR}^{{{self.phase2}}}$"
                    basictext+=f", ${self.brackDirCL}$u,v,w${self.brackDirCR}$$^{{{self.phase2}}}={self.brackDirCL}{uvwc[0]},{uvwc[1]},{uvwc[2]}{self.brackDirCR}^{{{self.phase2}}}$\n"
                    
                    #d2pcs=[]
                    # for CdCp,d2p in zip([self.Cp[:,:,self.varsel],self.Cd[:,:,self.varsel]],[hkl,uvw]):
                    #     d2pc=vector2miller(CdCp.dot(d2p))
                    #     if d2pc[2]<0:
                    #         d2pc=-1*d2pc 
                    #     if self.printcorrespascubicfamily:
                    #         d2pc=np.sort(np.abs(d2pc))[::-1]
                        
                    #     d2pcs.append(np.around(d2pc,decimals=2))
                    #print(d2pcs)
                    #basictext+=f"${self.brackPlaneCL}$h,k,l${self.brackPlaneCR}$$^{{{self.phase2}}}={self.brackPlaneCL}{d2pcs[0][0]},{d2pcs[0][1]},{d2pcs[0][2]}{self.brackPlaneCR}^{{{self.phase2}}}$"
                    #basictext+=f", ${self.brackDirCL}$u,v,w${self.brackDirCR}$$^{{{self.phase2}}}={self.brackDirCL}{d2pcs[1][0]},{d2pcs[1][1]},{d2pcs[1][2]}{self.brackDirCR}^{{{self.phase2}}}$\n"
                if self.scatterdata is not None:
                    for scdi,scattereqhkl in enumerate(self.scattereqhkl):
                        try:
                            scatteridx=scatteridxs.index(idx)
                            scatter1hkl=scattereqhkl[:,scatteridx]
                            #scatter2hkl=self.scatter2eqhkl[:,scatteridx]
                            scatter1hkl=np.sort(np.abs(scatter1hkl))[::-1]
                            try:
                                hpbpoint1.remove()
                            except:
                                pass
                        except:
                            scatter1hkl='No habit plane'  
                            #scatter2hkl='No habit plane'  
                        #basictext+=f", lenidx={len(scatteridxs)}, idx={idx}"
                        if self.ax not in self.ax2annot: 
                            basictext+=f"HBP{int(scdi+1)}={scatter1hkl}$^A$,".replace("[","(").replace("]",")")
                    basictext=basictext[:-1]
                    basictext+="\n"
                if self.showdatanames is not None:
                    if self.ax2annot is None or self.ax not in self.ax2annot:
                        for name,data in zip(self.showdatanames,self.showdata):                    
                            basictext+=f"{name}:{np.around(data[idx],decimals=3)}\n"
                    #print('ok')
                elif self.showdataasplotted and self.pickax is self.ax:
                    #basictext+=f"{self.showdatanames[0]}:{np.around(self.showdata[0][idx],decimals=3)}\n"
                    basictext+=f"{self.cbartitle}:{np.around(self.colmapdata[idx],decimals=3)}" 
        #except:
        #    basictext="Problem with data formater"
        basictext = r'{}'.format(basictext)      
        return basictext    


    def format_coord(self,x, y,**kwargs):
        """Format coordinates for display in plot toolbar."""
        self.setAttributes(**kwargs)
        #Meteric tensor ||[x,y,x]||^2=[uvw]^T*G_A*[uvw]
        LA=self.LPhase1
        self.GA = np.matmul(LA.T,LA)

        #Reciprocal Meteric tensor ||[x,y,x]||^2=d_hkl^2=[hkl]^T*Gr_A*[hkl]
        self.GrA = inv(self.GA)

        
        #Conversion to hkl, uvw
        basictext=self.format_annot(x, y)
        if basictext!=" ":
            if self.annot:    
                self.fig.texts[1].set_text(basictext)
                if self.withdraw:
                    self.fig.canvas.draw()
            self.annotationtext=basictext
        basictext="okook"
        return basictext

    def format_coord_test(self,x, y,**kwargs):
        """Test version of coordinate formatting for debugging."""

        return f"{x},{y}"

    def genPoris(self,**kwargs):
        """
            Generate pole figure orientation data for given crystal direction.
            
            Creates dense sampling of orientations where specified crystal direction
            aligns with projection direction, useful for inverse pole figures.
            
            Input:
                cr_dir: array (3,) - Crystal direction [u, v, w]
                symops: list - Crystal symmetry operations
                resol: int - HEALPix resolution (default: 4)
            
            Output:
                oris: array (N, 3, 3) - Orientation matrices
            
            Usage Example:
                >>> import numpy as np
                >>> from plotlib import plotter
                >>> 
                >>> # Generate orientations for [001] pole figure
                >>> p = plotter()
                >>> cubic_symops = [np.eye(3)]  # Simplified
                >>> oris_001 = p.genPoris([0, 0, 1], cubic_symops, resol=3)
            """
        if self.oris is not None:
            if self.ProjType=='equalarea':
                equalarea=True
                if self.oris2 is not None:
                    self.poris2 = equalarea_directions(self.oris2)
                self.poris = equalarea_directions(self.oris)
            else:
                equalarea=False
                if self.oris2 is not None:
                    self.poris2 = stereoprojection_directions(self.oris2)
                self.poris = stereoprojection_directions(self.oris)




    def generateSphericalHistSampleData(self,**kwargs):
        """
            Generate spherical histogram from orientation data.
            
            Bins orientations on sphere using equal-area binning scheme,
            useful for discrete texture representation.
            
            Input:
                oris: array (N, 3, 3) - Orientation matrices
                cr_dir: array (3,) - Crystal direction for projection
                bins: int - Number of bins (default: 128)
                histnorm: bool - Normalize histogram (default: True)
                symops: list - Symmetry operations
            
            Output:
                histogram: array - Binned density values
            
            Usage Example:
                >>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Generate histogram
                >>> orientations = R.random(1000).as_matrix()
                >>> 
                >>> p = plotter()
                >>> hist = p.generateSphericalHistSampleData(
                ...     orientations,
                ...     cr_dir=[0, 0, 1],
                ...     bins=64
                ... )
            """
        self.setAttributes(**kwargs)       
        
        #xy = self.sampleoris[0,:]**2 + self.sampleoris[1,:]**2
        #theta = np.arctan2(np.sqrt(xy), self.sampleoris[2,:])#elevation
        #phi=np.arctan2(self.sampleoris[1,:], self.sampleoris[0,:])
        theta,phi=xyz2spher(self.sampleoris.T)
        hist, phiv, thtv = np.histogram2d(phi,theta, bins=self.sphericalbins,range=[[-np.pi,np.pi], [0,np.pi]],weights=1/np.sin(theta))
        phig, thtg = np.meshgrid((phiv[:-1] + phiv[1:])/2.,
                    (thtv[:-1] + thtv[1:])/2.)
        self.oris=spher2xyz(thtg.flatten(),phig.flatten()).T
        Norm=(hist.flatten()*np.sin(thtg.flatten())*np.diff(phiv)[0]**2/4/np.pi).sum()
        
        self.colmapdata=hist.T.flatten()/Norm
        

    def generateSphericalKDESampleData(self,**kwargs):
        """
            Generate spherical kernel density estimate from orientation data.
            
            Computes density using Von Mises-Fisher distributions on the sphere,
            accounting for crystal symmetry and bandwidth optimization.
            
            Input:
                oris: array (N, 3, 3) - Orientation matrices
                cr_dir: array (3,) - Crystal direction for projection
                sdbandwidth: float - KDE bandwidth (default: 0.15)
                sdweights: array - Sample weights (optional)
                symops: list - Symmetry operations
            
            Output:
                density: array - Density values on sphere
            
            Usage Example:
                >>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Generate KDE from orientations
                >>> orientations = R.random(500).as_matrix()
                >>> 
                >>> p = plotter()
                >>> p.setAttributes(oris=orientations)
                >>> density = p.generateSphericalKDESampleData(
                ...     orientations,
                ...     cr_dir=[0, 0, 1],
                ...     sdbandwidth=0.15
                ... )
            """
        self.setAttributes(**kwargs)       
        
        #xy = self.sampleoris[0,:]**2 + self.sampleoris[1,:]**2
        #theta = np.arctan2(np.sqrt(xy), self.sampleoris[2,:])#elevation
        #phi=np.arctan2(self.sampleoris[1,:], self.sampleoris[0,:])
        theta,phi=xyz2spher(self.sampleoris.T)
        

        #calculation of the logarithm of kernel density function of inverse poles using corresponding spherical coordinates theta (elevation), phi
        #I do not know why SphericalKDE provides logarithm
        logkde = SphericalKDE(phi, theta, bandwidth=self.sdbandwidth,weights=self.sdweights)
        sigmahat = VonMises_std(phi, theta)
        weights = np.ones_like(phi)
        #print(1.06*sigmahat*len(weights)**-0.2)
        #print(len(weights))
        #print(sigmahat)
        #generation of orientation whithin whole space where we will evaluate logkde
        dangle=self.histdangle
        self.oris=genori(dangle=dangle,hemi='both', tol=1e-2, rot=np.eye(3), half='no')
        xy = self.oris[0,:]**2 + self.oris[1,:]**2
        thtg = np.arctan2(np.sqrt(xy), self.oris[2,:])#elevation==polar angle
        phig=np.arctan2(self.oris[1,:], self.oris[0,:])#azimuth angle
        #evaluation of kde for all orientations
        H=np.exp(logkde(phig, thtg))    
        
        #Normalization of H for MRD
        #normalizing data to multiples of random distribution
        #hist/(sum(hist)*pxarea)*numpixels*pxarea  
        #integral Hds=1
        #Norm=H.sum()/H.shape[0]
        Norm=(H*np.sin(thtg)*(dangle*np.pi/180)**2/4/np.pi).sum()
        self.colmapdata=H/Norm
        #print(self.colmapdata.sum()/H.shape[0])
        #H=H/H.sum()*H.shape[0]
        #check normalization
        #print(f'Norm={(self.colmapdata*np.sin(thtg)*(dangle*np.pi/180)**2/4/np.pi).sum()}')

        

    def getColormap(self,**kwargs):    
        """
            Generate orientation density colormap from sample data.
            
            Computes orientation density using spherical KDE or histograms,
            interpolates to grid, and prepares colormap data.
            
            Input:
                oris: array (N, 3, 3) - Orientation matrices
                nump: int - Grid resolution (default: 1001)
                colmapFromSampleData: bool - Generate from data (default: False)
                sdbandwidth: float - KDE bandwidth (default: 0.15)
                sdweights: array - Sample weights (optional)
            
            Output:
                colmapdata: array - Density values on grid
            
            Usage Example:
                >>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Generate colormap from orientations
                >>> orientations = R.random(1000).as_matrix()
                >>> 
                >>> p = plotter()
                >>> p.setAttributes(oris=orientations)
                >>> colormap_data = p.getColormap(nump=501)
            """
        if self.colmapFromSampleData:
            self.generateHistSampleData()
        self.setAttributes(**kwargs)
        if self.oris is not None:
            if self.ProjType=='equalarea':
                equalarea=True
                if self.oris2 is not None:
                    self.poris2 = equalarea_directions(self.oris2)
                self.poris = equalarea_directions(self.oris)
            else:
                equalarea=False
                if self.oris2 is not None:
                    self.poris2 = stereoprojection_directions(self.oris2)
                self.poris = stereoprojection_directions(self.oris)

        #Grid data
        titles=['gx','gy','gz','nummask','mask']

        if self.datadeviders is None:
            datadeviders=[[0,self.oris.shape[1]]]
        else:
            datadeviders=self.datadeviders
        if self.colormapdata is None:
            self.colormapdata=[]
        for datadevider in datadeviders: 
            #print('ok')
            gx,gy, gz, nummask, mask=genprojgrid(self.oris,#self.oris[:,datadevider[0]:datadevider[1]],
                                                    gdata=self.colmapdata[datadevider[0]:datadevider[1]],
                                                    nump=self.nump,proj=self.ProjType,method2='linear')
            #print(gx)
            #print(self.nump)
            if self.sphere=='half' and self.cut2half:
                
                mask[(gy<0)]=True
                nummask[(gy<0)]=0
                #nummask[(gz==np.nan)]=0
                idxs1=np.where(nummask.sum(axis=1)>1)[0]
                idxs2=np.where(nummask.sum(axis=0)>0)[0]
                #print(idxs1)
                #print(idxs2)
                self.nummask=nummask
                gz=gz[np.ix_(idxs1,idxs2)]
                gy=gy[np.ix_(idxs1,idxs2)]
                gx=gx[np.ix_(idxs1,idxs2)]
                #print('ok')
            if self.sphere=='triangle' and self.cut2triangle:
                
                oristri=genoritri(resolution=0.1,mesh="spherified_cube_edge")
                grid_x22,grid_y22, nummask2, mask2=genprojgrid(oristri,#self.oris[:,datadevider[0]:datadevider[1]],
                                                                nump=self.nump,proj=self.ProjType,minmax='full')
                
                nummask=nummask2
                mask=mask2
                gz.mask=mask2
                self.nummask=nummask2
                idxs1=np.where(self.nummask.sum(axis=1))[0]
                idxs2=np.where(self.nummask.sum(axis=0))[0]
                gz=gz[np.ix_(idxs1,idxs2)]
                gy=gy[np.ix_(idxs1,idxs2)]
                gx=gx[np.ix_(idxs1,idxs2)]
                #print(self.ProjType)
                #mask[(gy<0)]=True
                #nummask[(gy<0)]=1
                #mask[(gy<0)+(gx<0)+((np.divide(gy,gx,out=np.zeros_like(gy)+2,where=gx!=0))>1)]=True
                #nummask[(gy<0)+(gx<0)+((np.divide(gy,gx,out=np.zeros_like(gy)+2,where=gx!=0))>1)]=1
                #an=np.arctan(np.divide(gy,gx,out=np.zeros_like(gy),where=gx!=0))
                #an[gx==0]=np.pi/2
                #dirsxy1=np.cos(an)
                #dirsxy2=np.sin(an)
                
                #print(np.nanmin(alpha))
                #print(np.nanmax(alpha))
                #alpha=copy.deepcopy(gx)
                #alpha[gx!=0]=np.arcsin(gx[gx!=0]/dirsxy1[gx!=0]/2)*2
                #alpha[gx==0]=np.arcsin(gy[gx==0]/dirsxy2[gx==0]/2)*2
                #nummask[np.abs(alpha)>(np.pi/4)]=1
                #mask[np.abs(alpha)>(np.pi/4)]=True
            
        
            colormapdata={}
            for title in titles:
                exec(f'colormapdata["{title}"]={title}')  
            self.colormapdata.append(colormapdata)
            #print(colormapdata)
        

    def onclicactivate(self,**kwargs):
        """Activate mouse click event handling for the plot."""
        #print(self.fig)
        self.fig.canvas.mpl_connect('button_press_event', self.onclick2);


    def onclick(self,event):
        """Mouse click event handler for data selection."""
        #print('ok')
        x=event.xdata
        y=event.ydata
        self.pickax=event.inaxes
        basictext=self.format_annot(x, y)
        if basictext!=" ":
            if self.annot and self.pickax is self.ax:    
                self.fig.texts[1].set_text(basictext)
                self.fig.canvas.draw()
        if self.ax2annot is None or event.inaxes not in self.ax2annot:
            try:
                if self.oris2 is not None:
                    dp=np.sqrt((self.poris2[0,:]-x)**2+(self.poris2[1,:]-y)**2)
                else:
                    dp=np.sqrt((self.poris[0,:]-x)**2+(self.poris[1,:]-y)**2)
            except:
                dp=[10.]
            
            if min(dp)<0.05 and True:
                idx=np.where(dp==min(dp))[0][0]
            #if True:
                #for annot in self.ANNOTS:
                    #if self.ax2annot is None or self.ax in self.ax2annot:
                    #annot.xy = (0,0)
                    #annot.set_text(f'{ttt}:{min(dp)},idx:{idx},projlib:{self.scatterproj[0].shape}')
                    #annot.set_visible(True)
            if min(dp)<0.05 and True:
                idx=np.where(dp==min(dp))[0][0]
                        
                try:
                    #for annot in self.ANNOTS:
                    #    annot.set_visible(False)
                    scatx=[]
                    scaty=[]
                    showa=False
                    scatteridx=self.scatteridxs.index(idx)
                    for scattereqhkl,scatterproj,annot in zip(self.scattereqhkl,self.scatterproj,self.ANNOTS):
                        scatter1hkl=scattereqhkl[:,scatteridx]
                        if self.ax2annot is None or self.ax in self.ax2annot:
                            
                            if self.ax2annotcompres is not None and self.scattercolscale[scatteridx]<0:
                                if self.ax in self.ax2annotcompres:
                                    showannot=True
                                    gid=50000
                            if self.ax2annottension is not None and self.scattercolscale[scatteridx]>0:
                                if self.ax in self.ax2annottension:
                                    showannot=True
                                    gid=51000
                            if showannot:   
                                showa=True
                                scatx.append(scatterproj[0,scatteridx])
                                scaty.append(scatterproj[1,scatteridx])
                                #self.ax.plot(scatterproj[0,scatteridx],scatterproj[1,scatteridx],marker="o",markersize=8,markeredgecolor='k',markerfacecolor='k',gid=gid)
                                #self.fig.canvas.draw()
                                annot.xy = (scatterproj[0,scatteridx],scatterproj[1,scatteridx])
                                annot.set_text(f'Strain:{np.around(self.scattercolscale[scatteridx],decimals=1)},'+str(np.sort(np.abs(scatter1hkl))[::-1]).replace("[","(").replace("]",")")+'$^A$')
                                annot.set_visible(True)
                                #for scat in self.scat:
                                #    scat.set_visible(False)
                    if showa:
                        for curvei in self.ax.get_lines():
                            if curvei.get_gid()==gid:
                                curvei.remove()  
                        self.ax.plot(scatx,scaty,'o',markersize=8,markeredgecolor='k',markerfacecolor='k',gid=gid)
                        self.fig.canvas.draw()
                        
                    
                except:
                    pass
            
                if self.withdraw:
                    self.fig.canvas.draw()
                
            
        

    def onclick2(self,event):
        """Alternative click handler for different interaction mode."""
        plt.draw()
        #x=event.xdata
        #y=event.ydata
        #basictext=self.format_annot(x, y)
        self.fig.text(0.1,0.85,'ppppppppppppppppppppppp',fontsize=12)
        #if pppppp:
            #print('ok')
        
        #basictext="neco"
        #print("ook")
        #self.fig.texts[1].set_text(basictext)
        #if basictext!=" ":
        #    if self.annot:    
        #        self.fig.texts[1].set_text(basictext)
        #        self.fig.canvas.draw()
        plt.draw()

    def onclick3(self,event):
        """Third click handler variant for specialized interactions."""
        #print('ok')
        x=event.xdata
        y=event.ydata
        self.pickax=event.inaxes
        basictext=self.format_annot(x, y)
        if basictext!=" ":
            if self.annot and self.pickax is self.ax:    
                self.fig.texts[1].set_text(basictext)
                self.fig.canvas.draw()

    def onmove(self,event):  
        """Mouse move event handler for interactive data display."""
        if event.inaxes:
            #print(dir(event))
            self.fig.texts[1].set_text(self.format_coord(event.xdata,event.ydata))
            self.fig.texts[1].set_visible(True)
        else:
            self.fig.texts[1].set_visible(False)
        
        plt.draw()


    def onpress(self,event):
        """Keyboard press event handler for plot interaction."""
        sys.stdout.flush()
        #print(event.key)
        if event.key == ' ':
            for scat in self.scat:
                scat.set_visible(not scat.get_visible())
            #if self.withdraw:
            self.fig.canvas.draw()
        if event.key == 'alt':
            if self.ax2annot is not None:
                for ax in self.ax2annot:
                    for curvei in ax.get_lines():
                        if curvei.get_gid()==50000 or curvei.get_gid()==51000:
                            curvei.remove()
            for annot in self.ANNOTS:
                annot.set_visible(False)
            #if self.withdraw:
            self.fig.canvas.draw()


    def onpressActivate(self):
        """Activate keyboard event handling for the plot."""
        self.fig.canvas.mpl_connect('key_press_event', self.onpress)

    def plotColorbar(self,**kwargs):
        """
            Add colorbar to the plot with custom ticks and labels.
            
            Creates a horizontal or vertical colorbar with formatted tick labels,
            including support for inequality symbols and custom positioning.
            
            Input:
                cbartitle: str - Colorbar title (default: "")
                vmbar: list [min, max] - Value range (default: None, auto-detect)
                cbarh: float - Colorbar height (default: 0.04)
                cbarwfac: float - Width factor (default: 0.75)
                cbarhshift: float - Horizontal shift (default: -0.15)
                ticks: array - Tick positions (optional)
                ticklabels: list - Tick labels (optional)
            
            Output:
                None (adds colorbar to self.fig)
            
            Usage Example:
                >>> from plotlib import plotter
                >>> import numpy as np
                >>> 
                >>> # Basic colorbar
                >>> p = plotter()
                >>> p.plotProj(ProjType='equalarea')
                >>> orientations = np.random.randn(500, 3, 3)
                >>> p.setAttributes(oris=orientations)
                >>> p.plotColormap()
                >>> p.plotColorbar(cbartitle="MUD")
            """
        self.setAttributes(**kwargs)
        try:
            if  self.scatterplot:
                scbar=self.scat[0]
                plotbar=True
                vmcolbar=self.vmbar#scattercolscalevm
                vmcolbarticks=self.scattercolscaleticks
            elif self.scatterplotashist:
                scbar=self.histplot
                plotbar=True
                vmcolbar=None#self.vmbar#scattercolscalevm
                vmcolbarticks=None#self.scattercolscaleticks
                
            else:
                scbar=self.colormap
                plotbar=True
                vmcolbar=self.vmbar
                vmcolbarticks=self.ticks
            pos = self.ax.get_position()
            self.cbar_ax = self.fig.add_axes([pos.width*(1-self.cbarwfac)/2+pos.x0, pos.y0+self.cbarh+self.cbarhshift, pos.width*self.cbarwfac,self.cbarh])
            self.cbar = self.fig.colorbar(scbar, cax=self.cbar_ax, orientation='horizontal')  
            if self.labelpad is None:
                self.cbar.ax.set_xlabel(self.cbartitle)
            else:
                self.cbar.ax.set_xlabel(self.cbartitle,labelpad=self.labelpad)
            if not vmcolbar is None:
                self.cbar.ax.set_xlim(vmcolbar)
            if not vmcolbarticks is None:
                self.cbar.ax.set_xticks(vmcolbarticks)
            if self.ticklabels is not None:
                self.cbar.ax.set_xticklabels(self.ticklabels)
            self.cbar_ax.xaxis.set_ticks_position('bottom')
            self.cbar_ax.xaxis.set_label_position('bottom')
        except:
            pass


    def plotColormap(self,**kwargs):
        """
            Plot orientation density colormap with contours on stereographic projection.
            
            Creates density maps from orientation data using spherical KDE or histograms,
            with contour lines showing texture intensity.
            
            Input:
                nump: int - Number of points for grid (default: 1001)
                contourcol: str - Contour line color (default: 'k')
                plotmap: bool - Whether to plot colormap (default: True)
                colmapdata: array - Pre-computed colormap data (optional)
            
            Output:
                None (plots colormap on self.ax)
            
            Usage Example:
                >>> import numpy as np
                >>> from plotlib import plotter
                >>> 
                >>> # Create plotter with orientation data
                >>> p = plotter()
                >>> orientations = np.random.randn(1000, 3, 3)  # Random orientations
                >>> 
                >>> # Setup projection
                >>> p.plotProj(ProjType='equalarea', sphere='half')
                >>> 
                >>> # Set orientation data
                >>> p.setAttributes(oris=orientations)
                >>> 
                >>> # Plot colormap with contours
                >>> p.plotColormap(nump=501, contourcol='black')
                >>> 
                >>> # High-resolution colormap
                >>> p2 = plotter()
                >>> p2.setAttributes(oris=orientations, ProjType='equalarea')
                >>> p2.plotProj()
                >>> p2.plotColormap(nump=1001)
            """
        if self.colmapFromSampleData:
            self.generateHistSampleData()
        self.setAttributes(**kwargs)
        if self.oris is not None:
            if self.ProjType=='equalarea':
                equalarea=True
                if self.oris2 is not None:
                    self.poris2 = equalarea_directions(self.oris2)
                self.poris = equalarea_directions(self.oris)
            else:
                equalarea=False
                if self.oris2 is not None:
                    self.poris2 = stereoprojection_directions(self.oris2)
                self.poris = stereoprojection_directions(self.oris)
        if self.plotmap:
            #Grid data
            titles=['gx','gy','gz','nummask','mask']
            if not self.colmapdata is None and self.plotmap:
                #print(self.colmapdata)
                if self.colormapdata is None:
                    if self.datadeviders is None:
                        datadeviders=[[0,self.oris.shape[1]]]
                    else:
                        datadeviders=self.datadeviders
                    if self.colormapdata is None:
                        colormapdata=[]
                    if self.colormapdata is None:
                        for datadevider in datadeviders: 
                            #print('ok')
                            gx,gy, gz, nummask, mask=genprojgrid(self.oris,#self.oris[:,datadevider[0]:datadevider[1]],
                                                                    gdata=self.colmapdata[datadevider[0]:datadevider[1]],
                                                                    nump=self.nump,proj=self.ProjType,method2='linear')
                            #print(self.nump)
                            if self.sphere=='half' and self.cut2half:
                                mask[(gy<0)]=True
                                nummask[(gy<0)]=0
                                #nummask[(gz==np.nan)]=0
                                idxs1=np.where(nummask.sum(axis=1))[0]
                                idxs2=np.where(nummask.sum(axis=0))[0]
                                #print(idxs1)
                                #print(idxs2)
                                self.nummask=nummask
                                gz=gz[np.ix_(idxs1,idxs2)]
                                gy=gy[np.ix_(idxs1,idxs2)]
                                gx=gx[np.ix_(idxs1,idxs2)]
                            if self.sphere=='triangle' and self.cut2triangle:
                                
                                oristri=genoritri(resolution=0.1,mesh="spherified_cube_edge")
                                grid_x22,grid_y22, nummask2, mask2=genprojgrid(oristri,#self.oris[:,datadevider[0]:datadevider[1]],
                                                                                nump=self.nump,proj=self.ProjType,minmax='full')
                                nummask=nummask2
                                mask=mask2
                                gz.mask=mask2
                                self.nummask=nummask2
                                idxs1=np.where(self.nummask.sum(axis=1))[0]
                                idxs2=np.where(self.nummask.sum(axis=0))[0]
                                gz=gz[np.ix_(idxs1,idxs2)]
                                gy=gy[np.ix_(idxs1,idxs2)]
                                gx=gx[np.ix_(idxs1,idxs2)]
                                #print(self.ProjType)
                                #mask[(gy<0)]=True
                                #nummask[(gy<0)]=1
                                #mask[(gy<0)+(gx<0)+((np.divide(gy,gx,out=np.zeros_like(gy)+2,where=gx!=0))>1)]=True
                                #nummask[(gy<0)+(gx<0)+((np.divide(gy,gx,out=np.zeros_like(gy)+2,where=gx!=0))>1)]=1
                                #an=np.arctan(np.divide(gy,gx,out=np.zeros_like(gy),where=gx!=0))
                                #an[gx==0]=np.pi/2
                                #dirsxy1=np.cos(an)
                                #dirsxy2=np.sin(an)
                                
                                #print(np.nanmin(alpha))
                                #print(np.nanmax(alpha))
                                #alpha=copy.deepcopy(gx)
                                #alpha[gx!=0]=np.arcsin(gx[gx!=0]/dirsxy1[gx!=0]/2)*2
                                #alpha[gx==0]=np.arcsin(gy[gx==0]/dirsxy2[gx==0]/2)*2
                                #nummask[np.abs(alpha)>(np.pi/4)]=1
                                #mask[np.abs(alpha)>(np.pi/4)]=True

                            if self.colormapdata is None:
                                colormapdata.append({})
                                for title in titles:
                                    exec(f'colormapdata[-1]["{title}"]={title}')                               
                        self.colormapdata= colormapdata
                for colmapdat in self.colormapdata:
                    #print(colmapdat)
                    #for title in titles:
                        #exec(f'{title}=colmapdat["{title}"]')
                        
                    if not self.vm is None:
                        self.colormap=self.ax.pcolor(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],cmap=self.cmap,vmin=self.vm[0],vmax=self.vm[1],rasterized=True)
                    else:
                        #print('rast')
                        self.colormap=self.ax.pcolor(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],cmap=self.cmap,rasterized=True)
                    
                    
                    if self.levels is not None:
                        levels=self.levels
                    elif self.ticks is not None:
                        levels=self.ticks
                    else:
                        levels=9
                    #print(levels)
                    if not self.vm is None:    
                        try:
                            self.contours=self.ax.contour(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],levels=levels,colors=self.contourcol,vmin=self.vm[0],vmax=self.vm[1],
                                                        linewidths=self.linewidths,**kwargs)
                        except:
                            self.contours=self.ax.contour(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],levels=levels,colors=self.contourcol,vmin=self.vm[0],vmax=self.vm[1],
                                                        linewidths=self.linewidths)
                    else:
                        try:
                            self.contours=self.ax.contour(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],levels=levels,colors=self.contourcol,linewidths=self.linewidths,**kwargs)
                        except:
                            self.contours=self.ax.contour(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],levels=levels,colors=self.contourcol,linewidths=self.linewidths)

                    if self.contourlabel: 
                        self.ax.clabel(self.contours, fontsize=self.contourfontsize, inline=1)
            
        #self.ax.set_xlim([(1)*self.ax.get_xlim()[0]-0.05,1.05*self.ax.get_xlim()[1]])
        
        if abs(self.ax.get_xlim()[0])<0.01:
            self.ax.set_xlim([-.05+self.ax.get_xlim()[0],1.05*self.ax.get_xlim()[1]])
        else:
            self.ax.set_xlim([1.05*self.ax.get_xlim()[0],1.05*self.ax.get_xlim()[1]])
        
        if abs(self.ax.get_ylim()[0])<0.01:
            self.ax.set_ylim([-.05+self.ax.get_ylim()[0],1.05*self.ax.get_ylim()[1]])
        else:
            self.ax.set_ylim([1.05*self.ax.get_ylim()[0],1.05*self.ax.get_ylim()[1]])


    def plotColormaps(self,sel='all',**kwargs):
        """
            Plot multiple colormaps in multi-panel figure (typically 2x2).
            
            Creates comparative visualization of different pole figures,
            useful for showing texture evolution or multi-phase analysis.
            
            Input:
                Multiple attributes set via setAttributes for each panel
            
            Output:
                None (creates multi-panel figure)
            
            Usage Example:
                >>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Create 2x2 pole figure comparison
                >>> import matplotlib.pyplot as plt
                >>> fig, axes = plt.subplots(2, 2, figsize=(12, 12))
                >>> 
                >>> # Four different texture states
                >>> oris1 = R.random(500).as_matrix()
                >>> oris2 = R.random(500).as_matrix()
            """
        if sel=='all':
            idxs=range(len(self.colormapdata))
        elif type(sel)!=list:
            idxs=[sel]
        else:
            idxs=sel
        for idx in idxs:
            colmapdat = self.colormapdata[idx]
            #print(colmapdat)
            #for title in titles:
                #exec(f'{title}=colmapdat["{title}"]')
                
            if not self.vm is None:
                self.colormap=self.ax.pcolor(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],cmap=self.cmap,vmin=self.vm[0],vmax=self.vm[1],rasterized=True)
            else:
                #print('rast')
                self.colormap=self.ax.pcolor(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],cmap=self.cmap,rasterized=True)
            
            
            if self.levels is not None:
                levels=self.levels
            elif self.ticks is not None:
                levels=self.ticks
            else:
                levels=9
            #print(levels)
            if not self.vm is None:    
                try:
                    self.contours=self.ax.contour(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],levels=levels,colors=self.contourcol,vmin=self.vm[0],vmax=self.vm[1],
                                                linewidths=self.linewidths,**kwargs)
                except:
                    self.contours=self.ax.contour(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],levels=levels,colors=self.contourcol,vmin=self.vm[0],vmax=self.vm[1],
                                                linewidths=self.linewidths)
            else:
                try:
                    self.contours=self.ax.contour(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],levels=levels,colors=self.contourcol,linewidths=self.linewidths,**kwargs)
                except:
                    self.contours=self.ax.contour(colmapdat['gx'],colmapdat['gy'],colmapdat['gz'],levels=levels,colors=self.contourcol,linewidths=self.linewidths)

            if self.contourlabel: 
                self.ax.clabel(self.contours, fontsize=self.contourfontsize, inline=1)
            
        #self.ax.set_xlim([(1)*self.ax.get_xlim()[0]-0.05,1.05*self.ax.get_xlim()[1]])
        
        if abs(self.ax.get_xlim()[0])<0.01:
            self.ax.set_xlim([-.05+self.ax.get_xlim()[0],1.05*self.ax.get_xlim()[1]])
        else:
            self.ax.set_xlim([1.05*self.ax.get_xlim()[0],1.05*self.ax.get_xlim()[1]])
        
        if abs(self.ax.get_ylim()[0])<0.01:
            self.ax.set_ylim([-.05+self.ax.get_ylim()[0],1.05*self.ax.get_ylim()[1]])
        else:
            self.ax.set_ylim([1.05*self.ax.get_ylim()[0],1.05*self.ax.get_ylim()[1]])        
        

    def plotDirsNorms(self,**kwargs):
        """
            Plot crystallographic directions and plane normals on stereographic projection.
            
            Displays Miller indices with appropriate notation, supports family notation,
            and handles multiple phases with correspondence relationships.
            
            Input:
                dirs: list - Crystal directions [[u,v,w], ...] (optional)
                norms: list - Plane normals [[h,k,l], ...] (optional)
                printasfamily: bool - Use family notation {hkl} <uvw> (default: True)
                dhkl: bool - Use decimal Miller indices (default: False)
                Cd: array (3,3) - Direction transformation matrix (default: eye(3))
                Cp: array (3,3) - Plane transformation matrix (default: eye(3))
            
            Output:
                None (plots directions and normals on self.ax)
            
            Usage Example:
                >>> from plotlib import plotter
                >>> import numpy as np
                >>> 
                >>> # Basic pole figure with cubic directions
                >>> p = plotter()
                >>> p.plotProj(ProjType='equalarea', sphere='half')
                >>> 
                >>> # Plot low-index directions
                >>> dirs = [[1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,1,1]]
                >>> p.setAttributes(dirs=dirs)
                >>> p.plotDirsNorms()
                >>> 
                >>> # Plot with plane normals
                >>> norms = [[1,1,1], [1,1,0], [1,0,0]]
                >>> p.setAttributes(dirs=dirs, norms=norms)
                >>> p.plotDirsNorms()
            """
        self.setAttributes(**kwargs)
        if self.sphere=='full':
            self.textlim=-10    
        else:
            self.textlim=0 
        if self.printcorrespondentpoints:
            legbackg = [0.8,0.8,0.8]
            #legbackg = [1,1,1]
            legbackg2 = (0.8,0.8,0.8,0.)
            legsizecoeff=1.
            if self.correspondentpointsize!=self.scatterpointsize:
                legsizecoeff=1.4
                
            
            legend_elements = [Line2D([0], [0], color=legbackg2,marker='o', label=self.phase1,
                                markerfacecolor='k', markeredgecolor='w',markersize=self.scatterpointsize/70*10),
                                Line2D([0], [0], color=legbackg2,marker='o', label=self.phase2,
                                markerfacecolor='w', markeredgecolor='k',markersize=self.correspondentpointsize/self.scatterpointsize*10*legsizecoeff)]
            if not self.plot1phasedirs:
                legend_elements = [Line2D([0], [0], color=legbackg2,marker='o', label=self.phase2,
                                    markerfacecolor='w', markeredgecolor='k',markersize=self.correspondentpointsize/self.scatterpointsize*10*legsizecoeff)]
        if self.dirnormeq:      
            dirs,normals=gen_dirs_norms(self.LPhase1, self.LrPhase1,self.dirs,self.norms,R2Proj=self.R2Proj, recsymops=self.recsymops,symops=self.symops)
        else:
            dirs,normals=gen_dirs_norms(self.LPhase1, self.LrPhase1,self.dirs,self.norms,R2Proj=self.R2Proj, recsymops=[],symops=[])
        
        #for norm in normals:
        #    print(norm['hkl'])
        uvwkeys=[f"{int(d['uvw'][0])}{int(d['uvw'][1])}{int(d['uvw'][2])}" for d in dirs]
        hklkeys=[f"{int(d['hkl'][0])}{int(d['hkl'][1])}{int(d['hkl'][2])}" for d in normals]
        self.defaultdirtexthifts={key:[0,0] for key in uvwkeys}
        self.defaultnormtexthifts={key:[0,0] for key in hklkeys}
        self.scatterdirsfacecolors={key:self.scatterfacecolor for key in uvwkeys}
        self.scatternormsfacecolors={key:self.scatterfacecolor for key in hklkeys}
        self.scatterdirsedgecolors={key:self.scatteredgecolor for key in uvwkeys}
        self.scatternormsedgecolors={key:self.scatteredgecolor for key in hklkeys}
        #print(self.dirs)
        #print(self.scatterdirsfacecolors)
        self.scatterdirsfacecolors.update(self.dirsfacecolors)
        #print(self.scatterdirsfacecolors)
        self.scatterdirsedgecolors.update(self.dirsedgecolors)
        self.scatternormsfacecolors.update(self.normsfacecolors)
        self.scatternormsedgecolors.update(self.normsedgecolors)
        
        #print(self.dirsfacecolors)
        
        
        if len(self.dirtexthifts)!=0:
            for key in self.dirtexthifts.keys():
                try:
                    for i in [0,1]:
                        self.defaultdirtexthifts[key][i]+=self.dirtexthifts[key][i]
                except:
                    pass
        if len(self.normtexthifts)!=0:
            for key in self.normtexthifts.keys():
                try:
                    for i in [0,1]:
                        self.defaultnormtexthifts[key][i]+=self.normtexthifts[key][i]
                except:
                    pass
        #print(self.dirtexthifts)
        for key,facecolor,edgecolor,toplot,textshiftkeys,textshifts in zip(['uvw','hkl'],[self.scatterdirsfacecolors,self.scatternormsfacecolors],[self.scatterdirsedgecolors,self.scatternormsedgecolors],[dirs,normals],[uvwkeys,hklkeys],[self.defaultdirtexthifts,self.defaultnormtexthifts]):
            #print(toplot)
            if self.printasfamily and not self.sphere=='full':
                brackL=r'\{'
                brackR=r'\}'
                if key==r'uvw':
                    brackL='\\langle'
                    brackR='\\rangle'
            else:
                brackL='('
                brackR=')'
                if key=='uvw':
                    brackL='['
                    brackR=']'
            if self.printcorrespasfamily and not self.sphere=='full':
                brackCL=r'\{'
                brackCR=r'\}'
                if key=='uvw':
                    brackCL='\\langle'
                    brackCR='\\rangle'
            else:
                brackCL='('
                brackCR=')'
                if key=='uvw':
                    brackCL='['
                    brackCR=']'

            d2plot=[(d[self.ProjType][0:2],d[key],d['textshift'],textshiftkey) for d,textshiftkey in zip(toplot,textshiftkeys) if d['vector'][2]>=0 or np.abs(d['vector'][2])<1e-5]#np.array([d['equalarea'][0:2] for d in dirs if d['vector'][2]>=0 or np.abs(d['vector'][2])<1e-5])
            #print([d2p[0][:,0] for d2p in d2plot])
            #print(facecolor)
            facecolors=[facecolor[d2p[3]] for d2p in d2plot if d2p[0][1]>=self.textlim]
            edgecolors=[edgecolor[d2p[3]] for d2p in d2plot if d2p[0][1]>=self.textlim]
            #print(facecolors)
            #print(edgecolors)
            if self.plot1phasedirs:
                self.ax.scatter([d2p[0][0,0] for d2p in d2plot if d2p[0][1]>=self.textlim],[d2p[0][1,0] for d2p in d2plot if d2p[0][1]>=self.textlim],c=facecolors,s=self.scatterpointsize,edgecolors=edgecolors,linewidths=1, zorder=500000)
            C2plotx=[]
            C2ploty=[]
            COLS=[]
            for dpi, d2p in enumerate(d2plot):   
                if d2p[0][1]>=self.textlim:
                    if d2p[0][1]==0:
                        dy=self.dy1
                        dx=self.dx1
                    else:
                        dy=self.dy2
                        dx=-self.dx2   
                    #self.ax.scatter(d2p[0][0],d2p[0][1],c='k',s=70,edgecolors='w',alpha=1,linewidths=1, zorder=5000)
                    try:
                        textPh1=self.dirnames[dpi]
                    except:
                        if self.printascubicfamily:
                            d2pf=np.sort(np.abs(d2p[1]))[::-1]
                            textPh1=f'${brackL}{{{int(d2pf[0])}}}{{{int(d2pf[1])}}}{{{int(d2pf[2])}}}{brackR}^{{{self.phase1}}}$'
                        else:
                            textPh1=f'${brackL}{{{int(d2p[1][0])}}}{{{int(d2p[1][1])}}}{{{int(d2p[1][2])}}}{brackR}^{{{self.phase1}}}$'
                    if self.printcorrespondent and key in self.printcorrespondentuvwhkl:
                        if key == 'uvw':
                            #print(self.LPhase2)
                            CdCp = self.Cd[:,:,self.varsel]
                            LP2=self.LPhase2
                        else:
                            CdCp =self.Cp[:,:,self.varsel]
                            LP2=self.LrPhase2
                        d2pc=vector2miller(CdCp.dot(d2p[1]))
                        #rvd2pc=self.T[:,:,self.varsel].dot(LP2)
                        if key=='uvw':
                            cdirs=[d2pc]
                            #print(d2p[1])
                            #print(cdirs)
                            cnorms=[]
                            C2plot,Cnorms=gen_dirs_norms(self.LPhase2, self.LrPhase2,cdirs,cnorms,R2Proj=self.R2Proj.dot(self.T[:,:,self.varsel]))
                        else:
                            cdirs=[]
                            cnorms=[d2pc]
                            Cdirs,C2plot=gen_dirs_norms(self.LPhase2, self.LrPhase2,cdirs,cnorms,R2Proj=self.R2Proj.dot(self.T[:,:,self.varsel]))
                        #print(C2plot[0][self.ProjType])
                        #print(C2plot[0])
                        
                        C2plotx.append(C2plot[0][self.ProjType][0])
                        if self.sphere=='half':
                                C2ploty.append(abs(C2plot[0][self.ProjType][1]))
                        else:
                            C2ploty.append(C2plot[0][self.ProjType][1])
                        COLS.append(stereotriangle_colors(stereoprojection_intotriangle(np.array(d2pc))))
                        #self.ax.scatter(C2plot[0][self.ProjType][0],C2plot[0][self.ProjType][1],c='w',s=25,edgecolors='k',alpha=1,linewidths=1, zorder=5000)
                        if d2pc[2]<0:
                            d2pc=-1*d2pc 
                        if self.printcorrespascubicfamily and self.printcorrespasfamily:
                            d2pc=np.sort(np.abs(d2pc))[::-1]
                            
                        textPh2=f'${brackCL}{{{int(d2pc[0])}}}{{{int(d2pc[1])}}}{{{int(d2pc[2])}}}{brackCR}^{{{self.phase2}}}$'
                        
                        if self.printPhase2First:
                            text=f'{textPh2}{self.correspdelim}{textPh1}'
                            if not self.plot1phasedirs:
                                text=f'{textPh2}'
                        else:
                            text=f'{textPh1}{self.correspdelim}{textPh2}'
                    else:
                        text=f'{textPh1}'
                    if self.printcorrespondentpoints and key in self.printcorrespondentuvwhkl:
                        if self.printcorrespondentpointscolored:
                            self.ax.scatter(C2plotx,C2ploty,c=COLS,s=self.correspondentpointsize,edgecolors='k',alpha=1,linewidths=1, zorder=500000000)
                        else:
                            self.ax.scatter(C2plotx,C2ploty,c='w',s=self.correspondentpointsize,edgecolors='k',alpha=1,linewidths=1, zorder=5000000000)
                        if self.bbox_to_anchor is None:
                            self.legend=self.ax.legend(handles=legend_elements, loc=self.legend_loc,facecolor=legbackg,ncol=self.leg_ncol)
                        else:
                            #print(self.bbox_to_anchor)
                            self.legend=self.ax.legend(bbox_to_anchor=self.bbox_to_anchor, handles=legend_elements,facecolor=legbackg,ncol=self.leg_ncol)

                    text=text.replace('{-','\\overline{')#.replace('phase',phase).replace('phs2',phase2)
                    try:
                        text2apply=self.dirsnormstext[d2p[3]]
                    except:
                        text2apply=text
                    tt=self.ax.text(d2p[0][0]+dx+textshifts[d2p[3]][0],d2p[0][1]+dy+textshifts[d2p[3]][1],text2apply,color='k', zorder=600000)
                    tt.set_bbox(dict(boxstyle='square,pad=-0.',facecolor='white', alpha=self.dirsnormsalpha, edgecolor='None'))


    def plotHist():
        """
            Plot histogram of orientation data on stereographic projection.
            
            Creates binned representation of texture with discrete color levels,
            alternative to continuous colormap visualization.
            
            Input:
                oris: array (N, 3, 3) - Orientation matrices
                bins: int - Number of histogram bins (default: 128)
                histnorm: bool - Normalize histogram (default: True)
                histscale: list [min, max] - Value range (optional)
            
            Output:
                None (plots histogram on self.ax)
            
            Usage Example:
                >>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Basic histogram
                >>> orientations = R.random(1000).as_matrix()
                >>> 
                >>> p = plotter()
                >>> p.plotProj(ProjType='equalarea', sphere='half')
                >>> p.setAttributes(oris=orientations)
                >>> p.plotHist(bins=64)
            """
        
        return
        

    def plotScatter(self,**kwargs):  
        """
            Create scatter plot of crystallographic orientations or data points.
            
            Plots data points with color and size scaling based on additional parameters.
            Useful for showing orientation distributions with properties.
            
            Input:
                scatterdata: array (N, 3, 3) - Orientation matrices to plot
                scattercolscale: array (N,) - Values for color scaling (optional)
                scattersizescale: array (N,) - Values for size scaling (optional)
                scatterfacecolor: str/array - Point colors (default: 'k')
                scatteredgecolor: str - Edge color (default: 'w')
                scatterpointsize: float - Base point size (default: 70)
            
            Output:
                None (plots scatter points on self.ax)
            
            Usage Example:
                >>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Generate random orientations
                >>> N = 500
                >>> orientations = R.random(N).as_matrix()
                >>> 
                >>> # Create scatter plot
                >>> p = plotter()
                >>> p.plotProj(ProjType='equalarea', sphere='half')
                >>> p.setAttributes(scatterdata=orientations)
                >>> p.processScatterData()
                >>> p.plotScatter()
            """
        self.setAttributes(**kwargs)
        if self.scatterdata is None:
            self.scatterplot=False
        else:
            if self.scatteroris is None:
                self.scatteroris=self.scatterdata
            if self.scatteroris is None:
                self.processScatterData()
            if self.scatterzorder is None:
                self.scatterzorder=list(range(self.scatterproj[0].shape[1]))
            if self.scattercolscale is not None:
                c=self.scattercolscale[self.scatterzorder]
            else:
                c='k'
            if self.scattersizescale is not None:
                s=self.scattersizescale[self.scatterzorder]
            else:
                s=5
            self.scat=[]
            for scatterproj in self.scatterproj:
                if self.scattercolscalevm is not None:
                    self.scat.append(self.ax.scatter(scatterproj[0,self.scatterzorder],scatterproj[1,self.scatterzorder],s=s,c=c,linewidth=self.scatterlinewidth,edgecolors=self.scatteredgecolors, cmap=self.cmap,vmin=self.scattercolscalevm[0],vmax=self.scattercolscalevm[1]))
                else:
                    self.scat.append(self.ax.scatter(scatterproj[0,self.scatterzorder],scatterproj[1,self.scatterzorder],s=s,c=c,linewidth=self.scatterlinewidth, edgecolors=self.scatteredgecolors,cmap=self.cmap))
            

    def plotScatterAsHist(self,**kwargs):
        """
            Display scatter plot data as histogram representation.
            
            Converts scattered orientation points to binned histogram display,
            useful for comparing discrete and continuous representations.
            
            Input:
                scatterdata: array (N, 3, 3) - Orientation matrices
                bins: int - Number of bins (default: 128)
                histnorm: bool - Normalize (default: True)
            
            Output:
                None (plots histogram on self.ax)
            
            Usage Example:
                >>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Scatter as histogram
                >>> orientations = R.random(500).as_matrix()
                >>> 
                >>> p = plotter()
                >>> p.plotProj(ProjType='equalarea', sphere='half')
                >>> p.setAttributes(scatterdata=orientations)
                >>> p.plotScatterAsHist(bins=64)
            """
        self.setAttributes(**kwargs)
        histdata=np.empty((3,0))
        histdataweights=np.empty((0))
        for scatterproj in self.scatterproj:
            histdata=np.hstack((histdata,scatterproj))
            histdataweights=np.hstack((histdataweights,self.hbptrstrainnorm))
        self.histdata=histdata
        self.histdataweights=histdataweights
        if self.histdataweights is None:
            hist, xedges, yedges = np.histogram2d(self.histdata[1,:], self.histdata[0,:], bins=self.bins,range=[[-1, 1], [-1, 1]])
            #print('none')
        else:
            self.histdataweights=histdataweights
            hist, xedges, yedges = np.histogram2d(self.histdata[1,:], self.histdata[0,:], bins=self.bins,range=[[-1, 1], [-1, 1]],weights=self.histdataweights)
        
        if self.histscale=='sqrt':
            hist=hist**0.5
        elif self.histscale=='log':
            hist=np.log(hist)
        if self.histnorm:
            hist=hist/np.max(hist)
        if self.histlevels is None:
            self.histlevels=np.linspace(0, np.max(hist), 5)[1:]
        X, Y = np.meshgrid((xedges[:-1] + xedges[1:])/2.,
                            (yedges[:-1] + yedges[1:])/2.)
        #print('ok')
        self.histplot=self.ax.contourf(hist, extent=(-1, 1, -1, 1),levels=self.histlevels,zorder=40000)  
        

    def processScatterData(self,**kwargs):  
        """
            Process scatter plot data by transforming orientations to projection coordinates.
            
            Converts orientation matrices to stereographic projection coordinates,
            applies crystal symmetry, and prepares data for scatter plotting.
            
            Input:
                scatterdata: array (N, 3, 3) - Orientation matrices
                scatteroris: array (N, 3, 3) - Alternative orientation data (optional)
                symops: list of arrays - Symmetry operations (default: [eye(3)])
                R2Proj: array (3, 3) - Projection rotation matrix (default: eye(3))
            
            Output:
                None (sets self.scatterproj with processed coordinates)
            
            Usage Example:
                >>> import numpy as np
                >>> from plotlib import plotter
                >>> from scipy.spatial.transform import Rotation as R
                >>> 
                >>> # Process orientation data
                >>> orientations = R.random(100).as_matrix()
                >>> 
                >>> p = plotter()
                >>> p.setAttributes(scatterdata=orientations)
                >>> p.processScatterData()
                >>> print("Projection coords:", p.scatterproj.shape)
            """
        self.setAttributes(**kwargs)            
        if self.scatterdata is not None:
            if self.scatteroris is None:
                self.scatteroris=self.scatterdata
            self.scattereqhkl=[]
            self.scatterproj=[]
            for scatterdata in self.scatterdata:
                #print(scatterdata)
                if self.ProjType=='equalarea':
                    equalarea=True
                    if self.sphere=='triangle':
                        scatterproj,scattereq=equalarea_intotriangle(scatterdata,geteqdirs=True)
                        scattereq=scattereq['eqdirs']
                    else:
                        scatterproj=equalarea_directions(scatterdata)
                        scattereq=scatterdata
                else:   
                    equalarea=False
                    if self.sphere=='triangle':
                        scatterproj,scattereq=stereoprojection_intotriangle(scatterdata,geteqdirs=True)
                        scattereq=scattereq['eqdirs']
                    else:
                        scatterproj=stereoprojection_directions(scatterdata)
                        scattereq=scatterdata
                if self.scatterLr is None:
                    self.scattereqhkl.append(vectors2miller(self.LrPhase1.dot(scattereq)))
                else:
                    self.scattereqhkl.append(vectors2miller(self.scatterLr.dot(scattereq)))
                self.scatterproj.append(scatterproj)
            

    def scatterDataAnnot(self,**kwargs):
        """Annotate scatter plot data points with information."""
        self.setAttributes(**kwargs)
        if self.scatteridxs is None:
            self.scatteridxs=list(np.where(~np.isnan(self.colmapdata))[0])     
        else:
            scatteridxs=self.scatteridxs

        #self.scatteridxs=list(np.where(~np.isnan(self.colmapdata))[0])  
        self.ANNOTS=[]
        for scatterproj in self.scattereqhkl:
            self.ANNOTS.append(self.ax.annotate("", xy=(0,0), xytext=(5,5),textcoords="offset points",
                                bbox=dict(boxstyle="round", fc="w",alpha=0.5),zorder=1000000))
            self.ANNOTS[-1].set_visible(False)
        #if self.ax2annot is None:
            #self.ax2annot=[True]*len(self.ANNOTS)
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)




# =====================================
# Standalone Utility Functions
# =====================================

#End of added functions
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
    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])
    
    

def plot_lattice_plane(axl,PlanePoints,**kwargs):

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
    BasalPlane=False
    for idx in range(3):
        if (PlanePoints[idx,:]==0).all():
            BasalPlane=True
            break
    
    if BasalPlane:
        idxs=list(range(3))
        idxs.remove(idx)
        hull = ConvexHull(PlanePoints[idxs,:].T)
    
        axl.add_collection3d(Poly3DCollection([PlanePoints[:,hull.vertices].T],**kwargs))
    else:
        axl.plot_trisurf(PlanePoints[0,:],PlanePoints[1,:], PlanePoints[2,:],**kwargs)


def plot_lattice_boundaries(axl,LatticePointsNew,allPoints=None,polygon=False,tol=1e-1,**kwargs):

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
    if allPoints is None:
        allPoints=np.hstack([p for points in LatticePointsNew for p in points])
    if polygon:
        for idx in range(3):
            Xbound=copy.deepcopy(allPoints)
            for extrval in [np.min(Xbound[idx,:]),np.max(Xbound[idx,:])]:
                Xbound=copy.deepcopy(allPoints)
                #Xbound=Xbound[:,Xbound[idx,:]==extrval]
                Xbound=Xbound[:,np.abs(Xbound[idx,:]-extrval)<tol]
                #np.abs(vertices[0,:]-np.min(vertices[0,:]))<1e-1
                #Xbound[idx,:]=Xbound[idx,:]*0+extrval
                #axl.plot_trisurf(Xbound[0,:],Xbound[1,:], Xbound[2,:],\
                #                 alpha=0.5,color='r', linewidths=0., edgecolors='grey',linestyle='-',\
                #                 linewidth = 0.0, antialiased = True) 
                if Xbound.shape[1]>=3:
                    try:
                        hull = ConvexHull(Xbound[np.delete(range(3),idx,0),:].T)
                        axl.add_collection3d(Poly3DCollection([Xbound[:,hull.vertices].T], **kwargs))
                    except:
                        pass
    else:
        axl.plot_trisurf(allPoints[0,:],allPoints[1,:], allPoints[2,:],triangles=ConvexHull(allPoints.T).simplices,
                         **kwargs)

            


def plot_lattice3D(ax,VV,description,Parentlattice_points,Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2):

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
    xlim= np.array([-1.05,1.05])*np.sqrt(Product_uvw2xyz[:,0].dot(Product_uvw2xyz[:,0]))
    for point in Parentlattice_points:
        ax.plot(point[0],point[1],point[2],'r')
                
    for point in Productlattice_points:
        ax.plot(point[0],point[1],point[2],'b')
        
        
    Product_basal = np.matmul(Product_uvw_2_Parent_uvw_all_norm,np.matmul(Product_uvw2xyz,np.eye(3)))
    colors=['g','c','#800000']
        
    inc=-1;
    for v2 in Product_basal.T:
        inc+=1
        ax.plot([0,v2[0]],[0,v2[1]],[0,v2[2]],color=colors[inc],linewidth=linewidth)
    inc=0

    for v in VV.T:
        inc+=1
        #print(v)
        v2=1.5*np.sqrt(Product_uvw2xyz[:,0].dot(Product_uvw2xyz[:,0]))*v
        ax.plot([0,v2[0]],[0,v2[1]],[0,v2[2]],'k',linewidth=linewidth,linestyle='--')
        ax.text(v2[0],v2[1],v2[2],description.replace('{inc}','{'+str(inc)+'}'))       
        
    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
# equal aspect ratio
    ax.set_xlim3d(xlim)
    ax.set_ylim3d(xlim)
    ax.set_zlim3d(xlim)
    ax.set_xlim(xlim)
    ax.set_ylim(xlim)
    ax.set_zlim(xlim)
        
  

    set_aspect_equal_3d(ax)


def plot_latticefaces3D(ax,Parentlattices,linewidth=2,alpha=0.15,edgecolor='r',linestyle='-',facecolor=(1, 0, 0, 0.15)):
    """
    plot_latticefaces3D - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """
    
    for Lattice in Parentlattices:
        ax.add_collection3d(Poly3DCollection(Lattice, alpha=alpha,facecolors=facecolor, linewidths=linewidth, edgecolors=edgecolor,linestyle=linestyle))
        
    maxlimits=[]
    minlimits=[]
        
    for Lattice in Parentlattices:
        for face in Lattice:
            for point in face:
                if len(maxlimits)==0:
                    maxlimits=[point[0],point[1],point[2]]
                else:
                    for ii in range(0,3):
                        if maxlimits[ii]<point[ii]:
                            maxlimits[ii]=point[ii]
                if len(minlimits)==0:
                    minlimits=[point[0],point[1],point[2]]
                else:
                    for ii in range(0,3):
                        if minlimits[ii]>point[ii]:
                            minlimits[ii]=point[ii]
                    
                
    #ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
# equal aspect ratio
    ax.set_xlim3d([minlimits[0],maxlimits[0]])
    ax.set_ylim3d([minlimits[1],maxlimits[1]])
    ax.set_zlim3d([minlimits[2],maxlimits[2]])
    set_aspect_equal_3d(ax)


def plot_latticesfaces3D(ax,VV,description,Parentlattices,Productlattices,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2,alpha=0.15,xlim=[-2,2]):

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
    if xlim is None:
        xlim= np.array([-1.05,1.05])*max([np.sqrt(V.dot(V)) for V in Product_uvw2xyz.T])

    for Lattice in Parentlattices:
        ax.add_collection3d(Poly3DCollection(Lattice, alpha=alpha,facecolors=(1, 0, 0, alpha), linewidths=linewidth, edgecolors='r'))
    for Lattice in Productlattices:
        ax.add_collection3d(Poly3DCollection(Lattice, alpha=alpha,facecolors=(0, 0, 1, alpha), linewidths=linewidth, edgecolors='b'))
        
        
    Product_basal = np.matmul(Product_uvw_2_Parent_uvw_all_norm,np.matmul(Product_uvw2xyz,np.eye(3)))
    colors=['g','c','#800000']
        
    inc=-1;
    for v2 in Product_basal.T:
        inc+=1
        ax.plot([0,v2[0]],[0,v2[1]],[0,v2[2]],color=colors[inc],linewidth=linewidth)
    inc=0

    for v in VV.T:
        inc+=1
        #print(v)
        v2=1.5*np.sqrt(Product_uvw2xyz[:,0].dot(Product_uvw2xyz[:,0]))*v
        ax.plot([0,v2[0]],[0,v2[1]],[0,v2[2]],'k',linewidth=linewidth,linestyle='--')
        ax.text(v2[0],v2[1],v2[2],description.replace('{inc}','{'+str(inc)+'}'))       
        
    #ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
# equal aspect ratio
    ax.axis('auto')
    ax.set_xlim3d(xlim)
    ax.set_ylim3d(xlim)
    ax.set_zlim3d(xlim)
    ax.set_xlim(xlim)
    ax.set_ylim(xlim)
    ax.set_zlim(xlim)
        
  

    set_aspect_equal_3d(ax)





def plot_lattice2D(ax,VV,description,Parentlattice_points,Parent_lattice,\
                   Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2,xlim=None):

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
    if xlim is None:
        xlim= np.array([-1.05,1.05])*max([np.sqrt(V.dot(V)) for V in Product_uvw2xyz.T])
    shiftx=xlim[1];
    shifty=xlim[1]*0;
    facx=1.8
    facy=1.3
    shifts=[[facx*xlim[1],-facy*xlim[1]],[-facx*xlim[1],-facy*xlim[1]],[-facx*xlim[1],facy*xlim[1]]]
    pairs=[[0,1],[2,1],[0,2]]
    signs = [[1,1],[-1,1],[-1,1]]
    coordlength=xlim[1]*0.5
    coordabasalvecs = [[[1,0,0],[0,1,0]],[[0,0,1],[0,1,0]],[[0,0,1],[1,0,0]]]
    colors=['g','c','#800000']
    for shift,pair,sgn,vecs in zip(shifts,pairs,signs,coordabasalvecs):
#        for point,vecs in zip([[xlim[1]*0.5,0],[0,xlim[1]*0.5]],):
            #point=np.array(point);
        ax.plot(np.array([0,sgn[0]*coordlength])+2*shift[0],np.array([0,0])+2*shift[1],'k')
        ax.plot(np.array([0,0])+2*shift[0],np.array([0,sgn[1]*coordlength])+2*shift[1],'k')
#        ax.plot(np.array([0,sgn[0]*coordpoint[0]])+2*shift[0],np.array([0,sgn[1]*coordpoint[1]])+2*shift[1],'k')
        ax.text(2*shift[0],sgn[1]*coordlength+2*shift[1], dir2string(vecs[1], digits=0)+r'$_{'+Parent_lattice+'}$',fontsize=10)
        addshiftx=0.0
        addshifty=0.0
        if sgn[0]<0:
            addshiftx=-coordlength*1;
            addshifty=coordlength*0.1;
        ax.text(sgn[0]*coordlength+2*shift[0]+addshiftx,addshifty+2*shift[1], dir2string(vecs[0], digits=0)+r'$_{'+Parent_lattice+'}$',fontsize=10)
            
        for point in Parentlattice_points:
            point=np.array(point);
            ax.plot(sgn[0]*point[pair[0]]+shift[0],sgn[1]*point[pair[1]]+shift[1],'r')
                    
        for point in Productlattice_points:
            point=np.array(point);
            ax.plot(sgn[0]*point[pair[0]]+shift[0],sgn[1]*point[pair[1]]+shift[1],'b')
    
        inc=-1;
        Product_basal = np.matmul(Product_uvw_2_Parent_uvw_all_norm,np.matmul(Product_uvw2xyz,np.eye(3)))
        for v2 in Product_basal.T:
            inc+=1
            ax.plot(sgn[0]*np.array([0,v2[pair[0]]])+shift[0],sgn[1]*np.array([0,v2[pair[1]]])+shift[1],color=colors[inc],linewidth=linewidth)
    
        inc=0
        for v in VV.T:
            inc+=1
            v2=1.8*np.sqrt(Product_uvw2xyz[:,0].dot(Product_uvw2xyz[:,0]))*v
            ax.plot(sgn[0]*np.array([0,v2[pair[0]]])+shift[0],sgn[1]*np.array([0,v2[pair[1]]])+shift[1],'k',linewidth=linewidth,linestyle='--')
            ax.text(sgn[0]*v2[pair[0]]*1.1+shift[0],sgn[1]*v2[pair[1]]*1.1+shift[1],description.replace('{inc}','{'+str(inc)+'}'))       

    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
# equal aspect ratio
    ax.set_xlim([-facx*2*xlim[1],facx*2*xlim[1]])
    ax.set_ylim([-facx*2*xlim[1],facx*2*xlim[1]])


def plot_lattice_2Dprojection(ax,VV,description,Parentlattice_points,Parent_lattice,\
                   Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,normals,verticals, linewidth=2,xlim=None):

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
    if xlim is None:
        xlim= np.array([-1.05,1.05])*max([np.sqrt(V.dot(V)) for V in Product_uvw2xyz.T])
    shiftx=xlim[1];
    shifty=xlim[1]*0;
    facx=1.8
    facy=1.3
    shifts=[[facx*xlim[1],-facy*xlim[1]],[-1.5*facx*xlim[1],-facy*xlim[1]],[-1.5*facx*xlim[1],1.4*facy*xlim[1]]]
    pairs=[[0,1],[2,1],[0,2]]
    signs = [[1,1],[-1,1],[-1,1]]
#    normals = [[0.,0.,1.],[1.,0.,0.],[0.,-1.,0.]]
#    verticals = [[0.,1.,0.],[0.,1.,0.],[1.,0.,0.]]
#    normals = [[0.,0.,1.],[1.,0.,0.],[-1.,1.,0.]]
#    verticals = [[0.,1.,0.],[0.,1.,0.],[0.,0.,1.]]
    coordlength=xlim[1]*0.5
    coordabasalvecs = [[[1,0,0],[0,1,0]],[[0,0,1],[0,1,0]],[[0,0,1],[1,0,0]]]
    colors=['g','c','#800000']
    
    for shift,normal,vertical,vecs in zip(shifts,normals,verticals,coordabasalvecs):
        normal=np.array(normal);
        vertical=np.array(vertical);
        horizontal=np.cross(vertical,normal)
        
        ax.plot(np.array([0,coordlength])+2*shift[0],np.array([0,0])+2*shift[1],'k')
        ax.plot(np.array([0,0])+2*shift[0],np.array([0,coordlength])+2*shift[1],'k')
        ax.text(2*shift[0],coordlength+2*shift[1], dir2string(vertical, digits=1)+r'$_{'+Parent_lattice+'}$',fontsize=10)
        #print(vertical)
        addshiftx=0.0
        addshifty=0.0
#        if sgn[0]<0:
#            addshiftx=-coordlength*1;
#            addshifty=coordlength*0.1;
        ax.text(coordlength+2*shift[0]+addshiftx,addshifty+2*shift[1], dir2string(horizontal, digits=1)+r'$_{'+Parent_lattice+'}$',fontsize=10)
        ax.text(-coordlength+2*shift[0]+addshiftx,-0.5*coordlength+addshifty+2*shift[1], plane2string(normal, digits=1)+r'$_{'+Parent_lattice+'}$',fontsize=10)
#        print(normal)
#        print(plane2string(normal, digits=0))
        for point in Parentlattice_points:
            point=np.array(point);
            point_proj_x = horizontal.dot(point)
            point_proj_y = vertical.dot(point)
            ax.plot(point_proj_x+shift[0],point_proj_y+shift[1],'r')
                    
        for point in Productlattice_points:
            point=np.array(point);
            point_proj_x = horizontal.dot(point)
            point_proj_y = vertical.dot(point)
            ax.plot(point_proj_x+shift[0],point_proj_y+shift[1],'b')
            
        inc=-1;
        Product_basal = np.matmul(Product_uvw_2_Parent_uvw_all_norm,np.matmul(Product_uvw2xyz,np.eye(3)))
        for v2 in Product_basal.T:
            inc+=1
            point_proj_x = horizontal.dot(v2)
            point_proj_y = vertical.dot(v2)

            ax.plot(np.array([0,point_proj_x])+shift[0],np.array([0,point_proj_y])+shift[1],color=colors[inc],linewidth=linewidth)
    
        inc=0
        for v in VV.T:
            inc+=1
            v2=1.8*np.sqrt(Product_uvw2xyz[:,0].dot(Product_uvw2xyz[:,0]))*v
            point_proj_x = horizontal.dot(v2)
            point_proj_y = vertical.dot(v2)
            
            ax.plot(np.array([0,point_proj_x])+shift[0],np.array([0,point_proj_y])+shift[1],'k',linewidth=linewidth,linestyle='--')
            ax.text(point_proj_x*1.1+shift[0],point_proj_y*1.1+shift[1],description.replace('{inc}','{'+str(inc)+'}'))       

    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
# equal aspect ratio
    ax.set_xlim([-facx*2*xlim[1],facx*2*xlim[1]])
    ax.set_ylim([-facx*2*xlim[1],facx*2*xlim[1]])        
    

def plot_mohr_circles(mcircles,VV,DD,xyz2uvw,scale,xticks=None,yticks=None,ax=None,Parent_lattice='B2'):

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
    Return=True
    if ax==None:
        fig, ax = plt.subplots()
    else:
        fig=[]
        Return=False
#    fig, ax = plt.subplots()
#    ax.tick_params(
#        axis='both',
#        which='both',
#        bottom=False,
#        top=False,
#        left=False,
#        labelbottom=False,
#        labelleft=False)    
    phi=np.linspace(0,2*np.pi,1000)
    C13x = mcircles['C13']+mcircles['R13']*np.cos(phi)
    C13y = mcircles['R13']*np.sin(phi)
    C23x = mcircles['C23']+mcircles['R23']*np.cos(phi)
    C23y = mcircles['R23']*np.sin(phi)
    C12x = mcircles['C12']+mcircles['R12']*np.cos(phi)
    C12y = mcircles['R12']*np.sin(phi)
    
    ax.plot(C13x*scale,C13y*scale,'r')
    ax.plot(C12x*scale,C12y*scale,'g')
    ax.plot(C23x*scale,C23y*scale,'b')
    #spine placement data centered
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    
    ax.spines['left'].set_position(('data', 0.0))
    ax.spines['bottom'].set_position(('data', 0.0))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    xlim = [np.round(C13x.min()*scale*1.5),np.round(C13x.max()*scale*1.5)]
    ylim = [np.round(C13y.max()*scale*1.5),np.round(C13y.min()*scale*1.5)]
    
    if xticks is None:
        xlim2 = [np.round(C13x.min()*scale*1.5),np.round(C13x.max()*scale*1.5)]
        ax.set_xticks(np.round(np.linspace(xlim[0],xlim[1],10)))
    else:
        xlim2 = [xticks[0],xticks[-1]]
        ax.set_xticks(xticks)
    if yticks is None:
        ylim2 = [np.round(C13y.max()*scale*1.5),np.round(C13y.min()*scale*1.5)]
        ax.set_yticks(np.round(np.linspace(ylim[0],ylim[1],10)))
    else:
        ylim2 = [yticks[0],yticks[-1]]
        ax.set_yticks(yticks)
    ax.set_xlim(xlim2)
    ax.set_ylim(ylim2)
    
    #ax.set_yticks(np.linspace(9,-9,10))
    #ax.set_ylim([9,-10])
    ax.text(xlim[1]*0.65,max(ylim)/10*2,'Normal\nstrain, '+r'$\varepsilon$ [%]')
    ax.text(-max(ylim)/10*2,ylim[1]*0.75, 'Shear\nstrain, '+r'$\gamma$/2 [%]')
#    ax.set_xlim([-10,14])
#    ax.set_xticks(np.linspace(-10,12,12))
#    ax.set_ylim([9,-10])
#    ax.set_yticks(np.linspace(9,-9,10))
#    ax.set_ylim([9,-10])
#    ax.text(11,2,'Normal\nstrain, '+r'$\varepsilon$ [%]')
#    ax.text(-2,-9.5, 'Shear\nstrain, '+r'$\gamma$/2 [%]')
    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
  # equal aspect ratio
    #plt.show()
    
    #plot principal direction in mohr circles
    for i in [0,1,2]:
        ax.plot(DD[i]*scale,0,'ko')
        ax.text(DD[i]*scale+max(ylim)/100*2,-max(ylim)/100*2,r'$\varepsilon_'+str(i+1)+'$')
    textprincipaldirs=''
    for i in [0,1,2]:
        vd = xyz2fractional(xyz2uvw,VV[:,i])
        textprincipaldirs+=str(r'$\varepsilon_'+str(i+1)+'$~'+vec2string(vd)+'$_{'+Parent_lattice+'}$\n')
    ax.text(max(xlim)/2.2,max(ylim)/1.,textprincipaldirs)
    
    if Return:
        return fig,ax




def plot_planes_on_mohr_circle(ax,scale,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl, colors,text=False,Parent_lattice='B2'):

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
    #plot selected planes on 13 mohr cicle and triangle and wulffnet
    inc=-1;
    Upperhalftext = []
    Lowerhalftext=[]
    
    for ii in range(0,len(SelNormalsOnCircle['strainmag'])):

        if type(colors) is list:
            if len(colors)>=len(SelNormalsOnCircle['strainmag']):
                inc+=1
            else:
                colors=[colors[0]]
                inc=0
        else:
            colors=[colors]
            inc=0
            
        xx=[]
        yy=[]
        uvwtext=[]
        for i in [0,1]:
            xx.append(SelNormalsOnCircle['strainmag'][ii][i]*scale)
            yy.append(SelShearsOnCircle['strainmag'][ii][i]*scale)
            vd = xyz2fractional(Parent_xyz2hkl,SelNormalsOnCircle['inlatticedir'][ii][i])
            proj_dir=stereoprojection_intotriangle(SelNormalsOnCircle['inlatticedir'][ii][i])
            proj_dir2=stereoprojection_directions(SelNormalsOnCircle['inlatticedir'][ii][i])
    
            if yy[-1]>0:
                lowerhalftext=r'$\mathbf{\varepsilon}}$='+str(round((xx[-1]*10))/10)+','+r'$\mathbf{\gamma}$/2='+str(round((yy[-1]*10))/10)+',n='+r''+plane2string(vd)+r'$_{\mathbf{'+Parent_lattice+'}}$'
                Lowerhalftext.append(lowerhalftext)
                lowerhalftext=r'$\varepsilon$='+str(round((xx[-1]*10))/10)+r',$\gamma$/2='+str(round((yy[-1]*10))/10)+' ,n='+plane2string(vd)+r'$_{'+Parent_lattice+'}$'
                markerfacecolor=colors[inc]
                markersize=8
            else:
                markerfacecolor='None'
                markersize=12
                upperhalftext =r'$\mathbf{\varepsilon}$='+str(round((xx[-1]*10))/10)+','+r'$\mathbf{\gamma}$/2='+str(round((yy[-1]*10))/10)+',n='+r''+plane2string(vd)+r'$_{\mathbf{'+Parent_lattice+'}}$'
                Upperhalftext.append(upperhalftext)
                upperhalftext=r'$\varepsilon$='+str(round((xx[-1]*10))/10)+r',$\gamma$/2='+str(round((yy[-1]*10))/10)+r' ,n='+plane2string(vd)+'$_{'+Parent_lattice+'}$'
            ax.plot(xx[-1],yy[-1],'o',markerfacecolor=markerfacecolor,markeredgecolor=colors[inc])
        ax.plot(xx,yy,color=colors[inc])
        if text:
            idx = yy.index(max(yy))
            idx2 = yy.index(min(yy))
            ax.text(xx[idx],yy[idx]*1.1,lowerhalftext)
            ax.text(xx[idx2],yy[idx2]*1.05,upperhalftext)

    return Upperhalftext,Lowerhalftext


def plot_planes_on_stereotriangle(ax,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl,colors):

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
    inc=-1;
    for ii in range(0,len(SelNormalsOnCircle['strainmag'])):
        if type(colors) is list:
            if len(colors)>=len(SelNormalsOnCircle['strainmag']):
                inc+=1
            else:
                colors=[colors[0]]
                inc=0
        else:
            colors=[colors]
            inc=0
            
        xx=[]
        yy=[]
        for i in [0,1]:
            xx.append(SelNormalsOnCircle['strainmag'][ii][i])
            yy.append(SelShearsOnCircle['strainmag'][ii][i])
            proj_dir=stereoprojection_intotriangle(SelNormalsOnCircle['inlatticedir'][ii][i])
    
            if yy[-1]>0:
                markerfacecolor=colors[inc]
                markersize=8
            else:
                markerfacecolor='None'
                markersize=12
            ax.plot(proj_dir[0,:], proj_dir[1,:],'o',markerfacecolor=markerfacecolor,markeredgecolor=colors[inc],\
                     markeredgewidth=2,markersize=markersize)
#            ax4.plot(proj_dir2[0,:], proj_dir2[1,:],'o',markerfacecolor=markerfacecolor,markeredgecolor=colors[inc],\
#                     markeredgewidth=2,markersize=markersize)



def plot_planes_on_wulffnet(ax,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl,colors):

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
    inc=-1;
    for ii in range(0,len(SelNormalsOnCircle['strainmag'])):
        if type(colors) is list:
            if len(colors)>=len(SelNormalsOnCircle['strainmag']):
                inc+=1
            else:
                colors=[colors[0]]
                inc=0
        else:
            colors=[colors]
            inc=0
            
        xx=[]
        yy=[]
        for i in [0,1]:
            xx.append(SelNormalsOnCircle['strainmag'][ii][i])
            yy.append(SelShearsOnCircle['strainmag'][ii][i])
            proj_dir2=stereoprojection_directions(SelNormalsOnCircle['inlatticedir'][ii][i])
    
            if yy[-1]>0:
                markerfacecolor=colors[inc]
                markersize=8
            else:
                markerfacecolor='None'
                markersize=12
            ax.plot(proj_dir2[0,:], proj_dir2[1,:],'o',markerfacecolor=markerfacecolor,markeredgecolor=colors[inc],\
                     markeredgewidth=2,markersize=markersize)


def plot_princip_dir_on_stereotriangle(ax,VV,description,markersize=10,markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5):

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
    proj_dirs = stereoprojection_intotriangle(VV)
    
    for proj_dir,i in zip(proj_dirs.T,[0,1,2]):
        ax.plot(proj_dir[0], proj_dir[1], 'o',markerfacecolor=markerfacecolor,\
                 markeredgecolor=markeredgecolor,markeredgewidth=markeredgewidth,markersize=markersize)
        text=str(description.replace('{inc}','{'+str(i+1)+'}'))
        ax.text(proj_dir[0]-0.03, proj_dir[1]+0.01,text)


def plot_princip_dir_on_wulffnet(ax,VV,description,markersize=10,markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5):

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
    proj_dirs = stereoprojection_directions(VV)
    
    for proj_dir,i in zip(proj_dirs.T,[0,1,2]):
        ax.plot(proj_dir[0], proj_dir[1], 'o',markerfacecolor=markerfacecolor,\
                 markeredgecolor=markeredgecolor,markeredgewidth=markeredgewidth,markersize=markersize)
        text=str(description.replace('{inc}','{'+str(i+1)+'}'))
#        ax2.text(proj_dir[0]-0.03, proj_dir[1]+0.01,text)
        ax.text(proj_dir[0]-0.1, proj_dir[1]+0.03,text)
        

def plot_lattice(Points,LatticeVectors,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],Q=np.eye(3),shift=np.zeros(3),atoms=True,linewidth=1,move=np.zeros(3),normal=np.array([0,0,0]),halfspace='upper',s=200,plot=True):
    """
    plot_lattice - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """
    #colors=['r','b','g']
    #normal=np.array([1,1,1.5])
    
    halfscp=False
    
    pplot=True
    if not (normal==np.array([0,0,0])).all():
        halfscp=True
        normal=normal/np.sqrt(normal.dot(normal))
    if ax==None and plot:
        fig = plt.figure() 
        ax = fig.add_subplot(111, projection='3d', proj_type = 'ortho') 
    elif not plot:
        ax=[]
        fig=[]
    else:
        fig=[]
    eps=1e-1
    #S = 1
    if atoms:
        PointsNew=[]
        for Points1,col,ecol in zip(Points,colors,edgecolors):
            PointNew=[]
            for points in Points1:
                pointnew=[]
                points=np.matmul(np.eye(3),points)+np.repeat(np.array([shift]).T,points.shape[1],1)
                point_proj_n = normal.dot(points)
                if halfspace=='lower':
                    idxs = np.where((point_proj_n)<=(eps))[0]
                else:
                    idxs = np.where((point_proj_n)>(-eps))[0]
                #print(len(idxs))    
                if len(idxs)>0:
                    points=Q.dot(points)
                    PointNew.append(points[:,idxs])
                    if plot:
                        ax.scatter(points[0,idxs]+move[0], points[1,idxs]+move[1], points[2,idxs]+move[2], s = s,color=col,edgecolors=ecol,alpha=salpha,linewidths=2) 
                #plt.show()
            PointsNew.append(PointNew)
        #PointsNew.append(PointNew)
    LatticeVectorsNew=[]
    for LatticeVectors1 in LatticeVectors:     
        LatticeVectorNew=[]
        for vec in LatticeVectors1:
            pplot=True

            point1=vec[:,0]
            point2=vec[:,1]
            #print(vec)

            if halfscp:
                point1_proj_n = normal.dot(point1)
                point2_proj_n = normal.dot(point2)
                point_proj_n = normal.dot(vec)
                #print(point_proj_n)
                if halfspace=='lower':
                    idxs = np.where((point_proj_n)<=(eps))[0]
                else:
                    idxs = np.where((point_proj_n)>(-eps))[0]
                #print(len(idxs))    
                if len(idxs)==0:
                    pplot=False
                elif len(idxs)==1:
                    inp=plane_line_intersection(normal,shift,vec[:,0],vec[:,1])
                    vec=np.vstack((vec[:,idxs[0]],inp)).T
            
            if pplot:               
                vec=np.matmul(Q,vec)+np.repeat(np.array([shift]).T,vec.shape[1],1)+np.repeat(np.array([move]).T,vec.shape[1],1)
                LatticeVectorNew.append(vec)
#                ax.plot(vec[0,:]+move[0], vec[1,:]+move[1], vec[2,:]+move[2],color=gridcolor,alpha=lalpha,linewidth=linewidth)
                if plot:
                    ax.plot(vec[0,:], vec[1,:], vec[2,:],color=gridcolor,alpha=lalpha,linewidth=linewidth)
        LatticeVectorsNew.append(LatticeVectorNew)

def plane_line_intersection(n,V0,P0,P1):

    """
    Calculate intersection of plane and line in 3D.
    
    Finds point where line intersects plane, if it exists.
    
    Input:
        plane_point (array [3]): Point on plane
        plane_normal (array [3]): Plane normal vector
        line_point (array [3]): Point on line
        line_direction (array [3]): Line direction vector
    
    Output:
        numpy.ndarray [3] or None: Intersection point, or None if parallel
    
    Usage Example:
        >>> import numpy as np
        >>> 
        >>> # (001) plane at z=5
        >>> plane_pt = np.array([0, 0, 5])
        >>> plane_n = np.array([0, 0, 1])
        >>> 
        >>> # Line from origin in [111] direction
        >>> line_pt = np.array([0, 0, 0])
        >>> line_dir = np.array([1, 1, 1])
        >>> 
        >>> intersection = plane_line_intersection(plane_pt, plane_n, line_pt, line_dir)
        >>> print(f"Intersection point: {intersection}")
        >>> # Should be [5, 5, 5]
    
    Notes:
        - Returns None if line parallel to plane
        - Returns None if line in plane
        - Used in geometric calculations
    
    Formula:
        t = n·(p₀ - l₀) / (n·d)
        intersection = l₀ + t·d
    """
    # n: normal vector of the Plane 
    # V0: any point that belongs to the Plane 
    # P0: end point 1 of the segment P0P1
    # P1:  end point 2 of the segment P0P1
            
#    normal=normal/np.sqrt(normal.dot(normal))
#    PointsOut=[]    
#    for Points1 in LatticePoints:
#        PointsOut1=[]
#        for points in Points1:   
#            #point_proj_n = normal.dot(points/np.sqrt(np.sum(points**2,axis=0)))
#            #print(points[:,0]+shift)
#            
#            point_proj_n = normal.dot(points)
#            if side=='bottom':
#                idxs = np.where((point_proj_n-shift)<(eps))[0]
#            else:
#                idxs = np.where((point_proj_n-shift)>(-eps))[0]
#            #print(point_proj_n)
#            PointsOut1.append(points[:,idxs])
#        PointsOut.append(PointsOut1)
    #    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
    if plot:  
        ax.axis('auto')      
    
    if halfscp:
        if atoms:
            return fig,ax,LatticeVectorsNew,PointsNew
        else:
            return fig,ax,LatticeVectorsNew
    else:
        return fig,ax

def plot_lattice_proj(LatticeVectors,normalproj,verticalproj, ax=None, linewidth=2,color='b',eps=1e-1,Q=np.eye(3),Qprojr=np.eye(2),shift=np.zeros(3),shiftproj=0,move=np.zeros(3),normal=np.array([0,0,0]),shifthalfspace=np.zeros(3),
                      halfspace='upper',shiftplot=np.array([0,0]),out=False):
    """
    plot_lattice_proj - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """
    if not isinstance(normalproj, np.ndarray):
        normalproj=np.array(normalproj);
    normalproj=normalproj/np.sqrt(normalproj.dot(normalproj))
    if not isinstance(verticalproj, np.ndarray):
        verticalproj=np.array(verticalproj);
    verticalproj=verticalproj/np.sqrt(verticalproj.dot(verticalproj))
    halfscp=True
    pplot=True
    if not isinstance(normal, np.ndarray):
        normal=np.array(normal);
    if not (normal==np.array([0,0,0])).all():
        halfscp=True
        normal=normal/np.sqrt(normal.dot(normal))

    if ax==None:
        fig = plt.figure() 
        ax = fig.add_subplot(111) 
    else:
        fig=[]

    
    
    horizontalproj=np.cross(verticalproj,normalproj)
    pointOut=[]
    pointOutproj=[]
    for LatticeVectors1 in LatticeVectors:  
        pout=[]
        poutp=[]
        for vec in LatticeVectors1:
            #print(vec)
            pplot=True
            point_proj_n=normalproj.dot(vec)
            #print(abs((point_proj_n-shiftproj)))
            idxs = np.where(abs((point_proj_n-shiftproj))<=(eps))[0]
            #print(abs((point_proj_n-shiftproj)))
            #print(len(idxs))
            if len(idxs)==0:
                pplot=False
            #elif len(idxs)==1:
            #    inp=plane_line_intersection(normal,shifthalfspace,vec[:,0],vec[:,1])
            #    vec=np.vstack((vec[:,idxs[0]],inp)).T
            #    print(vec)

            if halfscp:
                point_proj_n = normal.dot(vec)
                #print(point_proj_n)
                if halfspace=='lower':
                    idxs = np.where(((point_proj_n))<=(eps))[0]
                else:
                    idxs = np.where((point_proj_n)>(-eps))[0]
                    #idxs = np.where(abs((point_proj_n))<=(eps))[0]
                #print(len(idxs))    
                if len(idxs)==0:
                    pplot=False
                elif len(idxs)==1:
                    inp=plane_line_intersection(normal,shifthalfspace,vec[:,0],vec[:,1])
                    vec=np.vstack((vec[:,idxs[0]],inp)).T
            
            if pplot:
                #print(vec)
                point = np.matmul(Q,vec)+np.repeat(np.array([shift]).T,vec.shape[1],1)
                point[0,:]=point[0,:]+move[0]
                point[1,:]=point[1,:]+move[1]
                point[2,:]=point[2,:]+move[2]
                
                point=np.array(point);
                point_proj_x = horizontalproj.dot(point)
                point_proj_y = verticalproj.dot(point)
                pv = Qprojr.dot(np.vstack((point_proj_x,point_proj_y)))
                #ax.plot(point_proj_x,point_proj_y,color=color,linewidth=linewidth)
                ax.plot(pv[0,:]+shiftplot[0],pv[1,:]+shiftplot[1],color=color,linewidth=linewidth)
                pout.append(point)
                poutp.append([pv[0,:]+shiftplot[0],pv[1,:]+shiftplot[1]])
        pointOut.append(pout)
        pointOutproj.append(poutp)
                
    ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
 
    #plt.show()       
    if out:
        return fig,ax,pointOut,pointOutproj,horizontalproj,verticalproj
    else:                 
        return fig,ax


def plot_points_proj(Points,normalproj,verticalproj, ax=None, marker="o",markersize=10, color='b',Q=np.eye(3),Qprojr=np.eye(2),shift=np.zeros(3),move=np.zeros(3)):
    """
    plot_points_proj - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """
    if not isinstance(normalproj, np.ndarray):
        normalproj=np.array(normalproj);
    normalproj=normalproj/np.sqrt(normalproj.dot(normalproj))
    if not isinstance(verticalproj, np.ndarray):
        verticalproj=np.array(verticalproj);
    verticalproj=verticalproj/np.sqrt(verticalproj.dot(verticalproj))

    if ax==None:
        fig = plt.figure() 
        ax = fig.add_subplot(111) 
    else:
        fig=[]

    eps=1e-1
    
    horizontalproj=np.cross(verticalproj,normalproj)
    points_proj=[]
    for point in Points:   
        #print(point)
        point = Q.dot(point)+shift + move
        
        point_proj_x = horizontalproj.dot(point)
        point_proj_y = verticalproj.dot(point)
        point_proj_n= normalproj.dot(point)
        pv = Qprojr.dot([point_proj_x,point_proj_y])
        points_proj.append([point_proj_x,point_proj_y,point_proj_n])
        #ax.plot(point_proj_x,point_proj_y,color=color,linewidth=linewidth)
        ax.plot(pv[0],pv[1],color=color,marker=marker,linestyle='',markersize=markersize)
                
    #ax.set_aspect('equal',adjustable='box')  # equal aspect ratio
 
    #plt.show()                        
    return fig,ax,points_proj


def plot_atomic_plane2D(LatticePoints,normal,vertical,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],plot=True,salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],linewidths=[1,1,1],markersizes=[200,200,200],
                      Q=np.eye(3),xlim=[],ylim=[],out=False,zorder=1):
    """
    plot_atomic_plane2D - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """
    if not isinstance(normal, np.ndarray):
        normal=np.array(normal);
    if not isinstance(vertical, np.ndarray):
        vertical=np.array(vertical);
    normal=normal/np.sqrt(normal.dot(normal))
    vertical=vertical/np.sqrt(vertical.dot(vertical))
    horizontal=np.cross(vertical,normal)
    if ax==None and plot:
        fig = plt.figure() 
        ax = fig.add_subplot(111) 
    else:
        fig=[]
    Pointsout=[]
    for Points1,col,ecol,linewidth,markersize in zip(LatticePoints,colors,edgecolors,linewidths,markersizes):
        pointout=[]
        for points in Points1:
            points2=Q.dot(points)
            point_proj_x = horizontal.dot(points2)
            point_proj_y = vertical.dot(points2)
            if len(xlim)>0:
                idxs = np.where(point_proj_x<xlim[0])[0]
                point_proj_x=np.delete(point_proj_x,idxs)
                point_proj_y=np.delete(point_proj_y,idxs)
                idxs = np.where(point_proj_x>xlim[1])[0]
                point_proj_x=np.delete(point_proj_x,idxs)
                point_proj_y=np.delete(point_proj_y,idxs)               
            if len(ylim)>0:
                idxs = np.where(point_proj_y<ylim[0])[0]
                point_proj_x=np.delete(point_proj_x,idxs)
                point_proj_y=np.delete(point_proj_y,idxs)
                idxs = np.where(point_proj_y>ylim[1])[0]
                point_proj_x=np.delete(point_proj_x,idxs)
                point_proj_y=np.delete(point_proj_y,idxs)  
            if plot:
                ax.scatter(point_proj_x, point_proj_y, color=col,edgecolors=ecol,alpha=salpha,linewidths=linewidth,s = markersize, zorder=zorder) 
            pointout.append([point_proj_x,point_proj_y])
        Pointsout.append(pointout)
    #plt.show()        
    ax.set_aspect('equal', 'datalim')
    if out:
        return fig,ax,Pointsout,horizontal,vertical
    else:
        return fig,ax

def plot_atomic_plane3D(LatticePoints,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],Q=np.eye(3)):
    """
    plot_atomic_plane3D - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """

    if ax==None:
        fig = plt.figure() 
        ax = fig.add_subplot(111, projection='3d', proj_type = 'ortho') 
    else:
        fig=[]

    for Points1,col,ecol in zip(LatticePoints,colors,edgecolors):
        for points in Points1:
            points=np.matmul(Q,points)
            ax.scatter(points[0,:], points[1,:], points[2,:], s = 200,color=col,edgecolors=ecol,alpha=salpha) 
    
    #plt.show()        
    ax.axis('auto')      
        
    return fig,ax
   
    

def plot_atomlattice2D(atoms_xyz_position,uvw2xyz,normal,vertical,S=1,R=np.eye(3),ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5]):
    """
    plot_atomlattice2D - Crystallographic function for materials analysis.
    
    See full documentation in extended modules for detailed usage.
    """
    #colors=['r','b','g']
    if not isinstance(normal, np.ndarray):
        normal=np.array(normal);
    if not isinstance(vertical, np.ndarray):
        vertical=np.array(vertical);
    normal=normal/np.sqrt(normal.dot(normal))
    vertical=vertical/np.sqrt(vertical.dot(vertical))
    horizontal=np.cross(vertical,normal)

    if ax==None:
        fig = plt.figure() 
        ax = fig.add_subplot(111) 
    else:
        fig=[]
    Points=[]
    Max=[]
    Min=[]
    #S = 2
    S_range = list(range(-1,S+1)) 
    S_range = list(range(-S,S+1)) 
    for atoms,col,ecol in zip(atoms_xyz_position,colors,edgecolors):
        for atom in atoms:
            #print('ok')
            triplets = list(itertools.product(S_range, repeat=3)) 
            triplets = np.array(triplets) 
            triplets = triplets.T
            points = uvw2xyz.dot(triplets)+np.repeat([atom],triplets.shape[1],axis=0).T
            Points.append(points)
            point_proj_x = horizontal.dot(points)
            point_proj_y = vertical.dot(points)
            
            #pn=points
            #sq = np.sqrt(np.sum(points**2,axis=0))
            #idxs = np.where(sq>0.0)[0]
            #pn[:,idxs] = pn[:,idxs]/sq[idxs]
            point_proj_n = normal.dot(points)
            #print(point_proj_n)
            idxs = np.where(abs(point_proj_n)<1e-3)[0]
#            ax3i.plot(point_proj_x,point_proj_y,'r',alpha=alpha[0])            
            #if abs(point_proj_n)<1e-5:
            ax.scatter(point_proj_x[idxs], point_proj_y[idxs], s = 200,color=col,edgecolors=ecol,alpha=salpha) 
            
#            else:
#                ax.scatter(point_proj_x, point_proj_y, s = 200,color=col,edgecolors=ecol,alpha=salpha) 

#            #plt.show()
    
#    S = S+1
#    S_range = list(range(-1,S+1)) 
#    for atom in atoms_xyz_position[0][0]*0.:
#        triplets = list(itertools.product(S_range, repeat=3)) 
#        triplets = np.array(triplets) 
#        triplets = triplets.T
#        points = uvw2xyz.dot(triplets)+np.repeat([atom],triplets.shape[1],axis=0).T
#    Max=[max(p) for p in points]
#    Min=[min(p) for p in points]
#    
#            
#    for point in points.T:        
#        for lattice_vec in uvw2xyz.T:
#            vec=np.array([lattice_vec*0,lattice_vec]).T+np.repeat([point],2,axis=0).T
#            sum1=np.sum((vec>np.repeat([Max],2,axis=0).T).astype(int),axis=0)
#            sum1=sum1+np.sum((vec<np.repeat([Min],2,axis=0).T).astype(int),axis=0)
#            idx = np.where(sum1>0)[0]
#            if not len(idx)>0:
#                vec=np.matmul(R,vec)
#                ax.plot(vec[0,:], vec[1,:], vec[2,:],color=gridcolor,alpha=lalpha)
#            
    #plt.show()        
    ax.set_aspect('equal', 'datalim')
    
    return fig,ax,Points,normal,vertical, horizontal



def plot_cut2D(ax,Lattice_points,normal,horizontal,vertical,col,alpha=1):

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
    IP=[]
    eps=1e-2
    for point in Lattice_points:
        point=np.array(point);
        point_proj_x = horizontal.dot(point)
        point_proj_y = vertical.dot(point)
        point_proj_n = normal.dot(point)
        if (point_proj_n<= eps).any():
            if not (point_proj_n<= eps).all():
                idxp0=np.where(point_proj_n> eps)[0][0]                                           
                point[:,idxp0]=plane_line_intersection(normal,np.array([0,0,0]),point[:,0],point[:,1])
                point_proj_x = horizontal.dot(point)
                point_proj_y = vertical.dot(point)
                IP.append([point_proj_x[idxp0],point_proj_y[idxp0]])
            ax.plot(point_proj_x,point_proj_y,col,alpha=alpha)
                
#        else:
#            ax.plot(point_proj_x,point_proj_y,col)
        r=[]
        theta =[]
        eps=1e-4
        IPnew=[]
        for ip in IP:
            #print(ip)
            x=ip[0]
            y=ip[1]
            R=np.sqrt(x*x + y*y)
            TH=np.arctan2(y, x)*180./np.pi
            if TH<0:
                TH=360-abs(TH)
            if len(r)>0 and R>eps:
                include=True
                if min(abs(r-R))<eps: 
                    idx = np.where(abs(r-R)==min(abs(r-R)))[0][0]
                    if abs(theta[idx]-TH)<eps:
                        include=False
                        
                if include:
                    r.append(R)
                    theta.append(TH)
                    IPnew.append(ip)
            else:
                if R>eps:
                    r.append(R)
                    theta.append(TH)
                    IPnew.append(ip)
        idxs=np.argsort(theta)
#        for IPi in range(0,len(idxs)-1):
#            ax.plot([IPnew[idxs[IPi]][0],IPnew[idxs[IPi+1]][0]],[IPnew[idxs[IPi]][1],IPnew[idxs[IPi+1]][1]],'r')
        IPnew=np.array(IPnew)
        IPnew=IPnew[idxs]
    return IPnew
#        if len(IPnew)>0:
#            ax.add_patch(Polygon(IPnew, color=col,closed=True,fill=False, hatch=hatch))


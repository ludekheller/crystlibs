# PLOTLIB - Quick Reference

**Module**: `plotlib.py`  
**Functions**: 57  
**Last Updated**: December 19, 2025

---

**dataAnnot**
```python
def dataAnnot(self,**kwargs):
```
*Annotate data points with crystallographic information*

**dataShow**
```python
def dataShow(self,**kwargs):
```
*Display data information in plot or console*

**figsave**
```python
def figsave(self, **kwargs):
```
*Save figure to file(s) with optional cropping*

**figsaveproc**
```python
def figsaveproc(self, fname, **kwargs):
```
*Process and save figure to file*

**format_annot**
```python
def format_annot(self,x, y,**kwargs):
```
*Format annotation text for data points*

**format_coord**
```python
def format_coord(self,x, y,**kwargs):
```
*Format coordinates for display in plot toolbar*

**format_coord_test**
```python
def format_coord_test(self,x, y,**kwargs):
```
*Test version of coordinate formatting for debugging*

**genPoris**
```python
def genPoris(self,**kwargs):
```
*Generate pole figure orientation data for given crystal direction*

**generateSphericalHistSampleData**
```python
def generateSphericalHistSampleData(self,**kwargs):
```
*Generate spherical histogram from orientation data*

**generateSphericalKDESampleData**
```python
def generateSphericalKDESampleData(self,**kwargs):
```
*Generate spherical kernel density estimate from orientation data*

**getColormap**
```python
def getColormap(self,**kwargs):
```
*Generate orientation density colormap from sample data*

**getFigparam**
```python
def getFigparam(self, fontsize=None, save=False, phase='A,', figsize=None, **kwargs):
```
*Get figure parameter dictionary for saving high-quality figures*

**getScales**
```python
def getScales(self, vmcbar, numticks=None, ticks=None, tickslabels=None, geq=False, leq=False, cmapbins=100, cmapbinsmult=None):
```
*Generate color scale parameters including ticks, labels, and colormap*

**get_cmap**
```python
def get_cmap(colors, nbins=1000, name='my_cmap'):
```
*Create a custom colormap from a list of colors with smooth gradients*

**get_colors**
```python
def get_colors(values, cmap, vmin=None, vmax=None, cmin=0, cmax=None, to255=True, nancolor=[0, 0, 0, 1]):
```
*Map data values to colors using a colormap*

**onclicactivate**
```python
def onclicactivate(self,**kwargs):
```
*Activate mouse click event handling for the plot*

**onclick**
```python
def onclick(self,event):
```
*Mouse click event handler for data selection*

**onclick2**
```python
def onclick2(self,event):
```
*Alternative click handler for different interaction mode*

**onclick3**
```python
def onclick3(self,event):
```
*Third click handler variant for specialized interactions*

**onmove**
```python
def onmove(self,event):
```
*Mouse move event handler for interactive data display*

**onpress**
```python
def onpress(self,event):
```
*Keyboard press event handler for plot interaction*

**onpressActivate**
```python
def onpressActivate(self):
```
*Activate keyboard event handling for the plot*

**plotColorbar**
```python
def plotColorbar(self,**kwargs):
```
*Add colorbar to the plot with custom ticks and labels*

**plotColormap**
```python
def plotColormap(self,**kwargs):
```
*Plot orientation density colormap with contours on stereographic projection*

**plotColormaps**
```python
def plotColormaps(self,sel='all',**kwargs):
```
*Plot multiple colormaps in multi-panel figure (typically 2x2)*

**plotDirsNorms**
```python
def plotDirsNorms(self,**kwargs):
```
*Plot crystallographic directions and plane normals on stereographic projection*

**plotHist**
```python
def plotHist():
```
*Plot histogram of orientation data on stereographic projection*

**plotProj**
```python
def plotProj(self, **kwargs):
```
*Create a stereographic projection (pole figure) with appropriate grid*

**plotScatter**
```python
def plotScatter(self,**kwargs):
```
*Create scatter plot of crystallographic orientations or data points*

**plotScatterAsHist**
```python
def plotScatterAsHist(self,**kwargs):
```
*Display scatter plot data as histogram representation*

**plot_atomic_plane2D**
```python
def plot_atomic_plane2D(LatticePoints,normal,vertical,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],plot=True,salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],linewidths=[1,1,1],markersizes=[200,200,200], Q=np.eye(3),xlim=[],ylim=[],out=False,zorder=1):
```
*plot_atomic_plane2D - Crystallographic function for materials analysis*

**plot_atomic_plane3D**
```python
def plot_atomic_plane3D(LatticePoints,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],Q=np.eye(3)):
```
*plot_atomic_plane3D - Crystallographic function for materials analysis*

**plot_atomlattice2D**
```python
def plot_atomlattice2D(atoms_xyz_position,uvw2xyz,normal,vertical,S=1,R=np.eye(3),ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5]):
```
*plot_atomlattice2D - Crystallographic function for materials analysis*

**plot_cut2D**
```python
def plot_cut2D(ax,Lattice_points,normal,horizontal,vertical,col,alpha=1):
```
*Plot 2D cross-section of 3D points*

**plot_lattice**
```python
def plot_lattice(Points,LatticeVectors,ax=None,colors=['r','b','g'],edgecolors=['r','b','g'],salpha=1.,lalpha=1.,gridcolor=[0.5,0.5,0.5],Q=np.eye(3),shift=np.zeros(3),atoms=True,linewidth=1,move=np.zeros(3),normal=np.array([0,0,0]),halfspace='upper',s=200,plot=True):
```
*plot_lattice - Crystallographic function for materials analysis*

**plot_lattice2D**
```python
def plot_lattice2D(ax,VV,description,Parentlattice_points,Parent_lattice,\ Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2,xlim=None):
```
*Plot 2D projection of lattice*

**plot_lattice3D**
```python
def plot_lattice3D(ax,VV,description,Parentlattice_points,Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2):
```
*Complete 3D lattice visualization*

**plot_lattice_2Dprojection**
```python
def plot_lattice_2Dprojection(ax,VV,description,Parentlattice_points,Parent_lattice,\ Productlattice_points,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,normals,verticals, linewidth=2,xlim=None):
```
*Project lattice along specific crystallographic direction*

**plot_lattice_boundaries**
```python
def plot_lattice_boundaries(axl,LatticePointsNew,allPoints=None,polygon=False,tol=1e-1,**kwargs):
```
*Plot unit cell boundaries in 3D*

**plot_lattice_plane**
```python
def plot_lattice_plane(axl,PlanePoints,**kwargs):
```
*Plot a crystallographic plane in 3D lattice*

**plot_lattice_proj**
```python
def plot_lattice_proj(LatticeVectors,normalproj,verticalproj, ax=None, linewidth=2,color='b',eps=1e-1,Q=np.eye(3),Qprojr=np.eye(2),shift=np.zeros(3),shiftproj=0,move=np.zeros(3),normal=np.array([0,0,0]),shifthalfspace=np.zeros(3), halfspace='upper',shiftplot=np.array([0,0]),out=False):
```
*plot_lattice_proj - Crystallographic function for materials analysis*

**plot_latticefaces3D**
```python
def plot_latticefaces3D(ax,Parentlattices,linewidth=2,alpha=0.15,edgecolor='r',linestyle='-',facecolor=(1, 0, 0, 0.15)):
```
*plot_latticefaces3D - Crystallographic function for materials analysis*

**plot_latticesfaces3D**
```python
def plot_latticesfaces3D(ax,VV,description,Parentlattices,Productlattices,Product_uvw_2_Parent_uvw_all_norm,Product_uvw2xyz,linewidth=2,alpha=0.15,xlim=[-2,2]):
```
*Plot two lattices with transparent faces*

**plot_mohr_circles**
```python
def plot_mohr_circles(mcircles,VV,DD,xyz2uvw,scale,xticks=None,yticks=None,ax=None,Parent_lattice='B2'):
```
*Plot all three Mohr's circles*

**plot_planes_on_mohr_circle**
```python
def plot_planes_on_mohr_circle(ax,scale,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl, colors,text=False,Parent_lattice='B2'):
```
*Plot specific crystallographic planes on Mohr's circle*

**plot_planes_on_stereotriangle**
```python
def plot_planes_on_stereotriangle(ax,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl,colors):
```
*Plot planes on stereographic triangle*

**plot_planes_on_wulffnet**
```python
def plot_planes_on_wulffnet(ax,SelNormalsOnCircle,SelShearsOnCircle,Parent_xyz2hkl,colors):
```
*Plot plane normals on Wulff net (equal-angle projection)*

**plot_points_proj**
```python
def plot_points_proj(Points,normalproj,verticalproj, ax=None, marker="o",markersize=10, color='b',Q=np.eye(3),Qprojr=np.eye(2),shift=np.zeros(3),move=np.zeros(3)):
```
*plot_points_proj - Crystallographic function for materials analysis*

**plot_princip_dir_on_stereotriangle**
```python
def plot_princip_dir_on_stereotriangle(ax,VV,description,markersize=10,markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5):
```
*Plot principal strain/stress directions on stereographic triangle*

**plot_princip_dir_on_wulffnet**
```python
def plot_princip_dir_on_wulffnet(ax,VV,description,markersize=10,markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5):
```
*Plot principal directions on Wulff net*

**plotcolmap**
```python
def plotcolmap(fname=None, withdraw=False):
```
*Plot a single colormap from a saved pickle file*

**plotcolmaps**
```python
def plotcolmaps(fname=None, withdraw=False):
```
*Plot multiple colormaps from a saved pickle file*

**processScatterData**
```python
def processScatterData(self,**kwargs):
```
*Process scatter plot data by transforming orientations to projection coordinates*

**scatterDataAnnot**
```python
def scatterDataAnnot(self,**kwargs):
```
*Annotate scatter plot data points with information*

**setAttributes**
```python
def setAttributes(self, **kwargs):
```
*Set or update plotter attributes dynamically*

**set_aspect_equal_3d**
```python
def set_aspect_equal_3d(ax):
```
*Fix equal aspect ratio bug for 3D matplotlib plots*

**shiftedColorMap**
```python
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
```
*Shift the center of a colormap to a specific data value*


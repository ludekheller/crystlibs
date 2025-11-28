from orix import plot
import copy
import numpy as np
from orilib import *
from crystallography_functions import *
from projlib import *
from crystals import Crystal
from numpy import sqrt
from crystlib import *
from crystallography_functions import *
from effective_elastic_constants_functions import *


class getPhases:
    """
    Main class for managing crystallographic phases and their transformations.
    
    Handles phase definitions, orientation relationships, deformation gradients,
    and twin systems for materials like NiTi shape memory alloys.
    """
    
    def __init__(self):
        """
        Initialize the phase management system.
        
        Sets up dictionaries for phases, orientation relationships, variants,
        deformation gradients, and elastic properties.
        """
        self.phases={}
        self.OR={}
        self.vardict={}
        self.getOR()
        self.austenite='A'
        self.martensite='M'
        self.defGrad={}
        self.CV={}
        self.twinSys={}
        self.el={}
        
    def setAttributes(self,**kwargs):
        """
        Set arbitrary attributes for the class instance.
        
        Args:
            **kwargs: Key-value pairs of attributes to set
        """
        self.__dict__.update(kwargs)
        
    def fromCif(self,file, phasename):
        """
        Load crystal structure from CIF file.
        
        Args:
            file (str): Path to CIF file
            phasename (str): Name to assign to this phase
            
        Extracts lattice parameters, symmetry operations, and builds
        lattice matrices from CIF data.
        """
        self.phases[phasename]={}
        self.phases[phasename]['cif'] = Crystal.from_cif(file)
        self.phases[phasename]['ciffile'] = file
        self.phases[phasename]['symops'] = [sym[0:3,0:3] for cc,sym in enumerate(self.phases[phasename]['cif'].symmetry_operations())]
        self.phases[phasename]['recsymops'] = [sym[0:3,0:3] for cc,sym in enumerate(self.phases[phasename]['cif'].reciprocal_symmetry_operations())]

        for idx, key in enumerate('a b c alpha beta gamma'.split()):
            if idx<3:
                self.phases[phasename][key] = self.phases[phasename]['cif'].lattice_parameters[idx]
            else:
                self.phases[phasename][key] = self.phases[phasename]['cif'].lattice_parameters[idx]*np.pi/180.
        a= self.phases[phasename]['a']
        b= self.phases[phasename]['b']
        c= self.phases[phasename]['c']
        alpha= self.phases[phasename]['alpha']
        beta= self.phases[phasename]['beta']
        gamma= self.phases[phasename]['gamma']
        self.fromLatticeParams(phasename,a,b,c,alpha,beta,gamma)

    def fromLatticeParams(self,phasename,a,b=None,c=None,alpha=np.pi/2,beta=np.pi/2,gamma=np.pi/2):
        """
        Build lattice matrices from lattice parameters.
        
        Args:
            phasename (str): Name of the phase
            a,b,c (float): Lattice parameters
            alpha,beta,gamma (float): Lattice angles in radians
            
        Constructs direct and reciprocal lattice matrices, metric tensors,
        and their inverses for crystallographic calculations.
        """
        if b is None:
            b=a
        if c is None:
            c=a
        L = np.zeros((3,3))
        L[:,0] = np.array([a,0.,0])
        L[:,1] = np.array([b*np.cos(gamma),b*np.sin(gamma),0])
        cx=c*np.cos(beta)
        cy=c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        cz=np.sqrt(c**2-cx**2-cy**2)
        L[:,2] = np.array([cx,cy,cz])
        L=np.round(L,decimals=10)
        #lattice matrix converting [uvw]->[xyz]
        self.phases[phasename]['L']=L
        #inverse lattice matrix converting [xyz]->[uvw]
        self.phases[phasename]['LI']=inv(L)
        #Meteric tensor ||[x,y,x]||^2=[uvw]^T*G_M*[uvw]
        self.phases[phasename]['G']=L.T.dot(L)
        #Reciprocal Meteric tensor ||[x,y,x]||^2=d_hkl^2=[hkl]^T*Gr_M*[hkl]
        self.phases[phasename]['Gr']=inv(self.phases[phasename]['G'])
        #reciprocal lattice matrix converting [hkl]->[xyz]
        self.phases[phasename]['Lr']=inv(L).T
        #inverse reciprocal lattice plane matrix converting [xyz]->[hkl]
        self.phases[phasename]['LrI']=inv(self.phases[phasename]['Lr'])

    def getOR(self, name=None, OR='NiTi'):
        """
        Get orientation relationship matrices between phases.
        
        Args:
            name (str, optional): Name for the OR
            OR (str): Type of orientation relationship ('NiTi')
            
        Sets up correspondence matrices for directions (Cd, CId) and 
        planes (Cp, CIp) mapping between martensite and austenite lattices.
        """
        if name is None:
            name=OR
        self.vardict[name] = {}
        if OR == 'NiTi':
            #Cd, CId, Cp, CIp = B19p_B2_lattice_correspondence(notation='Waitz')
            ORD = B19p_B2_lattice_correspondence(notation='Waitz')
            self.OR[name] = {}
            for idx, key in enumerate('Cd CId Cp CIp'.split()):
                self.OR[name][key] = ORD[idx]

    def getEl(self,name='NiTi'):
        """
        Get elastic constants for phases.
        
        Args:
            name (str): Name identifier for elastic properties
            
        Sets up stiffness and compliance matrices for austenite and martensite
        phases using literature values for NiTi.
        """
        if name is None:
            self.el={}
        elif name=='NiTi':
            #Martensite elastic constants Wagner https://www.sciencedirect.com/science/article/pii/S1359645408006083
            CM_11=223e3;
            CM_12=129e3;
            CM_13=99e3;
            CM_15=27e3;
            CM_22=241e3;
            CM_23=125e3;
            CM_25=-9e3;
            CM_33=200e3;
            CM_35=4e3;
            CM_44=76e3;
            CM_46=-4e3;
            CM_55=21e3;
            CM_66=77e3;

            constants = {}
            constants[self.martensite]={'11':CM_11,'12':CM_12,'13':CM_13,'15':CM_15,\
                       '22':CM_22,'23':CM_23,'25':CM_25,\
                        '33':CM_33,'35':CM_35,\
                       '44':CM_44,'46':CM_46,\
                       '55':CM_55,'66':CM_66}

            C11=169e3
            C12=141e3
            C44=33e3

            constants[self.austenite]={'11':C11,'22':C11,'33':C11,'12':C12,'13':C12,'23':C12,'44':C44,'55':C44,'66':C44}

            self.el[name]={}
            for phase in [self.austenite,self.martensite]:
                self.el[name][phase]={}
                self.el[name][phase]['constants']=constants[phase]
                self.el[name][phase]['C'] = stiffness_matrix(constants[phase])
                self.el[name][phase]['CT'] = stiffness_from_voight_notation2tensor(self.el[name][phase]['C'])
                self.el[name][phase]['S'] = inv(self.el[name][phase]['C'])
                self.el[name][phase]['ST'] = compliance_from_voight_notation2tensor(self.el[name][phase]['S'])

    def getDefGrad(self, name=None, OR='NiTi'):
        """
        Compute deformation gradients for stress-free transformation.
        
        Args:
            name (str, optional): Name identifier
            OR (str): Orientation relationship type
            
        Calculates transformation deformation gradients, principal strains,
        and related quantities for all martensite variants.
        """
        if name is None:
            name=OR
        defGrad = {}
        if OR == 'NiTi':
            F_AM, U_AM, Q_M, T_MA, T_AM=def_gradient_stressfree(self.OR[name]['Cd'],self.phases[self.austenite]['L'], self.phases[self.martensite]['L'], CId=self.OR[name]['CId'])
            StressfreeRefs={}
            for space in [self.austenite,self.martensite]:
                TransformationStrain=[]
                EigVals=[]
                EigVecs=[]
                epsilon_2=[]
                lambda_2=[]
                strain=[]
                for var in range(0,F_AM.shape[2]):
                    if space==self.martensite:
                        #in Martensite space
                        #def. grad in mart space
                        F_AM_M=T_AM[:,:,var].dot(F_AM[:,:,var].dot(T_AM[:,:,var].T))
                        F=T_AM[:,:,var].dot(F_AM[:,:,var].dot(T_AM[:,:,var].T))
                        U_AM_M = scipy.linalg.sqrtm(F_AM_M.T.dot(F_AM_M))
                        Q_M_M = F_AM_M.dot(inv(U_AM_M));

                    else:
                        F=F_AM[:,:,var]
                    #right Cauchy-Green deformation tensor https://www.comsol.com/multiphysics/analysis-of-deformation
                    C=F.T.dot(F)
                    #Get eigenvalues (squares of principal stretches) and eigenvectors
                    D,V = np.linalg.eig(C)
                    #matmul(C,V) - matmul(V,D*np.eye(3))
                    Idxs = np.argsort(D)
                    #Sorted eigenvalues and eigenvectors
                    EigVals.append(D[Idxs])
                    EigVecs.append(V[:,Idxs])
                    #Eigenstrains using Green-Lagrange strain tensor
                    Eps=0.5*(EigVals[-1]-1)
                    strain.append(0.5*(EigVals[-1]-1))
                    epsilon_2.append(Eps[1])
                    lambda_2.append(EigVals[-1][1])

                #Convert to array 12 x number of  orientations

                defGrad[space]={}
                volchange = np.linalg.det(F_AM[:,:,0])-1
                #print(volchange)
                keys=['F_AM', 'U_AM', 'Q_M', 'F_AM_M', 'U_AM_M', 'Q_M_M', 'T_MA', 'T_AM','EigVals','EigVecs','epsilon_2','lambda_2','strain','volchange']
                for key in keys:
                    try:
                        exec(f"defGrad['{space}']['{key}']={key}")
                    except:
                        pass
        self.defGrad[name]=defGrad
        keys = ['F_AM', 'U_AM', 'Q_M', 'T_MA', 'T_AM']
        keys2 = ['EigVals', 'EigVecs', 'epsilon_2', 'lambda_2', 'strain']
        self.CV[name]=[]
        for idx in range(0,defGrad[self.austenite][keys[0]].shape[2]):
            dd={}
            for key in keys:
                if key=='F_AM':
                    dd['F_a'] = defGrad[self.austenite][key][:,:,idx]
                    dd[key] = defGrad[self.austenite][key][:,:,idx]
                else:
                    dd[key] = defGrad[self.austenite][key][:,:,idx]
            for key in keys2:
                dd[key] = defGrad[self.austenite][key][idx]
            self.CV[name].append(dd)
        self.getStrains(self.CV[name])
        self.getTwinSys()
        self.getHBVs()
        self.getHPVII011()
        self.getStrains(self.hpvII011[name])

    def getTrStrain(self,phase=None,name='NiTi', oris=None):
        """
        Calculate transformation strain for orientations.
        
        Args:
            phase (str, optional): Phase to compute strain for
            name (str): Name identifier
            oris: Orientations to evaluate strain at
            
        Returns transformation strain magnitude for given orientations.
        """
        if phase is None:
            phase=self.martensite
        if oris is None:
            oris=genori(dangle=1,hemi='upper', half='upper')
        TransformationStrain=[]
        for var in range(self.defGrad['NiTi'][self.austenite]['F_AM'].shape[2]):
            if phase==self.martensite:
                F=self.defGrad[name][phase]['F_AM_M']
            else:
                F=self.defGrad[name][phase]['F_AM'][:,:,var]
            TransformationStrain.append(np.sqrt(np.sum(oris*(F.T.dot(F.dot(oris))),axis=0))-1)
        self.trStrain=(oris,np.array(TransformationStrain))

    def getStrains(self,F_a,key='F_a'):
        """
        Calculate various strain measures from deformation gradient.
        
        Args:
            F_a: List of deformation gradient dictionaries
            key (str): Key to access deformation gradient
            
        Computes engineering, Lagrangian, Biot, Hencky, and Almansi strains
        for each deformation state.
        """
        for FF in F_a:
            F=FF[key]
            for strain_meas in ['Engineering strain','Lagrangian finite strain tensor','Biot strain tensor','Hencky strain tensor','Almansi strain tensor']:
                C=F.T.dot(F)
                if strain_meas=='Engineering strain':
                    FF[strain_meas]=0.5*(F+F.T)-np.eye(3)
                elif strain_meas=='Lagrangian finite strain tensor':
                    FF[strain_meas]=0.5*(C-np.eye(3))
                elif strain_meas=='Biot strain tensor':
                    FF[strain_meas]=(scipy.linalg.sqrtm(C)-np.eye(3))
                elif strain_meas == 'Hencky strain tensor':
                    FF[strain_meas]=0.5*scipy.linalg.logm(C)
                elif strain_meas == 'Hencky strain tensor2':
                    FF[strain_meas]=scipy.linalg.logm(scipy.linalg.sqrtm(C))
                elif strain_meas == 'Almansi strain tensor':
                    FF[strain_meas]=0.5*(np.eye(3)-np.linalg.inv(C))

            W2=0.5*(F-F.T)
            rv=np.array([W2[2,1],-W2[2,0],W2[1,0]])#W=np.cross(np.eye(3),rv)
            #rotation angle
            omega=norm(rv)
            #rotation axis normalized
            rvn=rv/norm(rv)
            #Normalized skew matrix is np.cross(np.eye(3),rvn) formed with normilized elements of rv
            #https://www.brainm.com/software/pubs/math/Rotation_matrix.pdf
            #Rotation matrix=I*cos(omega)+sin(omega)*np.cross(rvn,np.eye(3)) + (1-cos(omega))*np.outer(rvn,rvn)
            R=np.eye(3)*np.cos(omega)+np.sin(omega)*np.cross(np.eye(3),rvn)+(1-np.cos(omega))*np.outer(rvn,rvn)
            FF['RotAngle']=omega
            FF['RotMat']=R

    def getTwinSys(self, name=None, OR='NiTi'):
        """
        Compute twinning systems for martensite variants.
        
        Args:
            name (str, optional): Name identifier
            OR (str): Orientation relationship type
            
        Calculates Type I and Type II twinning systems and their
        crystallographic parameters.
        """
        if name is None:
            name=OR
        if OR == 'NiTi':
            B2_symops = self.phases[self.austenite]['symops']
            B2_recsymops = self.phases[self.austenite]['recsymops']
            B19p_symops = self.phases[self.martensite]['symops']
            B19p_recsymops = self.phases[self.martensite]['recsymops']

            U_AM = self.defGrad[name][self.austenite]['U_AM']
            T_AM = self.defGrad[name][self.austenite]['T_AM']
            Q_AM = self.defGrad[name][self.austenite]['Q_M']
            LA = self.phases[self.austenite]['L']
            LrA = self.phases[self.austenite]['Lr']
            LM = self.phases[self.martensite]['L']
            LrM = self.phases[self.martensite]['Lr']
            CId = self.OR[name]['CId']
            CIp = self.OR[name]['CIp']

            self.twinSys[name]= niti_twinning(B2_symops,B2_recsymops,B19p_recsymops,B19p_symops,U_AM,LA,LrA,LM,LrM,CId, CIp,T_AM,miller='no',Qv=Q_AM)
            self.defGrad[name]['notation']=[vr[0] for vr in self.twinSys[name]['001']['Variant pairs 2'][0]]
            self.defGrad[name]['realnotation']=[idx for idx,vr in enumerate(self.twinSys[name]['001']['Variant pairs 2'][0])]
        self.getCVdict()

    def getHBVs(self, name='NiTi'):
        """
        Calculate habit plane variants for NiTi.
        
        Args:
            name (str): Name identifier
            
        Computes Type I and Type II habit plane variants using the
        twinned habit plane approach.
        """
        #get Habit plane variants for NiTi:
        self.hpvs={}
        self.hpvs[name]={}
        for twt in ['Type I','Type II']:
            self.hpvs[name][twt]=[]
            for twtidx in range(0,len(self.twinSys[name][twt]['eta1_a'])):
                for idx in range(0,len(self.twinSys[name][twt]['eta1_a'][twtidx])):
                    var_i = self.twinSys[name][twt]['Variant pairs'][twtidx][idx][0]
                    var_j = self.twinSys[name][twt]['Variant pairs'][twtidx][idx][1]
                    T_AM_i = self.twinSys[name][twt]['T_AM'][twtidx][idx][0]
                    T_AM_j = self.twinSys[name][twt]['T_AM'][twtidx][idx][1]
                    RotMat_i = self.twinSys[name][twt]['RotMat'][twtidx][idx][0]
                    RotMat_j = self.twinSys[name][twt]['RotMat'][twtidx][idx][1]
                    RotAngle_i = self.twinSys[name][twt]['RotAngle'][twtidx][idx][0]
                    RotAngle_j = self.twinSys[name][twt]['RotAngle'][twtidx][idx][1]
                    var_i2 = self.twinSys[name][twt]['Variant pairs 2'][twtidx][idx][0]
                    var_j2 = self.twinSys[name][twt]['Variant pairs 2'][twtidx][idx][1]
                    Ui = self.defGrad[name][self.austenite]['U_AM'][:,:,var_i]
                    Uj = self.defGrad[name][self.austenite]['U_AM'][:,:,var_j]
                    Fi = self.defGrad[name][self.austenite]['F_AM'][:,:,var_i]
                    Fj = self.defGrad[name][self.austenite]['F_AM'][:,:,var_j]
                    Qij = self.twinSys[name][twt]['Q_a'][twtidx][idx]
                    a1 = self.twinSys[name][twt]['a1_a'][twtidx][idx]
                    n1 = self.twinSys[name][twt]['n1_a'][twtidx][idx]
                    K1_m = self.twinSys[name][twt]['K1_m'][twtidx][idx]
                    eta1_m = self.twinSys[name][twt]['eta1_m'][twtidx][idx]
                    K2_m = self.twinSys[name][twt]['K2_m'][twtidx][idx]
                    eta2_m = self.twinSys[name][twt]['eta2_m'][twtidx][idx]
                    K1_a = self.twinSys[name][twt]['K1_a'][twtidx][idx]
                    eta1_a = self.twinSys[name][twt]['eta1_a'][twtidx][idx]
                    K2_a = self.twinSys[name][twt]['K2_a'][twtidx][idx]
                    eta2_a = self.twinSys[name][twt]['eta2_a'][twtidx][idx]

                    addondata={'var_i':var_i,'var_j':var_j,'var_i2':var_i2,'var_j2':var_j2,'Type':twt,'twtidx':twtidx,'idx':idx,
                              'K1_m':K1_m,'K2_m':K2_m,'K1_a':K1_a,'K2_a':K2_a,
                              'eta1_m':eta1_m,'eta2_m':eta2_m,'eta1_a':eta1_a,'eta2_a':eta2_a,'T_AM_i':T_AM_i,'T_AM_j':T_AM_j,
                              'RotMat_i':RotMat_i,'RotMat_j':RotMat_j,'RotAngle_i':RotAngle_i,'RotAngle_j':RotAngle_j}
                    #solution to :
                    #1. QiUi - Uj =a1 x n1
                    #2. Q_a(Lambda*QiUi+(1-Lambda)Uj)= I + b_a x m_a
                    #twinnedhabitplane(Ui,Uj,Qj,...) solves Q_a(Lambda*QjUj+(1-Lambda)Ui)= I + b_a x m_a, Lambda - Lj
                    self.hpvs[name][twt] = twinnedhabitplane(Ui, Uj, Qij, a1, n1, hbplanes=self.hpvs[name][twt], addondata=addondata, method='bhata')

    def getHPVII011(self, name='NiTi'):
        """
        Select commonly observed Type II <011> habit plane variants.
        
        Args:
            name (str): Name identifier
            
        Filters habit plane variants according to experimental observations
        from literature.
        """
        ### Selection of mostly observed habit plane variants Type II <011>_m
        #habit plane variants selected according to  https://www.sciencedirect.com/science/article/pii/0001616089900722
        #compare with s page 1876
        self.hpvII011 = {}
        self.hpvII011[name] = []
        self.vardict[name]['HPV']={}
        self.vardict[name]['CV2HPV']={}
        count=0
        for hb in self.hpvs[name]['Type II']:
            #print(hb['twtidx'])
            if hb['twtidx']==0 and count==1:
                ma=hb['m_a']
                Lj=hb['Lj']
                break
            else:
                count+=1
        count = -1
        #print(ma)
        #print(Lj)
        for hb in self.hpvs[name]['Type II']:
            #if hb['twtidx']==0 and (abs(abs(hb['m_a'])-abs(ma[0]))<1e-6).any() and (abs(hb['Lj']-Lj)<1e-6).any():
            if (abs(abs(hb['m_a'])-abs(ma[0]))<1e-6).any() and (abs(hb['Lj']-Lj)<1e-6).any():
                notation=hb['var_i2']+"(+)"
                if "'" in hb['var_j2']:
                    notation=notation.replace("(+)","(-)")
                #print(notation)
                count += 1
                self.hpvII011[name].append(hb)
                self.hpvII011[name][-1]['notation'] = notation
                self.vardict[name]['HPV'][notation] = count
                self.vardict[name]['CV2HPV'][(hb['var_i'],hb['var_j'])] = notation

    def getCVdict(self, name='NiTi'):
        """
        Create dictionary mapping variant notations to indices.
        
        Args:
            name (str): Name identifier
        """
        self.vardict[name]['CV'] = {}
        for notation, var in zip(self.twinSys[name]['001']['Variant pairs 2'][0],self.twinSys[name]['001']['Variant pairs'][0]):
            self.vardict[name]['CV'][notation[0]]=var[0]

    def printHBVII011TrStrain(self, name='NiTi'):
        """
        Print transformation strain statistics for Type II <011> habit planes.
        
        Args:
            name (str): Name identifier
            
        Calculates and displays min/max transformation strains for
        the selected habit plane variants.
        """
        ### Get transformation strain of mostly observed habit plane variants Type II <011>_m
        #Generate orientations in space
        oris=genori(dangle=1,hemi='upper')
        #Generate variable to plot in color
        grid_eps_tr=[]
        grid_eps_tr_dic={}
        for hb in self.hpvII011[name]:
            #grid_eps_tr.append(np.sqrt(np.sum(oris*(hbplanes['F_a'][var].T.dot(hbplanes['F_a'][var].dot(oris))),axis=0))-1)

            engeps=0.5*(hb['F_a'].T.dot(hb['F_a'])-np.eye(3))

            grid_eps_tr.append(np.sum(oris*engeps.dot(oris)*100,axis=0))
            grid_eps_tr_dic[hb['notation']]=np.sum(oris*engeps.dot(oris)*100,axis=0)
        EPStr=[]
        EPStr.append(np.max(np.array(grid_eps_tr),axis=0))
        EPStr.append(np.min(np.array(grid_eps_tr),axis=0))
        EPStr.append(np.min(np.abs(np.array(grid_eps_tr)),axis=0))
        EPStr=np.asarray(EPStr).T
        print("Min. maximum transformation strain: {:0.2f}".format(min(EPStr[:,0])))
        print("Max. maximum transformation strain: {:0.2f}".format(max(EPStr[:,0])))
        print("Min. minimum transformation strain: {:0.2f}".format(min(EPStr[:,1])))
        print("Max. minimum transformation strain: {:0.2f}".format(max(EPStr[:,1])))
        print("Min. abs minimum transformation strain: {:0.2f}".format(min(EPStr[:,2])))
        print("Max. abs minimum transformation strain: {:0.2f}".format(max(EPStr[:,2])))

    def printHBVII011(self, name='NiTi'):
        """
        Print information about Type II <011> habit plane variants.
        
        Args:
            name (str): Name identifier
        """
        print('Mostly observed habit plane variants Type II <011>_m' )
        for hb in self.hpvII011[name]:
            print('{}:{}-{},m:{},b:{},Lambda:{}'.format(hb['notation'],hb['var_i2'],hb['var_j2'],hb['m_a'],hb['b_a'],hb['Lj']))

    def printCV(self, name='NiTi'):
        """
        Print correspondence variant information in formatted table.
        
        Args:
            name (str): Name identifier
            
        Returns:
            pandas Styler: Formatted DataFrame with variant information
        """
        #print('Corresponding variants' )
        rows_list = []
        Cd = self.OR[name]['Cd']
        T_MA = self.defGrad[name][self.austenite]['T_MA']
        for idx,notation in enumerate(self.twinSys[name]['001']['Variant pairs 2'][0]):
            var=self.twinSys[name]['001']['Variant pairs'][0][idx][0]
            #print(var)
            #print(notation[0])
            pddict={}
            pddict[r'\\bf{Martensite Variant}    ']=var+1
            pddict['Prime notation']=notation[0]
            M=[]
            for col,dirm in zip(['\\bf{[100]}','\\bf{[010]}','\\bf{[001]}'],[[1,0,0],[0,1,0],[0,0,1]]):
                pddict[col+r'$_\mathrm{\mathbf{M}}$']=Cd[:,:,var].dot(dirm).astype(int)
                M.append(Cd[:,:,var].dot(dirm))
            if int(var/2+1)!=(var/2+1):
                M=T_MA[:,:,var]
                #print(M.dot(MI.T))
                r = R.from_matrix(M.dot(MI.T))
                pddict['180 deg rotation sym.']=r'Q[180 deg. around '+ str(r.as_rotvec()/np.abs(r.as_rotvec().max()))+r'$_\mathrm{A}$]'
            else:
                MI=T_MA[:,:,var]
                pddict['180 deg rotation sym.']=r"(001)$_\mathrm{M}$ twin between "+notation[0]+'-'+notation[0]+"'"
            rows_list.append(pddict)
        #print(rows_list)
        pd.set_option("display.max_column", None)
        pd.set_option("display.max_colwidth", None)
        pd.set_option('display.width', -1)
        pd.set_option('display.max_rows', None)
        pd.set_option('display.colheader_justify', 'center')
        df = pd.DataFrame(rows_list)
        df.style.hide().set_properties(**{'text-align': 'center'})
        df_styled=df.style.hide().set_properties(**{'text-align': 'center'}).format(formatter={(r'\\bf{[100]}$_\mathrm{\mathbf{M}}$'):lambda x:str(x)+r'$_\mathrm{{A}}$',
                                                                                          (r'\\bf{[010]}$_\mathrm{\mathbf{M}}$'):lambda x:str(x)+r'$_\mathrm{{A}}$',
                                                                                          (r'\\bf{[001]}$_\mathrm{\mathbf{M}}$'):lambda x:str(x)+r'$_\mathrm{{A}}$'}).apply(highlight_everyother)
        #print(df_styled.to_latex(convert_css=True))
        return df_styled

    def plotHBVII011(self, name='NiTi',ProjType='stereo'):
        """
        Plot habit plane variants on stereographic or equal-area projection.
        
        Args:
            name (str): Name identifier
            ProjType (str): Projection type ('stereo' or 'equalarea')
        """
        #ProjType=='equalarea'
        #plot shape strains of mostly observed habit plane variants Type II <011>_m
        if ProjType=='equalarea':
            fig,ax=schmidtnet(basedirs=False,facecolor='None')
        else:
            fig,ax=wulffnet(basedirs=False)

        #fig.suptitle(twt+' habit plane normals \n for '+str(HB[list(HB.keys())[1]])+' type plane')
        L = self.phases[self.austenite]['L']
        Lr = self.phases[self.austenite]['Lr']
        symops = self.phases[self.austenite]['symops']
        recsymops = self.phases[self.austenite]['recsymops']

        for hb in self.hpvII011[name]:
            hbp=1
            if Lr.dot(hb['m_a'])[2]<0:
                hbp=-1
            if ProjType=='equalarea':
                equalarea=True
                poris = equalarea_directions(hbp*hb['m_a'])
            else:
                equalarea=False
                poris = stereoprojection_directions(hbp*hb['m_a'])
            shape_strain=hbp*hb['b_a']
            shape_strain_point=hbp*hb['m_a']+shape_strain
            hbp2=1
            if shape_strain_point[2]<0:
                hbp2=-1
            if ProjType=='equalarea':
                equalarea=True
                poris2 = equalarea_directions(hbp2*shape_strain_point)[:,0]

            else:
                equalarea=False
                poris2 = stereoprojection_directions(hbp2*shape_strain_point)[:,0]

            ax.scatter(poris[0,:],poris[1,:],c='r',s=50,edgecolors='k',alpha=1,linewidths=1, zorder=100)

            shape_strain_proj=poris2-poris[:,0]
            shape_strain_proj=shape_strain_proj/np.linalg.norm(shape_strain_proj)*0.1
            ax.arrow(poris[0,0],poris[1,0],shape_strain_proj[0],shape_strain_proj[1],head_width=0.03, head_length=0.05, fc='k', ec='k', zorder=99)
            ax.text(poris[0,:]+0.02,poris[1,:]+0.02,hb['notation'],color='r', zorder=100)

        toPlot = [[1,0,0],[1,1,0],[1,1,1]]
        dirs = self.getEqUVW(self.austenite, toPlot,R2Proj=np.eye(3),hemisphere = "upper")
        normals = self.getEqHKL(self.austenite, toPlot, R2Proj=np.eye(3),hemisphere = "upper")
        d2plot=[(d[ProjType][0:2],d['uvw']) for d in dirs if d['vector'][2]>=0 or np.abs(d['vector'][2])<1e-5]#np.array([d['equalarea'][0:2] for d in dirs if d['vector'][2]>=0 or np.abs(d['vector'][2])<1e-5])
        for d2p in d2plot:
            ax.scatter(d2p[0][0],d2p[0][1],c='b',s=50,edgecolors='k',alpha=1,linewidths=1, zorder=100)
            ax.text(d2p[0][0]+0.02,d2p[0][1]+0.02,'[{:d}{:d}{:d}]'.format(int(d2p[1][0]),int(d2p[1][1]),int(d2p[1][2])),color='b', zorder=100)

        toPlot = [[1,1,0],[1,1,1]]
        normals = self.getEqHKL(self.austenite, toPlot, R2Proj=np.eye(3),hemisphere = "upper")
        p2plot=[(n[ProjType+' plane'][0:2],n['hkl']) for n in normals]
        for p2p in p2plot:
            if p2p[1]!= [0,0,1] and p2p[1]!= [0,0,-1]:
                ax.plot(p2p[0][0,:],p2p[0][1,:],color=[0.6,0.6,0.6],linewidth=0.5)

    def reduceCorresp_uvws(self,uvws,name='NiTi'):
        """
        Reduce martensite directions to unique austenite equivalents.
        
        Args:
            uvws: List of martensite direction vectors
            name (str): Name identifier
            
        Returns:
            tuple: (martensite_uvws, austenite_uvws, rotation_angles)
        """
        #uvws - martensite dirs
        austenite_uvws = []
        martensite_uvws = []
        rot = []
        austenite = self.austenite
        martensite = self.martensite
        for muvw in uvws:
            auvws=[]
            for Variant in range(0,self.OR['NiTi']['Cd'].shape[2]):
                auvw = list(self.OR['NiTi']['Cd'][:,:,Variant].dot(muvw))
                auvws.append(auvw)
                Productuvw_xyz=self.defGrad['NiTi'][austenite]['T_MA'][:,:,Variant].dot(self.phases[martensite]['L'].dot(muvw))
                Productuvw_xyz=Productuvw_xyz/np.linalg.norm(Productuvw_xyz)
                Parentuvw_xyz=self.phases[austenite]['L'].dot(auvw)
                Parentuvw_xyz=Parentuvw_xyz/np.linalg.norm(Parentuvw_xyz)
                roti= np.arccos(np.round(Productuvw_xyz.dot(Parentuvw_xyz),10))

            eq_els = []
            for m in self.phases[austenite]['symops']:
                idxs=np.where(abs(m)<1e-10)
                m[idxs[0],idxs[1]]=0.0
                eq_els.append(m.dot(auvws[0]))
            for auvwi in eq_els:
                if auvwi[0]>=0 and auvwi[1]>=0 and auvwi[2]>=0:
                    auvw=auvwi
                    break
            austenite_uvws.append(auvw)
            martensite_uvws.append(muvw)
            rot.append(roti)
        return martensite_uvws, austenite_uvws, rot

    def reduceCorresphkls(self,martensite_planes_hkl2,name='NiTi'):
        """
        Reduce martensite plane indices to unique austenite equivalents.
        
        Args:
            martensite_planes_hkl2: List of martensite plane indices
            name (str): Name identifier
            
        Returns:
            tuple: (martensite_planes, austenite_planes, rotation_angles)
        """
        austenite_planes_hkl = []
        martensite_planes_hkl = []
        rot = []
        austenite = self.austenite
        martensite = self.martensite
        for mhkl in martensite_planes_hkl2:
            ahkls=[]
            for Variant in range(0,self.OR['NiTi']['Cp'].shape[2]):
                ahkl = list(self.OR['NiTi']['Cp'][:,:,Variant].dot(mhkl))
                ahkls.append(ahkl)
                Producthkl_xyz=self.defGrad['NiTi'][austenite]['T_MA'][:,:,Variant].dot(self.phases[martensite]['Lr'].dot(mhkl))
                Producthkl_xyz=Producthkl_xyz/np.linalg.norm(Producthkl_xyz)
                Parenthkl_xyz=self.phases[austenite]['Lr'].dot(ahkl)
                Parenthkl_xyz=Parenthkl_xyz/np.linalg.norm(Parenthkl_xyz)
                roti= np.arccos(np.round(Producthkl_xyz.dot(Parenthkl_xyz),10))

            eq_els = []
            for m in self.phases[austenite]['recsymops']:
                idxs=np.where(abs(m)<1e-10)
                m[idxs[0],idxs[1]]=0.0
                eq_els.append(m.dot(ahkls[0]))
            for ahkli in eq_els:
                if ahkli[0]>=0 and ahkli[1]>=0 and ahkli[2]>=0:
                    ahkl=ahkli
                    break
            austenite_planes_hkl.append(ahkl)
            martensite_planes_hkl.append(mhkl)
            rot.append(roti)
        return martensite_planes_hkl, austenite_planes_hkl, rot

    def getUniqueCorresp(self,dir1, name='NiTi',inphase='M'):
        """
        Get unique corresponding directions between phases.
        
        Args:
            dir1: Input direction vector
            name (str): Name identifier
            inphase (str): Starting phase ('M' for martensite)
            
        Returns:
            array: Unique corresponding directions
        """
        #a - array Nx3 of dir/plane miller indexes indexes
        #dir1 - direction in phase 1
        if inphase=='M':
            #find martensite correspondence
            dirs2=vectors2miller(self.getUnique(np.array([v for v in (np.einsum('ijk,j->ki',self.OR[name]['CId'],dir1))])).T).T
        else:
            dirs2=vectors2miller(self.getUnique(np.array([v for v in (np.einsum('ijk,j->ki',self.OR[name]['Cd'],dir1))])).T).T
        return dirs2

    def getUnique(self,a):
        """
        Get unique vectors from array.
        
        Args:
            a: Array of vectors
            
        Returns:
            array: Unique vectors
        """
        #a - array Nx3 of dir/plane miller indexes indexes
        return np.unique(a.view(np.dtype((np.void, a.dtype.itemsize*a.shape[1])))).view(a.dtype).reshape(-1, a.shape[1])

    def generate_hkls(self,hklmax, phase,hkls=[]):
        """
        Generate unique HKL indices up to maximum index.
        
        Args:
            hklmax (int): Maximum Miller index
            phase (str): Phase name
            hkls: Existing HKL list to extend
            
        Returns:
            tuple: HKL generation results
        """
        sysms = self.phases[phase]['recsymops']
        hkls = generate_hkls(hklmax, sysms,hkls=hkls)
        hkls = (hkls[0], hkls[1], hkls[2], np.array([u for key in hkls[1].keys() for u in (hkls[1][key])]))
        return hkls

    def generate_uvws(self,uvwmax, phase,uvws=[]):
        """
        Generate unique UVW indices up to maximum index.
        
        Args:
            uvwmax (int): Maximum direction index
            phase (str): Phase name
            uvws: Existing UVW list to extend
            
        Returns:
            tuple: UVW generation results
        """
        sysms = self.phases[phase]['symops']
        uvws = generate_hkls(uvwmax,sysms, hkls=uvws)
        uvws = (uvws[0], uvws[1], uvws[2], np.array([u for key in uvws[1].keys() for u in (uvws[1][key])]))
        return  uvws

    def getEqHKL(self,phasename, hkls, R2Proj=np.eye(3),hemisphere = "upper", **kwargs):
        """
        Get equivalent HKL planes with projections.
        
        Args:
            phasename (str): Phase name
            hkls: List of HKL indices
            R2Proj: Rotation matrix for projection
            hemisphere (str): Hemisphere for projection ('upper' or 'lower')
            
        Returns:
            list: Equivalent HKL information with projections
        """
        #get all equivalent hkls, including stereo and equalarea projection into upper/lower hemisphere of a standard lattice coordinate system rotated by R2proj
        #for all hkls in the list of hkls
        L = self.phases[phasename]['L']
        Lr = self.phases[phasename]['Lr']
        symops = self.phases[phasename]['symops']
        recsymops = self.phases[phasename]['recsymops']
        eqhkls=[]
        for hkl in hkls:
                nv=R2Proj.dot(Lr.dot(hkl))
                nv/=np.linalg.norm(nv)
                isin=False
                for n in eqhkls:
                    if list(n['hkl']) == list(hkl):#np.linalg.norm(n['vector']-nvs)<1e-10:
                        isin=True
                        break
                if not isin:

                    eqhkls.append({'vector':nv,'hkl':hkl,'hklf':hkl,'label':str(hkl).replace('[','(').replace(']',')').replace(' ',''),
                                    'equalarea':equalarea_directions(nv),'stereo':stereoprojection_directions(nv),
                                   'equalarea plane':equalarea_planes(nv,hemisphere=hemisphere),
                                    'stereo plane':stereoprojection_planes(nv,hemisphere=hemisphere)})
                    if not recsymops is None:
                        for rs in recsymops:

                            hklsym=np.round(rs.dot(hkl))
                            isin=False
                            for n in eqhkls:
                                if True:#n['hklf']==hkl:
                                    if list(hklsym) == list(n['hkl']):#np.linalg.norm(n['vector']-nvs)<1e-10:
                                        isin=True
                                        break
                            if not isin:
                                #print(list(hklsym))
                                #print(hkl)
                                #print('------------------------------')
                                nvs=R2Proj.dot(rs.dot(Lr.dot(hkl)))
                                nvs/=np.linalg.norm(nvs)
                                eqhkls.append({'vector':nvs,'hkl':list(hklsym),'hklf':hkl,'label':str(hklsym).replace('[','(').replace(']',')').replace(' ',''),
                                                'equalarea':equalarea_directions(nvs),'stereo':stereoprojection_directions(nvs),
                                               'equalarea plane':equalarea_planes(nvs,hemisphere=hemisphere),
                                                'stereo plane':stereoprojection_planes(nvs,hemisphere=hemisphere)})
                #print('=============================================')
        return eqhkls

    def getEqUVW(self,phasename, uvws,R2Proj=np.eye(3),hemisphere = "upper", **kwargs):
        """
        Get equivalent UVW directions with projections.
        
        Args:
            phasename (str): Phase name
            uvws: List of UVW indices
            R2Proj: Rotation matrix for projection
            hemisphere (str): Hemisphere for projection
            
        Returns:
            list: Equivalent UVW information with projections
        """
        #get all equivalent uvws, including stereo and equalarea projection into upper/lower hemisphere of a standard lattice coordinate system rotated by R2proj
        #for all uvws in the list of uvws
        L = self.phases[phasename]['L']
        Lr = self.phases[phasename]['Lr']
        symops = self.phases[phasename]['symops']
        recsymops = self.phases[phasename]['recsymops']

        equvws=[]
        for uvw in uvws:
                dv=R2Proj.dot(L.dot(uvw))
                dv/=np.linalg.norm(dv)
                isin=False
                for d in equvws:
                    if list(d['uvw']) == list(uvw):#np.linalg.norm(n['vector']-nvs)<1e-10:
                        isin=True
                        break
                if not isin:
                    equvws.append({'vector':dv,'uvw':uvw,'uvwf':uvw,'label':str(uvw),
                                    'equalarea':equalarea_directions(dv),'stereo':stereoprojection_directions(dv)})
                    if not symops is None:
                        for rs in symops:
                            #dvs=rs.dot(dv)
                            uvwsym=np.round(rs.dot(uvw))
                            isin=False
                            for d in equvws:
                                if True:#d['uvwf']==uvw:
                                    if list(uvwsym) == list(d['uvw']):#np.linalg.norm(d['vector']-dvs)<1e-10:
                                        isin=True
                                        break
                            if not isin:
                                dvs=R2Proj.dot(rs.dot(L.dot(uvw)))
                                dvs/=np.linalg.norm(dvs)
                                equvws.append({'vector':dvs,'uvw':list(uvwsym),'uvwf':uvw,'label':str(uvwsym).replace('[','(').replace(']',')').replace(' ',''),
                                                'equalarea':equalarea_directions(dvs),'stereo':stereoprojection_directions(dvs)})

        return equvws

    def genSlipSys(self,phasename,SlipHKL,SlipUVW,mag=1):
        """
        Generate slip system deformation gradient and strains.
        
        Args:
            phasename (str): Phase name
            SlipHKL: Slip plane indices
            SlipUVW: Slip direction indices
            mag (float): Magnitude of slip
            
        Returns:
            dict: Slip system information including deformation gradient and strains
        """
        #generate deformation gradient and strains of a slip system
        L = self.phases[phasename]['L']
        Lr = self.phases[phasename]['Lr']
        symops = self.phases[phasename]['symops']
        recsymops = self.phases[phasename]['recsymops']

        SlipSi={}

        SlipPlane=Lr.dot(SlipHKL)
        SlipPlane/=norm(SlipPlane)
        SlipDir=L.dot(SlipUVW)
        SlipDir/=norm(SlipDir)
        SlipDir*=mag
        sm={}
        sm['SlipSystemFamily']={'hkl':SlipHKL,'uvw':SlipHKL}
        sm['SlipSystem']={'n':SlipPlane,'b':SlipDir,'hkl':SlipHKL,'uvw':SlipHKL}
        sm['SlipStrainNonSym']=np.outer(SlipDir,SlipPlane)
        sm['SlipDefGrad']=np.eye(3)+np.outer(SlipDir,SlipPlane)
        sm['SlipStrainSym']=0.5*(sm['SlipStrainNonSym']+sm['SlipStrainNonSym'].T)
        sm['SlipStrainAntiSym']=0.5*(sm['SlipStrainNonSym']-sm['SlipStrainNonSym'].T)

        omega,R = getRfromSkew(sm['SlipStrainAntiSym'])

        sm['RotAngle']=omega
        sm['RotMat']=R
        #SlipS.append(sm)
        SlipSi['{'+str(SlipHKL)+'}'+'<'+str(SlipUVW)+'>']=sm
        if False:
            v=np.array([0,0,1])
            for key in list(sm.keys()):
                print(f'{key}:\n{sm[key]}')
                exec(f'{key}=sm["{key}"]')
            RotMat = sm['RotMat']
            print('Analyzing the deformation gradient of the slip')
            print('Analyzing the meaning of the asymetric deformation gradient of the slip vs. its symmetrized version')
            print('The symmetrized version subtract a rotation from the asymetric deform gradient')
            print(f'Deformation gradient of the slip \n{SlipDefGrad}')
            SymDefGrad = (SlipStrainSym+np.eye(3))
            print(f'Symmetrized Deformation gradient of the slip constructed from the symmetric shear strain \n{SymDefGrad}')
            print(f'Rotation that have to be added to the symmetrized deformation gradient to get the same effect as that of asymmetric def. gradient:\n{RotMat}')
            v=np.array([0,0,1])
            print(f'Probe vector {v}')
            print(f'Effect of rotation matrix on the probe vector \n{RotMat.dot(v)}')
            print(f'Effect of asymmetric deformation gradient of the slip on the probe vector \n{SlipDefGrad.dot(v)}')
            print(f'Effect of symmetrized deformation gradient of the slip on the probe vector \n{SymDefGrad.dot(v)}')
            print(f'Effect of symmetrized deformation gradient of the slip on the probe vector when adding the rotation \n{RotMat.dot(SymDefGrad.dot(v))}~={SlipDefGrad.dot(v)}')
            print(f'Effect of asymmetric deformation gradient of the slip on the probe vector when removing the rotation \n{RotMat.T.dot(SlipDefGrad.dot(v))}~={SymDefGrad.dot(v)}')
        return SlipSi

    def genEqSlipSys(self,phasename,SlipHKLs,SlipUVWs, getIndependentOnly=True, tol=1e-10,mag=1):
        """
        Generate equivalent slip systems considering crystal symmetry.
        
        Args:
            phasename (str): Phase name
            SlipHKLs: List of slip plane indices
            SlipUVWs: List of slip direction indices
            getIndependentOnly (bool): Return only independent slip systems
            tol (float): Tolerance for orthogonality check
            mag (float): Slip magnitude
            
        Returns:
            dict: Equivalent slip systems
        """
        L = self.phases[phasename]['L']
        Lr = self.phases[phasename]['Lr']
        symops = self.phases[phasename]['symops']
        recsymops = self.phases[phasename]['recsymops']

        if type(SlipHKLs) == list:
            if type(SlipHKLs[0]) != list:
                SlipHKLs=[SlipHKLs]
        if type(SlipHKLs) == np.ndarray:
            if len(SlipHKLs.shape) == 1:
                SlipHKLs=np.array([SlipHKLs])
        if type(SlipUVWs) == list:
            if type(SlipUVWs[0]) != list:
                SlipUVWs=[SlipUVWs]
        if type(SlipUVWs) == np.ndarray:
            if len(SlipUVWs.shape) == 1:
                SlipUVWs=np.array([SlipUVWs])

        SlipSi={}
        for SlipPlaneHKL,SlipDirUVW in zip(SlipHKLs,SlipUVWs):
            SlipS=[]
            SlipSAll=[]
            SlipPlane=Lr.dot(SlipPlaneHKL)
            SlipPlane/=norm(SlipPlane)
            SlipDir=L.dot(SlipDirUVW)
            SlipDir/=norm(SlipDir)
            SlipDir*=mag
            for rsym in recsymops:
                splane=rsym.dot(SlipPlane)
                splane[abs(splane)<1e-10]=np.round(splane[abs(splane)<1e-10])
                splaneHKL=rsym.dot(SlipPlaneHKL)
                splaneHKL[abs(splaneHKL)<1e-10]=np.round(splaneHKL[abs(splaneHKL)<1e-10])
                for sym in symops:
                    sdir=sym.dot(SlipDir)
                    sdir[abs(sdir)<1e-10]=np.round(sdir[abs(sdir)<1e-10])
                    sdirUVW=sym.dot(SlipDirUVW)
                    sdirUVW[abs(sdirUVW)<1e-10]=np.round(sdirUVW[abs(sdirUVW)<1e-10])
                    if abs(splane.dot(sdir))<tol:
                        isin=False
                        for sm in SlipS:
                            ss=sm['SlipSystem']
                            if (ss['n']==splane).all() and (ss['b']==sdir).all():
                                isin=True
                            if (ss['hkl']==splaneHKL).all() and (ss['uvw']==sdirUVW).all():
                                isin=True
                            if (ss['n']==-1*splane).all() and (ss['b']==-1*sdir).all():
                                isin=True
                            if (ss['hkl']==-1*splaneHKL).all() and (ss['uvw']==-1*sdirUVW).all():
                                isin=True
                        if not isin:
                            sm={}
                            sm['SlipSystemFamily']={'hkl':SlipPlaneHKL,'uvw':SlipDirUVW}
                            sm['SlipSystem']={'n':splane,'b':sdir,'hkl':rsym.dot(SlipPlaneHKL),'uvw':sym.dot(SlipDirUVW)}
                            sm['SlipStrainNonSym']=np.outer(sdir/norm(sdir),splane/norm(splane))
                            sm['SlipDefGrad']=np.eye(3)+sm['SlipStrainNonSym']
                            sm['SlipStrainSym']=0.5*(sm['SlipStrainNonSym']+sm['SlipStrainNonSym'].T)
                            sm['SlipStrainAntiSym']=0.5*(sm['SlipStrainNonSym']-sm['SlipStrainNonSym'].T)
                            W2=sm['SlipStrainAntiSym']
                            rv=np.array([W2[2,1],-W2[2,0],W2[1,0]])
                            omega=norm(rv)
                            rvn=rv/norm(rv)
                            R=np.eye(3)*np.cos(omega)+np.sin(omega)*np.cross(np.eye(3),rvn)+(1-np.cos(omega))*np.outer(rvn,rvn)
                            sm['RotAngle']=omega
                            sm['RotMat']=R
                            SlipS.append(sm)
                        SlipSAll.append(sm)

            Am=np.empty((6,len(SlipS)))
            for i in range(0,len(SlipS)):
                Am[:,i]=SlipS[i]['SlipStrainSym'][np.triu_indices(SlipS[i]['SlipStrainSym'].shape[0],k=0)]

            _, inds=Matrix(Am).rref()
            IndependentSS=[SlipS[idx] for idx in inds]
            key='{'+str(SlipPlaneHKL)+'}'+'<'+str(SlipDirUVW)+'>'
            if getIndependentOnly:
                SlipSi[key]=IndependentSS
            else:
                SlipSi[key]=SlipS
            for sS in SlipSi[key]:
                sS['No. of Independent SS']=len(inds)

        return SlipSi

    def setPlotAttributes(self,defs={},name='NiTi'):
        """
        Set default plotting attributes for different crystal systems.
        
        Args:
            defs (dict): Custom attribute definitions
            name (str): Name identifier
        """
        defaults = {}
        defaults['dirnormeq']=False
        defaults['printasfamily']=False
        defaults['printcorrespondent']=False
        defaults['printcorrespasfamily']=False
        defaults['cbarhshift']=-0.15
        defaults['printcorrespascubicfamily']=False
        defaults['printcorrespasfamily']=False
        defaults['printcorrespondentpoints']=False
        defaults['varsel']=0
        defaults['stereomesh']=True
        defaults['R2Proj']=np.eye(3)
        defaults['ProjType'] = 'equalarea'
        defaults['cmap']='jet'
        defaults['cbarhshift']=-0.15
        defaults['norms']=[]
        self.plotAttributes={}
        for phase in self.phases.keys():
            self.plotAttributes[phase] = {}
            for key in defaults.keys():
                self.plotAttributes[phase][key]=defaults[key]
            x=self.phases[phase]['L'].dot([1,0,0])
            x=x/sqrt(x.dot(x))
            z=self.phases[phase]['Lr'].dot([0,0,1])
            z=z/sqrt(z.dot(z))
            y=np.cross(z,x)
            self.plotAttributes[phase]['R2Proj']=np.array([x,y,z])
            self.plotAttributes[phase]['phase1'] = phase
            for title in 'L Lr'.split():
                self.plotAttributes[phase][f'{title}Phase1'] = self.phases[phase][title]
            for title in 'symops recsymops'.split():
                self.plotAttributes[phase][f'{title}'] = self.phases[phase][title]
            if phase == self.austenite:
                phase2  = self.martensite
                Tkey = 'T_MA'
                ORkeys = 'CId CIp'
            else:
                phase2  = self.austenite
                Tkey = 'T_AM'
                ORkeys = 'Cd Cp'
            self.plotAttributes[phase]['phase2'] = phase2
            self.plotAttributes[phase]['T'] = self.defGrad['NiTi']['A'][Tkey]
            for title in 'L Lr'.split():
                self.plotAttributes[phase][f'{title}Phase2'] = self.phases[phase2][title]
            for key,title in zip('Cd Cp'.split(),ORkeys.split()):
                self.plotAttributes[phase][f'{key}'] = self.OR[name][title]
            system=self.phases[phase]['cif'].lattice_system.name
            if system == 'cubic':
                self.plotAttributes[phase]['sphere']='triangle'
                self.plotAttributes[phase]['dirs'] = [[0,0,1],[1,0,1],[1,1,1]]
                self.plotAttributes[phase]['cut2triangle'] = True
                dirtexthiftsi={}
                dirtexthiftsi['223']=[-0.08,0]
                dirtexthiftsi['112']=[-0.08,0]
                dirtexthiftsi['212']=[0.03,0]
                normtexthiftsi={}
            elif system == 'monoclinic':
                self.plotAttributes[phase]['sphere']='half'
                self.plotAttributes[phase]['dirs'] = [[1,0,0],[-1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,1,0],[-1,1,0],[1,1,1],[4,0,1],
                     [-1,1,1],[-1,0,2],[1,0,2],[-1,0,1],[1,0,1],[2,1,1],[-2,1,1],[0,1,3]]
                uvws = self.generate_uvws(1,self.martensite)
                uvws = uvws[3][(uvws[3][:,1]>=0)*(uvws[3][:,2]>=0),:]
                self.plotAttributes[phase]['dirs'] = uvws
                self.plotAttributes[phase]['norms'] = [[0,0,1]]
                self.plotAttributes[phase]['cut2half'] = True
                dirtexthiftsi={}
                normtexthiftsi={}
                dirtexthiftsi['-110']=[-0.3,0]
                dirtexthiftsi['-111']=[-0.3,0.]
                dirtexthiftsi['-101']=[-0.25,0.]
                dirtexthiftsi['-102']=[-0.1,0.0]
                dirtexthiftsi['001']=[-0.02,-0.0]
                dirtexthiftsi['101']=[0.05,0.0]
                dirtexthiftsi['-100']=[-0.2,0.0]
                dirtexthiftsi['-211']=[-0.35,0.05]
                dirtexthiftsi['011']=[-0.05,0.0]
                dirtexthiftsi['111']=[-0.1,0.0]
                dirtexthiftsi['211']=[-0.1,0.0]
                dirtexthiftsi['013']=[-0.35,0.25]
                dirtexthiftsi['401']=[-0.05,0.]

                normtexthiftsi['001']=[-0.02,0.2]
            else:
                dirtexthiftsi={}
                normtexthiftsi={}
                self.plotAttributes[phase]['sphere']='full'
                self.plotAttributes[phase]['dirs'] = [[0,0,1],[1,0,1],[1,1,1]]

            self.plotAttributes[phase]['dirtexthifts']=dirtexthiftsi
            self.plotAttributes[phase]['normtexthifts']=normtexthiftsi

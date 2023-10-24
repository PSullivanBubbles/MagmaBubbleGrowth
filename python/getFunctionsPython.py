#include <cmath>
#include <iostream>
#include <vector>
#include <valarray>
#include <string>

""""
/*function [SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
    PTt_fun] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
    EOSModel, PTtModel)
%This script returns the functions that are used by the bubble growth
% model (Numerical_Model.m)
% 

%Below are the switch-case statements that 'point' to the desired functions
%as defined by the user in the Bubble_Growth_Modelling.m script
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
*/
"""
import math
import scipy
import scipy.optimize
import numpy
import numpy.matlib

#quick def of sqrt
def sqrt(num):
    return pow(num,0.5)

#Return the target Solubility 

def SolFun(T, P, SolModel= "Liu 2005"):
    
    if(SolModel=="Ryan 2015"):
        # Regressed solubility function from Ryan et al. (2015)
        H2Oeq = numpy.array(92.3/T + 0.0287)
        return H2Oeq
        
    elif(SolModel == "Liu 2005"):
        #Solubility function from Liu et al. (2005)
        # convert pressure to MPa from Pa
        P = P*1e-6
        H2Oeq =numpy.array((354.94*(P**0.5)+9.623*P-1.5223*(P**1.5))/T + 0.0012439*(P**1.5))
        return H2Oeq
    
    print("ERROR: solubility model not defined \n")
    return 0*T

#Return the target Diffusivity 

def DiffFun(T,P,H2Ot,W,DiffModel= "Zhang 2010 Peralkaline"):

    if (DiffModel== "Zhang 2010 Metaluminous simple"):
        #Diffusivity function for metaluminous rhyolite from Zhang and Ni (2010)
        #(Equation 15 therein)
        P=P/1e9 #Convert pressure to GPa from Pa
        D = H2Ot*math.exp( -18.1 +(1.888*P)-((9699+3626*P)/T))
        return D
    
    elif (DiffModel== "Zhang 2010 Metaluminous"):
        # Diffusivity function for metaluminous rhyolite from Zhang and Ni (2010)
        #(Equations 7a, 13, 14 therein)
        P=P/1e9 # Convert pressure to GPa from Pa
        # convert from wt.% to mole fraction
        XH2Ot=(H2Ot/18.015)/((H2Ot/18.015)+(100-H2Ot)/W)
        # Compute the diffusivity of molecular water (Zhang & Ni 2010, equation 14)
        DH2Om = numpy.exp(-14.26 + 1.888*P - 37.26*XH2Ot - ((12939 + 3626*P - 75884*XH2Ot)/T))
        #Get the speciation constant (Zhang & Ni 2010, equation 7a)
        K = numpy.exp(1.876 - 3110/T)
        # Compute the total water diffusivity from Diffusivity of
        # molecular water and speciation constant (Zhang & Ni 2010, equation 13)
        D = DH2Om*(1-(0.5-XH2Ot)/sqrt(((4/K-1)*(XH2Ot-(XH2Ot*XH2Ot))+0.25)))
        return D
    
    elif (DiffModel== "Zhang 2010 Peralkaline"):
        # Diffusivity function for peralkaline rhyolite from Zhang and Ni (2010)
        # (Equations 7a, 13, 16 therein)
        P=P/1e9 # Convert from Pa to GPa
        # convert from wt.% to mole fraction
        XH2Ot=(H2Ot/18.015)/((H2Ot/18.015)+(100-H2Ot)/W) 
        # Compute the diffusivity of molecular water (Zhang & Ni 2010,
        #equation 16)
        DH2Om = numpy.exp(-12.79 - 27.87*XH2Ot - ((13939 + 1230*P - 60559*XH2Ot)/T))
        # Get the speciation constant (Zhang & Ni 2010, equation 7a)
        K = numpy.exp(1.876 - 3110/T)
        D = DH2Om*(1-(0.5-XH2Ot)/sqrt(((4./K-1)*(XH2Ot-(XH2Ot*XH2Ot))+0.25)))
        return D
      
    elif (DiffModel== "Constant"):
        # Diffusivity function for a constant diffusivty
        DiffVal = 1e-11 # Diffusivity value in m**2 s**-1
        D =DiffVal*numpy.ones(len(T))
        return D
    
#Return the target Viscosity function

def ViscFun(H2Ot,T,Composition,ViscModel="Giordano 2008"):

    if (ViscModel=="Giordano 2008"):
        outvalues=numpy.array(Giordano2008_Model(H2Ot, Composition))
        At = outvalues[:,0]
        Bt = outvalues[:,1]
        Ct = outvalues[:,2]
        V = 10**(At+(Bt/(T-Ct)))
        
        
    elif (ViscModel=="Hess and Dingwell 1996"):
        V = 10**(-3.545+0.833*numpy.log(H2Ot)+(9601-2368*numpy.log(H2Ot))/(T-(195.7+32.25*numpy.log(H2Ot))))
            
    elif (ViscModel=="Peralkaline Giordano 2000"):
        V = 10**(-5.9-0.286*numpy.log(H2Ot)+(10775.4-394.8*H2Ot)/(T-148.7+21.65*numpy.log(H2Ot)))
    return V

# Return the target functions for Equation of State

def m0_fun(R,P,T,EOSModel="Pitzer and Sterner"):
    V=(4*numpy.pi/3)*R**3 # Volume of the bubble (M**3)
    if(EOSModel=="Ideal Gas Law"):
        # Moles of water in the bubble (from volume and P)
        # Gas constant = 8.314 J mol**-1K-1
        # 1 Pa*m**3 = 1 joule
        n=P*V/(8.314*T)
        # 18.015 g/mol * n moles = mass (g)* 1kg / 1000g
        m0=18.015*n/1000
    elif(EOSModel=="Pitzer and Sterner"):
        # calculate initial mass of water in the bubble (kg)
        # Returns the density of the bubble in kg/m**3
        rho = density(P,T,coefficients())
        m0 = rho*V
    return m0

def pb_fun(m,T,R,EOSModel="Pitzer and Sterner"):
    Vbub = (4/3)*numpy.pi*R **3
    if(EOSModel=="Ideal Gas Law"):
        pb=(((m*1000)/18.015)*8314*T)/(Vbub)
    elif(EOSModel=="Pitzer and Sterner"):
                #%Rho must be in terms of mol/cm**3 for the pitzer & sterner equation of
        #%state. Therefore Volume: m**3 * (1e6 cm**3 / 1 m **3.
        #%and Mass: (kg * (1000g / 1 kg))*(1 mole / 18.015 g)
        rho = ((m*1000)/18.015)/(Vbub*1e6)
        
        #%Get the pitzer and sterner coefficients
        b = coefficients()
        a=numpy.zeros(10)
        for i in range(10):
            a[i]=b[i,0]*T**-4 + b[i,1]*T**-2 + b[i,2]*T**-1 + b[i,3] + b[i,4]*T + b[i,5]*T**2
        

        #%P [bars], R [cm**3*bar/K/mol] and T [K]

        #%NOTE Gas constant is in terms of [cm**3*bar/K/mol] as the equation of PS,
        #%94 has pressure in terms of bar and density in terms of mol cm**-3
        pb = (rho+a[0]*rho**2-rho**2*((a[2]+2*a[3]*rho+3*a[4]*rho**2+4*a[5]*rho**3)/((a[1]+a[2]*rho+a[3]*rho**2+a[4]*rho**3+a[5]*rho**4)**2))+a[6]*rho**2*numpy.exp(-a[7]*rho)+a[8]*rho**2*numpy.exp(-a[9]*rho))*(83.14472*T)

        #% Convert P from bars (equation P&S) to pascals (model)
        pb = pb/1e-5 # 1 pascal a 1e-5 bars
    
    return pb

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%The below functions are used for solving equation of
#%state via Pitzer & Sterner, 1994
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%The coefficients from Pitzer and Sterner, 1994
def coefficients():
    #% matrix of coefficients for eqs. in Pitzer & Sterner 1994
    b=numpy.zeros((10,6))
    b[0,2]=0.24657688e6 
    b[0,3]=0.51359951e2 
    b[1,2]=0.58638965e0 
    b[1,3]=-0.28646939e-2 
    b[1,4]=0.31375577e-4 
    b[2,2]=-0.62783840e1 
    b[2,3]=0.14791599e-1 
    b[2,4]=0.35779579e-3 
    b[2,5]=0.15432925e-7 
    b[3,3]=-0.42719875e0 
    b[3,4]=-0.16325155e-4 
    b[4,2]=0.56654978e4 
    b[4,3]=-0.16580167e2 
    b[4,4]=0.76560762e-1 
    b[5,3]=0.10917883e0 
    b[6,0]=0.38878656e13 
    b[6,1]=-0.13494878e9 
    b[6,2]=0.30916564e6 
    b[6,3]=0.75591105e1 
    b[7,2]=-0.65537898e5 
    b[7,3]=0.18810675e3 
    b[8,0]=-0.14182435e14 
    b[8,1]=0.18165390e9 
    b[8,2]=-0.19769068e6 
    b[8,3]=-0.23530318e2 
    b[9,2]=0.92093375e5 
    b[9,3]=0.12246777e3 
    return b

#%Pitzer and Sterner, 1994 gas density function, calls PS_myfun
def density(P,T,b):
    # convert P to bars from input (pascals)
    P = P*1e-5 # 1 pascal a 1e-5 bars
    a=numpy.zeros(10)
    for i in range(10):
        a[i]=b[i,0]*T**-4 + b[i,1]*T**-2 + b[i,2]*T**-1 +b[i,3] + b[i,4]*T + b[i,5]*T**2

    # PRT = P/RT where P [bars], R [cm**3*bar/K/mol] and T [K]
    PRT = P/(83.14472*T)
    # solve implicit equation for rho and convert to kg/m**3    
    PS_myfun = PS_myfunLam(a,PRT)
    rho = (scipy.optimize.root(PS_myfun,0.001,args=(), method='hybr')).x*18.01528*1000
    return rho

# the function from Pitzer & Sterner 1994, which takes the matrix of
# coefficients a and P/RT as arguments; rho is a first guess for the
# density [g/mol]
def PS_myfunLam(a,PRT):
    return lambda rho : (rho+a[0]*rho**2-rho**2*((a[2]+2*a[3]*rho+3*a[4]*rho**2+4*a[5]*rho**3)/((a[1]+a[2]*rho+a[3]*rho**2+a[4]*rho**3+a[5]*rho**4)**2))+a[6]*rho**2*math.exp(-a[7]*rho)+a[8]*rho**2*math.exp(-a[9]*rho)) - PRT
    

def PtFun(Pin,tin,t):
    if (t==0): return Pin[0]
    if (t>max(tin)): return Pin[-1]
    for tau in range(len(tin)):
        if (t>tin[tau]): pass
        else:
            lever = (t-tin[tau-1])/(tin[tau]-tin[tau-1])
            a = Pin[tau-1]
            b=(Pin[tau]-Pin[tau-1])
            outPres = a+lever*b
            return outPres      

    return outPres        

def TtFun(Tin,tin,t):
    ## Need to add option for Newton cooling later
    if (t==0): return Tin[0]
    if (t>max(tin)): return Tin[-1]
    for tau in range(len(tin)):
        if (t>tin[tau]): 
            pass
        else:
            lever = (t-tin[tau-1])/(tin[tau]-tin[tau-1])
            a = Tin[tau-1]
            b=(Tin[tau]-Tin[tau-1])
            outPres = a+lever*b
            return outPres   

    return outPres        


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%The below functions are used for solving viscosity
#%using Giordano 2008
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def Giordano2008_Model(H2Ot, Composition):
    """
    %This is code is from Girodano 2008 (see citation below) modified
    %for use by the Coumans et al., 2020 bubble growth model

    % SCRIPT grdmodel08 Nov 2007
    %
    % MATLAB script to compute silicate melt viscosity.
    %
    % Citation: Giordano D, Russell JK, & Dingwell DB (2008) 
    %  Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science
    %  Letters, v. 271, 123-134.
    %
    % ________________________________________________________________
    %INPUT: Chemical compositions of silicate melts (Wt. % oxides) as:
    % SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1
    % One line for each melt composition
    % _________________________________________________________________ 

    % ________________________________________________________________
    % OUTPUT: VFT Model parameters A, B and C for each melt composition
    %         to model temperature dependence: log n = A + B/[T(K) - C]
    % VFT_model.out contains: No. , A_value , B_value , C_value , Tg , F
    % VFT_curves.out contains: Temperature Array in Kelvins &
    %                          Log(n) values at T(K) values (1 line per melt) 
    % ________________________________________________________________

    % VFT Multicomponent-Model Coefficients
    % -----------------------------------------------------------------
    """
    AT  = -4.55
    bb  = numpy.array([159.56, -173.34, 72.13, 75.69, -38.98, -84.08, 141.54, -2.43, -0.91, 17.62])
    cc  = numpy.array([2.75, 15.72, 8.32, 10.2, -12.29, -99.54, 0.3 ])

    row = numpy.size(H2Ot,0)
    Comp_rep = numpy.matlib.repmat(Composition,row,1)
    for i in range(row):
        Comp_rep[i,10] = H2Ot[i]

    wtm=Comp_rep

    nx=numpy.size(wtm,0)
    nc=numpy.size(wtm,1)

    bb_rep = numpy.matlib.repmat(bb,nx,1)
    cc_rep = numpy.matlib.repmat(cc,nx,1)
    AT_rep = numpy.matlib.repmat(AT,nx,1)

    #% Function molefrac_grd: converts wt % oxide bais to mole % oxide basis
    #%[ncomps xmf_t] = molefrac_grd(wtm)
    [ncomps, xmf_t] = molepct_grd(wtm)
                
    #%outvalues=[];
    #% Load composition-basis matrix for multiplication against model-coefficients
    #% Result is two matrices bcf[nx by 10] and ccf[nx by 7]
    siti    =   xmf_t[:,0] + xmf_t[:,1] 
    tial    =   xmf_t[:,1]+xmf_t[:,2] 
    fmm     =   xmf_t[:,3] + xmf_t[:,4] + xmf_t[:,5] 
    nak     =   xmf_t[:,7] + xmf_t[:,8] 
    b1  =   siti
    b2  =   xmf_t[:,2] 
    b3  =   xmf_t[:,3] + xmf_t[:,4] + xmf_t[:,9] 
    b4  =   xmf_t[:,5] 
    b5  =   xmf_t[:,6] 

    #print(1+ xmf_t[1,10])

    b6  =   xmf_t[:,7] + xmf_t[:,10] + xmf_t[:,11] 
    b7  =   xmf_t[:,10] + xmf_t[:,11] + numpy.log(1+xmf_t[:,10]) 
    b12 =   siti*fmm 
    b13 =   (siti + xmf_t[:,2] + xmf_t[:,9])*( nak + xmf_t[:,10] )
    b14 =   xmf_t[:,2]*nak

    c1      =   xmf_t[:,0]
    c2      =   tial
    c3      =   fmm
    c4      =   xmf_t[:,6]
    c5      =   nak
    c6      =   numpy.log(1+xmf_t[:,10] + xmf_t[:,11])
    c11     =   xmf_t[:,2] + fmm + xmf_t[:,6] - xmf_t[:,9]
    c11     =   c11*(nak + xmf_t[:,10] + xmf_t[:,11])
    bcf      =   numpy.transpose(numpy.vstack((b1, b2, b3, b4, b5, b6, b7, b12, b13, b14)))
    ccf      =   numpy.transpose(numpy.vstack((c1, c2, c3, c4, c5, c6, c11)))   
        
    # for iz in range(nx):                      #% step through each composition 
     #   BT          = numpy.sum(bb_rep*bcf[iz,:],1)
      #  CT          = numpy.sum(cc_rep*ccf[iz,:],1)
       # TG          = BT/(12-AT) + CT
        #F           = BT/(TG*(1 - CT/TG)*(1 - CT/TG))
         #outvalues   =[outvalues ; iz AT BT CT TG F];
         #end
         
    
    
    #%Calculate the coefficients using matrix algebra instead of a loop

    BT          = numpy.array(numpy.sum(bb_rep*bcf,1))
    CT          = numpy.array(numpy.sum(cc_rep*ccf,1))
    BT = numpy.reshape(BT,(-1,1))
    CT = numpy.reshape(CT,(-1,1))
    TG          = BT/(12-AT_rep) + CT
    F           = BT/(TG*(1 - CT/TG)*(1 - CT/TG))


    outvalues   =numpy.hstack((AT_rep, BT, CT, TG, F))
    return outvalues

    #%Function for computing the mole percent


#function [nox xmf] = molepct_grd(wtm):
def molepct_grd(wtm):
    nr = numpy.size(wtm,0)
    nc = numpy.size(wtm,1)
    #% [n x]=MOLEPCT(X) Calculates mole percent oxides from wt % 
    #% 1 SiO2 2 TiO2  3 Al2O3  4 FeO  5 MnO  6 MgO  7CaO 8 Na2O  9 K2O  10 P2O5 11 H2O 12 F2O-1  
    #% Output: mole fractions of equivalent
    mw=[60.0843, 79.8658, 101.961276, 71.8444, 70.937449,40.3044,56.0774, 61.97894, 94.1960, 141.9446,18.01528, 18.9984]
    mp=[]
    xmf=[]

    #%Replicate the molecular weight row by the number of spatial nodes used
    #%(number of data rows)
    mw_rep = numpy.matlib.repmat(mw,nr,1)

    #%Perform the calculation using matrix algebra

    wtn = numpy.zeros((nr,nc))
    ZZ=(numpy.sum(wtm[:,0:9],1)+wtm[:,11])
    for i in range(12):
        if (i==10):
            wtn[:,i]=wtm[:,10]
        else:
            wtn[:,i]=(100-wtm[:,10])*wtm[:,i]/(numpy.sum(wtm[:,0:9],1)+wtm[:,11])

    #wtn = [wtm[:,0:9]*(100-wtm[:,10])/(numpy.sum(wtm[:,0:9],1)+wtm[:,11]), wtm[:,10], 0.5*wtm[:,11]*(100-wtm[:,10])/(numpy.sum(wtm[:,0:9],1)+wtm[:,11])]
    mp=wtn/mw_rep
    div=numpy.sum(mp,1)
    mpv= 100*(mp/div[:,None])
    xmf=mpv

    nc=numpy.size(xmf,1)
    nox =nc

    return [nox, xmf]
    

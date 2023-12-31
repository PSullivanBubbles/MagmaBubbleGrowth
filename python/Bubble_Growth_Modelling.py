
import numpy
import math
import scipy.interpolate
import getFunctionsPython as gFP
import Numerical_Model as NM
import matplotlib.pyplot as plt
from shapely.geometry import LineString


import time
start_time = time.time()


#%==========================================================================
#%Function to check if the selected models are in the range of validity, for
#%this current program it is only checking Diffusivity and Solubility
#%==========================================================================

def Checks(SolModel, DiffModel, t_nuc, t_f, T, P, H2Ot_0):
    solFlag = 1
    diffFlag = 1

    T_0 = T[0]
    T_f = T[-1]
    P_0 = P[0]
    P_f = P[-1]

#%Get the selected functions
#[SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
#    PTt_fun] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
#    EOSModel, PTtModel)

#%Check time parameters for errors 
    #if (len(t_nuc)==0): print('Enter a nucleation time (s), (e.g., t_0 = 0)')
    #if (len(t_f)==0): print('Enter a final time (s), (e.g., t_f = 2000)')
    if (t_f <= t_nuc): print('The final time is less than or equal to the initial time')

#%Check water parameter for errors
    #if (len(H2Ot_0)==0): print('Enter a water concentration, e.g., H2Ot_0 = 1.5')
    if (H2Ot_0 <= 0): print('Enter a positive water concentration')

    """
    NOT NEEDED with new PTt pathways

    #%If the final temperatures and pressures are empty, if so set them equal to
    #%the initial temperature and pressure
    if (print(T_f)==0), T_f = T_0; end
    if (size(P_f)==0), P_f = P_0; end
    %If the rates are empty, set equal to zero.
    if (size(dTdt)==0), dTdt = 0; end
    if (size(dPdt)==0), dPdt = 0; end


    %Check if the temperature and pressure rates are the right direction, if
    %not raise an error
    if ((T_f-T_0)/dTdt < 0), error('Temperature rate is the wrong direction, check signs'),end
    if ((P_f-P_0)/dPdt < 0), error('Pressure rate is the wrong direction, check signs'),end

    %Check if there are rates defined, if not raise an error
    if (abs((T_f-T_0)) > 0) && (dTdt == 0), error("A temperature range is defined, but a rate isn't"),end
    if (abs((P_f-P_0)) > 0) && (dPdt == 0), error("A Pressure range is defined, but a rate isn't"),end
    """
    print(' ')
    print('%%%Model warnings below%%%')
    print(' ')


#%Switch-cases to get the parameter bounds for the solubility models
    if (SolModel== 'Liu 2005'):
        Sol_P_min = 0.1*1e6# %Pa
        Sol_P_max = 500*1e6# %Pa
        Sol_T_min = 700 + 273.15# %Kelvin
        Sol_T_max = 1200 + 273.15# %kelvin
        
    elif (SolModel=='Ryan 2015'):
        Sol_P_min = 0.1*1e6# %Pa
        Sol_P_max = 0.1*1e6# %Pa
        Sol_T_min = 900 + 273.15# %Kelvin
        Sol_T_max = 1100 + 273.15# %kelvin
        
    else:
        print('Bounds for solubility model not defined')
        solFlag = 0# %flag to ignore bound check


#%Switch-cases to get the parameter bounds for the diffusion models
    if( DiffModel=='Zhang 2010 Metaluminous simple'): # %Zhang and Ni (2010) equation 15
        Diff_P_min = 0.1*1e6#; %Pa
        Diff_P_max = 1900*1e6#; %Pa
        Diff_T_min = 500 + 273.15#; %Kelvin
        Diff_T_max = 1200 + 273.15#; %kelvin
        Diff_H2O_min = 0#; %wt. %
        Diff_H2O_max = 2#; %wt. %
        
        
    elif( DiffModel== 'Zhang 2010 Metaluminous'):# %Zhang and Ni (2010) equations 7a, 13, 14
        Diff_P_min = 0.1*1e6#; %Pa
        Diff_P_max = 1900*1e6#; %Pa
        Diff_T_min = 676#; %Kelvin
        Diff_T_max = 1900#; %kelvin
        Diff_H2O_min = 0#; %wt. %
        Diff_H2O_max = 8#; %wt. %
        
    elif( DiffModel=='Zhang 2010 Peralkaline'):# %Zhang and Ni (2010) equation 7a, 13, 16
        Diff_P_min = 0.1*1e6#; %Pa
        Diff_P_max = 1900*1e6#; %Pa
        Diff_T_min = 676#; %Kelvin
        Diff_T_max = 1900#; %kelvin
        Diff_H2O_min = 0#; %wt. %
        Diff_H2O_max = 8#; %wt. %
        
    else:
        print('Bounds for diffusion model not defined')
        diffFlag = 0# %flag to ignore bound check

    if (solFlag ==1):
    #%Check if conditions are appropriate for solubility model, if outside of
    #%range raise a warning
        if (P_0 < Sol_P_min): print('Initial pressure below solubility model range')
        if (P_f > Sol_P_max): print('Final pressure above solubility model range')
        if (T_0 < Sol_T_min): print('Initial temperature below solubility model range')
        if (T_f > Sol_T_max): print('Final temperature above solubility model range')


    if (diffFlag==1):
        #%Check if conditions are appropriate for the diffusion model, if outside of
        #%range raise a warning
        if (P_0 < Diff_P_min): print('Initial pressure below diffusion model range')
        if (P_f > Diff_P_max): print('Final pressure above diffusion model range')
        if (T_0 < Diff_T_min): print('Initial temperature below diffusion model range')
        if (T_f > Diff_T_max): print('Final temperature above diffusion model range')
        H2O = [gFP.SolFun(T_f,P_0,SolModel),gFP.SolFun(T_f,P_f,SolModel),gFP.SolFun(T_0,P_0), gFP.SolFun(T_0,P_f,SolModel), H2Ot_0]
        H2O_min = min(H2O)
        H2O_max = max(H2O)
        if (H2O_min < Diff_H2O_min): print('H2O will potentially be below diffusion model range')
        if (H2O_max > Diff_H2O_max): print('H2O will potentially be above diffusion model range')
    print(' ')
    print('%%%Warnings or errors from Numerical Model%%%')
    print(' ')
    return


#%==========================================================================
#%Function to explore solubility conditions of the model inputs. The user
#%can use the visuals of this function to determine the nucleation time
#%delay for a given initial water or the solubility of the initial
#%conditions..etc.
#%==========================================================================

def SolubilityExplore(SolModel, T_0,t_T,P_0,t_P, R_0, SurfTens, H2Ot_0):

#%Get the selected functions
#[SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
#    PTt_fun] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
#    EOSModel, PTtModel);

    #%Get range of solubilities following P-T-t Pathway
    t = numpy.linspace(0,t_T[-1],1000)
    #PT = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t');
    P = gFP.PtFun(P_0,t_P,t)
    P = gFP.PtFun(T_0,t_T,t)

#%Determine the intersections between the P-T-t solubility and initial
#%water, this determines the nucleation delay. The function called below is
#%from the MATLAB file exchange authored by Douglas M. Schwarz, see the
#%function for more details.


#%Solubility of water at the bubble wall including Capilary pressure
    H2Oeq = gFP.SolFun(T, P + (2*SurfTens/R_0),SolModel)
    first_line = LineString(numpy.column_stack((t, numpy.ones(len(t))*H2Ot_0)))
    second_line = LineString(numpy.column_stack((t, H2Oeq)))

    [x0,y0] =  first_line.intersection(second_line)

    #%Plot
    plt.figure(1)
    #plt.hold
    plt.plot(t,H2Oeq, 'k-', 'linewidth', 2)
    plt.plot(x0, y0,'marker','p','markerfacecolor','r','markeredgecolor','r','markersize',20)
    plt.yline(H2Ot_0, 'b-', 'linewidth', 2)
    plt.grid(1)
    plt.box(1)

    plt.ylabel('Solubility H_2O_t, (wt. %)')
    plt.xlabel ('Time, \itt\rm (seconds)')

    #%Create a legend
    h=[]
    h[1] = plt.plot(numpy.NaN,numpy.NaN,'linestyle','-','color','k','linewidth',2)
    h[2] = plt.plot(numpy.NaN,numpy.NaN,'linestyle','-','color','b','linewidth',2)
    h[3] = plt.plot(numpy.NaN,numpy.NaN,'pr','markerfacecolor','r','markeredgecolor','r','markersize',10)

    lbl = {['Model solubility'],['Inital H_2O_t, (wt. %)'],['Saturation']}
    #%Set positions of the Legend
    plt.legend(h,lbl,'Location','northeast','FontName','Times New Roman','FontSize',10)


#%==========================================================================
#%Function to plot the output.
#%==========================================================================
def plotting(t, R, phi, P, T, x_out, H2Ot_all):

    #%Plot the results of the numerical model
    plt.figure(1)
    plt.subplot(2,2,1)
    #plt.hold on
    plt.plot(t, R*1e6)
    plt.grid(1)
    plt.box(1)
    plt.ylabel('Bubble radius, (micron)')
    plt.xlabel ('Time, (seconds)')

    plt.subplot(2,2,2)
    #plt.hold on
    plt.plot(t, phi)
    plt.grid(1)
    plt.box(1)
    plt.ylabel('Gas volume fraction, \phi')
    plt.xlabel ('Time,(seconds)')

    plt.subplot(2,2,3)
    #plt.hold on
    plt.plot(t, P*1e-6)
    plt.ylabel('Pressure, (MPa)')
    plt.xlabel ('Time, (seconds)')
    plt.grid(1)
    plt.box(1)

    plt.subplot(2,2,4)
    plt.plot(t, T)
    plt.ylabel('Temperature, (Kelvin)')
    plt.xlabel ('Time, (seconds)')
    plt.grid(1)
    plt.box(1)

    #plt.show()

    plt.figure(2)
    ClrMap = plt.jet(); #%Set colourmap to jet
    #%Set the colour values based on time
    colour_values = numpy.linspace(t[0],t[-1])
    #plt.hold on
    for i in range(numpy.size(x_out,1)):
        #Clr = [scipy.interpolate.interp1d(colour_values,ClrMap[:,0],t[i]), scipy.interpolate.interp1d(colour_values,ClrMap[:,1],t[i]), scipy.interpolate.interp1d(colour_values,ClrMap[:,2],t[i])]
        plt.plot(x_out[0:-2,i]*1e6,H2Ot_all[0:-2,i])

    Cmap = ClrMap
    #plt.colormap(plt.gca,Cmap)
    hcb = plt.colorbar
    #plt.caxis([colour_values[0]],colour_values[-1])
    #colorTitleHandle = plt.get(hcb,'Title')
    #titleString = {'time, \itt\rm (s)'}
    #set(colorTitleHandle ,'String',titleString,'FontName','Times New Roman','FontSize',12)
    #set(hcb,'YTick',colour_values[0]:colour_values[-1]/5:colour_values[-1])

    plt.grid(1)
    plt.box(1)

    plt.xlabel('Distance from bubble, (micron)')
    plt.ylabel('H_2O_t, (wt. %)')
    plt.show()




#%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
#%Choose which operation: Solubility explore or Model run

#Operation = 'Solubility Explore';
Operation = 'Run Model'


"""
%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
%Set the model properties
%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
"""

#%Melt composition
wt_dry = 0.01; #H2O wt. % that defines "dry"
#[SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1]
PhonoliteVAD79 = [57.15, 0.3, 21.34, 2.70, 0.14, 0.39, 3.26, 5.16, 9.46, 0.09, wt_dry, 0]
Krafla = [75.17, 0.22, 12.02, 3.13, 0.11, 0.09, 1.66, 4.58, 2.88, 0, wt_dry, 0]
PCD = [75.64, 0.08, 12.38, 0.94, 0.07, 0.03, 0.54, 4.17, 4.72, 0.01, wt_dry, 0]
ROGD = [76.51, 0.03, 12.56, 0.7, 0.07, 0.01, 0.25, 4.47, 4.24, 0, wt_dry, 0]

#%================================!!!
#%~~~Composition~~~
Composition = Krafla# %Tuffen and Castro (2009) locality AO
#%Composition = PhonoliteVAD79; %Iacono-Marziano et al. (2007)
#%Composition = PCD; %Mangan and Sisson (2000)
#%Composition = ROGD; %Mourtada-Bonnefoi and Laporte (2004)
#%================================!!!
#%================================!!!
#%~~~solubility model~~~
#%SolModel = 'Ryan 2015'# %Ryan et al. (2015)
SolModel = 'Liu 2005'# %Liu et al. (2005)
#%================================!!!
#%================================!!!
#%~~~Diffusion Model~~~
#DiffModel = 'Zhang 2010 Metaluminous simple'# %Zhang and Ni (2010) equation 15
DiffModel = 'Zhang 2010 Metaluminous'# %Zhang and Ni (2010) equations 7a, 13, 14
#DiffModel = 'Zhang 2010 Peralkaline'# %Zhang and Ni (2010) equation 7a, 13, 16
#DiffModel = 'Constant'# %A constant value that can be defined in getFunctions.m
#================================!!!
#%================================!!!
#%~~~Viscosity Model~~~
#ViscModel = 'Hess and Dingwell 1996'# %Hess and Dingwell (1996)
#ViscModel = 'Peralkaline Giordano 2000'# %Giordano et al. (2000)
ViscModel = 'Giordano 2008'# %Giordano et al. (2008)
#%================================!!!
#%================================!!!
#%~~~EOS model~~~
EOSModel = 'Pitzer and Sterner'# %Pitzer and Sterner (1994)
#EOSModel = 'Ideal Gas Law'
#%================================!!!
#%================================!!!


#%~~~P-T-t Profile~~~
#PTtModel = 'P: Isobaric, T: Isothermal'
#PTtModel = 'P: Isobaric, T: Polythermal-Dwell'
#%PTtModel = 'P: Polybaric-Dwell, T: Isothermal-quench'
#%================================!!!



#Constants used for calculations
SurfTens = 0.22; #Value for surface tension (N/m)
melt_Rho = 2350; #Melt density in kg/m^3

#Spatial parameters
Nb = 1e11; #Bubble number density (number per m^3)
Phi_0 = 1*1e-6; #Initial gas volume fraction (Phi)
#R_0 = Radius(Nb,Phi_0)# %R_0 (m) calculated from Phi and Nb
R_0 = 3e-6; #R_0 (m) set independently

#Finite difference parameters
Nodes = 1000; #Number of nodes in spatial discretization

#%Numerical tolerance:
#%[Absolute tolerance, relative tolerance], see:
#% https://www.mathworks.com/help/simbio/ref/absolutetolerance.html
#% https://www.mathworks.com/help/simulink/gui/relative-tolerance.html
#%For additional information
Numerical_Tolerance = [1e-5, 1e-5]



t_nuc = 0
t_f = 3000

T_0 = [1100+ 273.15, 1200 + 273.15 , 1100 + 273.15] 
t_T = [0, 1000, 2000]

P_0 = [ 20e6, 10e6, 20e6, 10e6]
t_P = [0, 500, 1000, 1500]

H2Ot_0 = 1.4

#%%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
#%Run the check function, solubility explore or numerical model
#%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!

#%Checks input parameters to detect if any are out of component model
#%calibration (from literature), or if there are errors in parameter
#%signs (e.g., negative temperature rate for positive temperature change)
Checks(SolModel, DiffModel, t_nuc, t_f,T_0, P_0, H2Ot_0)

#%The switch which either runs the parameter exploration or the model
if (Operation=='Solubility Explore'):
        #%Explore the solubility for the P-T-t pathway
    SolubilityExplore(SolModel, T_0,t_T,P_0,t_P, R_0, SurfTens, H2Ot_0)
        
elif (Operation=='Run Model'):
        
    #%Run the numerical model with the set parameters
    #tic
    [t, R, phi, P, T, x_out, H2Ot_all] =  NM.Numerical_Model_v2(Composition, SolModel, DiffModel, ViscModel, EOSModel, SurfTens, melt_Rho, Nodes, R_0, H2Ot_0, Nb, t_nuc, T_0, t_T, P_0, t_P, Numerical_Tolerance)
    #toc
    print("--- %s seconds ---" % (time.time() - start_time))
    #%Plot the output of the numerical model
    plotting(t, R, phi, P, T, x_out, H2Ot_all)
        
#end %End operation switch



#%
#%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
#%Additional functions used by this script
#%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!




def  Radius(Nb,Phi):
    meltvolume = 1./Nb
    gasvolume = Phi*meltvolume/(1 - Phi)
    return(gasvolume/((4/3)*math.pi()))**(1/3)

#function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)


def  mvmin(x):
#% Faster implementation of movmin(x,k) when k = 1.
    return min(x[0:-2],x[1:-1])
def mvmax(x):
#% Faster implementation of movmax(x,k) when k = 1.
    return max(x[0:-2],x[1:-1])


"""function [t, R, phi, P, T, x_out, H2Ot_all] =  Numerical_Model_v2(Composition, SolModel, DiffModel, ViscModel, EOSModel,...
    PTtModel, SurfTens, melt_Rho, Nodes, R_0, H2Ot_0, Nb, t_nuc, t_f, ...
    T_0, T_f, dTdt, P_0, P_f, dPdt, Numerical_Tolerance)
%This script defines the algorithm that solves the bubble-growth
%numerical experiment outlined in Bubble_Growth_Modelling.m
% 
% 
% See the user manual for more info.
% 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% If using or adapting this code, please cite Coumans et al., (2020):
% 
% Authors: J. P. Coumans, E. W. Llewellin*, F. B. Wadsworth
% M. C. S Humphreys, S. A.  Mathias, B. M. Yelverton, and J. E. Gardner
% 
% Title: An experimentally validated numerical model 
% for bubble growth in magma 
% 
% Journal: Journal of Volcanology and Geothermal Research (JVGR)
% 
% Year: 2020
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""
import numpy
import scipy.integrate
import math
import getFunctionsPython as gFP

def Numerical_Model_v2(Composition, SolModel, DiffModel, ViscModel, EOSModel, SurfTens, melt_Rho, Nodes, R_0, H2Ot_0, Nb, t_nuc, T_0, t_T, P_0, t_P, Numerical_Tolerance):

#%Get the user selected functions used in the numerical model
#[SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
#    PTt_fun] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
#    EOSModel, PTtModel);

#%Get the molar mass of anhydrous  melt on a single oxygen basis
    W = Mass_SingleOxygen(Composition)

#%Get the initial porisity of the shell model
    phi_0 = Porosity(Nb,R_0)

# distance to shell edge (m) (computed from R and phi)
    L = R_0*(phi_0**(-1/3) - 1); 

#Set the numerical duration bounds (from nucleation to end)
    t_f = 3000
    tspan= (t_nuc, t_f)

# Get the P and T at the defined t_0 using the P-T-t pathway
    #PT_t0 = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t_nuc)
    P_t0 = gFP.PtFun(P_0,t_P,t_nuc)
    #P_t0 = P[0]
    T_t0 = gFP.TtFun(T_0,t_T,t_nuc)
    #T_t0 = T[0]

#%Compute the initial mass of gas in the bubble
    m_0 = gFP.m0_fun(R_0, P_t0 +(2*SurfTens/R_0), T_t0, EOSModel)

#generate logarithmically spaced xB values
#scale to appropriate length, and shift so that xB
#starts at the bubble wall
#(equation A.6 from manuscript)
    xB=numpy.append( R_0*numpy.ones(1), numpy.logspace(-2,0,Nodes)*L + R_0,0)

#%block-centered grid - use midpoints of xB
#%(equation A.3 from manuscript)
    x=(xB[1:len(xB)]+xB[0:(len(xB)-1)])/2

#%uniform concentration; system in equilibrium Although this can be
#%changed to represent non-uniform conditions. Any initial concentration
#%profile can be set by the user.
    H2Ot = numpy.zeros(Nodes)+H2Ot_0

#%Y_0 is a vector containing quantities to be solved for by ode15s
    Y_0 = numpy.append(H2Ot,R_0*numpy.ones(1),0)

#%Declare the listener function
    #xoverFcn=@(t,X)eventDetection(t,X)

#%Set the options for the ODE solver
    #options = odeset('AbsTol',Numerical_Tolerance(1),'RelTol',Numerical_Tolerance(2),'Events',xoverFcn)

#%Create an anonymous function of the ODE to solve for passing to the
#%solver
    funsolve = lambda t, Y: MYodeFun(t,Y,x,xB,m_0,melt_Rho,T_0,t_T,P_0,t_P,H2Ot_0,R_0,W,SurfTens,SolModel, DiffModel, ViscModel, EOSModel,Composition,Nb)

#%Solve the PDE by method of lines using the ODE solver function ODE15s
    #[t,Y]=ode15s(funsolve,tspan,Y_0,options)
    sol = scipy.integrate.solve_ivp(funsolve,tspan,Y_0,method='BDF')

    t = sol.t
    Y = sol.y

#%Get the outputs
    [R, phi, P, T, x_out, H2Ot_all] = Outputs(Nodes,R_0,L,Y,t,Nb,T_0,t_T,P_0,t_P)


    return [t, R, phi, P, T, x_out, H2Ot_all]

#%==========================================================================
#%ODE function to be solved. See appendix A for an explanation and
#%pseudocode
#%==========================================================================


#function [dYdt,pb] =  MYodeFun(t,Y,x,xB,m_0,melt_Rho,T_0,P_0,H2Ot_0,R_0,W,SurfTens,SolFun,DiffFun,ViscFun,pb_fun,PTt_fun,Composition,Nb,t_f,T_f, dTdt, P_f, dPdt)

def MYodeFun(t,Y,x,xB,m_0,melt_Rho,T_0,t_T,P_0,t_P,H2Ot_0,R_0,W,SurfTens,SolModel, DiffModel, ViscModel, EOSModel,Composition,Nb):

#debugging print time
    #print(t)
    #print('\n')

#%extract individual concentrations
    nx = (len(Y) - 1)
    H2Ot = Y[0:-1]
    R = Y[-1]

#%Get current temperature and pressure
    #PT = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t)
    P = gFP.PtFun(P_0,t_P,t)
    T = gFP.TtFun(T_0,t_T,t)

#%Trapz is a function that uses the trapezoidal rule to numerically
#%integrate. Dim is the dimension to operate on
    dim = 1

#%Get mass of gas (Equation 7 from the main manuscript)
    i1=(1.0/100)*H2Ot_0*x**2
    I1=numpy.trapz((1.0/100)*H2Ot_0*x**2,x)
    I2=numpy.trapz((1.0/100)*H2Ot*x**2,x)
    m=m_0+4*math.pi*melt_Rho*(I1-I2)

#%Compute pressure of the gas  in the bubble from EOS
#%(section 2.3.2 of main manuscript)
    pb = gFP.pb_fun(m, T, R, EOSModel)

#%====boundary conditions====
#%Determine the solubility condition of water in the system based on gas
#%pressure in the bubble
#%(section 2.3.1 of main manuscript)
    H2Oeq = numpy.ones(1)*gFP.SolFun(T,pb,SolModel)

#%get diffusion coefficients at boundaries
#%Creates a vector where the first value is the boundary condition as
#%determined from the equilibrium solubility of the system.
#%(section 2.3.1 of main manuscript)
    H2O_BC= numpy.append(H2Oeq,(H2Ot[0:(len(H2Ot)-1)]+H2Ot[1:len(H2Ot)])/2,0)
    #H2O_BC= H2Oeq.append((H2Ot[1:(len(H2Ot)-1),:]+H2Ot[2:len(H2Ot),:])/2)
    DH2Ot = gFP.DiffFun(T,P,H2O_BC,W,DiffModel)

#%====Solve water diffusion==== (equation A.4 from manuscript)
#%molecular diffusion flux

    JH2Ot = -DH2Ot*numpy.diff(numpy.append(H2Oeq,H2Ot,0),1,0)/numpy.diff(numpy.append(numpy.ones(1)*xB[0],x,0),1,0)
#%Gradient of the diffusion flux.
    dJH2Odx = (1/(x**2))*numpy.diff([numpy.append((((x**3+R**3-R_0**3)**(4/3))/(x**2))*JH2Ot,numpy.zeros(1),0)],1,1)/numpy.diff(xB,1,0)

#%====solve hydrodynamic equation====
#%Compute the viscosity
    v = gFP.ViscFun(H2Ot,T,Composition,ViscModel)

#%Compute integrated viscosity (equation A.5 from manuscript)
    I3=numpy.trapz((v*x**2)/((R**3-R_0**3+x**3)**2),x)

#%Solve Rayleigh-Plesset equation (equation A.6 from manuscript)
    dRdt= ((pb-P-(2*(SurfTens)/R))/(12*R**2))/I3

#%return rhs of ode
    dYdt = numpy.real(numpy.append(-dJH2Odx,dRdt))

    return dYdt
##%==========================================================================
##%Functions to get outputs from the ODE solution
##%==========================================================================


#function [R, phi, P, T, x_out, H2Ot_all] = Outputs(Nodes,R_0,L,Y,t,Nb,PTt_fun,P_0, P_f, dPdt,T_0,T_f,dTdt)

def Outputs(Nodes,R_0,L,Y,t,Nb,T_0,t_T,P_0,t_P):

    Y = Y
    #%t = t';

    #%Get the bubble radius and phi
    R = Y[-1,:]
    phi = Porosity(Nb,R)

    #%Get the shell thickness for each R(t)
    L_all =  numpy.array(((R_0+L)**3+R**3-R_0**3)**(1/3)-R)
    #%Get the x position for each output shell thickness
    x_out=numpy.outer(numpy.logspace(-2,0,Nodes),L_all)

    #%Get all of the water profiles
    H2Ot_all = Y[1:len(Y),:]

    #%Get P-T-t history
    #PT = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t)

    P=numpy.zeros(len(t))
    T=numpy.zeros(len(t))
    for i in range(len(t)):
        P[i] = gFP.PtFun(P_0,t_P,t[i])
        T[i] = gFP.TtFun(T_0,t_T,t[i])

    return [R, phi, P, T, x_out, H2Ot_all]

#%==========================================================================
#%Functions that return constants which are used in the modelling
#%==========================================================================

#%This function returns the molar mass of the anhydrous melt on a single
#%oxygen basis
#function W = Mass_SingleOxygen(Composition)
def Mass_SingleOxygen(Composition):
    comp = Composition

    """%Convert composition matrix from Viscosity input to Shishkina format
    %!! SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1 !!
    %Convert To
    %!! SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 Cr2O3] !!"""

    X = numpy.zeros(12)
    X[0] = comp[0]
    X[2] = comp[1]
    X[2] = comp[2]
    X[3] = 0
    X[4] = comp[3]
    X[5] = comp[4]
    X[6] = comp[5]
    X[7] = comp[6]
    X[8] = comp[7]
    X[9] = comp[8]
    X[10] = comp[9]
    X[11] = 0


    Total_Mass = numpy.sum(X)

    #%Molar mass (g/mol) of individual elements
    mSi = 28.0855
    mTi = 47.867
    mAl = 26.981539
    mFe = 55.845
    mMn = 54.938044
    mMg = 24.305
    mCa = 40.078
    mNa = 22.989769
    mK = 39.0983
    mP = 30.973762
    mCr = 51.9961
    mO = 15.999

    #% [SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 Cr2O3]
    #%Molar mass (g/mol) of oxides
    OxideMolarMass = numpy.zeros(12)
    OxideMolarMass[0] =  (mSi+2*mO)
    OxideMolarMass[1] = (mTi+2*mO)
    OxideMolarMass[2] = (2*mAl+3*mO)
    OxideMolarMass[3] = (2*mFe+3*mO)
    OxideMolarMass[4] = (1*mFe+1*mO)
    OxideMolarMass[5] = (1*mMn+1*mO)
    OxideMolarMass[6] = (1*mMg+1*mO)
    OxideMolarMass[7] = (1*mCa+1*mO)
    OxideMolarMass[8] = (2*mNa+1*mO)
    OxideMolarMass[9] = (2*mK+1*mO)
    OxideMolarMass[10] = (2*mP+5*mO)
    OxideMolarMass[11] = (2*mCr+3*mO)

    #%Compute number of moles of element, and Cation Fraction
    numMolesOxygen = [2, 2, 3, 3, 1, 1, 1, 1, 1, 1, 5, 3]
    numMolesElement = [1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2]

    #%Compute the number of moles of each oxide
    Moles_Oxide = X/OxideMolarMass

    #%Compute moles of oxygen by stoichiometry
    Moles_Oxygen = Moles_Oxide*numMolesOxygen

    #%W_melt is the mass of anhydrous melt per mole of oxygen
    return Total_Mass/sum(Moles_Oxygen)

#%Compute the gas volume fraction from Nb and radius
def Porosity(Nb,R):
    meltvolume = 1/Nb
    gasvolume = (4/3)*math.pi*R**3
    test = gasvolume/(gasvolume+meltvolume)
    return gasvolume/(gasvolume+meltvolume)
"""
#%==========================================================================
#%ODE solver listener functions
#%==========================================================================
function [value,isterminal,direction] = eventDetection(t,X)
def eventDetection(t,X):
#%a. Define the timeout in seconds
    TimeOut = 1000000
    #% modified to allow for parallelisation
    timeFlag = 100000-TimeOut

    if (X[-1,:] >1e-14):
        sizeFlag = 0
    else:
        sizeFlag=1
        print('bubble too small')



    if ((sizeFlag==1)||(timeFlag<=0))
        value =1
    else
        value=0
    end


    isterminal = 1
    direction =0
"""
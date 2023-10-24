function [t, R, phi, P, T, x_out, H2Ot_all] =  Numerical_Model_v2(Composition, SolModel, DiffModel, ViscModel, EOSModel,...
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

%Get the user selected functions used in the numerical model
[SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
    PTt_fun] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
    EOSModel, PTtModel);

%Get the molar mass of anhydrous  melt on a single oxygen basis
W = Mass_SingleOxygen(Composition);

%Get the initial porisity of the shell model
phi_0 = Porosity (Nb,R_0);

%distance to shell edge (m) (computed from R and phi)
L = R_0*(phi_0^(-1/3) - 1); 

%Set the numerical duration bounds (from nucleation to end)
tspan= [t_nuc t_f];

%Get the P and T at the defined t_0 using the P-T-t pathway
PT_t0 = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t_nuc);
P_t0 = PT_t0(:,1);
T_t0 = PT_t0(:,2);

%Compute the initial mass of gas in the bubble
m_0 = m0_fun(R_0, P_t0 +(2*SurfTens/R_0), T_t0);

%generate logarithmically spaced xB values
%scale to appropriate length, and shift so that xB
%starts at the bubble wall
%(equation A.6 from manuscript)
xB=[0 logspace(-2,0,Nodes)*L]' + R_0;

%block-centered grid - use midpoints of xB
%(equation A.3 from manuscript)
x=(xB(2:end,1)+xB(1:end-1,1))/2;

%uniform concentration; system in equilibrium Although this can be
%changed to represent non-uniform conditions. Any initial concentration
%profile can be set by the user.
H2Ot = zeros(Nodes,1)+H2Ot_0;

%Y_0 is a vector containing quantities to be solved for by ode15s
Y_0 = [H2Ot;R_0];

%Declare the listener function
xoverFcn=@(t,X)eventDetection(t,X);

%Set the options for the ODE solver
options = odeset('AbsTol',Numerical_Tolerance(1),'RelTol',Numerical_Tolerance(2),'Events',xoverFcn);

%Create an anonymous function of the ODE to solve for passing to the
%solver
funsolve = @(t,Y)MYodeFun(t,Y,x,xB,m_0,melt_Rho,T_0,P_0,H2Ot_0,R_0,W,SurfTens,...
    SolFun,DiffFun,ViscFun,pb_fun,PTt_fun,Composition,Nb,t_f,...
    T_f, dTdt, P_f, dPdt);

%Solve the PDE by method of lines using the ODE solver function ODE15s
[t,Y]=ode15s(funsolve,tspan,Y_0,options);


%Get the outputs
[R, phi, P, T, x_out, H2Ot_all] = Outputs(Nodes,R_0,L,Y,t,Nb,PTt_fun,P_0, P_f, dPdt,T_0,T_f,dTdt);

%==========================================================================
%ODE function to be solved. See appendix A for an explanation and
%pseudocode
%==========================================================================
function [dYdt,pb] =  MYodeFun(t,Y,x,xB,m_0,melt_Rho,T_0,P_0,H2Ot_0,R_0,W,SurfTens,...
    SolFun,DiffFun,ViscFun,pb_fun,PTt_fun,Composition,Nb,t_f,...
    T_f, dTdt, P_f, dPdt)

%extract individual concentrations
nx = (size(Y,1) - 1);
H2Ot = Y(1:nx,:);
R = Y(end,:);

%Get current temperature and pressure
PT = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t);
P = PT(:,1);
T = PT(:,2);

%Trapz is a function that uses the trapezoidal rule to numerically
%integrate. Dim is the dimension to operate on
dim = 1;

%Get mass of gas (Equation 7 from the main manuscript)
I1=trapz(x,(1/100)*H2Ot_0.*x.^2,dim);
I2=trapz(x,(1/100)*H2Ot.*x.^2,dim);
m=m_0+4.*pi.*melt_Rho.*(I1-I2);

%Compute pressure of the gas  in the bubble from EOS
%(section 2.3.2 of main manuscript)
pb = pb_fun(m, T, R);

%====boundary conditions====
%Determine the solubility condition of water in the system based on gas
%pressure in the bubble
%(section 2.3.1 of main manuscript)
H2Oeq = SolFun(T,pb);

%get diffusion coefficients at boundaries
%Creates a vector where the first value is the boundary condition as
%determined from the equilibrium solubility of the system.
%(section 2.3.1 of main manuscript)
H2O_BC= [H2Oeq ; (H2Ot(1:end-1,:)+H2Ot(2:end,:))/2];
DH2Ot = DiffFun(H2O_BC,T,P,W);

%====Solve water diffusion==== (equation A.4 from manuscript)
%molecular diffusion flux
JH2Ot = -DH2Ot.*diff([H2Oeq;H2Ot],1,1)./diff([xB(1,:);x],1,1);
%Gradient of the diffusion flux.
dJH2Odx = (1./(x.^2)).*diff([(((x.^3+R.^3-R_0.^3).^(4/3))./(x.^2)).*JH2Ot;0],1,1)./diff(xB,1,1);

%====solve hydrodynamic equation====
%Compute the viscosity
v = ViscFun(H2Ot,T,Composition);

%Compute integrated viscosity (equation A.5 from manuscript)
I3=trapz(x, (v.*x.^2)./((R.^3-R_0.^3+x.^3).^2),dim);

%Solve Rayleigh-Plesset equation (equation A.6 from manuscript)
dRdt= ((pb-P-(2.*(SurfTens)./R))./(12.*R.^2))./I3;

%return rhs of ode
dYdt = real([-dJH2Odx;dRdt]);

%==========================================================================
%Functions to get outputs from the ODE solution
%==========================================================================
function [R, phi, P, T, x_out, H2Ot_all] = Outputs(Nodes,R_0,L,Y,t,Nb,PTt_fun,P_0, P_f, dPdt,T_0,T_f,dTdt)
Y = Y';
%t = t';

%Get the bubble radius and phi
R = Y(end,:);
phi = Porosity (Nb,R);

%Get the shell thickness for each R(t)
L_all =  ((R_0+L)^3+R.^3-R_0^3).^(1/3)-R;
%Get the x position for each output shell thickness
x_out=(logspace(-2,0,Nodes).*L_all')';

%Get all of the water profiles
H2Ot_all = Y(1:end-1,:);

%Get P-T-t history
PT = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t);
P = PT(:,1);
T = PT(:,2);

%==========================================================================
%Functions that return constants which are used in the modelling
%==========================================================================

%This function returns the molar mass of the anhydrous melt on a single
%oxygen basis
function W = Mass_SingleOxygen(Composition)

comp = Composition';

%Convert composition matrix from Viscosity input to Shishkina format
%!! SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1 !!
%Convert To
%!! SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 Cr2O3] !!

X = zeros(12,1);
X(1) = comp(1);
X(2) = comp(2);
X(3) = comp(3);
X(4) = 0;
X(5) = comp(4);
X(6) = comp(5);
X(7) = comp(6);
X(8) = comp(7);
X(9) = comp(8);
X(10) = comp(9);
X(11) = comp(10);
X(12) = 0;


Total_Mass = sum(X);

%Molar mass (g/mol) of individual elements
mSi = 28.0855;
mTi = 47.867;
mAl = 26.981539;
mFe = 55.845;
mMn = 54.938044;
mMg = 24.305;
mCa = 40.078;
mNa = 22.989769;
mK = 39.0983;
mP = 30.973762;
mCr = 51.9961;
mO = 15.999;

% [SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 Cr2O3]
%Molar mass (g/mol) of oxides
OxideMolarMass = zeros(12,1);
OxideMolarMass(1) =  (mSi+2*mO);
OxideMolarMass(2) = (mTi+2*mO);
OxideMolarMass(3) = (2*mAl+3*mO);
OxideMolarMass(4) = (2*mFe+3*mO);
OxideMolarMass(5) = (1*mFe+1*mO);
OxideMolarMass(6) = (1*mMn+1*mO);
OxideMolarMass(7) = (1*mMg+1*mO);
OxideMolarMass(8) = (1*mCa+1*mO);
OxideMolarMass(9) = (2*mNa+1*mO);
OxideMolarMass(10) = (2*mK+1*mO);
OxideMolarMass(11) = (2*mP+5*mO);
OxideMolarMass(12) = (2*mCr+3*mO);

%Compute number of moles of element, and Cation Fraction
numMolesOxygen = [2 2 3 3 1 1 1 1 1 1 5 3]';
numMolesElement = [1 1 2 2 1 1 1 1 2 2 2 2]';

%Compute the number of moles of each oxide
Moles_Oxide = X./OxideMolarMass;

%Compute moles of oxygen by stoichiometry
Moles_Oxygen = Moles_Oxide.*numMolesOxygen;

%W_melt is the mass of anhydrous melt per mole of oxygen
W = Total_Mass./sum(Moles_Oxygen);

%Compute the gas volume fraction from Nb and radius
function phi = Porosity (Nb,R)
meltvolume = 1/Nb;
gasvolume = (4/3).*pi()*R.^3;
phi = gasvolume./(gasvolume+meltvolume);

%==========================================================================
%ODE solver listener functions
%==========================================================================
function [value,isterminal,direction] = eventDetection(t,X)
%a. Define the timeout in seconds
TimeOut = 1000000;
% modified to allow for parallelisation
timeFlag = 100000-TimeOut;

if X(end,:) >1e-14
    sizeFlag = 0;
else
    sizeFlag=1;
    disp('bubble too small');
end


if ((sizeFlag==1)||(timeFlag<=0))
    value =1;
else
    value=0;
end


isterminal = 1;
direction =0;

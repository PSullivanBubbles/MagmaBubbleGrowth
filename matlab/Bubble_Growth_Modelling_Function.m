function [RScale,tscale]=Bubble_Growth_Modelling_Function(volScale)
%clear; %Clear all variables before running the script
%close all; %Close all figures before running the script

%This script runs the Bubble growth numerical model (Numerical_Model.m)
%using the functions defined in (getFunctions.m). The results are plotted
% as figures.
%
%
% See the user manual for more info. References to individual models
% can be found at the bottom of this script and in the user manual.
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
%

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%====================================================
%Gets the matlab filename
mfilename;
%Gets the folder (local computer) where the script is saved
mpath = strrep(which(mfilename),[mfilename '.m'],'');
mFolder = mpath;
%Adds the folder to the path so that local functions can be called
addpath(mFolder);
%===================================================


%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
%Choose which operation: Solubility explore or Model run
%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
%================================!!!
%Operation = 'Solubility Explore';
 Operation = 'Run Model';
%================================!!!

%
%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
%Set the model properties
%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!

%Melt composition
wt_dry = 0.01; %H2O wt. % that defines "dry"
%[SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1]
PhonoliteVAD79 = [57.15 0.3 21.34 2.70 0.14 0.39 3.26 5.16 9.46 0.09 wt_dry 0];
Krafla = [75.17 0.22 12.02 3.13 0.11 0.09 1.66 4.58 2.88 0 wt_dry 0];
PCD = [75.64 0.08 12.38 0.94 0.07 0.03 0.54 4.17 4.72 0.01 wt_dry 0];
ROGD = [76.51 0.03 12.56 0.7 0.07 0.01 0.25 4.47 4.24 0 wt_dry 0];

%================================!!!
%~~~Composition~~~
Composition = Krafla; %Tuffen and Castro (2009) locality AO
%Composition = PhonoliteVAD79; %Iacono-Marziano et al. (2007)
%Composition = PCD; %Mangan and Sisson (2000)
%Composition = ROGD; %Mourtada-Bonnefoi and Laporte (2004)
%================================!!!
%================================!!!
%~~~solubility model~~~
SolModel = 'Ryan 2015'; %Ryan et al. (2015)
%SolModel = 'Liu 2005'; %Liu et al. (2005)
%================================!!!
%================================!!!
%~~~Diffusion Model~~~
%DiffModel = 'Zhang 2010 Metaluminous simple'; %Zhang and Ni (2010) equation 15
DiffModel = 'Zhang 2010 Metaluminous'; %Zhang and Ni (2010) equations 7a, 13, 14
%DiffModel = 'Zhang 2010 Peralkaline'; %Zhang and Ni (2010) equation 7a, 13, 16
%DiffModel = 'Constant'; %A constant value that can be defined in getFunctions.m
%================================!!!
%================================!!!
%~~~Viscosity Model~~~
%ViscModel = 'Hess and Dingwell 1996'; %Hess and Dingwell (1996)
%ViscModel = 'Peralkaline Giordano 2000'; %Giordano et al. (2000)
ViscModel = 'Giordano 2008'; %Giordano et al. (2008)
%================================!!!
%================================!!!
%~~~EOS model~~~
EOSModel = 'Pitzer and Sterner'; %Pitzer and Sterner (1994)
%EOSModel = 'Ideal Gas Law';
%================================!!!
%================================!!!
%~~~P-T-t Profile~~~
%PTtModel = 'P: Isobaric, T: Isothermal';
PTtModel = 'P: Isobaric, T: Polythermal-Dwell';
%PTtModel = 'P: Polybaric-Dwell, T: Isothermal-quench';
%================================!!!

%Constants used for calculations
SurfTens = 0.22; %Value for surface tension (N/m)
melt_Rho = 2350; %Melt density in kg/m^3

%Spatial parameters
Nb = 1e11/volScale; %Bubble number density (number per m^3)
Phi_0 = 1*1e-6; %Initial gas volume fraction (Phi)
%R_0 = Radius(Nb,Phi_0); %R_0 (m) calculated from Phi and Nb
R_0 = 1e-5; %R_0 (m) set independently

%Finite difference parameters
Nodes = 2000; %Number of nodes in spatial discretization

%Numerical tolerance:
%[Absolute tolerance, relative tolerance], see:
% https://www.mathworks.com/help/simbio/ref/absolutetolerance.html
% https://www.mathworks.com/help/simulink/gui/relative-tolerance.html
%For additional information
Numerical_Tolerance = [1e-5, 1e-5];

%Run a switch to select model parameters depending on the PTt pathway
%chosen by the user. Note that these are examples, and the user should
%ensure they choose the correct consituent models (above) for the given
%input parameters
switch PTtModel
    case 'P: Isobaric, T: Isothermal'
        T_0 = 1100 + 273.15; %Initial temperature in K (C to K)
        T_f = []; %Final temperature in K
        dTdt = []; %Rate of temperature change in K/s
        P_0 = 20 * 1e6; %Initial pressure in Pa (MPa to Pa)
        P_f = []; %Final pressure in Pa (MPa to Pa)
        dPdt = []; %Rate of pressure change in Pa S (Mpa to Pa)
        H2Ot_0 = 1.4; %initial water concentration (wt. %)
        t_nuc = 0; %Nucleation time
        t_f = 3000; %Experiment duration (seconds)
        
    case 'P: Isobaric, T: Polythermal-Dwell'
        T_0 = 900 + 273.15; %Initial temperature in K
        T_f = 1000 + 273.15; %Final temperature in K
        dTdt = 0.5; %Rate of temperature change in K/s
        P_0 = 1e5; %Initial pressure in Pa (MPa to Pa)
        P_f = [];%Final pressure in Pa (MPa to Pa)
        dPdt = []; %Rate of pressure change in Pa S (Mpa to Pa)
        H2Ot_0 = 0.114; %initial water concentration (wt. %)
        t_nuc =0; %Nucleation time
        t_f = 200*3600; %Experiment duration (seconds)
        
    case 'P: Polybaric-Dwell, T: Isothermal-quench'
        T_0 = 850+ 273.15; %Initial temperature in K
        T_f = 500 + 273.15; %Final temperature in K
        dTdt = -50; %Rate of temperature change in K/s
        P_0 = 150 * 1e6; %Initial pressure in Pa (MPa to Pa)
        P_f = 30 * 1e6; %Final pressure in Pa (MPa to Pa)
        dPdt = -1 * 1e6; %Rate of pressure change in Pa S (Mpa to Pa)
        H2Ot_0 = 5.5; %initial water concentration (wt. %)
        t_nuc = 0; %Nucleation time
        %Experiment duration (seconds) note that the experimental duration
        %in this case is calculated as the sum of the depressurization
        %time, the quench time and an additional end of experiment time
        t_f = abs((P_f-P_0)/dPdt)...
            + abs((T_f-T_0)/dTdt) + 25;
        
end %End P-T-t pathway switch
%
%%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
%Run the check function, solubility explore or numerical model
%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!

%Checks input parameters to detect if any are out of component model
%calibration (from literature), or if there are errors in parameter
%signs (e.g., negative temperature rate for positive temperature change)
Checks(SolModel, DiffModel, ViscModel, EOSModel, PTtModel,...
    t_nuc, t_f, T_0, T_f, P_0, P_f, dPdt, dTdt, H2Ot_0, R_0, Nb, SurfTens);

%The switch which either runs the parameter exploration or the model
switch Operation
    case 'Solubility Explore'
        
        %Explore the solubility for the P-T-t pathway
        SolubilityExplore(SolModel, DiffModel, ViscModel, EOSModel, PTtModel,...
            P_0, P_f, dPdt, T_0, T_f, dTdt, t_f, R_0, SurfTens, H2Ot_0)
        
    case 'Run Model'
        
        %Run the numerical model with the set parameters
       % tic
        [t, R, phi, P, T, x_out, H2Ot_all] =  Numerical_Model_v2(Composition, SolModel, DiffModel,...
            ViscModel, EOSModel,PTtModel, SurfTens, melt_Rho, Nodes,...
            R_0, H2Ot_0, Nb, t_nuc, t_f, T_0, T_f, dTdt, P_0,...
            P_f, dPdt, Numerical_Tolerance);
       % toc
        
        %Plot the output of the numerical model
        %plotting(t, R, phi, P, T, x_out, H2Ot_all);
        
RScale={R};
phiScale = phi(end);
tscale={t};
end %End operation switch

%
%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
%Additional functions used by this script
%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!

%==========================================================================
%Function to check if the selected models are in the range of validity, for
%this current program it is only checking Diffusivity and Solubility
%==========================================================================
function Checks(SolModel, DiffModel, ViscModel, EOSModel, PTtModel,...
    t_nuc, t_f,T_0, T_f, P_0, P_f, dPdt, dTdt, H2Ot_0, R_0, Nb, SurfTens)
solFlag = 1;
diffFlag = 1;

%Get the selected functions
[SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
    PTt_fun] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
    EOSModel, PTtModel);

%Check time parameters for errors 
if (size(t_nuc)==0), error('Enter a nucleation time (s), (e.g., t_0 = 0)'),end
if (size(t_f)==0), error('Enter a final time (s), (e.g., t_f = 2000)'),end
if (t_f <= t_nuc), error('The final time is less than or equal to the initial time'), end

%Check water parameter for errors
if (size(H2Ot_0)==0), error('Enter a water concentration, e.g., H2Ot_0 = 1.5'),end
if (H2Ot_0 <= 0), error('Enter a positive water concentration'),end

%If the final temperatures and pressures are empty, if so set them equal to
%the initial temperature and pressure
if (size(T_f)==0), T_f = T_0; end
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

disp(' ');
disp('%%%Model warnings below%%%');
disp(' ');


%Switch-cases to get the parameter bounds for the solubility models
switch SolModel
    case 'Liu 2005'
        Sol_P_min = 0.1*1e6; %Pa
        Sol_P_max = 500*1e6; %Pa
        Sol_T_min = 700 + 273.15; %Kelvin
        Sol_T_max = 1200 + 273.15; %kelvin
        
    case 'Ryan 2015'
        Sol_P_min = 0.1*1e6; %Pa
        Sol_P_max = 0.1*1e6; %Pa
        Sol_T_min = 900 + 273.15; %Kelvin
        Sol_T_max = 1100 + 273.15; %kelvin
        
    otherwise
        warning('Bounds for solubility model not defined');
        solFlag = 0; %flag to ignore bound check
end

%Switch-cases to get the parameter bounds for the diffusion models
switch DiffModel
    case 'Zhang 2010 Metaluminous simple' %Zhang and Ni (2010) equation 15
        Diff_P_min = 0.1*1e6; %Pa
        Diff_P_max = 1900*1e6; %Pa
        Diff_T_min = 500 + 273.15; %Kelvin
        Diff_T_max = 1200 + 273.15; %kelvin
        Diff_H2O_min = 0; %wt. %
        Diff_H2O_max = 2; %wt. %
        
        
    case 'Zhang 2010 Metaluminous' %Zhang and Ni (2010) equations 7a, 13, 14
        Diff_P_min = 0.1*1e6; %Pa
        Diff_P_max = 1900*1e6; %Pa
        Diff_T_min = 676; %Kelvin
        Diff_T_max = 1900; %kelvin
        Diff_H2O_min = 0; %wt. %
        Diff_H2O_max = 8; %wt. %
        
    case 'Zhang 2010 Peralkaline' %Zhang and Ni (2010) equation 7a, 13, 16
        Diff_P_min = 0.1*1e6; %Pa
        Diff_P_max = 1900*1e6; %Pa
        Diff_T_min = 676; %Kelvin
        Diff_T_max = 1900; %kelvin
        Diff_H2O_min = 0; %wt. %
        Diff_H2O_max = 8; %wt. %
        
    otherwise
        warning('Bounds for diffusion model not defined');
        diffFlag = 0; %flag to ignore bound check
end
if solFlag ==1
    %Check if conditions are appropriate for solubility model, if outside of
    %range raise a warning
    if (P_0 < Sol_P_min), warning('Initial pressure below solubility model range'), end
    if (P_f > Sol_P_max), warning('Final pressure above solubility model range'),end
    if (T_0 < Sol_T_min), warning('Initial temperature below solubility model range'), end
    if (T_f > Sol_T_max), warning('Final temperature above solubility model range'),end
end

if diffFlag==1
    %Check if conditions are appropriate for the diffusion model, if outside of
    %range raise a warning
    if (P_0 < Diff_P_min), warning('Initial pressure below diffusion model range'), end
    if (P_f > Diff_P_max), warning('Final pressure above diffusion model range'),end
    if (T_0 < Diff_T_min), warning('Initial temperature below diffusion model range'), end
    if (T_f > Diff_T_max), warning('Final temperature above diffusion model range'),end
    H2O = [SolFun(T_f,P_0),SolFun(T_f,P_f),SolFun(T_0,P_0), SolFun(T_0,P_f), H2Ot_0];
    H2O_min = min(H2O);
    H2O_max = max(H2O);
    if (H2O_min < Diff_H2O_min), warning('H2O will potentially be below diffusion model range'),end
    if (H2O_max > Diff_H2O_max), warning('H2O will potentially be above diffusion model range'),end
end
%disp(' ');
%disp('%%%Warnings or errors from Numerical Model%%%');
%disp(' ');


%==========================================================================
%Function to explore solubility conditions of the model inputs. The user
%can use the visuals of this function to determine the nucleation time
%delay for a given initial water or the solubility of the initial
%conditions..etc.
%==========================================================================
function SolubilityExplore(SolModel, DiffModel, ViscModel, EOSModel, PTtModel,...
    P_0, P_f, dPdt, T_0, T_f, dTdt, t_f, R_0, SurfTens, H2Ot_0)

%Get the selected functions
[SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
    PTt_fun] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
    EOSModel, PTtModel);

%Get range of solubilities following P-T-t Pathway
t = linspace(0,t_f,1000);
PT = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t');
P = PT(:,1);
T = PT(:,2);

%Determine the intersections between the P-T-t solubility and initial
%water, this determines the nucleation delay. The function called below is
%from the MATLAB file exchange authored by Douglas M. Schwarz, see the
%function for more details.


%Solubility of water at the bubble wall including Capilary pressure
H2Oeq = SolFun(T, P + (2*SurfTens/R_0));
[x0,y0,iout,jout] =  intersections(t,ones(length(t),1)*H2Ot_0,t,H2Oeq,1);

%  Commented out to supress when used in voronoi model
% %Plot
% figure(1)
% hold on
% plot(t,H2Oeq, 'k-', 'linewidth', 2);
% plot(x0, y0,'marker','p','markerfacecolor','r','markeredgecolor','r','markersize',20);
% yline(H2Ot_0, 'b-', 'linewidth', 2);
% grid on
% box on
% 
% ylabel('Solubility H_2O_t, (wt. %)');
% xlabel ('Time, \itt\rm (seconds)');
% 
% %Create a legend
% h(1) = plot(NaN,NaN,'linestyle','-',...
%     'color','k','linewidth',2);
% h(2) = plot(NaN,NaN,'linestyle','-',...
%     'color','b','linewidth',2);
% h(3) = plot(NaN,NaN,'pr','markerfacecolor','r','markeredgecolor','r','markersize',10);
% 
% lbl = {['Model solubility'];['Inital H_2O_t, (wt. %)'];['Saturation']};
% %Set positions of the Legend
% legend(h,lbl,'Location','northeast','FontName','Times New Roman','FontSize',10);

%==========================================================================
%Function to plot the output.
%==========================================================================
%function plotting(t, R, phi, P, T, x_out, H2Ot_all)
%Plot the results of the numerical model
% figure(1)
% subplot(2,2,1)
% hold on
% plot(t, R*1e6,'k-','linewidth',2);
% grid on
% box on
% ylabel('Bubble radius, \itR\rm (\mum)');
% xlabel ('Time, \itt\rm (seconds)');
% 
% subplot(2,2,2)
% hold on
% plot(t, phi,'k-','linewidth',2);
% grid on
% box on
% ylabel('Gas volume fraction, \phi');
% xlabel ('Time, \itt\rm (seconds)');
% 
% subplot(2,2,3)
% hold on
% plot(t, P.*1e-6,'k-','linewidth',2);
% ylabel('Pressure, (MPa)');
% xlabel ('Time, \itt\rm (seconds)');
% grid on
% box on
% 
% subplot(2,2,4)
% plot(t, T,'k-', 'linewidth',2);
% ylabel('Temperature, \itT\rm  (Kelvin)');
% xlabel ('Time, \itt\rm (seconds)');
% grid on
% box on
% 
% figure(2)
% ClrMap = jet; %Set colourmap to jet
% %Set the colour values based on time
% colour_values = linspace(t(1),t(end),length(ClrMap));
% hold on
% for i = 1:1:size(x_out,2)
%     Clr = [interp1(colour_values,ClrMap(:,1),t(i))...
%         interp1(colour_values,ClrMap(:,2),t(i))...
%         interp1(colour_values,ClrMap(:,3),t(i))];
%     plot(x_out(:,i).*1e6,H2Ot_all(:,i),':','color',Clr, 'linewidth',2);
% end
% Cmap = ClrMap;
% colormap(gca,Cmap)
% hcb = colorbar;
% caxis([colour_values(1),colour_values(end)]);
% colorTitleHandle = get(hcb,'Title');
% titleString = {'time, \itt\rm (s)'};
% set(colorTitleHandle ,'String',titleString,'FontName','Times New Roman','FontSize',12);
% set(hcb,'YTick',colour_values(1):colour_values(end)/5:colour_values(end))
% 
% grid on
% box on
% 
% xlabel('Distance from bubble, (\mum)');
% ylabel('H_2O_t, (wt. %)');


function R = Radius(Nb,Phi)
meltvolume = 1./Nb;
gasvolume = Phi.*meltvolume./(1 - Phi);
R = nthroot(gasvolume./((4/3)*pi()),3);

function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)
%INTERSECTIONS Intersections of curves.
%   Computes the (x,y) locations where two curves intersect.  The curves
%   can be broken with NaNs or have vertical segments.
%
% Example:
%   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% where X1 and Y1 are equal-length vectors of at least two points and
% represent curve 1.  Similarly, X2 and Y2 represent curve 2.
% X0 and Y0 are column vectors containing the points at which the two
% curves intersect.
%
% ROBUST (optional) set to 1 or true means to use a slight variation of the
% algorithm that might return duplicates of some intersection points, and
% then remove those duplicates.  The default is true, but since the
% algorithm is slightly slower you can set it to false if you know that
% your curves don't intersect at any segment boundaries.  Also, the robust
% version properly handles parallel and overlapping segments.
%
% The algorithm can return two additional vectors that indicate which
% segment pairs contain intersections and where they are:
%
%   [X0,Y0,I,J] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% For each element of the vector I, I(k) = (segment number of (X1,Y1)) +
% (how far along this segment the intersection is).  For example, if I(k) =
% 45.25 then the intersection lies a quarter of the way between the line
% segment connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for
% the vector J and the segments in (X2,Y2).
%
% You can also get intersections of a curve with itself.  Simply pass in
% only one curve, i.e.,
%
%   [X0,Y0] = intersections(X1,Y1,ROBUST);
%
% where, as before, ROBUST is optional.
% Version: 2.0, 25 May 2017
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})
% Theory of operation:
%
% Given two line segments, L1 and L2,
%
%   L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
%   L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))
%
% we can write four equations with four unknowns and then solve them.  The
% four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of
% L1 and L2, t1 is the distance from the starting point of L1 to the
% intersection relative to the length of L1 and t2 is the distance from the
% starting point of L2 to the intersection relative to the length of L2.
%
% So, the four equations are
%
%    (x1(2) - x1(1))*t1 = x0 - x1(1)
%    (x2(2) - x2(1))*t2 = x0 - x2(1)
%    (y1(2) - y1(1))*t1 = y0 - y1(1)
%    (y2(2) - y2(1))*t2 = y0 - y2(1)
%
% Rearranging and writing in matrix form,
%
%  [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
%        0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
%   y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
%        0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]
%
% Let's call that A*T = B.  We can solve for T with T = A\B.
%
% Once we have our solution we just have to look at t1 and t2 to determine
% whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two
% line segments cross and we can include (x0,y0) in the output.
%
% In principle, we have to perform this computation on every pair of line
% segments in the input data.  This can be quite a large number of pairs so
% we will reduce it by doing a simple preliminary check to eliminate line
% segment pairs that could not possibly cross.  The check is to look at the
% smallest enclosing rectangles (with sides parallel to the axes) for each
% line segment pair and see if they overlap.  If they do then we have to
% compute t1 and t2 (via the A\B computation) to see if the line segments
% cross, but if they don't then the line segments cannot cross.  In a
% typical application, this technique will eliminate most of the potential
% line segment pairs.
% Input checks.
if verLessThan('matlab','7.13')
    error(nargchk(2,5,nargin)) %#ok<NCHKN>
else
    narginchk(2,5)
end
% Adjustments based on number of arguments.
switch nargin
    case 2
        robust = true;
        x2 = x1;
        y2 = y1;
        self_intersect = true;
    case 3
        robust = x2;
        x2 = x1;
        y2 = y1;
        self_intersect = true;
    case 4
        robust = true;
        self_intersect = false;
    case 5
        self_intersect = false;
end
% x1 and y1 must be vectors with same number of points (at least 2).
if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
        length(x1) ~= length(y1)
    error('X1 and Y1 must be equal-length vectors of at least 2 points.')
end
% x2 and y2 must be vectors with same number of points (at least 2).
if sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1 || ...
        length(x2) ~= length(y2)
    error('X2 and Y2 must be equal-length vectors of at least 2 points.')
end
% Force all inputs to be column vectors.
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);
% Compute number of line segments in each curve and some differences we'll
% need later.
n1 = length(x1) - 1;
n2 = length(x2) - 1;
xy1 = [x1 y1];
xy2 = [x2 y2];
dxy1 = diff(xy1);
dxy2 = diff(xy2);
% Determine the combinations of i and j where the rectangle enclosing the
% i'th line segment of curve 1 overlaps with the rectangle enclosing the
% j'th line segment of curve 2.
% Original method that works in old MATLAB versions, but is slower than
% using binary singleton expansion (explicit or implicit).
% [i,j] = find( ...
% 	repmat(mvmin(x1),1,n2) <= repmat(mvmax(x2).',n1,1) & ...
% 	repmat(mvmax(x1),1,n2) >= repmat(mvmin(x2).',n1,1) & ...
% 	repmat(mvmin(y1),1,n2) <= repmat(mvmax(y2).',n1,1) & ...
% 	repmat(mvmax(y1),1,n2) >= repmat(mvmin(y2).',n1,1));
% Select an algorithm based on MATLAB version and number of line
% segments in each curve.  We want to avoid forming large matrices for
% large numbers of line segments.  If the matrices are not too large,
% choose the best method available for the MATLAB version.
if n1 > 1000 || n2 > 1000 || verLessThan('matlab','7.4')
    % Determine which curve has the most line segments.
    if n1 >= n2
        % Curve 1 has more segments, loop over segments of curve 2.
        ijc = cell(1,n2);
        min_x1 = mvmin(x1);
        max_x1 = mvmax(x1);
        min_y1 = mvmin(y1);
        max_y1 = mvmax(y1);
        for k = 1:n2
            k1 = k + 1;
            ijc{k} = find( ...
                min_x1 <= max(x2(k),x2(k1)) & max_x1 >= min(x2(k),x2(k1)) & ...
                min_y1 <= max(y2(k),y2(k1)) & max_y1 >= min(y2(k),y2(k1)));
            ijc{k}(:,2) = k;
        end
        ij = vertcat(ijc{:});
        i = ij(:,1);
        j = ij(:,2);
    else
        % Curve 2 has more segments, loop over segments of curve 1.
        ijc = cell(1,n1);
        min_x2 = mvmin(x2);
        max_x2 = mvmax(x2);
        min_y2 = mvmin(y2);
        max_y2 = mvmax(y2);
        for k = 1:n1
            k1 = k + 1;
            ijc{k}(:,2) = find( ...
                min_x2 <= max(x1(k),x1(k1)) & max_x2 >= min(x1(k),x1(k1)) & ...
                min_y2 <= max(y1(k),y1(k1)) & max_y2 >= min(y1(k),y1(k1)));
            ijc{k}(:,1) = k;
        end
        ij = vertcat(ijc{:});
        i = ij(:,1);
        j = ij(:,2);
    end
    
elseif verLessThan('matlab','9.1')
    % Use bsxfun.
    [i,j] = find( ...
        bsxfun(@le,mvmin(x1),mvmax(x2).') & ...
        bsxfun(@ge,mvmax(x1),mvmin(x2).') & ...
        bsxfun(@le,mvmin(y1),mvmax(y2).') & ...
        bsxfun(@ge,mvmax(y1),mvmin(y2).'));
    
else
    % Use implicit expansion.
    [i,j] = find( ...
        mvmin(x1) <= mvmax(x2).' & mvmax(x1) >= mvmin(x2).' & ...
        mvmin(y1) <= mvmax(y2).' & mvmax(y1) >= mvmin(y2).');
    
end
% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take
% advantage of the fact that NaNs propagate through calculations, in
% particular subtraction (in the calculation of dxy1 and dxy2, which we
% need anyway) and addition.
% At the same time we can remove redundant combinations of i and j in the
% case of finding intersections of a line with itself.
if self_intersect
    remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
else
    remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
end
i(remove) = [];
j(remove) = [];
% Initialize matrices.  We'll put the T's and B's in matrices and use them
% one column at a time.  AA is a 3-D extension of A where we'll use one
% plane at a time.
n = length(i);
T = zeros(4,n);
AA = zeros(4,4,n);
AA([1 2],3,:) = -1;
AA([3 4],4,:) = -1;
AA([1 3],1,:) = dxy1(i,:).';
AA([2 4],2,:) = dxy2(j,:).';
B = -[x1(i) x2(j) y1(i) y2(j)].';
% Loop through possibilities.  Trap singularity warning and then use
% lastwarn to see if that plane of AA is near singular.  Process any such
% segment pairs to determine if they are colinear (overlap) or merely
% parallel.  That test consists of checking to see if one of the endpoints
% of the curve 2 segment lies on the curve 1 segment.  This is done by
% checking the cross product
%
%   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
%
% If this is close to zero then the segments overlap.
% If the robust option is false then we assume no two segment pairs are
% parallel and just go ahead and do the computation.  If A is ever singular
% a warning will appear.  This is faster and obviously you should use it
% only when you know you will never have overlapping or parallel segment
% pairs.
if robust
    overlap = false(n,1);
    warning_state = warning('off','MATLAB:singularMatrix');
    % Use try-catch to guarantee original warning state is restored.
    try
        lastwarn('')
        for k = 1:n
            T(:,k) = AA(:,:,k)\B(:,k);
            [unused,last_warn] = lastwarn; %#ok<ASGLU>
            lastwarn('')
            if strcmp(last_warn,'MATLAB:singularMatrix')
                % Force in_range(k) to be false.
                T(1,k) = NaN;
                % Determine if these segments overlap or are just parallel.
                overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
            end
        end
        warning(warning_state)
    catch err
        warning(warning_state)
        rethrow(err)
    end
    % Find where t1 and t2 are between 0 and 1 and return the corresponding
    % x0 and y0 values.
    in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
    % For overlapping segment pairs the algorithm will return an
    % intersection point that is at the center of the overlapping region.
    if any(overlap)
        ia = i(overlap);
        ja = j(overlap);
        % set x0 and y0 to middle of overlapping region.
        T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
            min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
        T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
            min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
        selected = in_range | overlap;
    else
        selected = in_range;
    end
    xy0 = T(3:4,selected).';
    
    % Remove duplicate intersection points.
    [xy0,index] = unique(xy0,'rows');
    x0 = xy0(:,1);
    y0 = xy0(:,2);
    
    % Compute how far along each line segment the intersections are.
    if nargout > 2
        sel_index = find(selected);
        sel = sel_index(index);
        iout = i(sel) + T(1,sel).';
        jout = j(sel) + T(2,sel).';
    end
else % non-robust option
    for k = 1:n
        [L,U] = lu(AA(:,:,k));
        T(:,k) = U\(L\B(:,k));
    end
    
    % Find where t1 and t2 are between 0 and 1 and return the corresponding
    % x0 and y0 values.
    in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
    x0 = T(3,in_range).';
    y0 = T(4,in_range).';
    
    % Compute how far along each line segment the intersections are.
    if nargout > 2
        iout = i(in_range) + T(1,in_range).';
        jout = j(in_range) + T(2,in_range).';
    end
end
% Plot the results (useful for debugging).

function y = mvmin(x)
% Faster implementation of movmin(x,k) when k = 1.
y = min(x(1:end-1),x(2:end));
function y = mvmax(x)
% Faster implementation of movmax(x,k) when k = 1.
y = max(x(1:end-1),x(2:end));

%==========================================================================
% References
%==========================================================================
% Giordano, D., Dingwell, D.B. and Romano, C., 2000. Viscosity of a
% Teide phonolite in the welding interval. J. Volcanol. Geotherm. Res., 103(1): 239-245.
%
% Giordano, D., Russell, J.K. and Dingwell, D.B., 2008. Viscosity of
% magmatic liquids: A model. Earth Planet. Sci. Lett., 271(1): 123-134.
%
% Hess, K.-U. and Dingwell, D.B., 1996. Viscosities of hydrous
% leucogranitic melts: A non-Arrhenian model. Am. Mineral., 81(9-10): 1297-1300.
%
% Iacono-Marziano, G., Schmidt, B.C. and Dolfi, D., 2007. Equilibrium and
% disequilibrium degassing of a phonolitic melt (Vesuvius AD 79 “white pumice”)
% simulated by decompression experiments. J. Volcanol. Geotherm. Res., 161(3): 151-164.
%
% Liu, Y., Zhang, Y. and Behrens, H., 2005. Solubility of H2O in rhyolitic
% melts at low pressures and a new empirical model for mixed H2O–CO2 solubility
% in rhyolitic melts. J. Volcanol. Geotherm. Res., 143(1): 219-235.
%
% Mangan, M. and Sisson, T., 2000. Delayed, disequilibrium degassing
% in rhyolite magma: decompression experiments and implications for
% explosive volcanism. Earth Planet. Sci. Lett., 183(3-4): 441-455.
%
% Mourtada-Bonnefoi, C.C. and Laporte, D., 2004. Kinetics of bubble
% nucleation in a rhyolitic melt: an experimental study of the effect of
% ascent rate. Earth Planet. Sci. Lett., 218(3-4): 521-537.
%
% Pitzer, K.S. and Sterner, S.M., 1994. Equations of state valid
% continuously from zero to extreme pressures for
% H2O and CO2. J. Chem. Phys., 101(4): 3111-3116.
%
% Ryan, A.G., Russell, J.K., Nichols, A.R., Hess, K.-U. and
% Porritt, L.A., 2015. Experiments and models on H2O retrograde
% solubility in volcanic systems. Am. Mineral., 100(4): 774-786.
%
% Tuffen, H. and Castro, J.M., 2009. The emplacement of an obsidian
% dyke through thin ice: Hrafntinnuhryggur, Krafla Iceland.
% J. Volcanol. Geotherm. Res., 185(4): 352-366.
%
% Zhang, Y. and Ni, H., 2010. Diffusion of H, C, and O components
% in silicate melts. Rev. Mineral. Geochem., 72(1): 171-225.


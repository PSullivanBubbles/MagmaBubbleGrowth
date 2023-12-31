#include <cmath>
#include <iostream>
#include <vector>
#include <valarray>
#include <string>

#include "Numerical_Model.h"
#include "petscerror.h"
  

void Checks(std::string SolModel, std::string DiffModel, double t_nuc, double t_f, std::valarray<double> T, std::valarray<double> P, double H2Ot_0){
    int solFlag = 1;
    int diffFlag = 1;


    return;
} 

    
void SolubilityExplore(std::string SolModel, std::valarray<double> T_0,std::valarray<double> t_T,std::valarray<double> P_0,std::valarray<double> t_P, double R_0, double SurfTens,double H2Ot_0){
/*
//#%Get the selected functions
//#[SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
//#    PTt_fun] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
//#    EOSModel, PTtModel);

//    #%Get range of solubilities following P-T-t Pathway
//    t = numpy.linspace(0,t_T[-1],1000)
//    #PT = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t);
//    P = gFP.PtFun(P_0,t_P,t)
//    P = gFP.PtFun(T_0,t_T,t)

//#%Determine the intersections between the P-T-t solubility and initial
//#%water, this determines the nucleation delay. The function called below is
//#%from the MATLAB file exchange authored by Douglas M. Schwarz, see the
//#%function for more details.


//#%Solubility of water at the bubble wall including Capilary pressure
    double H2Oeq = SolFun(T, P + (2*SurfTens/R_0),SolModel)
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
    plt.legend(h,lbl,'Location','northeast','FontName','Times New Roman','FontSize',10)*/
}

//#%==========================================================================
//#%Function to plot the output.
//#%==========================================================================
void plotting(std::valarray<double> t,std::valarray<double> R,std::valarray<double> phi,std::valarray<double> P,std::valarray<double> T,std::valarray<double> x_out,std::valarray<double> H2Ot_all){
/*
    #%Plot the results of the numerical model
    plt.figure(1)
    plt.subplot(2,2,1)
    #plt.hold on
    plt.plot(t, R*1e6)
    plt.grid(1)
    plt.box(1)
    plt.ylabel('Bubble radius, \itR\rm (\mum)')
    plt.xlabel ('Time, \itt\rm (seconds)')

    plt.subplot(2,2,2)
    #plt.hold on
    plt.plot(t, phi)
    plt.grid(1)
    plt.box(1)
    plt.ylabel('Gas volume fraction, \phi')
    plt.xlabel ('Time, \itt\rm (seconds)')

    plt.subplot(2,2,3)
    #plt.hold on
    plt.plot(t, P*1e-6)
    plt.ylabel('Pressure, (MPa)')
    plt.xlabel ('Time, \itt\rm (seconds)')
    plt.grid(1)
    plt.box(1)

    plt.subplot(2,2,4)
    plt.plot(t, T)
    plt.ylabel('Temperature, \itT\rm  (Kelvin)')
    plt.xlabel ('Time, \itt\rm (seconds)')
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
        plt.plot(x_out[:,i]*1e6,H2Ot_all[:,i])

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

    plt.xlabel('Distance from bubble, (\mum)')
    plt.ylabel('H_2O_t, (wt. %)')
    plt.show()

*/
}

//#%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
//#%Choose which operation: Solubility explore or Model run
     
int main(){
  
    std::string Operation = "Run Model";
    //"""
    //%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
    //%Set the model properties
    //%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
    //"""
        
    //#%Melt composition
    double wt_dry = 0.01;// #H2O wt. % that defines "dry"
    //#[SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1]
    std::valarray<double> PhonoliteVAD79, Krafla, PCD, ROGD;

    PhonoliteVAD79 = {57.15, 0.3, 21.34, 2.70, 0.14, 0.39, 3.26, 5.16, 9.46, 0.09, wt_dry, 0};
    Krafla = {75.17, 0.22, 12.02, 3.13, 0.11, 0.09, 1.66, 4.58, 2.88, 0, wt_dry, 0};
    PCD = {75.64, 0.08, 12.38, 0.94, 0.07, 0.03, 0.54, 4.17, 4.72, 0.01, wt_dry, 0};
    ROGD = {76.51, 0.03, 12.56, 0.7, 0.07, 0.01, 0.25, 4.47, 4.24, 0, wt_dry, 0};
     
    //#%================================!!!
    //#%~~~Composition~~~
    std::valarray<double> Composition = Krafla; //# %Tuffen and Castro (2009) locality AO
    //#%Composition = PhonoliteVAD79; %Iacono-Marziano et al. (2007)
    //#%Composition = PCD; %Mangan and Sisson (2000)
    //#%Composition = ROGD; %Mourtada-Bonnefoi and Laporte (2004)
    //#%================================!!!
    //#%================================!!!
    //#%~~~solubility model~~~  
    std::string SolModel = "Ryan 2015"; //# %Ryan et al. (2015)
    //std::string SolModel = "Liu 2005"; //# %Liu et al. (2005)  
    //#%================================!!!
    //#%================================!!!
    //#%~~~Diffusion Model~~~
    //#DiffModel = 'Zhang 2010 Metaluminous simple'# %Zhang and Ni (2010) equation 15
    std::string DiffModel = "Zhang 2010 Metaluminous"; //# %Zhang and Ni (2010) equations 7a, 13, 14
    //#DiffModel = 'Zhang 2010 Peralkaline'# %Zhang and Ni (2010) equation 7a, 13, 16
    //#DiffModel = 'Constant'# %A constant value that can be defined in getFunctions.m
    //#================================!!!
    //#%================================!!!
    //#%~~~Viscosity Model~~~
    //#ViscModel = 'Hess and Dingwell 1996'# %Hess and Dingwell (1996)   
    //#ViscModel = 'Peralkaline Giordano 2000'# %Giordano et al. (2000)
    std::string ViscModel = "Giordano 2008"; // # %Giordano et al. (2008)
    //#%================================!!!
    //#%================================!!!
    //#%~~~EOS model~~~ 
    std::string EOSModel;
    EOSModel = "Pitzer and Sterner"; //# %Pitzer and Sterner (1994)
    //#EOSModemalloc_consolidate(): invalid chunk sizel = 'Ideal Gas Law'
    //#%================================!!!
    //#%================================!!!  
                       
    //#Constants used for calculations
    double SurfTens = 0.22; //#Value for surface tension (N/m)
    double melt_Rho = 2350; //#Melt density in kg/m^3
         
    //#Spatial parameters
    double Nb = 1e11; //#Bubble number density (number per m^3)
    double Phi_0 = 1*1e-6; //#Initial gas volume fraction (Phi)
    //#R_0 = Radius(Nb,Phi_0); //# %R_0 (m) calculated from Phi and Nb
    double R_0 = 3e-5; //#R_0 (m) set independently
   
    //#Finite difference parameters
    int Nodes = 500; //#Number of nodes in spatial discretization
                 
    //#%Numerical tolerance:         
    //#%[Absolute tolerance, relative tolerance], see:
    //#% https://www.mathworks.com/help/simbio/ref/absolutetolerance.html
    //#% https://www.mathworks.com/help/simulink/gui/relative-tolerance.html
    //#%For additional information
    std::valarray<double> Numerical_Tolerance = {1e-5, 1e-5};
         
        
    double t_nuc = 0;     
    double t_f = 3000; 
 
    std::valarray<double> T_0 = {1100+ 273.15};
    std::valarray<double> t_T = {0};
  
    std::valarray<double> P_0 = {1*1e5};
    std::valarray<double> t_P = {0};
       
    double H2Ot_0 = 0.1;  
             
    //#%%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
    //#%Run the check function, solubility explore or numerical model
    //#%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!

    //#%Checks input parameters to detect if any are out of component model
    //#%calibration (from literature), or if there are errors in parameter
    //#%signs (e.g., negative temperature rate for positive temperature change)
    Checks(SolModel, DiffModel, t_nuc, t_f,T_0, P_0, H2Ot_0);
        
    //#%The switch which either runs the parameter exploration or the model
    if (Operation=="Solubility Explore") 
            //#%Explore the solubility for the P-T-t pathway
        SolubilityExplore(SolModel, T_0,t_T,P_0,t_P, R_0, SurfTens, H2Ot_0);
                 
    else if (Operation=="Run Model"){
                                   
        //#%Run the numerical model with the set parameters
        //#tic
        std::valarray<double> t;
        std::valarray<double> R;
        std::valarray<double> phi;
        std::valarray<double> P;
        std::valarray<double> T;
        std::valarray<std::valarray<double>> x_out;
        std::valarray<std::valarray<double>> H2Ot_all;
 
         std::cout<<"Starting Numerical Model from main \n\n";
        int argc =0;
        char **args =NULL;
        PetscInitialize(&argc, &args, (char *)0, NULL); 
        CHKERRQ(Numerical_Model_v2(Composition, SolModel, DiffModel, ViscModel, EOSModel, SurfTens, melt_Rho, Nodes, R_0, H2Ot_0, Nb, t_nuc, T_0, t_T, P_0, t_P, Numerical_Tolerance,
        t,R,phi, P, T,x_out,H2Ot_all));

        //#toc
        //std::cout<<"--- "<<(time.time() - start_time)<<" seconds --- \n" ;
        //#%Plot the output of the numerical model
        //plotting(t, R, phi, P, T, x_out, H2Ot_all);

        std::cout<<"Final Radius - "<<R[R.size()-1]<<"\n";

        for(int i =0; i<t.size();i++)
            std::cout<<t[i]<<"\t"<<R[i]<<"\r\n";

    }
    
   
    return 0;
}

//#%
//#%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
//#%Additional functions used by this script
//#%+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!




double  Radius(double Nb,double Phi){
    double meltvolume = 1./Nb;
    double gasvolume = Phi*meltvolume/(1 - Phi);
    return pow((gasvolume/((4.0/3.0)*PI)),(1.0/3.0));
    }

//#function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)



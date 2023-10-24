#include <cmath>
#include <iostream>
#include <vector>
#include <valarray>
#include <string>

#include "boost_1_82_0/boost/math/tools/roots.hpp"
 
namespace tools = boost::math::tools;

const double PI  =3.141592653589793238463;



// DECLARE ALL FUNCTIONS TO AVOID ISSUES WITH THEIR LOCATIONS

std::valarray<std::valarray<double>> molepct_grd(std::valarray<std::valarray<double>> wtm);
std::valarray<std::valarray<double>> Giordano2008_Model(std::valarray<double> H2Ot, std::valarray<double> Composition);
double PS_myfun(double rho, double a[10], double PRT);
double density(double P, double T);
std::valarray<double> Giordano_2008_visc(std::valarray<double> H2Ot, double T, std::valarray<double> Composition);
double TtFun(std::valarray<double> Tin,std::valarray<double> tin,double t);
double PtFun(std::valarray<double> Pin,std::valarray<double> tin,double t);
double pb_fun(double m, double T, double R,std::string EOSModel);
double m0_fun(double R, double P, double T,std::string EOSModel);
std::valarray<double> ViscFun(std::valarray<double> H2Ot,double T,std::valarray<double> Composition, std::string ViscModel);
std::valarray<double> DiffFun(double T,double P,std::valarray<double> H2Ot, double W, std::string DiffModel);
double SolFun(double T,double P, std::string SolModel);



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

//Return the target Solubility 

double SolFun(double T,double P, std::string SolModel){

    double H2Oeq;
    if(SolModel=="Ryan 2015"){
        // Regressed solubility function from Ryan et al. (2015)
        H2Oeq = 92.3/T + 0.0287;
        return H2Oeq;
    }    
    else if(SolModel == "Liu 2005"){
        //Solubility function from Liu et al. (2005)
        // convert pressure to MPa from Pa
        P = P*1e-6; 
        H2Oeq =((354.94*sqrt(P)+9.623*P-1.5223*(pow(P,1.5)))/T) + 0.0012439*(pow(P,1.5));
        return H2Oeq;
    }
    return 0;
}


//Return the target Diffusivity 

std::valarray<double> DiffFun(double T,double P,std::valarray<double> H2Ot, double W, std::string DiffModel){

    std::valarray<double> D;

    if (DiffModel== "Zhang 2010 Metaluminous simple"){
        //%Diffusivity function for metaluminous rhyolite from Zhang and Ni (2010)
        //%(Equation 15 therein)
        P=P/1e9; //Convert pressure to GPa from Pa
        D = H2Ot*exp( -18.1 +(1.888*P)-((9699+3626*P)/T));
        return D;
    }
    else if (DiffModel== "Zhang 2010 Metaluminous"){
        //Diffusivity function for metaluminous rhyolite from Zhang and Ni (2010)
        //(Equations 7a, 13, 14 therein)
        P=P/1e9; //Convert pressure to GPa from Pa
        //convert from wt.% to mole fraction
        std::valarray<double> XH2Ot=(H2Ot/18.015)/((H2Ot/18.015)+(100-H2Ot)/W); 
        //Compute the diffusivity of molecular water (Zhang & Ni 2010, equation 14)
        std::valarray<double> DH2Om = exp(-14.26 + 1.888*P - 37.26*XH2Ot - ((12939 + 3626*P - 75884*XH2Ot)/T));
        //%Get the speciation constant (Zhang & Ni 2010, equation 7a)
        double K = exp(1.876 - 3110/T);
        //%Compute the total water diffusivity from Diffusivity of
        //%molecular water and speciation constant (Zhang & Ni 2010, equation 13)
        D = DH2Om*(1-(0.5-XH2Ot)/sqrt(((4.0/K-1)*(XH2Ot-(XH2Ot*XH2Ot))+0.25)));
        return D;
    }
    else if (DiffModel== "Zhang 2010 Peralkaline"){
        //Diffusivity function for peralkaline rhyolite from Zhang and Ni (2010)
        //(Equations 7a, 13, 16 therein)
        P=P/1e9; //Convert from Pa to GPa
        //convert from wt.% to mole fraction
        std::valarray<double> XH2Ot=(H2Ot/18.015)/((H2Ot/18.015)+(100-H2Ot)/W); 
        //Compute the diffusivity of molecular water (Zhang & Ni 2010,
        //equation 16)
        std::valarray<double> DH2Om = exp(-12.79 - 27.87*XH2Ot - ((13939 + 1230.*P - 60559.*XH2Ot)/T));
        //%Get the speciation constant (Zhang & Ni 2010, equation 7a)
        double K = exp(1.876 - 3110./T);
        D = DH2Om*(1-(0.5-XH2Ot)/sqrt(((4.0/K-1)*(XH2Ot-(XH2Ot*XH2Ot))+0.25)));
        return D;
    }  
    else if (DiffModel== "Constant"){
        //Diffusivity function for a constant diffusivty
        double DiffVal = 1e-11; // Diffusivity value in m^2 s^-1
        D =DiffVal*(T/T);
        return D;
    }

return D;

}

//%==========================================================================
//%Return the target Viscosity function
//%==========================================================================

std::valarray<double> ViscFun(std::valarray<double> H2Ot,double T,std::valarray<double> Composition, std::string ViscModel){

    if (ViscModel=="Giordano 2008"){
        
        return Giordano_2008_visc(H2Ot,T,Composition);

    }


    else if (ViscModel=="Hess and Dingwell 1996"){
        return pow(10,(-3.545+0.833*log(H2Ot)+(9601-2368*log(H2Ot))/(T-(195.7+32.25*log(H2Ot)))));
    }

    else if (ViscModel=="Peralkaline Giordano 2000"){
        return pow(10,(-5.9-0.286*log(H2Ot)+(10775.4-394.8*H2Ot)/(T-148.7+21.65*log(H2Ot))));
    }
    return H2Ot-H2Ot;
}


//%==========================================================================
//%Return the target functions for Equation of State
//%==========================================================================


double m0_fun(double R, double P, double T,std::string EOSModel){
    double V=(4*PI/3)*pow(R,3); // # Volume of the bubble (M**3)
    if(EOSModel=="Ideal Gas Law"){
        // Moles of water in the bubble (from volume and P)
        // Gas constant = 8.314 J mol**-1K-1
        // 1 Pa*m**3 = 1 joule
        double n=P*V/(8.314*T);
        //# 18.015 g/mol * n moles = mass (g)* 1kg / 1000g
        return 18.015*n/1000;
    }
    else if(EOSModel=="Pitzer and Sterner"){
        // calculate initial mass of water in the bubble (kg)
        // Returns the density of the bubble in kg/m**3
        double rho = density(P,T);
        return rho*V;}
    return 0;
}

double pb_fun(double m, double T, double R,std::string EOSModel){
    if(EOSModel=="Ideal Gas Law")
        return (((m*1000)/18.015)*8.314*T)/((4.0*PI/3.0)*pow(R,3));
    else if (EOSModel=="Pitzer and Sterner"){
        
        //%Volume of the bubble in m^3
        double Vbub = (4.0/3.0)*PI*pow(R,3);

        //%Rho must be in terms of mol/cm^3 for the pitzer & sterner equation of
        //%state. Therefore Volume: m^3 * (1e6 cm^3 / 1 m ^3).
        //%and Mass: (kg * (1000g / 1 kg))*(1 mole / 18.015 g)
        double rho = ((m*1000)/18.015)/(Vbub*1e6);

        //%Get the pitzer and sterner coefficients
        static double b[10][6] = {0};
        b[0][2]=0.24657688e6; 
        b[0][3]=0.51359951e2 ;
        b[1][2]=0.58638965e0 ;
        b[1][3]=-0.28646939e-2; 
        b[1][4]=0.31375577e-4 ;
        b[2][2]=-0.62783840e1 ;
        b[2][3]=0.14791599e-1 ;
        b[2][4]=0.35779579e-3 ;
        b[2][5]=0.15432925e-7 ;
        b[3][3]=-0.42719875e0 ;
        b[3][4]=-0.16325155e-4 ;
        b[4][2]=0.56654978e4 ;
        b[4][3]=-0.16580167e2 ;
        b[4][4]=0.76560762e-1 ;
        b[5][3]=0.10917883e0 ;
        b[6][0]=0.38878656e13 ;
        b[6][1]=-0.13494878e9 ;
        b[6][2]=0.30916564e6 ;
        b[6][3]=0.75591105e1 ;
        b[7][2]=-0.65537898e5 ;
        b[7][3]=0.18810675e3 ;
        b[8][0]=-0.14182435e14; 
        b[8][1]=0.18165390e9 ;
        b[8][2]=-0.19769068e6 ;
        b[8][3]=-0.23530318e2 ;
        b[9][2]=0.92093375e5 ;
        b[9][3]=0.12246777e3 ;

        double a[10];
        for (int i=0; i<10; i++)
            a[i]=b[i][0]*pow(T,-4) + b[i][1]*pow(T,-2) + b[i][2]*pow(T,-1) + b[i][3] + b[i][4]*T + b[i][5]*pow(T,2);
        

        //%P [bars], R [cm^3*bar/K/mol] and T [K]

        //%NOTE Gas constant is in terms of [cm^3*bar/K/mol] as the equation of PS,
        //%94 has pressure in terms of bar and density in terms of mol cm^-3
        double pb = (rho+a[0]*pow(rho,2)-pow(rho,2)*((a[2]+2*a[3]*rho+3*a[4]*pow(rho,2)+4*a[5]*pow(rho,3))/(pow(a[1]+a[2]*rho+a[3]*pow(rho,2)+a[4]*pow(rho,3)+a[5]*pow(rho,4),2)))+a[6]*pow(rho,2)*exp(-a[7]*rho)+a[8]*pow(rho,2)*exp(-a[9]*rho))*(83.14472*T);

        //% Convert P from bars (equation P&S) to pascals (model)
        return  pb/1e-5; //%1 pascal a 1e-5 bars
        }
        return 0;
    }


//%==========================================================================
//%Return the pressure-temperature-time functions
//%==========================================================================
double PtFun(std::valarray<double> Pin,std::valarray<double> tin,double t){
    if (t==0) return Pin[0];
    if (t>tin.max()) return Pin[Pin.size()-1];
    for (int tau =0; tau< tin.size(); tau++){
        if (t>tin[tau]);
        else{
            double lever = (t-tin[tau-1])/(tin[tau]-tin[tau-1]);
            double a = Pin[tau-1];
            double b=(Pin[tau]-Pin[tau-1]);
            double outPres = a+lever*b;
            return outPres;}      
    }
    return 0;
}

double TtFun(std::valarray<double> Tin,std::valarray<double> tin,double t){
    //Need to add Newton cooling option
    if (t==0) return Tin[0];
    if (t>tin.max()) return Tin[Tin.size()-1];
    for (int tau =0; tau< tin.size(); tau++){
        if (t>tin[tau]);
        else{
            double lever = (t-tin[tau-1])/(tin[tau]-tin[tau-1]);
            double a = Tin[tau-1];
            double b=(Tin[tau]-Tin[tau-1]);
            double outTemp = a+lever*b;
            return outTemp;}      
    }
    return 0;
} 

/*%==========================================================================
%==========================================================================
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%Below are the functions that can be called by the model to perform
%Calculations. There are functions for Viscosity, Solubility, Diffusivity,
%The inital mass of gas (equation of state), and pressure in the bubble
%(equation of state)
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%==========================================================================
%==========================================================================




%==========================================================================
%Viscosity functions
%==========================================================================*/


std::valarray<double> Giordano_2008_visc(std::valarray<double> H2Ot, double T, std::valarray<double> Composition){
    //Get the Giordano2008 VFT parameters from their script
    //Giordano2008_Model.mat
    std::valarray<std::valarray<double>> outvalues=Giordano2008_Model(H2Ot, Composition);
    std::valarray<double> At, Bt, Ct = std::valarray<double>(H2Ot.size());
    At = outvalues[0];
    Bt = outvalues[1];
    Ct = outvalues[2];
    return pow(10,(At+(Bt/(T-Ct))));}



/*%==========================================================================
%==========================================================================
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%Below are the additional functions that are used by the various model
%functions for performing calculations
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%==========================================================================
%==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The below functions are used for solving equation of
%state via Pitzer & Sterner, 1994
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


//%Pitzer and Sterner, 1994 gas density function, calls PS_myfun
double density(double P, double T){
    double rho =0;
    static double b[10][6] = {0};
    b[0][2]=0.24657688e6; 
    b[0][3]=0.51359951e2 ;
    b[1][2]=0.58638965e0 ;
    b[1][3]=-0.28646939e-2; 
    b[1][4]=0.31375577e-4 ;
    b[2][2]=-0.62783840e1 ;
    b[2][3]=0.14791599e-1 ;
    b[2][4]=0.35779579e-3 ;
    b[2][5]=0.15432925e-7 ;
    b[3][3]=-0.42719875e0 ;
    b[3][4]=-0.16325155e-4 ;
    b[4][2]=0.56654978e4 ;
    b[4][3]=-0.16580167e2 ;
    b[4][4]=0.76560762e-1 ;
    b[5][3]=0.10917883e0 ;
    b[6][0]=0.38878656e13 ;
    b[6][1]=-0.13494878e9 ;
    b[6][2]=0.30916564e6 ;
    b[6][3]=0.75591105e1 ;
    b[7][2]=-0.65537898e5 ;
    b[7][3]=0.18810675e3 ;
    b[8][0]=-0.14182435e14; 
    b[8][1]=0.18165390e9 ;
    b[8][2]=-0.19769068e6 ;
    b[8][3]=-0.23530318e2 ;
    b[9][2]=0.92093375e5 ;
    b[9][3]=0.12246777e3 ;


//% convert P to bars from input (pascals)
    P = P*1e-5; //%1 pascal a 1e-5 bars

    double a[10];
    for (int i=0; i<10; i++){
        a[i]=b[i][0]*pow(T,-4) + b[i][1]*pow(T,-2) + b[i][2]*pow(T,-1) + b[i][3] + b[i][4]*T + b[i][5]*pow(T,2);
    }
    //% PRT = P/RT where P [bars], R [cm^3*bar/K/mol] and T [K]
    double PRT = P/(83.14472*T);
    //% solve implicit equation for rho and convert to kg/m^3

    std::pair<double, double> res ;
    double lower = 0;
    double upper = 100;
    res =  tools::bisect([a, PRT](double rho){return (rho+a[0]*pow(rho,2)-pow(rho,2)*((a[2]+2*a[3]*rho+3*a[4]*pow(rho,2)+4*a[5]*pow(rho,3))/pow((a[1]+a[2]*rho+a[3]*pow(rho,2)+a[4]*pow(rho,3)+a[5]*pow(rho,4)),2))+a[6]*pow(rho,2)*exp(-a[7]*rho)+a[8]*pow(rho,2)*exp(-a[9]*rho)) - PRT;}, lower, upper, [](double l, double r){return abs(l-r) < 1e-8;});
    rho = (res.first +res.second)/2.0;
    //fzero(@PS_myfun,0.001,[],a,PRT)*18.01528*1000;

    return rho*18.01528*1000;
}


//% the function from Pitzer & Sterner 1994, which takes the matrix of
//% coefficients a and P/RT as arguments; rho is a first guess for the
//% density [g/mol]
double PS_myfun(double rho, double a[10], double PRT){
return (rho+a[0]*pow(rho,2)-pow(rho,2)*((a[2]+2*a[3]*rho+3*a[4]*pow(rho,2)+4*a[5]*pow(rho,3))/pow((a[1]+a[2]*rho+a[3]*pow(rho,2)+a[4]*pow(rho,3)+a[5]*pow(rho,4)),2))+a[6]*pow(rho,2)*exp(-a[7]*rho)+a[8]*pow(rho,2)*exp(-a[9]*rho)) - PRT;}



std::valarray<std::valarray<double>> Giordano2008_Model(std::valarray<double> H2Ot, std::valarray<double> Composition){

    double AT      = -4.55;
    std::valarray<double> bb  = {159.56, -173.34, 72.13, 75.69, -38.98, -84.08, 141.54, -2.43, -0.91, 17.62};
    std::valarray<double> cc  = {2.75, 15.72, 8.32, 10.2, -12.29, -99.54, 0.3 };

    int row = H2Ot.size();
    std::valarray<std::valarray<double>> Comp_rep = std::valarray<std::valarray<double>>(Composition,row);
    
    for (int i = 0; i<row; i++)
        Comp_rep[i][10]=H2Ot[i];

    std::valarray<std::valarray<double>> wtm=Comp_rep;

    int nc = 12;
    int nx = row;

    std::valarray<std::valarray<double>> bb_rep = std::valarray<std::valarray<double>>(bb,nx);
    std::valarray<std::valarray<double>> cc_rep = std::valarray<std::valarray<double>>(cc,nx);
    std::valarray<double> AT_rep = std::valarray<double>(AT,nx);

    //% Function molefrac_grd: converts wt % oxide bais to mole % oxide basis
    //%[ncomps xmf_t] = molefrac_grd(wtm);
    std::valarray<std::valarray<double>> xmf_t = molepct_grd(wtm);
                
    //%outvalues=[];
    //% Load composition-basis matrix for multiplication against model-coefficients
    //% Result is two matrices bcf[nx by 10] and ccf[nx by 7]
    std::valarray<double> b1=std::valarray<double>(0.0,nx),b2=std::valarray<double>(0.0,nx),b3=std::valarray<double>(0.0,nx),b4=std::valarray<double>(0.0,nx),b5=std::valarray<double>(0.0,nx),b6=std::valarray<double>(0.0,nx),b7=std::valarray<double>(0.0,nx),b12=std::valarray<double>(0.0,nx),b13=std::valarray<double>(0.0,nx),b14 = std::valarray<double>(0.0,nx);
    std::valarray<double> c1= std::valarray<double>(0.0,nx),c2=std::valarray<double>(0.0,nx),c3=std::valarray<double>(0.0,nx),c4=std::valarray<double>(0.0,nx),c5=std::valarray<double>(0.0,nx),c6=std::valarray<double>(0.0,nx),c11 = std::valarray<double>(0.0,nx);
    std::valarray<double> siti= std::valarray<double>(0.0,nx), tial= std::valarray<double>(0.0,nx),fmm= std::valarray<double>(0.0,nx),nak = std::valarray<double>(0.0,nx);

    for (int i = 0; i<nx; i++){
    siti[i]    =   xmf_t[i][0] + xmf_t[i][1] ;
    tial[i]    =   xmf_t[i][1]+xmf_t[i][2] ;
    fmm[i]     =   xmf_t[i][3] + xmf_t[i][4] + xmf_t[i][5] ;
    nak[i]     =   xmf_t[i][7] + xmf_t[i][8] ;
    b1[i]  =   siti[i];
    b2[i]  =   xmf_t[i][2] ;
    b3[i]  =   xmf_t[i][3] + xmf_t[i][4] + xmf_t[i][9] ;
    b4[i]  =   xmf_t[i][5] ;
    b5[i]  =   xmf_t[i][6] ;

    //#print(1+ xmf_t[1,10])

    b6[i]  =   xmf_t[i][7] + xmf_t[i][10] + xmf_t[i][11] ;
    b7[i]  =   xmf_t[i][10] + xmf_t[i][11] + log(1+xmf_t[i][10]) ;
    b12[i] =   siti[i]*fmm[i] ;
    b13[i] =   (siti[i] + xmf_t[i][2] + xmf_t[i][9])*( nak[i] + xmf_t[i][10] );
    b14[i] =   xmf_t[i][2]*nak[i];

    c1[i]      =   xmf_t[i][0];
    c2[i]      =   tial[i];
    c3[i]      =   fmm[i];
    c4[i]      =   xmf_t[i][6];
    c5[i]      =   nak[i];
    c6[i]      =   log(1+xmf_t[i][10] + xmf_t[i][11]);
    c11[i]     =   xmf_t[i][2] + fmm[i] + xmf_t[i][6] - xmf_t[i][9];
    c11[i]     =   c11[i]*(nak[i] + xmf_t[i][10] + xmf_t[i][11]);

    }
    std::valarray<std::valarray<double>> bcf      =   bb_rep;
    std::valarray<std::valarray<double>> ccf      =   cc_rep;

    for (int i = 0; i<bb_rep.size(); i++){
        bcf[i]={b1[i], b2[i], b3[i], b4[i], b5[i], b6[i], b7[i], b12[i], b13[i], b14[i]};
        ccf[i]={c1[i], c2[i], c3[i], c4[i], c5[i], c6[i], c11[i]}; 
        
    }
        
    
    std::valarray<std::valarray<double>> bt,ct;
    bt = bb_rep*bcf;    
    ct = cc_rep*ccf;

    std::valarray<double> BT= std::valarray<double>(nx),CT= std::valarray<double>(nx), TG= std::valarray<double>(nx),F = std::valarray<double>(nx);
    for (int i = 0; i<nx; i++){
        BT[i]          = bt[i].sum();
        CT[i]          = ct[i].sum();
    }   
    TG          = BT/(12-AT_rep) + CT;
    F           = BT/(TG*(1 - CT/TG)*(1 - CT/TG));

    std::valarray<std::valarray<double>> outvalues = {AT_rep, BT, CT, TG, F};
    return outvalues;

}





//%Function for computing the mole percent
std::valarray<std::valarray<double>> molepct_grd(std::valarray<std::valarray<double>> wtm){
    int nr =  wtm.size();

    double wtmSum = wtm[0].sum();
    //std::valarray<double> wtmSumV = std::valarray<double>(wtm[0].sum(),nr);
    //% [n x]=MOLEPCT(X) Calculates mole percent oxides from wt % 
    //% 1 SiO2 2 TiO2  3 Al2O3  4 FeO  5 MnO  6 MgO  7CaO 8 Na2O  9 K2O  10 P2O5 11 H2O 12 F2O-1  
    //% Output: mole fractions of equivalent
    std::valarray<double> mw={60.0843, 79.8658, 101.961276, 71.8444, 70.937449,40.3044,56.0774, 61.97894, 94.1960, 141.9446,18.01528, 18.9984}; // Molecular weights of different oxides
    std::valarray<std::valarray<double>> mp, xmf,mpv;

    //%Replicate the molecular weight row by the number of spatial nodes used
    //%(number of data rows)
    std::valarray<std::valarray<double>> mw_rep = std::valarray<std::valarray<double>>(mw,nr); //reproduce molecular weights for each diffusion node
    std::valarray<std::valarray<double>> wtn = mw_rep;
    for (int j=0; j<nr; j++){
        for (int i =0; i<12; i++){
            if (i==10)
                wtn[j][i]=wtm[j][10];
            else
                wtn[j][i]=(100-wtm[j][10])*wtm[j][i]/(wtmSum);
        }
    }

    mp=wtn/mw_rep;
    std::valarray<double> div = std::valarray<double>(nr);
    for (int i=0; i<nr; i++){
        div[i]=mp[i].sum()/100.0;
    }

    mpv= (mp/div);
    xmf=mpv;

    return xmf;
    }


#include <cmath>
#include <iostream>
#include <vector>
#include <valarray>
#include <string>


#include "boost/numeric/odeint.hpp"
#include "getFunctions_v2.h"

typedef std::valarray<double> state_type;
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;


struct push_back_state_and_time
{
    std::vector< vector_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< vector_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const vector_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );

        std::cout<<"t = "<<t<<"\t r = "<<x[x.size()-1]<<"\n";
    }
};

std::vector<double> outputTime={};
std::vector<vector_type> outputY={};

std::valarray<double> xG,xBG,T_0G,t_TG,P_0G,t_PG,CompositionG;
double m_0G,melt_RhoG,H2Ot_0G,R_0G,WG,SurfTensG,NbG;
std::string SolModelG,EOSModelG,DiffModelG,ViscModelG;

std::valarray<std::valarray<double>> outer(std::valarray<double> first,std::valarray<double>second);
std::valarray<double> diff(std::valarray<double> input);
//double trapz(auto function, std::valarray<double> x);
std::valarray<double> logspace(double lower, double upper, int nodes);
std::valarray<double> Porosity(double Nb,std::valarray<double> R);
double Porosity(double Nb,double R);
double Mass_SingleOxygen(std::valarray<double> Composition);
void Outputs(double Nodes,double R_0,double L,std::valarray<std::valarray<double>> Y,std::valarray<double> t,double Nb,std::valarray<double> T_0,std::valarray<double> t_T,std::valarray<double> P_0,std::valarray<double> t_P, 
            std::valarray<double> &R, std::valarray<double> &phi, std::valarray<double> &P, std::valarray<double> &T, std::valarray<std::valarray<double>> &x_out ,std::valarray<std::valarray<double>> &H2Ot_all);
void MYodeFun(const vector_type  &Y, vector_type  &dYdt, double t,std::valarray<double> x,std::valarray<double> xB,double m_0,double melt_Rho,std::valarray<double> T_0,std::valarray<double> t_T,std::valarray<double> P_0,std::valarray<double> t_P,double H2Ot_0,double R_0,double W,double SurfTens,std::string SolModel,std::string DiffModel,std::string ViscModel,std::string EOSModel,std::valarray<double> Composition,double Nb);
void myObserver(const std::valarray<double>  &x, const double t);
//void ODEFunction(TS ts, PetscReal t, Vec U, Vec U_t, Vec r, void *ctx);  
void JacobianFunction(const vector_type &x, matrix_type &J, const double &t, vector_type &dxdt);

void setValues(std::valarray<double> x_in,std::valarray<double> xB_in,double m_0_in,double melt_Rho_in,std::valarray<double> T_0_in,std::valarray<double> t_T_in,std::valarray<double> P_0_in,std::valarray<double> t_P_in,double H2Ot_0_in,double R_0_in,double W_in,double SurfTens_in,std::string SolModel_in,std::string DiffModel_in,std::string ViscModel_in,std::string EOSModel_in,std::valarray<double> Composition_in,double Nb_in){
    xG=x_in;
    xBG=xB_in;
    T_0G=T_0_in;
    t_TG=t_T_in;
    P_0G=P_0_in;
    t_PG=t_P_in;
    CompositionG=Composition_in;

    m_0G=m_0_in;
    melt_RhoG=melt_Rho_in;
    H2Ot_0G=H2Ot_0_in;
    R_0G=R_0_in;
    WG=W_in;
    SurfTensG=SurfTens_in;
    NbG=Nb_in;

    SolModelG=SolModel_in;
    EOSModelG=EOSModel_in;
    DiffModelG=DiffModel_in;
    ViscModelG=ViscModel_in;
}


void solveSys( const vector_type &Y, vector_type &dYdt, const double t)
{
    MYodeFun(Y,dYdt,t,xG,xBG,m_0G,melt_RhoG,T_0G,t_TG,P_0G,t_PG,H2Ot_0G,R_0G,WG,SurfTensG,SolModelG, DiffModelG, ViscModelG, EOSModelG,CompositionG,NbG);
}




void Numerical_Model_v2(std::valarray<double> Composition,std::string SolModel,std::string DiffModel,std::string ViscModel,std::string EOSModel,double SurfTens, double melt_Rho, int Nodes, double R_0,double H2Ot_0,double Nb,double t_nuc,std::valarray<double> T_0,std::valarray<double> t_T,std::valarray<double> P_0,std::valarray<double> t_P,std::valarray<double> Numerical_Tolerance,
std::valarray<double> &t, std::valarray<double> &R, std::valarray<double> &phi, std::valarray<double> &P, std::valarray<double> &T, std::valarray<std::valarray<double>> &x_out ,std::valarray<std::valarray<double>> &H2Ot_all){

    std::cout<<"Starting Numerical Model \n\n";

//#%Get the molar mass of anhydrous  melt on a single oxygen basis
    double W = Mass_SingleOxygen(Composition);

//#%Get the initial porisity of the shell model
    double phi_0 = Porosity(Nb,R_0);

//# distance to shell edge (m) (computed from R and phi)
    double L = R_0*(pow(phi_0,(-1.0/3.0)) - 1); 

//#Set the numerical duration bounds (from nucleation to end)
    double t_f = 3000;
    //tspan = (t_nuc, t_f)

//# Get the P and T at the defined t_0 using the P-T-t pathway
    //#PT_t0 = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t_nuc)
    //P_t0 = gFP.PtFun(P_0,t_P,t_nuc)
    double P_t0 = P_0[0];
    //T_t0 = gFP.TtFun(T_0,t_T,t_nuc)
    double T_t0 = T_0[0];
    

//#%Compute the initial mass of gas in the bubble
    double m_0 = m0_fun(R_0, P_t0 +(2*SurfTens/R_0), T_t0, EOSModel);

//#generate logarithmically spaced xB values
//#scale to appropriate length, and shift so that xB
//#starts at the bubble wall
//#(equation A.6 from manuscript)
    std::valarray<double> xB = std::valarray<double>(Nodes+1);
    xB[0]=R_0;
    std::valarray<double> spacing = logspace(-2,0,Nodes)*L+R_0;
    for (int i = 0; i<Nodes; i++){
        xB[i+1]=spacing[i];
    }

//#%block-centered grid - use midpoints of xB
//#%(equation A.3 from manuscript)
    std::valarray<double> x = std::valarray<double>(Nodes);
    for (int i =0 ; i<Nodes; i++)
        x[i]=(xB[i+1]+xB[i])/2;


//#%uniform concentration; system in equilibrium Although this can be
//#%changed to represent non-uniform conditions. Any initial concentration
//#%profile can be set by the user.
    std::valarray<double> H2Ot = std::valarray<double>(H2Ot_0,Nodes);


//#%Y_0 is a vector containing quantities to be solved for by ode15s
    vector_type Y_0 = vector_type(Nodes+1);
    for (int i=0;i<H2Ot.size();i++){
        Y_0[i]=H2Ot_0;
    }
    Y_0[H2Ot.size()]=R_0;
    std::cout<<"Initial radius "<<R_0<<"\n";   

    setValues(x,xB,m_0,melt_Rho,T_0,t_T,P_0,t_P,H2Ot_0,R_0,W,SurfTens,SolModel, DiffModel, ViscModel, EOSModel,Composition,Nb);

    std::cout<<"Start solving \n\n";
    double initial_time = 0;
    double end_time = t_f;
    double dt = 1e-4;


    typedef boost::numeric::odeint::rosenbrock4<double> stepper_type;
    typedef boost::numeric::odeint::rosenbrock4_controller< stepper_type > controlled_stepper_type;
    typedef boost::numeric::odeint::rosenbrock4_dense_output< controlled_stepper_type > dense_output_type;

    std::vector<vector_type> x_vec;
    std::vector<double> times;
    double error = 1e-4;

    size_t steps = boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::rosenbrock4_controller<stepper_type>(error,error),make_pair(solveSys,JacobianFunction), Y_0, 0.0,t_f, dt, push_back_state_and_time( x_vec , times ) );



    // Process or save the solution as needed

    // Clean up


    t = std::valarray<double>(outputTime.size());
    std::valarray<std::valarray<double>> Y = std::valarray<std::valarray<double>>(std::valarray<double>(Y_0.size()),t.size());
for( size_t i=0; i<=steps; i++ )
{
    std::cout << times[i] << '\t' << x_vec[i][0] << '\t' << x_vec[i][1] << '\n';
}

    Outputs(Nodes,R_0,L,Y,t,Nb,T_0,t_T,P_0,t_P,R, phi, P, T, x_out, H2Ot_all);
    return;
}

void JacobianFunction(const vector_type &u, matrix_type &K, const double &t, vector_type &udot) {
    // Extract the array of state values
//std::cout<<"Start Jacobian \n";
    double shift = 1e-5;
    int N = u.size();

//std::cout<<"Vectors are length"<<N<<"\n";
    matrix_type J = K;
    

    vector_type X = u;
    vector_type Y = u;
    vector_type dXdt = udot;
    vector_type dYdt = udot;

//std::cout<<"Jacobian calc 1 at t="<<t<<"\n";
    MYodeFun(Y,dYdt,t,xG,xBG,m_0G,melt_RhoG,T_0G,t_TG,P_0G,t_PG,H2Ot_0G,R_0G,WG,SurfTensG,SolModelG, DiffModelG, ViscModelG, EOSModelG,CompositionG,NbG);

    for (int i =0; i<N-1; i++){
        if(i%3==0)
            X[i]=X[i]+shift;
    }
//std::cout<<"Jacobian calc 2 at t="<<t<<"\n";
    MYodeFun(X,dXdt,t,xG,xBG,m_0G,melt_RhoG,T_0G,t_TG,P_0G,t_PG,H2Ot_0G,R_0G,WG,SurfTensG,SolModelG, DiffModelG, ViscModelG, EOSModelG,CompositionG,NbG);

    J(0,0)=(dXdt[0]-dYdt[0])/shift;
    J(0,1)=(dXdt[1]-dYdt[1])/shift;

    for (int i =3; i<N-1; i=i+3){
        J(i,i)= (dXdt[i]-dYdt[i])/shift;
        J(i,i-1)= (dXdt[i-1]-Y[i-1])/shift;
        J(i,i+1)= (dXdt[i+1]-Y[i+1])/shift;
    }
    X=Y;
    for (int i =0; i<N-1; i++){
            X[i]=Y[i]+shift*(i%3==1);
    }
//std::cout<<"Jacobian calc 3 at t="<<t<<"\n";    
MYodeFun(X,dXdt,t,xG,xBG,m_0G,melt_RhoG,T_0G,t_TG,P_0G,t_PG,H2Ot_0G,R_0G,WG,SurfTensG,SolModelG, DiffModelG, ViscModelG, EOSModelG,CompositionG,NbG);

    for (int i =1; i<N-1; i=i+3){
        J(i,i)= (dXdt[i]-dYdt[i])/shift;
        J(i,i-1)= (dXdt[i-1]-Y[i-1])/shift;
        J(i,i+1)= (dXdt[i+1]-Y[i+1])/shift;
    }
    X=Y;
    for (int i =1; i<N-1; i++){
            X[i]=Y[i]+shift*(i%3==2);
    }
//std::cout<<"Jacobian calc 4 at t="<<t<<"\n";
MYodeFun(X,dXdt,t,xG,xBG,m_0G,melt_RhoG,T_0G,t_TG,P_0G,t_PG,H2Ot_0G,R_0G,WG,SurfTensG,SolModelG, DiffModelG, ViscModelG, EOSModelG,CompositionG,NbG);

    for (int i =2; i<N-1; i=i+3){
        J(i,i)= (dXdt[i]-dYdt[i])/shift;
        J(i,i-1)= (dXdt[i-1]-Y[i-1])/shift;
        J(i,i+1)= (dXdt[i+1]-Y[i+1])/shift;
    }

    X=Y;

    shift = 1e-6;
    X[N-1]=Y[N-1]+shift;
//std::cout<<"Jacobian calc 5 at t="<<t<<"\n";
    MYodeFun(X,dXdt,t,xG,xBG,m_0G,melt_RhoG,T_0G,t_TG,P_0G,t_PG,H2Ot_0G,R_0G,WG,SurfTensG,SolModelG, DiffModelG, ViscModelG, EOSModelG,CompositionG,NbG);


    for (int i =0; i<N; i++)
    J(i,N-1)=(dXdt[i]-Y[i])/shift;
//std::cout<<"Calculated Jacobian \n";
    K=J;
    udot = dYdt;
    //std::cout<<"Done Jacobian \n";
    return;

}

//#function [dYdt,pb] =  MYodeFun(t,Y,x,xB,m_0,melt_Rho,T_0,P_0,H2Ot_0,R_0,W,SurfTens,SolFun,DiffFun,ViscFun,pb_fun,PTt_fun,Composition,Nb,t_f,T_f, dTdt, P_f, dPdt)

void MYodeFun(const vector_type &X, vector_type  &dXdt, double t,std::valarray<double> x,std::valarray<double> xB,double m_0,double melt_Rho,std::valarray<double> T_0,std::valarray<double> t_T,std::valarray<double> P_0,std::valarray<double> t_P,double H2Ot_0,double R_0,double W,double SurfTens,std::string SolModel,std::string DiffModel,std::string ViscModel,std::string EOSModel,std::valarray<double> Composition,double Nb){


    std::valarray<double> Y= std::valarray<double>(0.0,X.size());
    std::valarray<double> dYdt= std::valarray<double>(0.0,X.size());
    for (int i = 0; i<X.size();i++){
        Y[i]=X[i];
        //std::cout<<"Y is "<<Y[i]<<"\n";
    }

//#%extract individual concentrations
    int nx = (Y.size() - 1);
    std::valarray<double> H2Ot = std::valarray<double>(0.0,nx);
    for (int i = 0 ; i< nx; i++){
            H2Ot[i]=Y[i];
    }
    double R = Y[nx];

//#%Get current temperature and pressure
   // #PT = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t)
    double P = PtFun(P_0,t_P,t);
    double T = TtFun(T_0,t_T,t);

//#%Trapz is a function that uses the trapezoidal rule to numerically
//#%integrate. Dim is the dimension to operate on
    int dim = 1;

//#%Get mass of gas (Equation 7 from the main manuscript)
    std::valarray<double> i1=(1.0/100.0)*H2Ot_0*pow(x,2);
    auto I1fun = [H2Ot_0,x ](int i){return (1.0/100)*H2Ot_0*pow(x[i],2);};
        double I1=0;
    for (int i =1; i<x.size(); i++){
        I1=I1+ (I1fun(i)+I1fun(i-1))*((x[i]-x[i-1]))*0.5;
    }


    auto I2fun = [H2Ot, x](int i)->double {return (1.0/100)*H2Ot[i]*pow(x[i],2);};
    double I2=0;
    for (int i =1; i<x.size(); i++){
        I2=I2+ (I2fun(i)+I2fun(i-1))*((x[i]-x[i-1]))*0.5;
    }


    double m=m_0+4*PI*melt_Rho*(I1-I2);
   // std::cout<<"\n\n m is "<<m<<"\n\n";

//#%Compute pressure of the gas  in the bubble from EOS
//#%(section 2.3.2 of main manuscript)
    double pb = pb_fun(m, T, R, EOSModel);

//#%====boundary conditions====
//#%Determine the solubility condition of water in the system based on gas
//#%pressure in the bubble
//#%(section 2.3.1 of main manuscript)
    double H2Oeq = SolFun(T,pb,SolModel);

//#%get diffusion coefficients at boundaries
//#%Creates a vector where the first value is the boundary condition as
//#%determined from the equilibrium solubility of the system.
//#%(section 2.3.1 of main manuscript)
    std::valarray<double> H2O_BC= std::valarray<double>(H2Ot.size());
    H2O_BC[0]=H2Oeq;
    for (int i=1; i<H2Ot.size(); i++){
        H2O_BC[i]= (H2Ot[i]+H2Ot[i-1])/2.0;
    }

    std::valarray<double> DH2Ot = DiffFun(T,P,H2O_BC,W,DiffModel);

//#%====Solve water diffusion==== (equation A.4 from manuscript)
//#%molecular diffusion flux

    std::valarray<double> diffTop =std::valarray<double>(0.0,H2Ot.size()+1), diffBottom =std::valarray<double>(0.0,H2Ot.size()+1);
    diffTop[0] = H2Oeq;
    diffBottom[0]=xB[0];
    for (int i =1; i<diffTop.size(); i++){
        diffTop[i]=H2Ot[i-1];
        diffBottom[i]=x[i-1];
    }

    std::valarray<double> JH2Ot = -DH2Ot*diff(diffTop)/diff(diffBottom);


//#%Gradient of the diffusion flux.

    std::valarray<double> dJdiff = std::valarray<double>(JH2Ot.size()+1);
    for (int i=0; i<JH2Ot.size(); i++)
        dJdiff[i]=(((pow((pow(x,3)+pow(R,3)-pow(R_0,3)),(4.0/3.0)))/(pow(x,2)))*JH2Ot)[i];
    dJdiff[JH2Ot.size()]=0;

    std::valarray<double> dJH2Odx = (1.0/(pow(x,2.0)))*diff(dJdiff)/diff(xB);

//#%====solve hydrodynamic equation====
//#%Compute the viscosity
    std::valarray<double> v = ViscFun(H2Ot,T,Composition,ViscModel);

//#%Compute integrated viscosity (equation A.5 from manuscript)
    auto I3fun = [v, R, R_0,x](int i){return (v[i]*pow(x[i],2))/(pow((pow(R,3)-pow(R_0,3)+pow(x[i],3)),2));};
    //double I3=trapz(I3fun,x);

    double I3=0;



    for (int i =1; i<x.size(); i++){
        I3+= (I3fun(i)+I3fun(i-1))*((x[i]-x[i-1]))*0.5;
    }

//#%Solve Rayleigh-Plesset equation (equation A.6 from manuscript)
    double dRdt= ((pb-P-(2*(SurfTens)/R))/(12*pow(R,2)))/I3;

//#%return rhs of ode
    dYdt = std::valarray<double>(dJH2Odx.size()+1);
    for (int i =0; i< dJH2Odx.size(); i++){
        dYdt[i] = -dJH2Odx[i];
        dXdt[i] = -dJH2Odx[i];
    }
    dYdt[dJH2Odx.size()]=dRdt;
    dXdt[dJH2Odx.size()]=dRdt;


    //std::cout<<"pb is "<<pb<<"\n"; 
    //std::cout<<"dR/dt is "<<dRdt<<"\n"; 

}




//##%==========================================================================
//##%Functions to get outputs from the ODE solution
//##%==========================================================================


//#function [R, phi, P, T, x_out, H2Ot_all] = Outputs(Nodes,R_0,L,Y,t,Nb,PTt_fun,P_0, P_f, dPdt,T_0,T_f,dTdt)

void Outputs(double Nodes,double R_0,double L,std::valarray<std::valarray<double>> Y,std::valarray<double> t,double Nb,std::valarray<double> T_0,std::valarray<double> t_T,std::valarray<double> P_0,std::valarray<double> t_P, 
            std::valarray<double> &R, std::valarray<double> &phi, std::valarray<double> &P, std::valarray<double> &T, std::valarray<std::valarray<double>> &x_out ,std::valarray<std::valarray<double>> &H2Ot_all){
    //pass ouputs by reference - convert fuction to void;
    //{R, phi, P, T, x_out, H2Ot_all}
    //#%Get the bubble radius and phi
    R = std::valarray<double>(t.size());
    for (int i = 0 ; i< t.size(); i++){
        R[i] = Y[i][Y[0].size()-1];
    }
    Y[Y.size() -1];
    phi = Porosity(Nb,R);

    //#%Get the shell thickness for each R(t)
    std::valarray<double> L_all = (pow((pow((R_0+L),3)+pow(R,3)-pow(R_0,3)),(1/3))-R);
    //#%Get the x position for each output shell thickness
    x_out=outer(logspace(-2,0,Nodes),L_all);

    //#%Get all of the water profiles

    H2Ot_all = std::valarray<std::valarray<double>>(std::valarray<double>(Y.size()-1),t.size());
    for (int i =0; i<t.size(); i++){
        for (int j=0; j<Y.size()-1; j++)
            H2Ot_all[i][j]=Y[i][j];
    }

    //#%Get P-T-t history
    //#PT = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t)

    P=std::valarray<double>(t.size());
    T=std::valarray<double>(t.size());
    for (int i=0; i<t.size();i++){ 
        P[i] = PtFun(P_0,t_P,t[i]);
        T[i] = TtFun(T_0,t_T,t[i]);
    }
    //Function redefined as a void to update values passed by reference
    //std::valarray<std::valarray<double>> output = {R, phi, P, T, x_out, H2Ot_all};
    //return output;
    return;
}

//#%==========================================================================
//#%Functions that return constants which are used in the modelling
//#%==========================================================================

// #%This function returns the molar mass of the anhydrous melt on a single
// #%oxygen basis
// #function W = Mass_SingleOxygen(Composition)

double Mass_SingleOxygen(std::valarray<double> Composition){
    std::valarray<double> comp = Composition;

    //"""%Convert composition matrix from Viscosity input to Shishkina format
    //%!! SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1 !!
    //%Convert To
    //%!! SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 Cr2O3] !!"""

    std::valarray<double> X = std::valarray<double>(12);
    X[0] = comp[0];
    X[2] = comp[1];
    X[2] = comp[2];
    X[3] = 0;
    X[4] = comp[3];
    X[5] = comp[4];
    X[6] = comp[5];
    X[7] = comp[6];
    X[8] = comp[7];
    X[9] = comp[8];
    X[10] = comp[9];
    X[11] = 0;


    double Total_Mass = X.sum();

    //#%Molar mass (g/mol) of individual elements
    double mSi = 28.0855;
    double mTi = 47.867;
    double mAl = 26.981539;
    double mFe = 55.845;
    double mMn = 54.938044;
    double mMg = 24.305;
    double mCa = 40.078;
    double mNa = 22.989769;
    double mK = 39.0983;
    double mP = 30.973762;
    double mCr = 51.9961;
    double mO = 15.999;

    //#% [SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 Cr2O3]
    //#%Molar mass (g/mol) of oxides
    std::valarray<double> OxideMolarMass = std::valarray<double>(12);
    OxideMolarMass[0] =  (mSi+2*mO);
    OxideMolarMass[1] = (mTi+2*mO);
    OxideMolarMass[2] = (2*mAl+3*mO);
    OxideMolarMass[3] = (2*mFe+3*mO);
    OxideMolarMass[4] = (1*mFe+1*mO);
    OxideMolarMass[5] = (1*mMn+1*mO);
    OxideMolarMass[6] = (1*mMg+1*mO);
    OxideMolarMass[7] = (1*mCa+1*mO);
    OxideMolarMass[8] = (2*mNa+1*mO);
    OxideMolarMass[9] = (2*mK+1*mO);
    OxideMolarMass[10] = (2*mP+5*mO);
    OxideMolarMass[11] = (2*mCr+3*mO);

    //#%Compute number of moles of element, and Cation Fraction
    std::valarray<double> numMolesOxygen = {2.0, 2.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 3.0};
    std::valarray<double> numMolesElement = {1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0};

    //#%Compute the number of moles of each oxide
    std::valarray<double> Moles_Oxide = X/OxideMolarMass;

    //#%Compute moles of oxygen by stoichiometry
    std::valarray<double> Moles_Oxygen = Moles_Oxide*numMolesOxygen;

    //#%W_melt is the mass of anhydrous melt per mole of oxygen
    return Total_Mass/(Moles_Oxygen.sum());
}

//#%Compute the gas volume fraction from Nb and radius
double Porosity(double Nb,double R){
    double meltvolume = 1.0/Nb;
    double gasvolume = (4.0/3.0)*PI*pow(R,3);
    return gasvolume/(gasvolume+meltvolume);
}

std::valarray<double> Porosity(double Nb,std::valarray<double> R){
    double meltvolume = 1.0/Nb;
    std::valarray<double> gasvolume = (4.0/3.0)*PI*pow(R,3);
    return gasvolume/(gasvolume+meltvolume);
}



std::valarray<double> logspace(double lower, double upper, int nodes){
    std::valarray<double> output = std::valarray<double>(nodes);
    std::valarray<double> regular = std::valarray<double>(nodes);
    for( int i = 0 ; i< nodes; i++){
        regular[i] = lower + (upper-lower)*i/(nodes-1);
    }
    output = pow(10,regular);
    return output;

}


//return first difference 
std::valarray<double> diff(std::valarray<double> input){
    std::valarray<double> output = std::valarray<double>(input.size()-1);
    for (int i = 0; i<output.size(); i++){
        output[i]=input[i+1] -input[i];
    }
    return output;
}

std::valarray<std::valarray<double>> outer(std::valarray<double> first,std::valarray<double>second){
    std::valarray<std::valarray<double>> output = std::valarray<std::valarray<double>>(std::valarray<double>(0.0,second.size()),first.size());
    for (int i =0; i<first.size(); i++){
        for(int j = 0 ; j< second.size(); j++){
            output[i][j]=first[i]*second[j];
        }
    }
    return output;
}

void myObserver(const vector_type &x, const double t){
    outputTime.push_back(t);
    outputY.push_back(x);

    std::cout<<"\n"<<t<<"\t"<<x[x.size()-1]<<"\n";
}

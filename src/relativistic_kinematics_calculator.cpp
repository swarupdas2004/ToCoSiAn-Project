#include <iostream>
#include <cmath>
#include <iomanip>
//#include <numbers>__update to C++(20+)
using namespace std;

const double pi=acos(-1.0);
//const double PI=numbers::pi;

int main(){
    double px,py,pz,m;
    double p,E,eta,theta,phi,pt;
    int c;
    c=0;

    cout<<"Enter the mass of the particle (in MeV) : ";
    cin>>m;
    if(m<0){cerr<<"Error! Mass can't be negative.";}

    cout<<"Enter its momentum in x,y and z-direction (in MeV) : ";
    cin>>px>>py>>pz;

    p=sqrt(px*px+py*py+pz*pz);
    E=sqrt(p*p+m*m);
    pt=sqrt(px*px+py*py);

    phi=atan(py/px);
    if(p>1e-10 && (p-abs(pt))>1e-10)   {theta=acos(pz/p);
        eta=-log(tan(theta/2.0));}
    else    {c=1;}
    cout<<"RESULTS-----------\n";
    cout<<"======Cartesian======\n";
    cout<<"3-Momentum components (x,y,z in order) = "<<px<<" MeV"<<py<<" MeV"<<pz<<" MeV\n";
    cout<<"Invariant mass = "<<m<<" MeV\n";
    cout<<"Energy = "<<E<<" MeV\n";
    cout<<"3-Momentum Magnitude = "<<p<<" MeV\n";
    cout<<"======Spherical Polar======\n";
    cout<<"Polar Angle (theta) = "<<theta*180/pi<<" degrees\n";
    cout<<"Azimuthal Angle (phi) = "<<phi*180/pi<<" degrees\n";
    cout<<"Tranverse Momentum = "<<pt<<" MeV\n";
    if(c==0){cout<<setprecision(4)<<"Pseudorapidity = "<<eta<<"(upto 4 decimal places)\n";}
    else    {cout<<"Pseudorapidity : UNDEFINED.\n";}

}

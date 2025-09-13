#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <iomanip>

#include <Security/Security.h>
#include <cstdint>
#include <cstring>

using namespace std;

class Particle{
    private:
    int m_pid;
    double m_px,m_py,m_pz,m_pt,m_eta;

    public:
    Particle(int pid, double px, double py, double pz): m_pid(pid),m_px(px),m_py(py),m_pz(pz){}
    int getpid() const {return m_pid;}
    double getpx() const {return m_px;}
    double getpy() const {return m_py;}
    double getpz() const {return m_pz;}
    double getpt() const {return m_px*m_px + m_py*m_py;}
};

//use -framework Security
double distrib(double min,double max)
{
    uint64_t r;
    if(SecRandomCopyBytes(kSecRandomDefault, sizeof(r), (uint8_t*)&r) != errSecSuccess){return 0;}
    double rand = (double)r / (double)UINT64_MAX;
    return min+(max-min)*rand;
}

double distrib1(double min, double max)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(min,max);
    return dist(gen);
}

int main()
{
    double min=-1500.0;double max=1500.0;
    vector<Particle> finalstate;

    double p1x=distrib1(min,max);double p1y=distrib1(min,max);double p1z=distrib1(min,max);
    Particle pion1(211,p1x,p1y,p1z);
    double p2x=distrib1(min,max);double p2y=distrib1(min,max);double p2z=distrib1(min,max);
    Particle pion2(211,p2x,p2y,p2z);

    double p3x=-(p1x+p2x);double p3y=-(p1y+p2y);double p3z=-(p1z+p2z);
    Particle pion3(-211,p3x,p3y,p3z);

    finalstate.push_back(pion1);
    finalstate.push_back(pion2);
    finalstate.push_back(pion3);

    ofstream outFile("/Users/swarupdas/icloud/Documents/ToCoSiAn/data/event.dat");
    if(outFile.is_open()){
        for (const auto& p: finalstate){
            outFile<<p.getpid()<<" "<<p.getpx()<<" "<<p.getpy()<<" "<<p.getpz()<<"\n";
        }
        outFile.close();
        cout<<"Successfully wrote one event to event.dat"<<endl;
    }
    else {cerr<<"Warning! Can't open file event.dat\n"; return -1;}

    cout<<"Generated Particle after the event--------\n";
    for (const auto& p: finalstate){
    cout<<"Particle PID : "<<p.getpid()<<", Transverse Momenta = "<<p.getpt()<<" MeV\n";}
    
    return 0;
}
    

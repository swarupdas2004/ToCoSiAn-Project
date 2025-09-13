#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <iomanip>

using namespace std;

class Particle{
    private:
    int m_pid;
    double m_px,m_py,m_pz,m_pt;

    public:
    Particle(int pid, double px, double py, double pz): m_pid(pid),m_px(px),m_py(py),m_pz(pz){}
    int getpid() const {return m_pid;}
    double getpx() const {return m_px;}
    double getpy() const {return m_py;}
    double getpz() const {return m_pz;}
    double getpt() const {return sqrt(m_px*m_px + m_py*m_py);}
    double getp() const {return sqrt((m_px*m_px + m_py*m_py + m_pz*m_pz));}
    double geteta() const {return -log(tan(acos(m_pz/sqrt((m_px*m_px + m_py*m_py + m_pz*m_pz)))/2.0));}
};

double distrib(double min, double max)
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

    double p0x=distrib(min,max);double p0y=distrib(min,max);double p0z=distrib(min,max);
    Particle img(-9999,p0x,p0y,p0z);
    finalstate.push_back(img);


    double p1x=distrib(-min,max);double p1y=distrib(min,max);double p1z=distrib(min,max);
    Particle pion1(211,p1x,p1y,p1z);
    double p2x=distrib(min,max);double p2y=distrib(min,max);double p2z=distrib(min,max);
    Particle pion2(211,p2x,p2y,p2z);

    double p3x=distrib(min,max);double p3y=distrib(min,max);double p3z=distrib(min,max);
    Particle pion3(-211,p3x,p3y,p3z);
    double p4x=p0x-(p1x+p2x+p3x);double p4y=p0y-(p1y+p2y+p3y);double p4z=p0z-(p1z+p2z+p3z);
    Particle pion4(-211,p4x,p4y,p4z);

    finalstate.push_back(pion1);
    finalstate.push_back(pion2);
    finalstate.push_back(pion3);
    finalstate.push_back(pion4);

    ofstream outFile("/Users/swarupdas/icloud/Documents/ToCoSiAn/data/event1.dat");
    if(outFile.is_open()){
        for (const auto& p: finalstate){
            outFile<<p.getpid()<<setw(10)<<p.getpx()<<setw(10)<<p.getpy()<<setw(10)<<p.getpz()<<"\n";
        }
        outFile.close();
        cout<<"Successfully wrote one event to event1.dat"<<endl;
    }
    else {cerr<<"Warning! Can't open file event1.dat\n"; return -1;}

    cout<<"Generated Particle after the event--------\n";
    for (const auto& p: finalstate){
        if (p.getp()>1e-10 && (p.getp()-abs(p.getpt()))>1e-10){
            cout<<"Particle PID : "<<p.getpid()<<", Total Momentum = "<<p.getp()<<"MeV, Transverse Momentum = "<<p.getpt()<<" MeV, Pseudorapidity = "<<p.geteta()<<"\n";
        }
        else{
            cout<<"Particle PID : "<<p.getpid()<<", Total Momentum = "<<p.getp()<<"MeV, Transverse Momentum = "<<p.getpt()<<" MeV, Pseudorapidity : Undefined\n";
        }
    }
    return 0;
}
    

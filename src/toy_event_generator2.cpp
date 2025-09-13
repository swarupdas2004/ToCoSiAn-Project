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

    for (int i=1; i<=1000000;i++){

    double p1x=distrib(min,max);double p1y=distrib(min,max);double p1z=distrib(min,max);
    Particle pion1(211,p1x,p1y,p1z);
    double p2x=distrib(min,max);double p2y=distrib(min,max);double p2z=distrib(min,max);
    Particle pion2(211,p2x,p2y,p2z);

    double p3x=-(p1x+p2x);double p3y=-(p1y+p2y);double p3z=-(p1z+p2z);
    Particle pion3(-211,p3x,p3y,p3z);

    finalstate.push_back(pion1);
    finalstate.push_back(pion2);
    finalstate.push_back(pion3);

    }

    ofstream outFile("/Users/swarupdas/icloud/Documents/ToCoSiAn/data/event2.dat");
    if(outFile.is_open()){
        for (const auto& p: finalstate){
            outFile<<p.getpid()<<setw(15)<<p.getpx()<<setw(15)<<p.getpy()<<setw(15)<<p.getpz()<<"\n";
        }
        outFile.close();
        cout<<"Successfully wrote one event to event2.dat"<<endl;
    }
    else {cerr<<"Warning! Can't open file event.dat\n"; return -1;}

    
    return 0;
}
    

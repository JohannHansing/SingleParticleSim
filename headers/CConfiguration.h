#ifndef CCONFIGURATION_H_
#define CCONFIGURATION_H_

#include <string>
#include <vector>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include "CPolymers.h"
#include "CRod.h"


using namespace std;

class CConfiguration {
    /*Class where all the configuration variables such as potRange etc. and also most functions for the
     * simulation are stored
     */
private:
    //SCALING
    double _timestep;         //This is the RESCALED timestep! timestep = dt * kT / (frictionCoeffcient * particlesize)
    double _mu_sto;
    double _pradius;     //particle size is most relevant for scaling! (default = 1)
    double _boxsize;          // ALWAYS define boxsize through particlesize due to scaling!
    double _epsilon;
	double _hpi_u;
	double _hpi_k;


    //EXPONENTIAL Potential
    double _potRange;         // Avoid values like 10/3 = 3.33333... ALWAYS define the Range of the exp potential through boxsize due to scaling!
    double _potStrength;      // rescaled U_0 for exponential Potential
    double _rodDistance;
    CPolymers _poly;


    //bool Parameters
    bool _potMod;             // true if the modified exponential potential version is used that is not 3*U_0 at the intersections.
    bool _LJPot;              // if true, then LJ interaction is calculated for steric hindrance
    bool _ranU;
    bool _hpi;
    bool _ranRod;

    //COUNTERS AND INIT VALUES
    int _boxnumberXYZ[3];           //counter to calculate the actual position of the particle
    unsigned int _wallcrossings[3]; //counts how many times the particle has crossed the [entry, opposite, side] wall, while
                                            //travelling through the lattice
    int _entryside[3];            //records through which side and in which direction the particle last entered a box. For the face of the cube in the x-z
                                            //plane at y=0, it is entryside[0] = 1, in the x-y plane at z=L it is entryside[1] = -1!
    double _resetpos;
    double _startpos[3];          //Stores where the particle starting position was. This is needed to calculate the mean square displacement
    double _prevpos[3];           //Stores previous particle position before particle is moved.
    vector<vector<vector<int> > > _posHistoM;

    int _min, _max;        // parameters for determining up to which order neighboring rods are considered for the potential

    //Particle parameters
    double _ppos[3];    //initialize particle position (DO IT LIKE resetpos FOR MOVEPARTICLEFOLLOW/RESET)
    double _upot;
    double _f_mob[3];   //store mobility and stochastic force
    double _f_sto[3];


    boost::mt19937 *m_igen;                      //generate instance of random number generator "twister".
    double zerotoone(){
        boost::uniform_01<boost::mt19937&> dist(*m_igen);
        return dist();
    }
    
    // NEW RandomMesh Stuff
    /*TODO :
     * Implement random number rods depending on n_rods.
     * Test if I can use C++11 on sheldon. If so, use:
     * for(auto& s: _rodvec[ortho[oa]]){
     *    do smth with s;
     *} 
     * instead of  
     * (int i=0;i<_rodvec[ortho[oa]].size();i++){
     *     do smth with _rodvec[ortho[oa]][i]
     * }
     * for looping over arrys/vectors
     */
    std::vector<vector<CRod>> _rodvec; // vector to store polymer rods in cell, one vector stores polymers that are parallel to the same axis
    
public:
    void initRodsVec(double n_rods=1){
        double xipos, xjpos;
        _rodvec.resize(3);
        int Nrods[3]; // number of rods in certain plane, i.e. parallel to a certain axis.
        
        for (int i=0;i<3;i++){
            //TODO if zerotoone()> ... MAKE Nrods[i] random
            Nrods[i] =  9 * n_rods;// 3x3 cells in one plane -> 9*n_rods
        }
        for (int axis=0;axis<3;axis++){//axis 0 is x axis.
            for (int i=0; i<Nrods[axis];i++){
                xipos = (zerotoone() * 3 - 1) *_boxsize ;
                xjpos = (zerotoone() * 3 - 1) *_boxsize;
                CRod newRod = CRod(axis, xipos, xjpos );
                _rodvec[axis].push_back(newRod);
            }
        }
    }
    void updateRodsVec(int crossaxis,int exitmarker){//exitmarker is -1 for negative direction, or 1 for positive
        //delete all polymers orthogonal to crossaxis, that are outside the box now
        //update other polymer positions
        int ortho[2] = {1,2};
        if (crossaxis == 1){
            ortho[0]=2; 
            ortho[1]=0;
        }
        else if (crossaxis == 2){
            ortho[0]=0; 
            ortho[1]=1;
        }
        // shift positions of rods
        for (int oa=0;oa<2;oa++){
            for (int i=0;i<_rodvec[ortho[oa]].size();i++){
                //shift rod positions parallel to crossaxis. ortho[oa] is direction that the shifted rods are parallel to.
                _rodvec[ortho[oa]][i].coord[crossaxis] -= exitmarker * _boxsize;
                if (abs(_rodvec[ortho[oa]][i].coord[crossaxis] - _boxsize/2.  ) > 1.5*_boxsize){// TODO CHECK THIS AGAIN!
                    //in direction parallel to crossaxis, choose new position in side cell 
                    _rodvec[ortho[oa]][i].coord[crossaxis] = (zerotoone()  + exitmarker) * _boxsize;
                    int ortho2 = 3 - (ortho[oa] * crossaxis);
                    // in direction orthogonal to both ortho[oa] and crossaxis
                    _rodvec[ortho[oa]][i].coord[ortho2] = (zerotoone() * 3 -1) * _boxsize;// TODO CHECK THIS AGAIN! Draw a sketch
                }
            }
        }
    }
    
    void prinRodPos(int axis){
        for (int irod=0;irod<_rodvec[axis].size();irod++){
            double rx = _rodvec[axis][irod].coord[0];
            double ry =_rodvec[axis][irod].coord[1];     
            double rz =_rodvec[axis][irod].coord[2];
            cout << ",[" << rx << "," << ry << "," << rz << "]";
        }
        cout << "]," << endl;
    }



private:
    void setRanNumberGen(double seed);
    void countWallCrossing(int crossaxis, int exitmarker);
    void calculateExpHPI(const double r, double& U, double& Fr);
    void calculateExpPotential(const double r, double& U, double& Fr);
    void calculateExpPotentialMOD(const double r, double& U, double& Fr, int index);
    void modifyPot(double& U, double& Fr, double dist);
    void calcLJPot(const double r, double &U, double &dU);
    void initPosHisto();
    




public:
    CConfiguration();
    CConfiguration(
            double timestep,  double potRange,  double potStrength,  double boxsize, double rodDistance, const bool potMod, double psize,
            const bool posHisto, const bool steric, const bool ranU,  bool hpi, double hpi_u, double hpi_k, bool ranRods);
    void resetParameters(double timestep, double potRange, double potStrength, double boxsize);
    void updateStartpos();
    void resetposition();
    void makeStep();
    void checkBoxCrossing();
    void calcStochasticForces();
    void calcMobilityForces();
    void calc_1RODONLY_MobilityForces();
    void calc_ZRODONLY_MobilityForces();
    void calc_YZRODONLY_MobilityForces();
    void saveXYZTraj(string name, const int& move, string a_w);
    void positionHistogram(double block, double possq, double pposprev, int l, int posHisto[]);

    double getPosVariance();
    double get1DPosVariance(int dim);
    double getUpot(){ return _upot; }
    double getDisplacement();
    unsigned int getwallcrossings(int i){ return _wallcrossings[i]; }
    bool checkFirstPassage(double mfpPos, int dim);
    bool testOverlap();
    void moveBack();
    void addHistoValue();
    void printHistoMatrix(string folder);
    std::vector<double> getppos(){ // returns pointer to current particle position array
    	std::vector<double> pos (3);
    	for (int i = 0; i < 3; i++){
            pos[i] = _ppos[i] + _boxsize * _boxnumberXYZ[i];
    	}
    	return pos;
    }





};



#endif /* CCONFIGURATION_H_ */

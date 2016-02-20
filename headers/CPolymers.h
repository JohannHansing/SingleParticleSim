/*
 * CPolymers.h
 *
 *  Created on: Aug 9, 2013
 *      Author: jh
 */

#ifndef CPOLYMERS_H_
#define CPOLYMERS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <time.h>
#include "boost/random.hpp"

using namespace std;



class CPolymers {
private:
    vector<vector<vector<int> > > _sign;
    boost::mt19937 *_igen;         // produces randomness out of thin air
                                        // see pseudo-random number generators
    int ran_sign();
    void printSign();
    void setRanNumberGen(double seed);

public:
    int axis; // rod parallel to axis 0 = x, axis 1 = y, etc.
    double planeCoord[3]; // Coordinate of rod in 2D plane orthogonal to axis. the planCoordinate parallel to axis is always 0. (see initiaion)
    
    CPolymers(int ax, double xi, double xj);
    void shiftPolySign(const int , const int);
    bool samesign(const int direction, const int plane_axis, const int edgenumber);
    int get_sign(const int plane_axis, const int edgenumber);
	void printSignFile(string a_w);

};



#endif /* CPOLYMERS_H_ */

/**
 * @file main.cpp
 * @brief main program file for the bilateral filter program
 * @author Julie Digne julie.digne@liris.cnrs.fr
 * @date 2014/08/14
 * @copyright This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */


#include <iostream>
#include <sstream>
#include <ctime>
#include <vector>

#include "Octree.h"
#include "OctreeIterator.h"
#include "BilateralFilter.h"
#include "utilities.h"
#include "FileIO.h"

#include "types.h"
#include "cmdLine.h"

using namespace std;

/** @brief display the usage of the program
 */
void usage()
{
    std::cout<<"USAGE"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"bilateralfilter INPUT OUTPUT"
             <<"-r <radius> -n <normal_radius> -N <Niter> -p"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"INPUT  (mandatory) file containing the points (oriented or not)."<<std::endl;
    std::cout<<std::endl;
    std::cout<<"OUTPUT (mandatory) output_file."<<std::endl;
    std::cout<<std::endl;
    std::cout<<"-r     (recommended) neighborhood radius"<<std::endl;
    std::cout<<"       will be translated into a distance weight."<<std::endl;
    std::cout<<std::endl;
    std::cout<<"-n     (recommended) normal neighborhood radius"<<std::endl;
    std::cout<<"       will be translated into a normal distance weight."<<std::endl;
    std::cout<<std::endl;
    std::cout<<"-N     Number of filter iterations (default: 1)."<<std::endl;
    std::cout<<std::endl;
    std::cout<<"-p     (recommended) will perform the computation in parallel."<<std::endl;
    std::cout<<std::endl;
    std::cout<<"-R     use this option if the point cloud is oriented"<<std::endl;
    std::cout<<std::endl;
}

/**
 * @brief main function for the ball pivoting reconstruction
 * @param argc
 * @param argv
 * @return 1 if the program exited successfully
 */
int main(int argc, char **argv)
{
    CmdLine cmd;
    //handling command line options
    string infile, outfile;
    double radius = -1;
    double normal_radius = -1;
    int parallel_flag = -1;
    int niter = 1;

    cmd.add(make_option('r', radius));
    cmd.add(make_option('n',normal_radius));
    cmd.add(make_option('N', niter));
    cmd.add(make_switch('p'));
    cmd.add(make_switch('R'));


    try {
        cmd.process(argc, argv);
    } catch(std::string str) {
        std::cerr << "Error: " << str << std::endl<<std::endl;
        usage();
        return 1;
    }

    if(argc < 3)
    {
        usage();
        return EXIT_FAILURE;
    }

    infile = argv[1];
    outfile = argv[2];

    std::cout<<"Input file:  "<<argv[1]<<std::endl;
    std::cout<<"Output file: "<<argv[2]<<std::endl;


    if(cmd.used('p')) parallel_flag = 1;

    if( (!cmd.used('N')) || (niter <= 0) )
    {
        niter = 1;
        std::cout<<"No number of iterations given or nonpositive number of";
        std::cout<<" iterations, using default: "<<niter<<std::endl;
    }
    if( ( (! cmd.used('r')) || (radius < 0) ) && (! cmd.used('d')))
    {
        std::cout<<"No radius given or positive radius and no depth given!"<<std::endl;
        std::cout<<"Estimating a radius ... "<<std::flush;
        radius = FileIO::estimateRadius(infile.c_str());
        std::cout<<radius<<std::endl;;
    }
    if( (! cmd.used('n'))||( normal_radius < 0))
    {
        std::cout<<"No normal radius given; setting it equal to the radius: "<<radius<<std::endl;
        normal_radius = radius;
    }

    time_t start, end, globstart, globend;

    Octree octree;

    std::time(&globstart);

    bool ok;
    if(cmd.used('R'))
        ok = FileIO::readAndSortOrientedPoints(infile.c_str(),octree,radius); 
    else
        ok = FileIO::readAndSortUnorientedPoints(infile.c_str(),octree,radius); 

    if( !ok )
    {
        std::cerr<<"Pb opening the file; exiting."<<std::endl;
        return EXIT_FAILURE;
    }
    std::time(&end);

    std::cout<<"Parameters:"<<std::endl;
    std::cout<<"Radius: "<<radius<<std::endl;
    std::cout<<"Number of iterations "<<niter<<"."<<std::endl;
    std::cout<<"Normal radius: "<<normal_radius<<"."<<std::endl;
    std::cout<<"Octree with depth "<<octree.getDepth()<<" created."<<std::endl;
    std::cout<<"Octree contains "<<octree.getNpoints()<<" points."<<std::endl;
    std::cout<<"The bounding box size is "<<octree.getSize()<<std::endl;
    std::cout<<"Reading and sorting the points took "<<difftime(end,globstart)
             <<" s."<<std::endl;

    std::cout<<"Octree statistics"<<std::endl;
    octree.printOctreeStat();


    OctreeIterator iterator(&octree);

    if(radius<0)
        radius = octree.getSmallestCellSize();

    //creating the bilateral filter
    std::time(&start);
    BilateralFilter bilateralfilter(&octree, radius, normal_radius, niter);

    if(parallel_flag == 1)
        bilateralfilter.parallelApplyBilateralFilter();
    else
        bilateralfilter.applyBilateralFilter();

    std::time(&end);
    std::cout<<"Applying "<<niter<<" bilateral filter iteration(s) took "
             <<difftime(end,start) <<" s."<<std::endl;

    if(! FileIO::savePoints(outfile.c_str(), octree,
               bilateralfilter.getSetIndex(), cmd.used('R')) )
    {
      std::cerr<<"Trouble saving the result file, aborting"<<std::endl;
      return EXIT_FAILURE;
    }
    std::time(&globend);
    
    std::cout<<"Total computation time (including I/O): "
             <<difftime(globend,globstart) <<" s."<<std::endl;

    return EXIT_SUCCESS;
}

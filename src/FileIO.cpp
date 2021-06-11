/** @file FileIO.cpp
 * @brief file defining methods to read points from a file
 * @author Julie Digne julie.digne@liris.cnrs.fr
 * @date 2012/10/22
 * @copyright This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "FileIO.h"
#include<iostream>
#include<fstream>
#include<istream>
#include<ostream>
#include<sstream>
#include<algorithm>
#include<deque>
#include<limits>

using namespace std;

FileIO::FileIO()
{
}

FileIO::~FileIO()
{
}

double FileIO::estimateRadius(const char* filename)
{
    ifstream in;
    in.open(filename);

    if(!in)
    {
        std::cerr<<"File "<<filename<<" could not be opened"<<std::endl;
        return false;
    }

    string line;

    getline(in, line);

    double x,y,z;
    in >> x >> y >> z;

    double xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = xmax = x;
    ymin = ymax = y;
    zmin = zmax = z;

    unsigned int npts = 1;
    while(getline(in,line))
    {
        in >> x >> y >> z;
        xmin = std::min(x,xmin);
        xmax = std::max(x,xmax);
        ymin = std::min(y,ymin);
        ymax = std::max(y,ymax);
        zmin = std::min(z,zmin);
        zmax = std::max(z,zmax);
        ++npts;
    }
    in.close();

    double lx = xmax - xmin;
    double ly = ymax - ymin;
    double lz = zmax - zmin;

    double size = std::max(lx,ly);
    size = std::max(lz,size);

    return std::sqrt(20.0/(double)npts)*size;
}

bool FileIO::readAndSortUnorientedPoints(const char* filename,
                               Octree& octree,
                               double min_radius)
{
    ifstream in;
    in.open(filename);

    if(!in)
    {
        std::cerr<<"File "<<filename<<" could not be opened"<<std::endl;
        return false;
    }

    string firstline;
    getline(in, firstline);

    istringstream line_in(firstline);
    string word;
    int nword = 0;
    while (line_in >> word)
        nword++;

    if(nword < 3)
    {
        cerr<< "File must contain at least 3 doubles per line (x y z)"<<endl;
        return false;
    }

    size_t nprop = nword - 3;
    std::cout<<"Unoriented points with "<<nprop<<" properties"<<std::endl;

    in.clear() ;
    in.seekg(0, ios::beg);//starting back from the beginning

    double x,y,z;
    in >> x >> y >> z;

    deque<Sample> input_vertices;
    Sample s(x,y,z);
    s.setIndex(0);
    input_vertices.push_back(s);

    octree.addProperty(nprop);

    for(size_t t = 0 ; t < nprop ; ++t)
    {
      double val;
      in >> val;
      octree.setProperty(0, t, val);
    }

    double xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = xmax = x;
    ymin = ymax = y;
    zmin = zmax = z;

    size_t index = 1;
    while( in >> x >> y >> z)
    {
      Sample temp(x,y,z);
      temp.setIndex(index);
      input_vertices.push_back(temp);
      xmin = std::min(x, xmin);
      xmax = std::max(x, xmax);
      ymin = std::min(y, ymin);
      ymax = std::max(y, ymax);
      zmin = std::min(z, zmin);
      zmax = std::max(z, zmax);
    
    //reading properties
      octree.addProperty(nprop);
      for(size_t t = 0; t < nprop; ++t)
      {
        double val;
        in >> val;
        octree.setProperty(index, t, val);
      }
      ++index;
    }
    in.close();

    std::cout<<input_vertices.size()<<" points read"<<std::endl;

    double lx = xmax - xmin;
    double ly = ymax - ymin;
    double lz = zmax - zmin;

    double size = std::max(lx,ly);
    size = std::max(lz,size);
    size = 1.1 * size;//loose bounding box around the object

    double margin;
    if(min_radius > 0)
    {
        unsigned int depth = (unsigned int)ceil( log2( size / (min_radius) ));
        //adapting the bouding box size so that the smallest cell size is 
        //exactly 2*min_radius
        double adapted_size = pow2(depth) * min_radius;
        //margins of the bouding box around the object
        margin = 0.5 * (adapted_size - size); 
        size = adapted_size;
        octree.setDepth(depth);
    }
    else
    {
        margin = 0.05 * size; //margins of the bouding box around the object
    }

    double ox = xmin - margin;
    double oy = ymin - margin;
    double oz = zmin - margin;
    Point origin(ox,oy,oz);

    octree.initialize(origin, size);

    //add the points to the set with index 0 (initial set)
    octree.addInitialPoints(input_vertices.begin(), input_vertices.end());

    return true;
}

bool FileIO::readAndSortOrientedPoints(const char* filename,
                               Octree& octree,
                               double min_radius)
{
    ifstream in;
    in.open(filename);

    if(!in)
    {
        std::cerr<<"File "<<filename<<" could not be opened"<<std::endl;
        return false;
    }

    string firstline;
    getline(in, firstline);

    istringstream line_in(firstline);
    string word;
    int nword = 0;
    while (line_in>> word)
        nword++;

    if( nword < 6)
    {
        cerr<<"each point must be given by at least 6 values (position + normal) : "
        <<"x y z nx ny nz"<<endl;
        return false;
    }
    
    size_t nprop = nword - 6;
    std::cout<<"Oriented points with "<<nprop<<" properties"<<std::endl;

    in.clear() ;
    in.seekg(0, ios::beg);//starting back from the beginning

    double x,y,z,nx,ny,nz;
    in >> x >> y >> z >> nx >> ny >> nz;
    Sample s(x, y, z, nx, ny, nz);
    s.setIndex(0);
    deque<Sample> input_vertices;
    input_vertices.push_back(s);

    octree.addProperty(nprop);
    for(size_t t = 0 ; t < nprop ; ++t)
    {
      double val;
      in >> val;
      octree.setProperty(0, t, val);
    }

    double xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = xmax = x;
    ymin = ymax = y;
    zmin = zmax = z;

    size_t index = 1;
    while( in >> x >> y >> z >> nx >> ny >> nz)
    {
        Sample temp(x,y,z,nx,ny,nz);
        temp.setIndex(index);
        input_vertices.push_back(temp);
        xmin = std::min(x,xmin);
        xmax = std::max(x,xmax);
        ymin = std::min(y,ymin);
        ymax = std::max(y,ymax);
        zmin = std::min(z,zmin);
        zmax = std::max(z,zmax);

        octree.addProperty(nprop);
        for(size_t t = 0 ; t < nprop ; ++t)
        {
          double val;
          in >> val;
          octree.setProperty(index, t, val);
        }

        ++index;
    }
    in.close();

    std::cout<<input_vertices.size()<<" points read"<<std::endl;

    double lx = xmax - xmin;
    double ly = ymax - ymin;
    double lz = zmax - zmin;

    double size = std::max(lx,ly);
    size = std::max(lz,size);
    size = 1.1 * size;//loose bounding box around the object

    double margin;
    if(min_radius > 0)
    {
        unsigned int depth = (unsigned int)ceil( log2( size / (min_radius) ));
        //adapting the bouding box size so that the smallest cell size is 
        //exactly 2*min_radius
        double adapted_size = pow2(depth) * min_radius;
        //margins of the bouding box around the object
        margin = 0.5 * (adapted_size - size); 
        size = adapted_size;
        octree.setDepth(depth);
    }
    else
    {
        margin = 0.05 * size; //margins of the bouding box around the object
    }

    double ox = xmin - margin;
    double oy = ymin - margin;
    double oz = zmin - margin;
    Point origin(ox,oy,oz);

    octree.initialize(origin, size);

    //add the points to the set with index 0 (initial set)
    octree.addInitialPoints(input_vertices.begin(), input_vertices.end());

    return true;
}

bool FileIO::savePoints(const char* filename, Octree& octree,
                        unsigned int index, bool isoriented)
{ 
    ofstream out;
    out.open(filename);
    out.precision( numeric_limits<double>::digits10 + 1);

    if(!out)
        return false;

    OctreeNode *node = octree.getRoot();
    saveContent(node, octree, index, out, isoriented);

    out.close();

    return true;
}


void FileIO::saveContent(OctreeNode* node, Octree &octree, unsigned int index,
                         ofstream& f, bool isoriented)
{
    if(node->getDepth() != 0)
    {
        for(int i = 0; i < 8 ;i++)
            if(node->getChild(i) != NULL)
                saveContent(node->getChild(i), octree, index, f, isoriented);
    }
    else if(node->getNpts(index) != 0)
    {
        Sample_deque::const_iterator iter;
        std::vector<double>::iterator pi;
        for(iter = node->points_begin(index);
            iter != node->points_end(index);++iter)
        {
            const Sample &s = *iter;
            size_t sindex = s.index();
            if(isoriented)
                f << s;
            else
                f << s.x() << "\t" << s.y() << "\t" << s.z();
            if(octree.getNproperties() > 0)
            {
                f<<"\t";
                for(size_t pindex = 0; pindex < octree.getNproperties(); ++pindex)
                {
                    double val = octree.getProperty(sindex, pindex);
                    f<<val<<"\t";
                }
            }
            f<<endl;
        }
    }
}



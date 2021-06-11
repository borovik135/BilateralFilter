/**
 * @file BilateralFilter.h
 * @brief methods for applying the bilateral filter
 * @author Julie Digne
 * @date 2014/08/14
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

#ifndef BILATERAL_FILTER_H
#define BILATERAL_FILTER_H

#include<cstdlib>
#include <list>
#include <map>
#include "Point.h"
#include "Octree.h"
#include "OctreeNode.h"
#include "OctreeIterator.h"
#include <cmath>
#include <cassert>
/**
 * @class TBilateralFilter
 * @brief class providing access to the bilateral filter 
 */
template<class T>
class TBilateralFilter
{
    private ://class members

        /**octree containing the points to mesh
         * */
        TOctree<T> *m_octree;

        /**iterator over the octree*/
        TOctreeIterator<T> *m_iterator;

        /**radius for the neighborhood
         */
        double m_radius;
        
        /**radius for the normal
         */
        double m_normal_radius;

        /**number of iterations
         */
        unsigned int m_niter;

        /**index of the set to process
         */
        unsigned int m_setIndex;

        /**gaussian weight on the distance*/
        double m_weight;

        /**gaussian weight on the normal distance*/
        double m_normal_weight;

        /**min number of neighbors so a point is not removed*/
        static const size_t MIN_NEIGHBORS=5;

    public ://constructors+destructors
        /**Default Constructor
         */
        TBilateralFilter<T>();

        /**Constructor
         * @param octree the octree containing the point set to process
         * @param radius neighborhood radius (also sets the weight on the euclidean distance)
         * @param normal_radius sets the weight on the normal distance
         * @param niter number of iterations to perform
         */
        TBilateralFilter<T>(TOctree<T> *octree, double radius,
                            double normal_radius, unsigned int niter);

        /**Destructor
         * */
        ~TBilateralFilter<T>();

    public : //accessors
        /**get index of the set
         *     @return index of the set that is being meshed
         */
        unsigned int getSetIndex() const;

        /**set index of the set 
         *     @param index of the set that is being meshed
         */
        void setSetIndex(unsigned int index);

        /**get the next index to save the next iteration result
         * set 0 is the initial set, its filtering should be put
         * in set 1, set 0 will never be destroyed
         * The result of filtering set 1 is put in set 2
         * The result of filtering set 2 is put in set 1
         * Therefore at each step only the last iteration result
         * and the initial set are remembered 
         * @return next index
         */
        unsigned int getNextSetIndex() const;


    public : //filter methods

        /** apply the bilateral filter iterations to the point set
         */
        void applyBilateralFilter();

        /**apply the bilateral filter iterations in parallel
         */
        void parallelApplyBilateralFilter();

    private : //auxiliary methods for applying the bilateral filter

        /**apply the bilateral filter to a given cell
         * @param cell cell to process
         * @param nextindex index of the set to save the iteration result to
         */
        void applyBilateralFilter(TOctreeNode<T> *cell, unsigned int nextindex);

        /**apply the bilateral filter to a given cell with a
         * parent at the right scale
         * @param cell cell to process
         * @param parent parent of the cell whose size corresponds to the 
         * processing radius
         * @param nextindex index of the set to save the iteration result to
         */
        void applyBilateralFilter(TOctreeNode<T> *cell, TOctreeNode<T> *parent,
                                            unsigned int nextindex);

        /**apply the bilateral filter to a given point
         * @param p point to apply the bilateral filter to
         * @param parent cell containing the point and whose size
         * corresponds to the processing radius
         * @param nextindex index of the set to save the iteration result to
         */
         void applyBilateralFilter(T &p, TOctreeNode<T> *parent,
                              unsigned int nextindex);

        /**perform local PCA: the weighted barycenter and normal to the local 
         * regression plane are computed
         * The algorithm for decomposing 3x3 real symmetric matrices is
         * described in [Smith, 1961]
         * @param neighbors list of neighbors
         * @param distances distances from the neighbors to the center point 
         * @param bar barycenter of the distance-weighted neighbors
         * @param nx estimated normal (x coeff)
         * @param ny estimated normal (y coeff)
         * @param nz estimated normal (z coeff)
         */
         void performLocalPCA(const std::list<T*> &neighbors,
                              const std::list<double> &distances, Point &bar,
                              double &nx, double &ny, double &nz) const;
};

template<class T>
TBilateralFilter<T>::TBilateralFilter()
{
    m_octree = NULL;
    m_iterator = NULL;
    m_radius = 0.0;
    m_niter = 0;
    m_weight = 0.0;
    m_setIndex = 0;
}

template<class T>
TBilateralFilter<T>::TBilateralFilter(TOctree<T>* octree,
                            double radius,
                            double normal_radius,
                            unsigned int niter)
{
    m_octree = octree;
    m_iterator = new TOctreeIterator<T>(octree);
    m_radius = radius;
    m_normal_radius = normal_radius;
    m_niter = niter;
    m_iterator->setR(m_radius);
    m_weight = - 4.5/(radius*radius);//= -1/(2*(1/3radius)^2)
    m_normal_weight = -4.5/(m_normal_radius * m_normal_radius); 
    m_setIndex = m_iterator->getSetIndex();
}


template<class T>
TBilateralFilter<T>::~TBilateralFilter()
{
    delete m_iterator;
}


template<class T>
unsigned int TBilateralFilter<T>::getSetIndex() const
{
    return m_setIndex;
}

template<class T>
void TBilateralFilter<T>::setSetIndex(unsigned int index)
{
    m_setIndex = index;
    m_iterator->setSetIndex(m_setIndex);
}

template<class T>
unsigned int TBilateralFilter<T>::getNextSetIndex() const
{
  return (m_setIndex == 1 ? 2:1);
}


    template<class T>
void TBilateralFilter<T>::parallelApplyBilateralFilter()
{
    std::cout<<"Filtering with radius "<<m_radius<<"."<<std::endl;
    TOctreeNode<T> *root = m_octree->getRoot();
    unsigned int depth = m_iterator->getDepth();

    const double d = 2.1 * m_radius;

    //look for the smallest depth at which the cell size is above 2.1radius
    depth = (unsigned int)(m_octree->getDepth() - floor( log2(
                    m_octree->getSize() / d)));

    if(depth > m_octree->getDepth() )
        depth = m_octree->getDepth();
    //depth should be at least 1 to ensure that the same cell set is not
    //accessed at the same time
    else if(depth < 1)
        depth = 1;

    std::cout<<"Processing depth "<< depth <<" ; size "
        << m_octree->getSize()/(double)pow2(m_octree->getDepth()-depth)
        <<" ; dilatation radius "<<d<<std::endl;

    typename TOctree<T>::OctreeNode_collection node_collection;

    m_octree->getNodes(depth, root, node_collection);

    for(unsigned int n = 0; n < m_niter ; ++n)
    {
        unsigned int nextindex = getNextSetIndex();

        for(unsigned int i = 0; i < 8; ++i)
        {

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
            for(int j = 0; j < (int)node_collection[i].size(); ++j)
            {
                TOctreeNode<T> *node = node_collection[i][j];
                applyBilateralFilter(node, nextindex);
            }
        }

        if(m_setIndex > 0)
        {
            for(unsigned int i = 0; i < 8; ++i)
            {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
                for(int j = 0; j < (int)node_collection[i].size(); ++j)
                {
                    TOctreeNode<T> *node = node_collection[i][j];
                    node->clearSet(m_setIndex);
                }
            }
        }
        setSetIndex(nextindex);
        std::cout<<"Iteration "<<(n+1)<<" done."<<std::endl;
    }
}


    template<class T>
void TBilateralFilter<T>::applyBilateralFilter()
{
    std::cout<<"Filtering with radius "<<m_radius<<"."<<std::endl;
    for(unsigned int i=0; i < m_niter; ++i)
    {
        unsigned int nextindex = getNextSetIndex();
        TOctreeNode<T> *root = m_octree->getRoot();
        applyBilateralFilter(root, nextindex);
        if(m_setIndex > 0)
            m_octree->clearSet(m_setIndex);
        setSetIndex(nextindex);
        std::cout<<"Iteration "<<(i+1)<<" done."<<std::endl;
    }
}


    template<class T>
void TBilateralFilter<T>::applyBilateralFilter(TOctreeNode<T>* cell,
        const unsigned int index)
{
    if(cell->getDepth() > m_iterator->getDepth())
    {
        for(unsigned int i = 0; i < 8 ; ++i)
        {
            if(cell->getChild(i) != NULL)
                applyBilateralFilter(cell->getChild(i), index);
        }
    }
    else if(cell->getDepth() == m_iterator->getDepth())
    {
        applyBilateralFilter(cell, cell, index);
    }
}


    template<class T>
void TBilateralFilter<T>::applyBilateralFilter(TOctreeNode<T>* cell,
        TOctreeNode<T>* parent,
        const unsigned int index)
{
    if(cell->getDepth() == 0)
    {
        typename std::deque<T>::iterator pi;
        for(pi = cell->points_begin(m_setIndex);
                pi != cell->points_end(m_setIndex); ++pi)
        {
            applyBilateralFilter(*pi, parent, index);
        }
    }
    else
    {
        for(unsigned int i = 0; i < 8 ; ++i)
        {
            if(cell->getChild(i) != NULL)
                applyBilateralFilter(cell->getChild(i), parent, index);
        }
    }
}


template<class T>
void TBilateralFilter<T>::applyBilateralFilter(T& p, TOctreeNode<T>* parent,
                                     const unsigned int index)
{
    std::list< T* > neighbors;
    std::list<double> distances;
    m_iterator->getNeighbors(p, parent, neighbors, distances);
 
    if(neighbors.size()<MIN_NEIGHBORS)
        return;//not enough neighbors

    Point barycenter;
    double nx,ny,nz;
    performLocalPCA(neighbors, distances, barycenter, nx, ny, nz);

    if(nx * p.nx() + ny * p.ny() + nz * p.nz() < 0)
    {
        nx = -nx;
        ny = -ny;
        nz = -nz;
    }

    double diffx,diffy,diffz,h,wc,ws;
    double sum = 0.0;
    double normalizer = 0.0;
    std::list<double>::iterator itd = distances.begin();  
    typename std::list<T*>::iterator it;
    for(it = neighbors.begin(); it != neighbors.end(); ++it, ++itd)
    {
      diffx = (*it)->x() - p.x();
      diffy = (*it)->y() - p.y();
      diffz = (*it)->z() - p.z();
      h = nx * diffx + ny * diffy + nz * diffz;

      wc = exp( h*h * m_normal_weight);
      ws = exp((*itd) * m_weight);
      double w = wc * ws;

      sum = sum + w * h;
      normalizer = normalizer + w;
    }
    if(normalizer < 1e-15)
    {
      sum = 0.0;
      return;//suppress the point
    }
    else
      sum = sum / normalizer;

    T pnew( p.x() + sum * nx, p.y() + sum * ny, p.z() + sum * nz,
            p.nx(), p.ny(), p.nz());
    pnew.setIndex(p.index());

    m_octree->addPoint(pnew, index);
}


template<class T>
void TBilateralFilter<T>::performLocalPCA(const std::list< T* >& neighbors,
        const std::list< double >& distances,
        Point& bar,
        double& nx, double& ny, double& nz) const
{
    double xo,yo,zo;
    xo = yo = zo = 0.0;
    double a00, a01, a02, a11, a12, a22;
    a00 = a01 = a02 = a11 = a12 = a22 = 0.0;
    double sum_w = 0.0;

    typename std::list<T*>::const_iterator ni;
    typename std::list<double>::const_iterator di;


    for(ni = neighbors.begin(), di = distances.begin() ;
            ni != neighbors.end() ; ++ni, ++di)
    {
        const T* s = *ni;
        const double w = exp((*di)*m_weight);

        double temp = sum_w + w;
        double dx = s->x() - xo;
        double dy = s->y() - yo;
        double dz = s->z() - zo;
        double ratio = w / temp;
        double rx = dx * ratio;
        double ry = dy * ratio;
        double rz = dz * ratio;

        xo += rx;
        yo += ry;
        zo += rz;
        a00 += sum_w * dx * rx;
        a01 += sum_w * dx * ry;
        a02 += sum_w * dx * rz;
        a11 += sum_w * dy * ry;
        a12 += sum_w * dy * rz;
        a22 += sum_w * dz * rz;

        sum_w = temp;
    }

    double t = neighbors.size();
    t = t/(sum_w * (t - 1.0));
    a00 = a00 * t;
    a01 = a01 * t;
    a02 = a02 * t;
    a11 = a11 * t;
    a12 = a12 * t;
    a22 = a22 * t;

    //computing the least eigenvalue of the covariance matrix
    //Algorithm described in [Smith, 1961] for decomposing 3x3
    //real symmetric matrices
    double q = (a00 + a11 + a22)/3.0;
    double pq = (a00 - q)*(a00 - q) + (a11 - q) * (a11 - q)
        + (a22 - q) * (a22 - q)
        + 2 * (a01 * a01 + a02 * a02 + a12 * a12);
    pq = sqrt(pq / 6.0);
    double mpq = 1./pq;
    mpq = pow(mpq,3);

    double detB = mpq * (  (a00 - q)*( (a11 - q) * (a22 - q) - a12 * a12)
            - a01 * ( a01 * (a22 -q) - a12 * a02)
            + a02 * ( a01 * a12 - ( a11 -q) * a02)); 

    double r = 0.5 * detB;

    double phi = 0.0;
    if (r <= -1) 
        phi = M_PI / 3.0;
    else if (r >= 1.0)
        phi = 0.0;
    else
        phi = acos(r) / 3.0;

    double eig = q + 2.0 * pq * cos(phi + M_PI * (2.0/3.0));

    //computing the corresponding eigenvector
    nx =   a01 * a12 - a02 * (a11 - eig);
    ny =   a01 * a02 - a12 * (a00 - eig);
    nz =   (a00 - eig) * (a11 -eig) - a01*a01;

    normalize(nx,ny,nz);

    bar = Point(xo, yo, zo);
}


#endif

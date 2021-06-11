/**
 *@file types.h
 * @brief defines useful types
 * @author Julie Digne
 * @date 2012/10/25
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


#ifndef TYPES_H
#define TYPES_H

#include<list>
#include<vector>
#include<set>
#include<map>
#include<deque>
#include "Point.h"

class Sample;


typedef std::deque<Sample> Sample_deque;
typedef std::deque<Sample*> Sample_star_deque;

typedef std::map<double, Sample*> Neighbor_star_map;
typedef std::list<Sample*> Neighbor_star_list;

typedef std::list<double> Distance_list;
typedef std::list<Sample*> Sample_star_list;

#include "Octree.h"
typedef TOctree<Sample> Octree;

#include "OctreeNode.h"
typedef TOctreeNode<Sample> OctreeNode;

#include "OctreeIterator.h"
typedef TOctreeIterator<Sample> OctreeIterator;


#include "BilateralFilter.h"
typedef TBilateralFilter<Sample> BilateralFilter;


#endif
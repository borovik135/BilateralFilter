/**
 * @file Sample.cpp
 * @brief implementation of the sample methods declared in Sample.h
 * @author Julie Digne
 * @date 2012-10-08
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

#include "Sample.h"
#include "utilities.h"
#include <cstdio>
#include <iostream>

Sample::Sample() : Point()
{
    m_nx=m_ny=m_nz=0.0;
    m_index = -1;
}

Sample::Sample(double x, double y, double z, double nx, double ny, double nz)
                : Point(x,y,z)
{
    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
    m_index = -1;
}

Sample::Sample(double x, double y, double z)
                : Point(x,y,z)
{
    m_nx = 0;
    m_ny = 0;
    m_nz = 0;
    m_index = -1;
}

Sample::~Sample()
{
    m_nx=m_ny=m_nz=0.0;
    m_index = -1;
}

int Sample::index() const
{
    return m_index;
}

void Sample::setIndex(int index)
{
    m_index = index;
}


double Sample::nx() const
{
    return m_nx;
}

double Sample::ny() const
{
    return m_ny;
}

double Sample::nz() const
{
    return m_nz;
}



std::ostream& operator << (std::ostream& out, const Sample& v)
{
    out << v.x() << "\t" << v.y() << "\t" << v.z()
        << "\t" << v.nx() << "\t" << v.ny() << "\t" << v.nz();
    return out;
}




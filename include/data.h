/*
 *
 * tsp.cpp
 *
 * Traveling salesmen problem (TSP).
 *
 * Created by Isaac Balster on 26/01/2023.
 *
 */

#ifndef _TSP_DATA_H_
#define _TSP_DATA_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>

class ATSPDataC
{
public :
    /// Members.
    int _size;
    std::vector<std::vector<int>> _distances;

    /// Constructor(s).
    ATSPDataC(const std::string & filename, const bool & silent = true)
    {
        std::ifstream inf(filename);
        if (!inf)
        {
            _size = -1;
            return;
        }

        std::string temp;
        do
        {
            inf >> temp >> std::ws;
            if (temp == "DIMENSION:") inf >> _size >> std::ws;
        } while (temp != "EDGE_WEIGHT_SECTION");

        _distances.resize(_size,std::vector<int>(_size));

        for (size_t i = 0; i < _size; ++i)
        {
            for (size_t j = 0; j < _size; ++j)
            {
                inf >> _distances[i][j] >> std::ws;
            }
        }
        inf.close();

        if (!silent)
        {
            std::cout << "Size : " << _size << std::endl;
            std::cout << "Distances :" << std::endl;
            for (size_t i = 0; i < _size; ++i)
            {
                for (size_t j = 0; j < _size; ++j)
                {
                    std::cout << _distances[i][j] << " ";
                }
                std::cout << std::endl;
            }
        }
    }

    /// Destructor(s).
    ~ATSPDataC() = default;
};

#endif /// _TSP_DATA_H_

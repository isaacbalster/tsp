/*
 *
 * tsp.cpp
 *
 * Traveling salesmen problem (TSP).
 *
 * Created by Isaac Balster on 26/01/2023.
 *
 */

#include "data.h"
#include "model.h"

int main(int argc, char *argv[])
{
    /// Passing instance path to instantiate data class ATSPDataC.
    ATSPDataC data(argv[1]);

    /// Calling the desired model having the data class as single parameter.
    //mtzModel(data);
    lazyModel(data, "SEC", true);

    return 0;
}
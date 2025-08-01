/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of turbinesFoam, which is based on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fv::penalisedSource

Description
    Class implements the penalised momentum source for immersed boundary representation
    of a moving body.

SourceFiles
    penalisedSource.C

\*---------------------------------------------------------------------------*/

#ifndef solidMasker_H
#define solidMasker_H

#include "fvOption.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvCFD.H"
#include "OFstream.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "searchableSurfaces.H"
#include "Function1.H"
// #include "../penalisedSource.H"

// Forward declare penalisedSource (we only need to know it exists)
namespace Foam {
namespace fv {
    class penalisedSource;
}
}

namespace Foam
{
namespace fv
{
class solidMasker
{
protected:
    const fvMesh & 	mesh_;

    const Time& runTime_;
        
public:
    solidMasker
    (
        const fvMesh& mesh
    );

    virtual ~solidMasker();

    virtual void updateMask(penalisedSource& source) = 0;
};


} // End namespace fv
} // End namespace Foam


#endif

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

#ifndef cellCentresolidMasker_H
#define cellCentresolidMasker_H

#include "solidMasker.H"
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


namespace Foam {
namespace fv {
    class penalisedSource;
}
}


namespace Foam
{
namespace fv
{
class cellCentreSolidMasker 
:
    public solidMasker
{
protected:

    pointField searchPoints_;

    pointField transformedSearchPoints_;

public:
    cellCentreSolidMasker
    (
        const fvMesh& mesh
    );

    virtual ~cellCentreSolidMasker();

    virtual void updateMask(penalisedSource& source) override;
};


} // End namespace fv
} // End namespace Foam


#endif

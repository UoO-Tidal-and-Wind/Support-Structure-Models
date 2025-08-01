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

\*---------------------------------------------------------------------------*/

#include "fvMatrices.H"
#include "meshSearch.H"
#include "meshTools.H"
#include "unitConversion.H"


#include "cellCentreSolidMasker.H"
#include "../penalisedSource.H"



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::cellCentreSolidMasker::cellCentreSolidMasker
(
    const fvMesh& mesh
)
:
    solidMasker(mesh)
{
    searchPoints_ = mesh.C();
    transformedSearchPoints_ = searchPoints_;
}

void Foam::fv::cellCentreSolidMasker::updateMask(penalisedSource& source)
{
    if (source.isMoving())
    {

        vector rotationalDOF = source.getRotationalDOF();
        vector translationalDOF = source.getTranslationalDOF();
        vector centreOfRotation = source.getCentreOfRotation();

        tensor R0 = Rx(-degToRad(rotationalDOF[0]));
        tensor R1 = Rx(-degToRad(rotationalDOF[1]));
        tensor R2 = Rx(-degToRad(rotationalDOF[2]));
        tensor R = R2 & R1 & R0;

        // translate points
        vector shiftedCentre = centreOfRotation -translationalDOF;

        forAll(searchPoints_, pti)
        {
            transformedSearchPoints_[pti] = 
                (R & (searchPoints_[pti] - translationalDOF - shiftedCentre))
                 + shiftedCentre;
        }
    }

    List<List<volumeType>> volTypes;
    source.findVolumeType(transformedSearchPoints_, volTypes);

    volScalarField& solidMask = source.getSolidMask();
    solidMask = scalar(0.0);

    forAll(volTypes, listi)
    {
        forAll(transformedSearchPoints_, pti)
        {
            if (volTypes[listi][pti] == volumeType::INSIDE)
            {
                solidMask[pti] = 1.0;
            }
            else if (volTypes[listi][pti] == volumeType::OUTSIDE)
            {
                solidMask[pti] = Foam::max(0.0, solidMask[pti]);
            }
            else if (volTypes[listi][pti] == volumeType::MIXED)
            {
                solidMask[pti] = Foam::max(1.0, solidMask[pti]);
            }
            else if (volTypes[listi][pti] == volumeType::UNKNOWN)
            {
                Info << "Cell " << pti << " has unknown volume type." << endl;
            }
        }
    }
}



// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::cellCentreSolidMasker::~cellCentreSolidMasker()
{}


// ************************************************************************* //

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


#include "cellPointsSolidMasker.H"
#include "../penalisedSource.H"



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::cellPointsSolidMasker::cellPointsSolidMasker
(
    const fvMesh& mesh
)
:
    solidMasker(mesh)
{
    searchPoints_ = mesh.points();
    transformedSearchPoints_ = searchPoints_;
}

void Foam::fv::cellPointsSolidMasker::updateMask(penalisedSource& source)
{
    if (source.isMoving())
    {

        vector rotationalDOF = source.getRotationalDOF();
        vector translationalDOF = source.getTranslationalDOF();
        vector centreOfRotation = source.getCentreOfRotation();

        tensor R0 = Rx(-degToRad(rotationalDOF[0]));
        tensor R1 = Ry(-degToRad(rotationalDOF[1]));
        tensor R2 = Rz(-degToRad(rotationalDOF[2]));
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

    const faceList& ff = mesh_.faces();
    // for each cell
    forAll(mesh_.C(), celli)
    {   
        const cell & cc = mesh_.cells()[celli];
        labelList pLabels(cc.labels(ff));
        forAll(volTypes, listi)
        {
            int nVerticesIn = 0;
            int nVerticesOut = 0;
            forAll(pLabels,pLabelsi)
            {   
                if (volTypes[listi][pLabels[pLabelsi]] == volumeType::INSIDE)
                {
                    nVerticesIn++;
                }
                else if (volTypes[listi][pLabels[pLabelsi]] == volumeType::OUTSIDE)
                {
                    nVerticesOut++;
                }
                else if (volTypes[listi][pLabels[pLabelsi]] == volumeType::MIXED)
                {
                    nVerticesIn++;
                }
                else if (volTypes[listi][pLabels[pLabelsi]] == volumeType::UNKNOWN)
                {
                    Info << "Point " << pLabelsi << " has unknown volume type." << endl;
                }
            }

            scalar fraction = double(nVerticesIn) / double(nVerticesIn + nVerticesOut);
            solidMask[celli] = Foam::max(fraction, solidMask[celli]);

        }
    }
}



// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::cellPointsSolidMasker::~cellPointsSolidMasker()
{}


// ************************************************************************* //

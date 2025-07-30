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

#include "penalisedSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "OFstream.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "meshTools.H"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(penalisedSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        penalisedSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::fv::penalisedSource::read(const dictionary& dict)
{
    if (option::read(dict)) 
    {
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        coeffs_.lookup("penalisationFactor") >> penalisationFactor_;
        moving_ = coeffs_.lookupOrDefault<bool>("moving", false);
        coeffs_.lookup("baseVelocity") >> baseVelocity_;
        coeffs_.lookup("centreOfRotation") >> centreOfRotation_;

        dictionary geometryDict;
        geometryDict = dict.subDict("geometry");

        surfacesPtr_.reset(new searchableSurfaces(
            IOobject
            (
                "abc",                      // dummy name
                mesh_.time().constant(),     // instance
                //mesh.time().findInstance("triSurface", word::null),// instance
                "triSurface",               // local
                mesh_.time(),                // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            geometryDict,
            true
        ));
        
        return true;
    }
    else
    {
        return false;
    }
}

void Foam::fv::penalisedSource::writeData(Ostream& os) const
{
}

void Foam::fv::penalisedSource::updateSolidMask()
{
    // forAll(mesh_.C(), celli)
    // {
        // const point& pt = mesh_.C()[celli];
        

        // // // first check if cell is within the bounding box
        // if (boundingBox_.contains(pt))
        // {
        //     // point is inside bouding box so check if inside surface
            
        // }    
        // else
        // {
        //     // point is outside bounding box so set mask to -1
        //     solidMask_[celli] = -1.0;
        //     continue;
        // }
    // }

    forAll(surfacesPtr_->names(), geomi)
    {
        const searchableSurface& s = (*surfacesPtr_)[geomi];
        List<volumeType> volTypes;

        // perform based on cell centres (not perfect but will do for now)
        s.getVolumeType(mesh_.C(), volTypes);

        forAll(mesh_.C(), pti)
        {
            if (volTypes[pti] == volumeType::INSIDE)
            {
                solidMask_[pti] = 1.0;
            }
            else if (volTypes[pti] == volumeType::OUTSIDE)
            {
                // solidMask_[pti] = -1.0;
            }
            else if (volTypes[pti] == volumeType::MIXED)
            {
                // solidMask_[pti] = 0.0; 
            }
            else if (volTypes[pti] == volumeType::UNKNOWN)
            {
                Info << "Point " << pti << " has unkown volume type." << endl;
            }
        }

    }


}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::penalisedSource::penalisedSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    runTime_(mesh_.time()),
    searchEngine_(mesh_),
    bodyForce_
    (
        IOobject
        (
            "force." + name,
            mesh_.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForce",dimForce/dimVolume/dimDensity,vector::zero)
    ),
    solidMask_
    (
        IOobject
        (
            "solidMask." + name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("cellFinder",dimless,-1.0)
    )
{
    read(dict);
    
    Info << "Surfaces found:" << endl;
    forAll(surfacesPtr_->names(), geomi)
    {
        Info << "Name:    " << surfacesPtr_->names()[geomi] << endl;
        Info << "Number of points:    " << (*surfacesPtr_)[geomi].size() << endl;
    }
    Info << endl;
    
    boundingBox_ = surfacesPtr_->bounds();

    updateSolidMask();
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::penalisedSource::~penalisedSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::penalisedSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // zero the body forces
    bodyForce_ *= 0.0;

    // check if the source is moving
    if (moving_)
    {
        updateSolidMask();
    }
}



// ************************************************************************* //

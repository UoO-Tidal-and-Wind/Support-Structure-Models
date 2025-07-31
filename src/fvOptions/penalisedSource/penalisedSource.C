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
    // This function could be abstracted into a class to allow for
    // more advanced methods later on...
    
    solidMask_ = scalar(0.0);
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
                // Info << "Point: " << pti << " Volume Type: " << "INSIDE" << endl;
            }
            else if (volTypes[pti] == volumeType::OUTSIDE)
            {
                solidMask_[pti] = Foam::max(0.0, solidMask_[pti]);
            }
            else if (volTypes[pti] == volumeType::MIXED)
            {
                // bounday points
                solidMask_[pti] = Foam::max(1.0, solidMask_[pti]); 
            }
            else if (volTypes[pti] == volumeType::UNKNOWN)
            {
                Info << "Point " << pti << " has unkown volume type." << endl;
            }
        }

    }
}

void Foam::fv::penalisedSource::updateBodyVelocity()
{
    // bodyVelocity_ = dimensionedVector(baseVelocity_, dimVelocity);
    bodyVelocity_ = dimensionedVector("bodyVelocity",dimVelocity,baseVelocity_);
}

void Foam::fv::penalisedSource::updateBodyForce()
{
    // bodyForceLHSCoeff_ = -penalisationFactor_ * solidMask_;
    // bodyForceLHSCoeff_ *= mesh_.V();
    forAll(mesh_.C(), celli)
    {
        bodyForceLHSCoeff_[celli] = -penalisationFactor_ * solidMask_[celli] / mesh_.V()[celli];
        bodyForceRHS_[celli] = -penalisationFactor_ * solidMask_[celli] * bodyVelocity_[celli] / mesh_.V()[celli];
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
            "bodyForce." + name,
            mesh_.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForce",dimForce/dimVolume/dimDensity,vector::zero)
    ),
    bodyForceLHS_
    (
        IOobject
        (
            "bodyForceLHS." + name,
            mesh_.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForceLHS",dimForce/dimVolume/dimDensity,vector::zero)
    ),
    bodyForceRHS_
    (
        IOobject
        (
            "bodyForceRHS." + name,
            mesh_.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForceRHS",dimForce/dimVolume/dimDensity,vector::zero)
    ),
    bodyForceLHSCoeff_
    (
        IOobject
        (
            "bodyForceLHSCoeff." + name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("bodyForceLHSCoeff",dimless/dimTime,0.0)
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
        dimensionedScalar("solidMask",dimless,-1.0)
    ),
    bodyVelocity_
    (
        IOobject
        (
            "bodyVelocity." + name,
            mesh_.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyVelocity",dimVelocity,vector::zero)
    )
{
    read(dict);
    
    Info << "Surfaces found:" << endl;
    forAll(surfacesPtr_->names(), geomi)
    {
        Info << "Name:    " << surfacesPtr_->names()[geomi] << endl;
        Info << "Number of points:    " << (*surfacesPtr_)[geomi].size() << endl;
        Info << "Has volume Type:    " << (*surfacesPtr_)[geomi].hasVolumeType() << endl;
    }
    Info << endl;
    
    boundingBox_ = surfacesPtr_->bounds();

    updateSolidMask();
    updateBodyVelocity();
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
    // check if the source is moving
    if (moving_)
    {
        updateSolidMask();
        updateBodyVelocity();
        updateBodyForce();
    }

    // add term to the LHS of the momentum equation
    const volVectorField& U = eqn.psi();
    eqn -= fvm::Sp(bodyForceLHSCoeff_,U);
    eqn -= bodyForceRHS_;

    bodyForceLHS_ = bodyForceLHSCoeff_ * U;
    bodyForce_ = bodyForceLHS_ + bodyForceRHS_;
}



// ************************************************************************* //

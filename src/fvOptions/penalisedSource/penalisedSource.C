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
#include "solidMasker/solidMasker.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "OFstream.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "meshTools.H"
#include "unitConversion.H"
#include "tensor.H"
#include "clockValue.H"

#include "solidMasker/cellCentreSolidMasker.H"
#include "solidMasker/cellPointsSolidMasker.H"

#define TIMING_MSG(msg) if (showTiming_) Info << msg


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
        showTiming_ = coeffs_.lookupOrDefault<bool>("showTiming", false);
        coeffs_.lookup("referenceDensity") >> referenceDensity_;
        
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
        

        word maskerType = word(coeffs_.lookup("masker"));
        if (maskerType == "cellCentre")
        {
            solidMaskerPtr_.reset(new cellCentreSolidMasker(mesh_));
        }
        else if (maskerType == "cellPoints")
        {
            solidMaskerPtr_.reset(new cellPointsSolidMasker(mesh_));
        }
        else
        {
            FatalErrorInFunction
            << "Invalid masker selection: " 
            << maskerType << abort(FatalIOError);
        }
        

        if (moving_)
        {
            dictionary DOFDict;
            DOFDict = dict.subDict("degreesOfFreedom");
            
            DOFDict.lookup("centreOfRotation") >> centreOfRotation_;
            rotationalDOFFuncPtr_.reset(
                Function1<vector>::New("rotationalDOF", DOFDict)
            );
            translationalDOFFuncPtr_.reset(
                Function1<vector>::New("translationalDOF", DOFDict)
            );
            rotationalDOFVelocityFuncPtr_.reset(
                Function1<vector>::New("rotationalDOFVelocity", DOFDict)
            );
            translationalDOFVelocityFuncPtr_.reset(
                Function1<vector>::New("translationalDOFVelocity", DOFDict)
            );
        }
        
        return true;
    }
    else
    {
        return false;
    }
}

void Foam::fv::penalisedSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}

void Foam::fv::penalisedSource::updateSolidMask()
{

    solidMaskerPtr_->updateMask(*this);
}

void Foam::fv::penalisedSource::updateBodyVelocity()
{
    if (moving_)
    {
        forAll(mesh_.C(), celli)
        {
            vector transU = translationalDOFVelocity_; // same for all points
            if (solidMask_[celli] > 0.0)
            {
                // dimensionedVector tmp = dimensionedVector("bodyVelocity",dimVelocity,baseVelocity_);
                vector rotU = rotationalDOFVelocity_ ^ (mesh_.C()[celli] - centreOfRotation_);
                bodyVelocity_[celli] = baseVelocity_ + rotU + transU;
            }
            else
            {
                bodyVelocity_[celli] = vector::zero;
            }
        }
    }
    else
    {
        forAll(mesh_.C(), celli)
        {
            if (solidMask_[celli] > 0.0)
            {
                bodyVelocity_[celli] = baseVelocity_;
            }
            else
            {
                bodyVelocity_[celli] = vector::zero;
            }
        }
    }
}

void Foam::fv::penalisedSource::updateBodyForce()
{
    forAll(mesh_.C(), celli)
    {
        bodyForceLHSCoeff_[celli] = -penalisationFactor_ * solidMask_[celli] / mesh_.V()[celli];
        bodyForceRHS_[celli] = -penalisationFactor_ * solidMask_[celli] * bodyVelocity_[celli] / mesh_.V()[celli];
    }
}


vector Foam::fv::penalisedSource::getRotationalDOF() const
{
    return rotationalDOF_;
}

vector Foam::fv::penalisedSource::getTranslationalDOF() const
{
    return translationalDOF_;
}

vector Foam::fv::penalisedSource::getCentreOfRotation() const
{
    return centreOfRotation_;
}

bool Foam::fv::penalisedSource::isMoving() const
{
    return moving_;
}

void Foam::fv::penalisedSource::updateDOFs()
{
    if (moving_)
    {
        rotationalDOF_ = rotationalDOFFuncPtr_->value(runTime_.value());
        translationalDOF_ = translationalDOFFuncPtr_->value(runTime_.value());
        rotationalDOFVelocity_ = rotationalDOFVelocityFuncPtr_->value(runTime_.value());
        translationalDOFVelocity_ = translationalDOFVelocityFuncPtr_->value(runTime_.value());
    }
}

void Foam::fv::penalisedSource::findVolumeType
(
    const pointField& searchPoints,
    List<List<volumeType>>& volTypes
) const
{
    volTypes.setSize(surfacesPtr_->size());

    forAll(surfacesPtr_->names(), geomi)
    {
        const searchableSurface& s = (*surfacesPtr_)[geomi];
        s.getVolumeType(searchPoints, volTypes[geomi]);
    }
}


volScalarField& Foam::fv::penalisedSource::getSolidMask()
{
    return solidMask_;
}


void Foam::fv::penalisedSource::createOutputFile()
{
    fileName dir;
    if (Pstream::parRun())
    {
        dir = runTime_.path()/"../postProcessing/IBMForce"
            / runTime_.timeName();
    }
    else
    {
        dir = runTime_.path()/"postProcessing/IBMForce"
            / runTime_.timeName();
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }
    forceOutputFile_ = new OFstream(dir/name_);
    *forceOutputFile_ << "time\tbodyForceX\tbodyForceY\tbodyForceZ\t" << endl;
}

void Foam::fv::penalisedSource::writeOutput()
{
    vector totForce = vector(0,0,0);
    forAll(bodyForce_, cellI)
    {
        totForce += bodyForce_[cellI] * mesh_.V()[cellI] * referenceDensity_;
    }
    reduce(totForce, sumOp<vector>());

    // vector totForce = returnReduce(sum(bodyForce_).value(),sumOp<vector>()) * referenceDensity_;
    *forceOutputFile_ << runTime_.value()  << "\t"
             << totForce[0] << "\t" << totForce[1] 
             << "\t" << totForce[2] << "\t" << endl;
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
    ),
    rotationalDOF_(vector::zero),
    translationalDOF_(vector::zero),
    rotationalDOFVelocity_(vector::zero),
    translationalDOFVelocity_(vector::zero)
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
    
    updateDOFs();
    updateSolidMask();
    updateBodyVelocity();

    createOutputFile();
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
    clockTime timing;
    TIMING_MSG("Timing:\n");

    if (moving_)
    {
        TIMING_MSG("   updateDOFs"); 
        updateDOFs();
        TIMING_MSG(": Done (" << timing.timeIncrement() << "s)\n");
        
        TIMING_MSG("   updateSolidMask");
        updateSolidMask();
        TIMING_MSG(": Done (" << timing.timeIncrement() << "s)\n");

        TIMING_MSG("   updateBodyVelocity");
        updateBodyVelocity();
        TIMING_MSG(": Done (" << timing.timeIncrement() << "s)\n");
    }

    TIMING_MSG("   updateBodyForce"); 
    updateBodyForce();
    TIMING_MSG(": Done (" << timing.timeIncrement() << "s)\n");

    // add term to the LHS of the momentum equation
    const volVectorField& U = eqn.psi();
    eqn -= fvm::Sp(bodyForceLHSCoeff_, U);
    eqn -= bodyForceRHS_;

    bodyForceLHS_ = bodyForceLHSCoeff_ * U;
    bodyForce_ = bodyForceLHS_ + bodyForceRHS_;

    TIMING_MSG("   writeOutput"); 
    writeOutput();
    TIMING_MSG(": Done (" << timing.timeIncrement() << "s)\n");
    TIMING_MSG("addSup: Done (" << timing.elapsedTime() << "s)\n");
}



// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "IBMForces.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"

#define TIMING_MSG(msg) if (showTiming_) Info << msg


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(IBMForces, 0);
    addToRunTimeSelectionTable(functionObject, IBMForces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //



void Foam::functionObjects::IBMForces::reset()
{
    for (vector& v : sumIBMForces_)
    {
        v = vector::zero;
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::IBMForces::IBMForces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    runTime_(mesh_.time()),
    showTiming_(false),
    referenceDensity_(VGREAT),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        createOutputFiles();
        Log << endl;
    }
}


Foam::functionObjects::IBMForces::IBMForces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name),
    runTime_(mesh_.time()),
    showTiming_(false),
    referenceDensity_(VGREAT),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        createOutputFiles();
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::IBMForces::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    initialised_ = false;

    Info<< type() << " " << name() << ":" << endl;

    dict.lookup("referenceDensity") >> referenceDensity_;
    Info << "reference density: " << referenceDensity_ << endl;


    // read list of IBM names
    dict.lookup("IBMNames") >> IBMObjectNames_;
    // dict.readIfPresent<bool>("showTiming", showTiming_);
    showTiming_ = dict.getOrDefault<bool>("showTiming", false);
    
    // sums
    sumIBMForces_ = vectorList(IBMObjectNames_.size(), vector::zero);
    sumIBMInertForces_ = vectorList(IBMObjectNames_.size(), vector::zero);
    
    // pointers to fields
    IBMForceFieldPtrs_.resize(IBMObjectNames_.size());
    IBMInertForceFieldPtrs_.resize(IBMObjectNames_.size());
    
    // output files
    forceFiles_.resize(IBMObjectNames_.size());

    forAll(IBMObjectNames_, nameI)
    {
        const word name = IBMObjectNames_[nameI];
        Info << "object name: " << name << endl;

        const word bodyForceName = "bodyForce." + name;
        if (foundObject<volVectorField>(bodyForceName))
        {
            Info << "\tfound: " << bodyForceName << endl;
            IBMForceFieldPtrs_[nameI] = &lookupObject<volVectorField>(bodyForceName);
        }
        else
        {
            Info << "could not find the field:" << bodyForceName << endl;
            FatalErrorInFunction
            << "Could not find field: " 
            << bodyForceName << abort(FatalIOError);
        }
        
        const word inertBodyForceName = "bodyInertForce." + name;
        if (foundObject<volVectorField>(inertBodyForceName))
        {
            Info << "\tfound: " << inertBodyForceName << endl;
            IBMInertForceFieldPtrs_[nameI] = &lookupObject<volVectorField>(inertBodyForceName);
        }
        else
        {
            Info << "could not find the field:" << inertBodyForceName << endl;
            FatalErrorInFunction
            << "Could not find field: " 
            << inertBodyForceName << abort(FatalIOError);
        }
        
        reset();
    }
    
    return true;
}


bool Foam::functionObjects::IBMForces::execute()
{
    Log << type() << " " << name() << " write:" << nl;

    calculateIBMForces();
    return true;
}

void Foam::functionObjects::IBMForces::calculateIBMForces()
{
    forAll(IBMForceFieldPtrs_, fieldI)
    {
        const volVectorField& forceField = *IBMForceFieldPtrs_[fieldI];
        const volVectorField& inertForceField = *IBMInertForceFieldPtrs_[fieldI];

        vector totForce = vector::zero;
        vector totInertForce = vector::zero;
        forAll(forceField, cellI)
        {   
            totForce += (forceField[cellI]) *
                            mesh_.V()[cellI] * referenceDensity_;
            totInertForce += (inertForceField[cellI]) *
                            mesh_.V()[cellI] * referenceDensity_;
        } 
        
        reduce(totForce, sumOp<vector>());
        reduce(totInertForce, sumOp<vector>());

        totForce -= totInertForce;

        sumIBMForces_[fieldI] = totForce;
        sumIBMInertForces_[fieldI] = totInertForce;
    }
}

bool Foam::functionObjects::IBMForces::write()
{
    if (writeToFile())
    {
        Log << "    writing force files." << endl;
        writeToOutputFiles();
    }
    Log << endl;

    return true;
}

void Foam::functionObjects::IBMForces::createOutputFiles()
{
    fileName dirForce;
    if (Pstream::parRun())
    {
        dirForce = runTime_.path()/"../postProcessing/IBMForce"
            / runTime_.timeName();
    }
    else
    {
        dirForce = runTime_.path()/"postProcessing/IBMForce"
            / runTime_.timeName();
    }

    if (not isDir(dirForce))
    {
        mkDir(dirForce);
    }

    forAll(IBMObjectNames_, nameI)
    {
        forceFiles_[nameI] = new OFstream(dirForce/IBMObjectNames_[nameI]);
        *forceFiles_[nameI] << "time\tbodyForceX\tbodyForceY\tbodyForceZ\tbodyInertForceX\tbodyInertForceY\tbodyInertForceZ\t" << endl;
    }
}

void Foam::functionObjects::IBMForces::writeToOutputFiles()
{
    forAll(IBMObjectNames_, nameI)
    {
        *forceFiles_[nameI] << runTime_.value() << "\t"
                            << sumIBMForces_[nameI][0] << "\t"
                            << sumIBMForces_[nameI][1] << "\t"
                            << sumIBMForces_[nameI][2] << "\t"
                            << sumIBMInertForces_[nameI][0] << "\t"
                            << sumIBMInertForces_[nameI][1] << "\t"
                            << sumIBMInertForces_[nameI][2] << "\t"
                            << endl;

    }
}


// ************************************************************************* //

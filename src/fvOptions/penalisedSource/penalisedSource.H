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

#ifndef penalisedSource_H
#define penalisedSource_H

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declare penalisedSource
namespace Foam {
namespace fv {
    class solidMasker;
}
}


namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class penalisedSource Declaration
\*---------------------------------------------------------------------------*/

class penalisedSource
:
    public option
{

protected:

    // Protected data
        //- Pointers to runTime
        const Time& runTime_;

        //- Output file stream
        OFstream* forceOutputFile_;

        //- Penalisation factor (lambda)
        scalar  penalisationFactor_;

        //- Boolean to indicate if source is moving
        bool moving_;

        //- the base velocity of the solid
        vector baseVelocity_;

        //- Centre of rotation for transformations
        vector  centreOfRotation_;

        //- Mesh cell search engine
        meshSearch searchEngine_;

        //- Body force to apply to field
        volVectorField bodyForce_;

        //- Body force applied to LHS of momentum equations
        volVectorField bodyForceLHS_;

        //- Forcing applied to the RHS of momentum equations
        volVectorField bodyForceRHS_;
        
        //- Forcing applied to the LHS of momentum equations
        // (bodyForceLHS_ = bodyForceLHScoeff_ * U)
        volScalarField bodyForceLHSCoeff_;

        //- Mask for the solid cells
        volScalarField solidMask_;

        //- Field of body velocities
        volVectorField bodyVelocity_;

        // searchableSurfaces solidGeometry_;
        autoPtr<searchableSurfaces> surfacesPtr_;

        // - Function1 for rotational DOF
        autoPtr<Function1<vector>> rotationalDOFFuncPtr_;

        // - Function1 for translational DOF
        autoPtr<Function1<vector>> translationalDOFFuncPtr_;

        //- Function1 for rotational DOF velocity
        autoPtr<Function1<vector>> rotationalDOFVelocityFuncPtr_;

        //- Function1 for translational DOF velocity
        autoPtr<Function1<vector>> translationalDOFVelocityFuncPtr_;

        //- Rotational DOF (roll, pitch, yaw)
        vector rotationalDOF_;
        
        //- Translational DOF
        vector translationalDOF_;
        
        //- rotational DOF velocity 
        vector rotationalDOFVelocity_;

        //- translational DOF velocity
        vector translationalDOFVelocity_;

        //- Pointer to solidMasker
        autoPtr<solidMasker> solidMaskerPtr_;

        //- Flag to show function timings
        bool showTiming_;

    // Protected Member Functions

        //- Move the solid mask
        void updateSolidMask();

        //- Move the solid mask
        void updateBodyVelocity();

        //- Move the solid mask
        void updateBodyForce();

        void updateDOFs();

        //- Read dictionary
        bool read(const dictionary& dict);

        //- Write data to output stream
        virtual void writeData(Ostream& os) const;

        //- Create directory for saving outputs to
        void createOutputFile();

        //- Write to output file
        virtual void writeOutput();

public:

    //- Runtime type information
    TypeName("penalisedSource");

    // Selectors

    //- Return a reference to the selected fvOption model
    static autoPtr<penalisedSource> New
    (
        const word& name,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Constructor
    penalisedSource
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~penalisedSource();


    // Member functions

        // Access
            //- Return rotational DOF
            vector getRotationalDOF() const;

            //- Return translational DOF
            vector getTranslationalDOF() const;

            //- Return centre of rotation
            vector getCentreOfRotation() const;

            //- Return if the source is moving
            bool isMoving() const;

            //- For a given point field, store the corresponding volume type
            // in a List of Lists
            void findVolumeType(
                const pointField& searchPoints,
                List<List<volumeType>>& volTypes
            ) const;

            //- Return reference to solid mask
            volScalarField& getSolidMask();

        // Source term addition

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

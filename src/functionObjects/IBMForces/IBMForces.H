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

Class
    Foam::functionObjects::IBMForces

Group
    grpIBMForcesFunctionObjects

Description
    Computes IBMForces and moments over a given list of patches by integrating
    pressure and viscous IBMForces and moments, and optionally resistance IBMForces
    and moments from porous zones.

    IBMForces and moments are output in their total and constituent components:
    - total IBMForces and moments
    - pressure contributions
    - viscous contributions
    - porous resistance contributions (optional)

    IBMForces and moments can be computed and output in:
    - the global Cartesian coordinate system (default)
    - a user-defined Cartesian coordinate system

    Operands:
    \table
      Operand       | Type           | Location
      input         | -              | -
      output file   | dat            | postProcessing/\<FO\>/\<time\>/\<file\>s
      output field  | volVectorField | \<time\>/\<outField\>s
    \endtable

    where \c \<file\>s:
    \verbatim
      force.dat        | IBMForces
      moment.dat       | Moments
    \endverbatim

    where \c \<outField\>s:
    \verbatim
      <namePrefix>:force   | Force field
      <namePrefix>:moment  | Moment field
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    <namePrefix>
    {
        // Mandatory entries
        type                IBMForces;
        libs                (IBMForces);
        patches             (<wordRes>);

        // Optional entries
        directForceDensity  <bool>;
        porosity            <bool>;
        writeFields         <bool>;
        useNamePrefix       <bool>;

        // Conditional mandatory entries

            // if directForceDensity == true
            fD              <word>;


            // Cartesian coordinate system specification when
            // evaluating IBMForces and moments, either of the below

            // Define the centre of rotation
            // with implicit directions e1=(1 0 0) and e3=(0 0 1)
            CofR                (0 0 0); // Centre of rotation

            // Define local coordinate system by origin + axes
            origin              (0 0 0);
            e1                  (1 0 0);
            e3                  (0 0 1); // (e1, e2) or (e2, e3) or (e3, e1)

            // General coordinate system specification (always cartesian)
            coordinateSystem
            {
                origin          (0 0 0);
                rotation
                {
                    type        axes;
                    e3          (0 0 1);
                    e1          (1 0 0); // (e1, e2) or (e2, e3) or (e3, e1)
                }
            }

        // Conditional optional entries

            // if directForceDensity == false
            p               <word>;
            U               <word>;
            rho             <word>;
            rhoInf          <scalar>; // enabled if rho=rhoInf
            pRef            <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description               | Type | Reqd    | Deflt
      type       | Type name: IBMForces         | word | yes     | -
      libs       | Library name: IBMForces      | word | yes     | -
      patches    | Names of operand patches  | wordRes | yes  | -
      directForceDensity | Flag to directly supply force density <!--
                 -->                         | bool | no      | false
      porosity   | Flag to include porosity contributions | bool | no | false
      writeFields | Flag to write force and moment fields | bool | no | false
      useNamePrefix | Flag to include prefix for field names | bool | no | false
      coordinateSystem | Coordinate system specifier | dictionary | cndtnl | -
      CofR    | Centre of rotation          | vector | cndtnl   | -
      origin  | Origin of coordinate system | vector | cndtnl   | -
      e3      | e3 coordinate axis          | vector | cndtnl   | -
      e1      | e1 coordinate axis          | vector | cndtnl   | -
      fD      | Name of force density field | word   | cndtnl   | -
      p       | Name of pressure field      | word   | cndtnl   | p
      U       | Name of velocity field      | word   | cndtnl   | U
      rho     | Name of density field       | word   | cndtnl   | rho
      rhoInf  | Value of reference density  | scalar | cndtnl   | -
      pRef    | Value of reference pressure | scalar | cndtnl   | 0
    \endtable

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link writeFile.H \endlink
      - \link coordinateSystem.H \endlink

Note
  - For incompressible cases, set \c rho to \c rhoInf.
    You will then be required to provide a \c rhoInf
    value corresponding to the constant freestream density.
  - \c writeControl and \c writeInterval entries of function
    object do control when to output force and moment files and fields.
  - If a \c coordinateSystem entry exists, it is taken in favour of \c CofR.

SourceFiles
    IBMForces.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_IBMForces_H
#define Foam_functionObjects_IBMForces_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "coordinateSystem.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class IBMForces Declaration
\*---------------------------------------------------------------------------*/

class IBMForces
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data
        
        //- Pointers to runTime
        const Time& runTime_;
        
        // Fields

            //- Sum of patch pressure IBMForces
            List<vector> sumIBMForces_;

            List<vector> sumIBMInertForces_;

            List<const volVectorField*> IBMForceFieldPtrs_;

            List<const volVectorField*> IBMInertForceFieldPtrs_;

        // File streams

            List<OFstream*> forceFiles_;

        // Read from dictionary

            //- Reference density needed for incompressible calculations
            scalar referenceDensity_;

            //- List of IBM object names to calcualte force of
            wordList IBMObjectNames_;

            //- Flag of initialisation (internal)
            bool initialised_;


    // Protected Member Functions

        //- Initialise containers and fields
        void initialise();

        //- Reset containers and fields
        void reset();

    // Evaluation
        
        void calculateIBMForces();

    // I-O
        
        void createOutputFiles();

        void writeToOutputFiles();


public:

    //- Runtime type information
    TypeName("IBMForces");


    // Constructors

        //- Construct from Time and dictionary
        IBMForces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        IBMForces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );

        //- No copy construct
        IBMForces(const IBMForces&) = delete;

        //- No copy assignment
        void operator=(const IBMForces&) = delete;


    //- Destructor
    virtual ~IBMForces() = default;


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary& dict);

        //- Execute the function object
        virtual bool execute();

        //- Write to data files/fields and to streams
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

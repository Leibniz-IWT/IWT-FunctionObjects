/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::functionObjects::velocityAverageBox

Description
    Calculates the Average Velocity between two heights in the chosen direction.

Usage
    \table
        Property | Description                   | Required    | Default value
        type     | type name: velocityAverageBox | yes         |
        field    | the alpha field name          | no          | U
        axis     | the directional axis          | no          | x
        min      | the minimum height            | yes         |
        max      | the maximum height            | yes         |
    \endtable

Note
    Writing field 'velocityAverageBox' is done by default, but it can be overridden
    by defining an empty \c objects list. For details see writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    velocityAverageBox.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_velocityAverageBox_H
#define functionObjects_velocityAverageBox_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "volFieldsFwd.H"
//#include "HashSet.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class velocityAverageBox Declaration
\*---------------------------------------------------------------------------*/

class velocityAverageBox
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{

protected:

    // Protected data

        //- Runtime Variables
        word field_;
        word axis_;
        scalar min_;
        scalar max_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the surface area
        dimensionedScalar calcAverageVelocity
        (
            const word& field,
            const word& axis,
            const scalar& min,
            const scalar& max
        );


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        velocityAverageBox(const velocityAverageBox&);

        //- Disallow default bitwise assignment
        void operator=(const velocityAverageBox&);


public:

    //- Runtime type information
    TypeName("velocityAverageBox");


    // Constructors

        //- Construct from Time and dictionary
        velocityAverageBox
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~velocityAverageBox();


    // Member Functions

        //- Read the wallShearStress data
        virtual bool read(const dictionary&);

        //- Calculate the wall shear-stress
        virtual bool execute();

        //- Write the wall shear-stress
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

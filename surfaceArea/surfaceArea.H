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
    Foam::functionObjects::surfaceArea

Description
    Calculates the free Surface Area in a Two-Phase Simulation.

Usage
    \table
        Property | Description               | Required    | Default value
        type     | type name: wallShearStress | yes        |
        patches  | list of patches to process | no         | all wall patches
    \endtable

Note
    Writing field 'wallShearStress' is done by default, but it can be overridden
    by defining an empty \c objects list. For details see writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    surfaceArea.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceArea_H
#define functionObjects_surfaceArea_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
//#include "volFieldsFwd.H"
//#include "HashSet.H"
#include "sampledIsoSurfaceCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class wallShearStress Declaration
\*---------------------------------------------------------------------------*/

class surfaceArea
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{

protected:

    // Protected data

        //- Runtime Variables
        sampledIsoSurfaceCell isoSurf_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the surface area
        scalar calcSurfaceArea
        (
            const sampledIsoSurfaceCell& isoSurf
        );


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        surfaceArea(const surfaceArea&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceArea&);


public:

    //- Runtime type information
    TypeName("surfaceArea");


    // Constructors

        //- Construct from Time and dictionary
        surfaceArea
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~surfaceArea();


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

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

\*---------------------------------------------------------------------------*/

#include "velocityAverageBox.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(velocityAverageBox, 0);
    addToRunTimeSelectionTable(functionObject, velocityAverageBox, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::velocityAverageBox::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Velocity Average in Box");
    writeCommented(file(), "Time");
    writeTabbed(file(), "AverageVelocity");
    file() << endl;
}


Foam::dimensionedScalar
Foam::functionObjects::velocityAverageBox::calcAverageVelocity
(
    const word& field,
    const word& axis,
    const scalar& min,
    const scalar& max
)
{
    volVectorField velocity
    (
        IOobject
        (
            "velocity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("velocity",dimensionSet(0,1,-1,0,0,0,0),{0,0,0})
    );
    
    if (obr_.foundObject<volVectorField>(field_))
    {
        Info<< "    getting " << field_ << " from object registry" << nl;
        velocity = obr_.lookupObject<volVectorField>(field_);
    }
    else
    {
        Info<< "    reading " << field_ << " from file" << nl;
        IOobject velocityHeader
        (
            field_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ
        );
        
        if (velocityHeader.typeHeaderOk<volVectorField>(true))
        {
            velocity = volVectorField(velocityHeader, mesh_);
        }
    }
    
    volScalarField positions
    (
        IOobject
        (
            "positions",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("positions",dimensionSet(0,1,0,0,0,0,0), 0.0)
    );
    
    scalar ax = vector::X;
    if (axis == "y" or axis == "Y")
    {
        Info<< "    setting direction: Y-axis" << nl << endl;
        ax = vector::Y;
    }
    else if (axis == "z" or axis == "Z")
    {
        Info<< "    setting direction: Z-axis" << nl << endl;
        ax = vector::Z;
    }
    else
    {
        Info<< "    setting direction: X-axis" << nl << endl;
    }
    
    forAll(positions,iter)
    {
        positions[iter]=mesh_.C().internalField()[iter].component(ax);
    }

    dimensionedScalar average = 0.0;
    int count = 0;
    forAll(positions,iter)
    {
        if (positions[iter] <= max && positions[iter] >= min)
        {
            count+=1;
            average+=mag(velocity[iter]);
        }
    }
    
    average=average/count;
    
    return average;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::velocityAverageBox::velocityAverageBox
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    field_("U"),
    axis_("x"),
    min_(readScalar(dict.lookup("min"))),
    max_(readScalar(dict.lookup("max")))
{
    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::velocityAverageBox::~velocityAverageBox()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::velocityAverageBox::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);
    dict.readIfPresent("field", field_);
    dict.readIfPresent("axis", axis_);
    
    return true;
}


bool Foam::functionObjects::velocityAverageBox::execute()
{
    //Info<< "    updating isoSurfaceCell" << nl << endl;

    return true;
}


bool Foam::functionObjects::velocityAverageBox::write()
{
    Log << type() << " " << name() << " write:" << nl;

    logFiles::write();

    dimensionedScalar average = calcAverageVelocity(field_, axis_, min_, max_);
    reduce(average, sumOp<dimensionedScalar>());
    
    if (Pstream::master())
    {
        file() << mesh_.time().timeName() << tab << average.value() << endl;
    }

    return true;
}


// ************************************************************************* //

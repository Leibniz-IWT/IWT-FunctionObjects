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

#include "phaseVolumeBox.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(phaseVolumeBox, 0);
    addToRunTimeSelectionTable(functionObject, phaseVolumeBox, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::phaseVolumeBox::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Phase Volume in Box");
    writeCommented(file(), "Time");
    writeTabbed(file(), "Volume");
    file() << endl;
}


Foam::dimensionedScalar
Foam::functionObjects::phaseVolumeBox::calcPhaseVolume
(
    const word& field,
    const word& axis,
    const scalar& min,
    const scalar& max
)
{
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alpha",dimensionSet(0,0,0,0,0,0,0),0)
    );
    
    if (obr_.foundObject<volScalarField>(field_))
    {
        Info<< "    getting " << field_ << " from object registry" << nl;
        alpha = obr_.lookupObject<volScalarField>(field_);
    }
    else
    {
        Info<< "    reading " << field_ << " from file" << nl;
        IOobject alphaHeader
        (
            field_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ
        );
        
        if (alphaHeader.typeHeaderOk<volScalarField>(true))
        {
            alpha = volScalarField(alphaHeader, mesh_);
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

    dimensionedScalar volume = 0.0;
    forAll(positions,iter)
    {
        if (positions[iter] <= max && positions[iter] >= min)
        {
            volume+=mesh_.V()[iter]*alpha[iter];
        }
    }
    
    return volume;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::phaseVolumeBox::phaseVolumeBox
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    field_("alpha.water"),
    axis_("x"),
    min_(readScalar(dict.lookup("min"))),
    max_(readScalar(dict.lookup("max")))
{
    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::phaseVolumeBox::~phaseVolumeBox()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::phaseVolumeBox::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);
    dict.readIfPresent("field", field_);
    dict.readIfPresent("axis", axis_);
    
    return true;
}


bool Foam::functionObjects::phaseVolumeBox::execute()
{
    //Info<< "    updating isoSurfaceCell" << nl << endl;

    return true;
}


bool Foam::functionObjects::phaseVolumeBox::write()
{
    Log << type() << " " << name() << " write:" << nl;

    logFiles::write();

    dimensionedScalar volume = calcPhaseVolume(field_, axis_, min_, max_);
    reduce(volume, sumOp<dimensionedScalar>());
    
    if (Pstream::master())
    {
        file() << mesh_.time().timeName() << tab << volume.value() << endl;
    }

    //Log << mesh_.time().value() << ", "  << volume << endl;

    //Log << endl;

    return true;
}


// ************************************************************************* //

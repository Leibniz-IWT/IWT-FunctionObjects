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

#include "surfaceArea.H"
//#include "volFields.H"
//#include "surfaceFields.H"
//#include "turbulentTransportModel.H"
//#include "turbulentFluidThermoModel.H"
//#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "sampledIsoSurfaceCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(surfaceArea, 0);
    addToRunTimeSelectionTable(functionObject, surfaceArea, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::surfaceArea::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Free Surface Area");
    writeCommented(file(), "Time");
    writeTabbed(file(), "surfaceArea");
    file() << endl;
}


Foam::scalar
Foam::functionObjects::surfaceArea::calcSurfaceArea
(
    const sampledIsoSurfaceCell& isoSurf
)
{
    scalar surfaceArea = isoSurf.area();
    
    //Info << mesh_.time().value() << ", "  << surfaceArea << endl;

    return surfaceArea;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::surfaceArea::surfaceArea
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    isoSurf_("isoSurface", mesh_,  dict.subDict("isoSurfaceCell"))
{
    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::surfaceArea::~surfaceArea()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::surfaceArea::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);
    
    return true;
}


bool Foam::functionObjects::surfaceArea::execute()
{
    Info<< "    updating isoSurfaceCell" << nl << endl;
    
    isoSurf_.update();

    return true;
}


bool Foam::functionObjects::surfaceArea::write()
{
    Log << type() << " " << name() << " write:" << nl;

    logFiles::write();

    scalar surfaceArea = calcSurfaceArea(isoSurf_);
    if (Pstream::master())
    {
        file() << mesh_.time().timeName() << tab << surfaceArea << endl;
    }

    //Log << mesh_.time().value() << ", "  << surfaceArea << endl;

    //Log << endl;

    return true;
}


// ************************************************************************* //

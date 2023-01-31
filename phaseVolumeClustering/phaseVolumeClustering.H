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
    Foam::functionObjects::phaseVolumeClustering

Description
    Calculates the Volume of a Phase between two heights in the chosen direction in a Two-Phase Simulation.

Usage
    \table
        Property            | Description                      | Required    | Default value
        type                | type name: phaseVolumeClustering | yes         |
        cellDistance        | the distance between two cells   | yes         |
        field               | the field name                   | no          | alpha.water
        newClusterThreshold | the threshold to start a cluster | no          | 1
        acceptanceTolerance | the threshold to accept a cell   | no          | 0.01
        cutoff              | activate cutoff                  | no          | false
        cutoffAxis          | the cutoff axis                  | no          | x
        cutoffValue         | the cutoff value                 | no          | 0
        cutoffAboveValue    | cutoff above or below value      | no          | true
    \endtable

Note
    Writing field 'phaseVolumeClustering' is done by default, but it can be overridden
    by defining an empty \c objects list. For details see writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    phaseVolumeClustering.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_phaseVolumeClustering_H
#define functionObjects_phaseVolumeClustering_H

#include "fvMeshFunctionObject.H"
//#include "logFiles.H"
//#include "writeLocalObjects.H"
#include "writeFile.H"
#include "writer.H"
#include "volFieldsFwd.H"
//#include "HashSet.H"
#include "dimensionedScalar.H"
//#include "scalar.H"
//#include "vector.H"
//#include "curveTools.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class phaseVolumeClustering Declaration
\*---------------------------------------------------------------------------*/

class phaseVolumeClustering
:
    public fvMeshFunctionObject//,
    //public logFiles,
    //public writeLocalObjects
{

protected:

    // Protected data

        //- Runtime Variables
        writeFile file_;
        word field_;
        scalar newClusterThreshold_;
        scalar acceptanceTolerance_;
        scalar cellDistance_;
        bool cutoff_;
        word cutoffAxis_;
        scalar cutoffValue_;
        bool cutoffAboveValue_;


    // Protected Member Functions

        //- File header information
        //virtual void writeFileHeader(const label i);
        
        void writeClusterFile
        (
            const writeFile& file,
            const word& field,
            const std::vector< Tuple2<dimensionedScalar, int> >& processedClusters
        );

        //- Calculate the surface area
        std::vector< Tuple2<dimensionedScalar, int> > calcClusters
        (
            const word& field,
            const scalar& newClusterThreshold,
            const scalar& acceptanceTolerance,
            const scalar& cellDistance,
            const bool& cutoff,
            const word& cutoffAxis,
            const scalar& cutoffValue,
            const bool& cutoffAboveValue
        );
        
        Foam::Tuple2<dimensionedScalar, int> processCluster
        (
            const std::vector<int>& cluster,
            const volScalarField& alpha
        );


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        phaseVolumeClustering(const phaseVolumeClustering&);

        //- Disallow default bitwise assignment
        void operator=(const phaseVolumeClustering&);


public:

    //- Runtime type information
    TypeName("phaseVolumeClustering");


    // Constructors

        //- Construct from Time and dictionary
        phaseVolumeClustering
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~phaseVolumeClustering();


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

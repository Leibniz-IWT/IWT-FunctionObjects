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

#include "phaseVolumeClustering.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(phaseVolumeClustering, 0);
    addToRunTimeSelectionTable(functionObject, phaseVolumeClustering, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

/*
void Foam::functionObjects::phaseVolumeClustering::writeFileHeader(const label i)
{
    // Add headers to output data
    //writeHeader(file_, "Volume Cell Clusters");
    //writeCommented(file_, "NoCells");
    //writeTabbed(file_, "Volume");
    //writeTabbed(file_, "SphereEquivalentDiameter");
    //file_ << endl;
}
*/
void Foam::functionObjects::phaseVolumeClustering::writeClusterFile
(
    const writeFile& file,
    const word& field,
    const std::vector< Tuple2<dimensionedScalar, int> >& processedClusters
)
{
    fileName outputPath = file.baseTimeDir();
    mkDir(outputPath);
    
    OFstream outputFile
    (
        outputPath/"clusters.dat"
    );
    
    Log << "    Writing Clusters of " << field << " to " << outputFile.name() << endl;
    
    outputFile << "# Volume Cell Clusters" << nl;
    outputFile << "# Volume" << tab << "NoCells" << tab << "SphereEquivalentDiameter" << nl;
    
    for (unsigned int i=0;i<processedClusters.size();i++)
    {
        dimensionedScalar volume = processedClusters[i].first();
        int NoCells = processedClusters[i].second();
        
        reduce(volume, sumOp<dimensionedScalar>());
        reduce(NoCells, sumOp<int>());
        
        dimensionedScalar sphereEquD = cbrt(6*volume/Foam::constant::mathematical::pi);
        
        outputFile << volume.value() << tab << NoCells << tab << sphereEquD.value() << endl;
    }
}

std::vector< Foam::Tuple2<Foam::dimensionedScalar, int> >
Foam::functionObjects::phaseVolumeClustering::calcClusters
(
    const word& field,
    const scalar& newClusterThreshold,
    const scalar& acceptanceTolerance,
    const scalar& cellDistance,
    const bool& cutoff,
    const word& cutoffAxis,
    const scalar& cutoffValue,
    const bool& cutoffAboveValue
)
{
    Info<< "    initiating Clustering" << nl;
    scalar clusterDistance = 1.05*sqrt((2*cellDistance*cellDistance)+cellDistance*cellDistance);
    
    //Info<< "    clusterDistance = " << clusterDistance << nl;
    
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
    
    scalar ax = vector::X;
    if (cutoffAxis == "y" or cutoffAxis == "Y")
    {
        Info<< "    setting direction: Y-axis" << nl;
        ax = vector::Y;
    }
    else if (cutoffAxis == "z" or cutoffAxis == "Z")
    {
        Info<< "    setting direction: Z-axis" << nl;
        ax = vector::Z;
    }
    else
    {
        Info<< "    setting direction: X-axis" << nl;
    }
    
    Info<< "    calculating Cell Clusters" << nl;
    std::vector< std::vector<int> > clusters;
    std::vector< std::vector<int> > killedClusters;
    
    std::vector<int> relevantCells;
    forAll(alpha,iter)
    {
        if (alpha[iter] >= acceptanceTolerance)
        {
            relevantCells.push_back(iter);
        }
    }
    
    Info<< "        relevant Cells: " << relevantCells.size() << nl;
    
    for (unsigned int iter=0;iter<relevantCells.size();iter++)
    {
        //Info<< "    currentCellExt: " << relevantCells[iter] << nl;
        if (alpha[relevantCells[iter]] >= newClusterThreshold)
        {
            bool inCluster = false;
            for (unsigned int i=0;i<clusters.size();i++)
            {
                for (unsigned int j=0;j<clusters[i].size();j++)
                {
                    if (relevantCells[iter] == clusters[i][j])
                    {
                        inCluster = true;
                    }
                }
            }
            
            for (unsigned int i=0;i<killedClusters.size();i++)
            {
                for (unsigned int j=0;j<killedClusters[i].size();j++)
                {
                    if (relevantCells[iter] == killedClusters[i][j])
                    {
                        inCluster = true;
                    }
                }
            }
            
            if (not inCluster)
            {
                bool allowedNewCluster = true;
                if (cutoff)
                {
                    if (cutoffAboveValue)
                    {
                        if (mesh_.C().internalField()[relevantCells[iter]].component(ax) >= cutoffValue)
                        {
                            allowedNewCluster = false;
                        }
                    }
                    else
                    {
                        if (mesh_.C().internalField()[relevantCells[iter]].component(ax) <= cutoffValue)
                        {
                            allowedNewCluster = false;
                        }
                    }
                }
                
                if (allowedNewCluster)
                {
                    std::vector<int> newCluster = { relevantCells[iter] };
                    bool killCluster = false;
                    bool newCell = true;
                    
                    while (newCell)
                    {
                        newCell = false;
                        for (unsigned int i=0;i<relevantCells.size();i++)
                        {
                            //Info<< "    currentCellInt: " << relevantCells[i] << nl;
                            for (unsigned int j=0;j<newCluster.size();j++)
                            {
                                vector pos1 = mesh_.C().internalField()[relevantCells[i]];
                                vector pos2 = mesh_.C().internalField()[newCluster[j]];
                                scalar dist = sqrt(pow((pos1.x()-pos2.x()),2)+pow((pos1.y()-pos2.y()),2)+pow((pos1.z()-pos2.z()),2));
                                if (dist <= clusterDistance)
                                {
                                    bool allowedNewCell = true;
                                    if (cutoff)
                                    {
                                        if (cutoffAboveValue)
                                        {
                                            if (mesh_.C().internalField()[relevantCells[i]].component(ax) >= cutoffValue)
                                            {
                                                allowedNewCell = false;
                                            }
                                        }
                                        else
                                        {
                                            if (mesh_.C().internalField()[relevantCells[i]].component(ax) <= cutoffValue)
                                            {
                                                allowedNewCell = false;
                                            }
                                        }
                                    }
                                    
                                    if (allowedNewCell)
                                    {
                                        bool inCurrentCluster = false;
                                        for (unsigned int k=0;k<newCluster.size();k++)
                                        {
                                            if (relevantCells[i] == newCluster[k])
                                            {
                                                inCurrentCluster = true;
                                            }
                                        }
                                        
                                        if (not inCurrentCluster)
                                        {
                                            //Info<< "    addingCellToCluster: " << relevantCells[i] << nl;
                                            newCluster.push_back(relevantCells[i]);
                                            //Info<< "    newClusterSize: " << newCluster.size() << nl;
                                            newCell = true;
                                        }
                                    }
                                    else
                                    {
                                        killCluster = true;
                                        break;
                                    }
                                }
                            }
                            if (killCluster)
                            {
                                break;
                            }
                        }
                        if (killCluster)
                        {
                            break;
                        }
                    }
                    if (not killCluster)
                    {
                        clusters.push_back(newCluster);
                        Info<< "        Cluster: " << clusters.size() << ", NoCells: " << newCluster.size() << nl;
                    }
                    else
                    {
                        killedClusters.push_back(newCluster);
                        //Info<< "        Killed Cluster: " << killedClusters.size() << ", NoCells: " << newCluster.size() << nl;
                    }
                }
            }
        } 
    }
    
    Info<< "    processing Cell Clusters" << nl;
    std::vector< Tuple2<dimensionedScalar, int> > processedClusters;
    
    for (unsigned int i=0;i<clusters.size();i++)
    {
        Tuple2<dimensionedScalar, int> clusterProperties;
        
        clusterProperties = processCluster(clusters[i], alpha);
        
        processedClusters.push_back(clusterProperties);
    }
    
    return processedClusters;
}

Foam::Tuple2<Foam::dimensionedScalar, int>
Foam::functionObjects::phaseVolumeClustering::processCluster
(
    const std::vector<int>& cluster,
    const volScalarField& alpha
)
{
    dimensionedScalar volume = 0.0;
    
    for (unsigned int i=0;i<cluster.size();i++)
    {
        volume+=mesh_.V()[cluster[i]]*alpha[cluster[i]];
    }
    
    Tuple2<dimensionedScalar, int> clusterProperties(volume, cluster.size());
    
    return clusterProperties; 
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::phaseVolumeClustering::phaseVolumeClustering
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    //logFiles(obr_, name),
    //writeLocalObjects(obr_, log),
    file_(obr_, name),
    field_("alpha.water"),
    newClusterThreshold_(1),
    acceptanceTolerance_(0.1),
    cellDistance_(readScalar(dict.lookup("cellDistance"))),
    cutoff_(false),
    cutoffAxis_("x"),
    cutoffValue_(0),
    cutoffAboveValue_(true)
{
    read(dict);
    //resetName(typeName);
    //resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::phaseVolumeClustering::~phaseVolumeClustering()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::phaseVolumeClustering::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    //writeLocalObjects::read(dict);
    dict.readIfPresent("field", field_);
    dict.readIfPresent("newClusterThreshold", newClusterThreshold_);
    dict.readIfPresent("acceptanceTolerance", acceptanceTolerance_);
    dict.readIfPresent("cutoff", cutoff_);
    dict.readIfPresent("cutoffAxis", cutoffAxis_);
    dict.readIfPresent("cutoffValue", cutoffValue_);
    dict.readIfPresent("cutoffAboveValue", cutoffAboveValue_);
    
    return true;
}


bool Foam::functionObjects::phaseVolumeClustering::execute()
{
    //Info<< "    updating isoSurfaceCell" << nl << endl;

    return true;
}


bool Foam::functionObjects::phaseVolumeClustering::write()
{
    Log << type() << " " << name() << " write:" << nl;
    
    std::vector< Tuple2<dimensionedScalar, int> > processedClusters = calcClusters(field_, newClusterThreshold_, acceptanceTolerance_, cellDistance_, cutoff_, cutoffAxis_, cutoffValue_, cutoffAboveValue_);
    
    if (Pstream::master())
    {
        writeClusterFile(file_, field_, processedClusters);
    }
    
    return true;
}


// ************************************************************************* //

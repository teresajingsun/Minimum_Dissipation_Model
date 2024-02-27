/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "QR.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  namespace incompressible
  {
    namespace LESModels
    {

      // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

      defineTypeNameAndDebug(QR, 0);
      addToRunTimeSelectionTable(LESModel, QR, dictionary);

      // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

      void QR::updateSubGridScaleFields(const volTensorField& gradU)
      {
	nuSgs_ = ck_*delta()*sqrt(k());
	nuSgs_.correctBoundaryConditions();
      }


      // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

      QR::QR
      (
       const volVectorField& U,
       const surfaceScalarField& phi,
       transportModel& transport,
       const word& turbulenceModelName,
       const word& modelName
       )
	:
	LESModel(modelName, U, phi, transport, turbulenceModelName),
	GenEddyVisc(U, phi, transport),

	ck_
	(
	 dimensioned<scalar>::lookupOrAddToDict
	 (
	  "ck",
	  coeffDict_,
	  0.094
	  )
	 )
      {
	updateSubGridScaleFields(fvc::grad(U));

	printCoeffs();
      }


      // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

      //- Return SGS kinetic energy
      //  calculated from the given velocity gradient
      tmp<volScalarField> QR::k(const tmp<volTensorField>& gradU) const
      {
	volSymmTensorField D(symm(gradU));
	volScalarField r(max(-det(D),dimensionedScalar("r",dimensionSet(0,0,-3,0,0,0,0),0.0)));
	volScalarField q(max((1.0/2.0)*dev(D) && D,dimensionedScalar("q",dimensionSet(0,0,-2,0,0,0,0),1e-15)));
	return sqr(delta()*mag(r)/q);//(2.0*ck_/ce_)*sqr(delta())*magSqr(dev(symm(gradU)));
      }

      void QR::correct(const tmp<volTensorField>& gradU)
      {
	GenEddyVisc::correct(gradU);
	updateSubGridScaleFields(gradU());
      }


      bool QR::read()
      {
	if (GenEddyVisc::read())
	  {
	    ck_.readIfPresent(coeffDict());

	    return true;
	  }
	else
	  {
	    return false;
	  }
      }


      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace LESModels
  } // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

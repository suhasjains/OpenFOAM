//Reinitialisation of Level Set function

psi.correctBoundaryConditions();


volScalarField psiNew=psi;

//pseudo time step
scalar deltaTau(deltaX/3);


//No of Cells along x,y,z
scalar x(readScalar(transportProperties.lookup("cellsAlongX")));
scalar y(readScalar(transportProperties.lookup("cellsAlongY")));
scalar z(readScalar(transportProperties.lookup("cellsAlongZ")));



//No of cells in a plane and the domain
scalar xy=x*y;
scalar xyz=x*y*z;


dimensionedScalar psiBoundary;


Info<< "delta tau = " << deltaTau << nl << endl;

       
for (int reinit=0; reinit<nReinit; reinit++)
{                                  
  const faceList & ff = mesh.faces();
	const pointField & pp = mesh.points();

	//const volVectorField& C = mesh.C();



	//Smeared sign function
	volScalarField sign=psi/sqrt(pow(psi,2)+pow(epi,2));

	




	
	forAll ( mesh.cells(),i)
 {
    	const cell & cc = mesh.cells()[i];
    	labelList pLabels(cc.labels(ff));
    	pointField pLocal(pLabels.size(), vector::zero);

    	forAll (pLabels, pointi)
           pLocal[pointi] = pp[pLabels[pointi]];
     

	//calculating deltax,deltay,deltaz 
    	scalar xDim = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
	scalar yDim = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
	scalar zDim = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));

	

	//const unallocLabelList& neigh = mesh.cellCells()[i];

	//Info << " Cell: " << i << " Neighbours: " << neigh << " Center of Cell: " << C[i];
	

	//Declaration
	scalar l=0;
	scalar r=0;
	scalar s=0;
	scalar n=0;
	scalar t=0;
	scalar b=0;
	scalar iniPsi=0;
	scalar G=0;
	scalar ll=0;
	scalar rr=0;
	scalar ss=0;
	scalar nn=0;
	scalar left=0;
	scalar right=0;
	scalar north=0;
	scalar south=0;
	scalar lOne=0;
	scalar nOne=0;
	scalar rOne=0;
	scalar sOne=0;
	scalar lTwo=0;
	scalar rTwo=0;
	scalar nTwo=0;
	scalar sTwo=0;
	scalar Gamma1=0;
	scalar Gamma2=0;
	scalar rBeta1=0;
	scalar rBeta2=0;
	scalar lBeta1=0;
	scalar lBeta2=0;
	scalar nBeta1=0;
	scalar nBeta2=0;
	scalar sBeta1=0;
	scalar sBeta2=0;
	scalar rOmega1=0;
	scalar rOmega2=0;
	scalar lOmega1=0;
	scalar lOmega2=0;	
	scalar nOmega1=0;
	scalar nOmega2=0;
	scalar sOmega1=0;
	scalar sOmega2=0;
	scalar rWeight1=0;
	scalar rWeight2=0;
	scalar lWeight1=0;
	scalar lWeight2=0;
	scalar nWeight1=0;
	scalar nWeight2=0;
	scalar sWeight1=0;
	scalar sWeight2=0;
	

	//psi at 1st Pseudo time step
	if(reinit==0)   iniPsi=psi[i];



	if(z>1)
  {


	// get faces of a cell
	const cell& faces = mesh.cells()[i];

	// loop over faces
	forAll(faces, facei)
    {
    	label faceI = faces[facei];

   	//find patch
    	label patchI = psi.mesh().boundaryMesh().whichPatch(faceI);

   	// if facei is in patch: patchI < -1
  	if (patchI > -1)
     {
      	faceI = psi.mesh().boundaryMesh()[patchI].whichFace(faceI);
  	psiBoundary = psi.boundaryField()[patchI][faceI];
	
	//Info <<" size:  "<< size; 
     }
    }

	
	
	//Differential Coefficients
	if((i+xy)<xyz)	t=(psi[i+xy]-psi[i])/zDim;  //top differential

	if((i-xy)>=0)	b=(psi[i]-psi[i-xy])/zDim;  //bottom differential

	if((i-1)>=0)	l=(psi[i]-psi[i-1])/xDim;   //left differential
	
	if((i+1)<xyz)  r=(psi[i+1]-psi[i])/xDim;   //right differential

	if((i-x)>=0) 	s=(psi[i]-psi[i-x])/yDim;   //south differential

	if((i+x)<xyz) 	n=(psi[i+x]-psi[i])/yDim;   //north differential



	//Godunov Coefficients
	if(iniPsi>0)       G=Foam::sqrt(max(pow(max(l,0),2),pow(min(r,0),2))+max(pow(max(s,0),2),pow(min(n,0),2))+max(pow(max(b,0),2),pow(min(t,0),2)))-1;

	else if(iniPsi<0)  G=Foam::sqrt(max(pow(min(l,0),2),pow(max(r,0),2))+max(pow(min(s,0),2),pow(max(n,0),2))+max(pow(min(b,0),2),pow(max(t,0),2)))-1;

	else if(iniPsi==0) G=0;
  }



	

	else if(z==1)
    {

/*
	// get faces of a cell
	const cell& faces = mesh.cells()[i];

	// loop over faces
	forAll(faces, faceIndex)
    {
    	label faceI = faces[faceIndex];

   	//find patch
    	label patchI = psi.mesh().boundaryMesh().whichPatch(faceI);

   	// if facei is in patch: patchI < -1
  	if (patchI > -1)
        {
      	   //faceI = psi.mesh().boundaryMesh()[patchI].whichFace(faceI);
	
	  if (!isA<emptyFvPatch>(psi.boundaryField()[patchI]))
	    Info <<" value "  <<    psi.boundaryField()[patchI];
	
  	   //psiBoundary = psi.boundaryField()[patchI][faceI];

	}
	
	

     
    }
*/




	//Differential Coefficients
	
	//left
	if(((i-1)>0 )&&(i/x==0))
	{
	        l = 0;  
		ll= 0;
	}

	if(((i-2)>0 )&&((i-1)/x==0))
	{
	        l=(psi[i]-psi[i-1])/xDim;  
		ll=0;
	}

	else   
	{
		l=(psi[i]-psi[i-1])/xDim;
		ll=(psi[i-1]-psi[i-2])/xDim;
	} 
	


	//right
	if(((i+1)<xy)&&((i+1)/x==0))
	{
	        r = 0;  
		rr= 0;
	}

	if(((i+2)<xy)&&((i+2)/x==0))
	{
	        r=(psi[i+1]-psi[i])/xDim;  
		rr=0;
	}

	else   
	{
		r=(psi[i+1]-psi[i])/xDim;
		rr=(psi[i+2]-psi[i+1])/xDim;
	} 



	//south
	if((i-x)>=x)
	{
	        s=(psi[i]-psi[i-x])/yDim;
		ss=(psi[i-x]-psi[i-(2*x)])/yDim;
	}

	if((i-x)<x)
	{
	        s=(psi[i]-psi[i-x])/yDim;  
		ss=0;
	}

	if(i<x)   
	{
		s=0;
		ss=0;
	} 

	

	//north
	if((i+x)>=(xy-x))
	{
	        n=(psi[i+x]-psi[i])/yDim;  
		nn=0;
	}

	if((i+x)>=xy)
	{
	        n=0;  
		nn=0;
	}

	else   
	{
		n=(psi[i+x]-psi[i])/yDim;
		nn=(psi[i+(2*x)]-psi[i+x])/yDim;
	} 

	


	Gamma1=2/3;
	Gamma2=1/3;  
	

	rBeta1=pow((rr-r),2);
	rBeta2=pow((r-l) ,2);
	lBeta1=pow((ll-l),2);
	lBeta2=pow((l-r) ,2);
	nBeta1=pow((nn-n),2);
	nBeta2=pow((n-s) ,2);
	sBeta1=pow((ss-s),2);
	sBeta2=pow((s-n) ,2);

	rOmega1=Gamma1/pow((1e-7+rBeta1),2);
	rOmega2=Gamma2/pow((1e-7+rBeta2),2);	//1e-06
	lOmega1=Gamma1/pow((1e-7+lBeta1),2);
	lOmega2=Gamma2/pow((1e-7+lBeta2),2);
	nOmega1=Gamma1/pow((1e-7+nBeta1),2);
	nOmega2=Gamma2/pow((1e-7+nBeta2),2);
	sOmega1=Gamma1/pow((1e-7+sBeta1),2);
	sOmega2=Gamma2/pow((1e-7+sBeta2),2);

	rWeight1=rOmega1/(rOmega1+rOmega2);
	rWeight2=rOmega2/(rOmega1+rOmega2);
	lWeight1=lOmega1/(lOmega1+lOmega2);
	lWeight2=lOmega2/(lOmega1+lOmega2);
	nWeight1=nOmega1/(nOmega1+nOmega2);
	nWeight2=nOmega2/(nOmega1+nOmega2);
	sWeight1=sOmega1/(sOmega1+sOmega2);
	sWeight2=sOmega2/(sOmega1+sOmega2);

	rOne=0.5*r+0.5*rr;
	rTwo=1.5*r-0.5*l;
	lOne=0.5*l+0.5*ll;
	lTwo=1.5*l-0.5*r;
	nOne=0.5*n+0.5*nn;
	nTwo=1.5*n-0.5*s;
	sOne=0.5*s+0.5*ss;
	sTwo=1.5*s-0.5*n;

	
	right=rWeight1*rOne+rWeight2*rTwo;
	left =lWeight1*lOne+lWeight2*lTwo;
	north=nWeight1*nOne+nWeight2*nTwo;
	south=sWeight1*sOne+sWeight2*sTwo;
	


	//Godunov Coefficients
	
	if(iniPsi>0)       G = Foam::sqrt(max(pow(max(l,0),2),pow(min(r,0),2))+max(pow(max(s,0),2),pow(min(n,0),2)))-1;  

	else if(iniPsi<0)  G = Foam::sqrt(max(pow(min(l,0),2),pow(max(r,0),2))+max(pow(min(s,0),2),pow(max(n,0),2)))-1;
	
	else if(iniPsi==0) G = 0;

    }	
	

	
	//Reinitilisation step
	psi[i] = psi[i] - deltaTau * sign[i] * G; 
    	
 }

	//psi=psiNew;

	psi.correctBoundaryConditions();

	
}

interface.correct();

unsigned int nparams = 10;//number of dimensions

unsigned int nstat = 20000;//number of tries
TVectorD last(nparams);
TVectorD curr(nparams);
//should also initialise last here

int jobid = 0;
TRandom3* rnd = new TRandom3(359895+9836*jobid);

Bool_t accepted;
Bool_t inchain;
unsigned int index = 0;


unsigned int ntested = 0;
unsigned int naccepted = 0;
TMatrixDSym identity(nparams);
identity.UnitMatrix();

unsigned int nlast = 200;
std::vector<bool> lastaccepted;

TMatrixDSym S1(nparams);
S1.Zero();
S1 = identity;

TMatrixDSym SNminusone(nparams);
SNminusone = S1;
TMatrixDSym SN(nparams);
SN.Zero();

double minllh = getllh(last);
//last.negllh = minllh;

double alphastar = 0.1; //forced acceptance rate

for (unsigned int i=0; i<nstat; i++)
    {

      //use values of last for current, then vary
      curr = last;

      TVectorD WN(nparams);
      for (unsigned int j=0; j<WN.GetNrows(); j++)
	WN[j] = rnd->Gaus(0.0, 1.0);//could also use students t distribution

      TVectorD SW(SNminusone*WN);
      curr += SW;

      //get neg log lh, can also save this for last and current event to be more optimal
      //negllh = getllh(curr);

      accepted = false;
      double r = rnd->Rndm();
      double alpha = std::min(1.0, exp(getllh(last)-getllh(curr));//can cache lhs...
      double acceptrate = double(naccepted)/double(ntested);
      if (r < alpha)
	{
	  accepted = true;
	  naccepted++;
	  //success, candidate filled below
	  last = curr;
	}
      else
	{
	  accepted = false;
	  inchain = false;
	  //first save rejected candidate to file
	  //TODO IMPLEMENT
	  //t->Fill();
	  //then reset to last candidate
	  curr = last;
	}
      inchain = true;

      //Write candidate to file here
      //fill output file (this adds the previous candidate again if the current one failed)
      //TODO IMPLEMENT
      //t->Fill();

      lastaccepted.insert(lastaccepted.begin(), accepted);
      double lastratio = 0.0;
      if (lastaccepted.size() > nlast)
	{
	  //lastaccepted.push_back(accepted);
	  unsigned int nlastaccepted = 0;
	  for (unsigned int j=0; j<lastaccepted.size(); j++)
	    if (lastaccepted.at(j))
	      nlastaccepted++;
	  lastratio = double(nlastaccepted)/double(lastaccepted.size());
	  lastaccepted.pop_back();
	}
      //        if (accepted)
      //  	std::cout << "set " << i << " accepted, rate " << std::fixed << std::setprecision(3) << acceptrate << " last " << nlast << " " << lastratio << std::endl;
      //        else
      //  	std::cout << "set " << i << " rejected, rate " << std::fixed << std::setprecision(3) << acceptrate << " last " << nlast << " " << lastratio << std::endl;

      ntested++;

      //update S matrix!
      TMatrixDSym SNminusoneT(SNminusone);
      SNminusoneT.T();
      double etan = std::min(1.0, nparams*pow(double(i), -2.0/3.0));
      TMatrixDSym WNWNT(nparams);
      WNWNT.Zero();
      for (unsigned int row = 0; row < WNWNT.GetNrows(); row++)
	for (unsigned int col = 0; col < WNWNT.GetNcols(); col++)
	  WNWNT[row][col] = WN[row]*WN[col]/WN.Norm2Sqr();

      TMatrixDSym SNSNT(identity + WNWNT*etan*(alpha-alphastar));
      SNSNT = SNSNT.Similarity(SNminusone);

      //SNSNT = (SNminusone*identity*SNminusoneT);
      TDecompChol chol(SNSNT);
      bool success = chol.Decompose();
      assert(success);
      TMatrixD SNT = chol.GetU();
      TMatrixD SN(SNT);
      SN.T();

      for (unsigned int row = 0; row < SN.GetNrows(); row++)
	for (unsigned int col = 0; col < SN.GetNcols(); col++)
	  SNminusone(row,col) = SN(row,col);
      index++;
    }
  std::cout << ntested << " points were tested" << std::endl;
  std::cout << naccepted << " points were accepted" << std::endl;
  std::cout << "The accept fraction is " << double(naccepted)/double(ntested) << std::endl;

// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


/**
   HypoTestInverterResult class: holds the array of hypothesis test results and compute a confidence interval.
   Based on the RatioFinder code available in the RooStatsCms package developed by Gregory Schott and Danilo Piparo
   Ported and adapted to RooStats by Gregory Schott
   Some contributions to this class have been written by Matthias Wolf (error estimation)
**/


// include header file of this class 
#include "RooStats/HypoTestInverterResult.h"

#include "RooStats/HybridResult.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooMsgService.h"
#include "RooGlobalFunc.h"

#include "TF1.h"
#include "TGraphErrors.h"
#include <cmath>
#include "Math/BrentRootFinder.h"
#include "Math/WrappedFunction.h"

#include "TCanvas.h"
#include "RooStats/SamplingDistPlot.h"

#include <algorithm>

ClassImp(RooStats::HypoTestInverterResult)

using namespace RooStats;
using namespace RooFit;
using namespace std;


// initialize static value 
double HypoTestInverterResult::fgAsymptoticMaxSigma      = 5; 


HypoTestInverterResult::HypoTestInverterResult(const char * name ) :
   SimpleInterval(name),
   fUseCLs(false),
   fIsTwoSided(false),
   fInterpolateLowerLimit(true),
   fInterpolateUpperLimit(true),
   fFittedLowerLimit(false),
   fFittedUpperLimit(false),
   fInterpolOption(kLinear),
   fLowerLimitError(-1),
   fUpperLimitError(-1)
{
  // default constructor
   fLowerLimit = TMath::QuietNaN();
   fUpperLimit = TMath::QuietNaN();
}


HypoTestInverterResult::HypoTestInverterResult( const char* name,
						const RooRealVar& scannedVariable,
						double cl ) :
   SimpleInterval(name,scannedVariable,TMath::QuietNaN(),TMath::QuietNaN(),cl), 
   fUseCLs(false),
   fIsTwoSided(false),
   fInterpolateLowerLimit(true),
   fInterpolateUpperLimit(true),
   fFittedLowerLimit(false),
   fFittedUpperLimit(false),
   fInterpolOption(kLinear),
   fLowerLimitError(-1),
   fUpperLimitError(-1)
{
   // constructor 
   fYObjects.SetOwner();
   // put a cloned copy of scanned variable to set in the interval
   // to avoid I/O problem of the Result class - 
   // make the set owning the cloned copy (use clone istead of Clone to not copying all links)
   fParameters.removeAll();
   fParameters.takeOwnership();
   fParameters.addOwned(*((RooRealVar *) scannedVariable.clone(scannedVariable.GetName()) ));
}


HypoTestInverterResult::~HypoTestInverterResult()
{
   // destructor
   // no need to delete explictly the objects in the TList since the TList owns the objects
}


bool HypoTestInverterResult::Add( const HypoTestInverterResult& otherResult   )
{
   // Merge this HypoTestInverterResult with another
   // HypoTestInverterResult passed as argument
   // The merge is done by combining the HypoTestResult when the same point value exist in both results.
   // If results exist at different points these are added in the new result
   // NOTE: Merging of the expected p-values obtained with pseudo-data.
   //  When expected p-values exist in the result (i.e. when rebuild option is used when getting the expected 
   // limit distribuiton in the HYpoTestInverter) then the expected p-values are also merged. This is equivalent 
   // at merging the pseudo-data. However there can be an inconsistency if the expected p-values have been 
   // obtained with different toys. In this case the merge is done but a warning message is printed.


   int nThis = ArraySize();
   int nOther = otherResult.ArraySize();
   if (nOther == 0) return true;
   if (nOther != otherResult.fYObjects.GetSize() ) return false; 
   if (nThis != fYObjects.GetSize() ) return false; 

   // cannot merge in case of inconsistent memebrs
   if (fExpPValues.GetSize() > 0 && fExpPValues.GetSize() != nThis ) return false;
   if (otherResult.fExpPValues.GetSize() > 0 && otherResult.fExpPValues.GetSize() != nOther ) return false;

   oocoutI(this,Eval) << "HypoTestInverterResult::Add - merging result from " << otherResult.GetName()  
                                << " in " << GetName() << std::endl;

   bool addExpPValues = (fExpPValues.GetSize() == 0 && otherResult.fExpPValues.GetSize() > 0);
   bool mergeExpPValues = (fExpPValues.GetSize() > 0 && otherResult.fExpPValues.GetSize() > 0);

   if (addExpPValues || mergeExpPValues)
      oocoutI(this,Eval) << "HypoTestInverterResult::Add - merging also the expected p-values from pseudo-data" << std::endl;


   // case current result is empty
   // just make a simple copy of the other result
   if (nThis == 0) {
      fXValues = otherResult.fXValues;
      for (int i = 0; i < nOther; ++i) 
         fYObjects.Add( otherResult.fYObjects.At(i)->Clone() );
      for (int i = 0; i <  fExpPValues.GetSize() ; ++i)
         fExpPValues.Add( otherResult.fExpPValues.At(i)->Clone() );
   }
   // now do teh real mergemerge combining point with same value or adding extra ones 
   else { 
      for (int i = 0; i < nOther; ++i) {
         double otherVal = otherResult.fXValues[i];
         HypoTestResult * otherHTR = (HypoTestResult*) otherResult.fYObjects.At(i);
         if (otherHTR == 0) continue;
         bool sameXFound = false;
         for (int j = 0; j < nThis; ++j) {
            double thisVal = fXValues[j];
            
            // if same value merge the result 
            if ( (std::abs(otherVal) < 1  && TMath::AreEqualAbs(otherVal, thisVal,1.E-12) ) || 
                 (std::abs(otherVal) >= 1 && TMath::AreEqualRel(otherVal, thisVal,1.E-12) ) ) {
               HypoTestResult * thisHTR = (HypoTestResult*) fYObjects.At(j);               
               thisHTR->Append(otherHTR);
               sameXFound = true;
               if (mergeExpPValues) { 
                  ((SamplingDistribution*) fExpPValues.At(j))->Add( (SamplingDistribution*)otherResult.fExpPValues.At(i) );
                  // check if same toys have been used for the test statistic distribution
                  int thisNToys = (thisHTR->GetNullDistribution() ) ? thisHTR->GetNullDistribution()->GetSize() : 0;
                  int otherNToys = (otherHTR->GetNullDistribution() ) ? otherHTR->GetNullDistribution()->GetSize() : 0;
                  if (thisNToys != otherNToys ) 
                     oocoutW(this,Eval) << "HypoTestInverterResult::Add expexted p values have been generated with different toys " << thisNToys << " , " << otherNToys << std::endl;                  
               }
               break;
            }
         }
         if (!sameXFound) { 
            // add the new result 
            fYObjects.Add(otherHTR->Clone() );
            fXValues.push_back( otherVal);
         }
         // add in any case also when same x found
         if (addExpPValues)  
            fExpPValues.Add( otherResult.fExpPValues.At(i)->Clone() );
         
         
      }
   }
   
   if (ArraySize() > nThis) 
      oocoutI(this,Eval) << "HypoTestInverterResult::Add  - new number of points is " << fXValues.size()
                         << std::endl;
   else 
      oocoutI(this,Eval) << "HypoTestInverterResult::Add  - new toys/point is " 
                         <<  ((HypoTestResult*) fYObjects.At(0))->GetNullDistribution()->GetSize() 
                         << std::endl;
      
   return true;
}

bool HypoTestInverterResult::Add (Double_t x, const HypoTestResult & res) 
{
   // Add a single point result (an HypoTestResult)
   int i= FindIndex(x);
   if (i<0) {
      fXValues.push_back(x);
      fYObjects.Add(res.Clone());
   } else {
      HypoTestResult* r= GetResult(i);
      if (!r) return false;
      r->Append(&res);
   }
   return true;
}


 
double HypoTestInverterResult::GetXValue( int index ) const
{
 // function to return the value of the parameter of interest for the i^th entry in the results
  if ( index >= ArraySize() || index<0 ) {
    oocoutE(this,InputArguments) << "Problem: You are asking for an impossible array index value\n";
    return -999;
  }

  return fXValues[index];
}

double HypoTestInverterResult::GetYValue( int index ) const
{
// function to return the value of the confidence level for the i^th entry in the results
  if ( index >= ArraySize() || index<0 ) {
    oocoutE(this,InputArguments) << "Problem: You are asking for an impossible array index value\n";
    return -999;
  }

  if (fUseCLs) 
    return ((HypoTestResult*)fYObjects.At(index))->CLs();
  else 
     return ((HypoTestResult*)fYObjects.At(index))->CLsplusb();  // CLs+b
}

double HypoTestInverterResult::GetYError( int index ) const
{
// function to return the estimated error on the value of the confidence level for the i^th entry in the results
  if ( index >= ArraySize() || index<0 ) {
    oocoutE(this,InputArguments) << "Problem: You are asking for an impossible array index value\n";
    return -999;
  }

  if (fUseCLs) 
    return ((HypoTestResult*)fYObjects.At(index))->CLsError();
  else 
    return ((HypoTestResult*)fYObjects.At(index))->CLsplusbError();
}

double HypoTestInverterResult::CLb( int index ) const
{
  // function to return the observed CLb value  for the i-th entry
  if ( index >= ArraySize() || index<0 ) {
    oocoutE(this,InputArguments) << "Problem: You are asking for an impossible array index value\n";
    return -999;
  }
  return ((HypoTestResult*)fYObjects.At(index))->CLb();  
}

double HypoTestInverterResult::CLsplusb( int index ) const
{
  // function to return the observed CLs+b value  for the i-th entry
  if ( index >= ArraySize() || index<0 ) {
    oocoutE(this,InputArguments) << "Problem: You are asking for an impossible array index value\n";
    return -999;
  }
  return ((HypoTestResult*)fYObjects.At(index))->CLsplusb();  
}
double HypoTestInverterResult::CLs( int index ) const
{
  // function to return the observed CLs value  for the i-th entry
  if ( index >= ArraySize() || index<0 ) {
    oocoutE(this,InputArguments) << "Problem: You are asking for an impossible array index value\n";
    return -999;
  }
  return ((HypoTestResult*)fYObjects.At(index))->CLs();  
}

double HypoTestInverterResult::CLbError( int index ) const
{
  // function to return the error on the observed CLb value  for the i-th entry
  if ( index >= ArraySize() || index<0 ) {
    oocoutE(this,InputArguments) << "Problem: You are asking for an impossible array index value\n";
    return -999;
  }
  return ((HypoTestResult*)fYObjects.At(index))->CLbError();  
}

double HypoTestInverterResult::CLsplusbError( int index ) const
{
  // function to return the error on the observed CLs+b value  for the i-th entry
  if ( index >= ArraySize() || index<0 ) {
    oocoutE(this,InputArguments) << "Problem: You are asking for an impossible array index value\n";
    return -999;
  }
  return ((HypoTestResult*)fYObjects.At(index))->CLsplusbError();  
}
double HypoTestInverterResult::CLsError( int index ) const
{
   // function to return the error on the observed CLs value  for the i-th entry
   if ( index >= ArraySize() || index<0 ) {
      oocoutE(this,InputArguments) << "Problem: You are asking for an impossible array index value\n";
      return -999;
   }
   return ((HypoTestResult*)fYObjects.At(index))->CLsError();  
}


HypoTestResult* HypoTestInverterResult::GetResult( int index ) const
{
   // get the HypoTestResult object at the given index point
   if ( index >= ArraySize() || index<0 ) {
      oocoutE(this,InputArguments) << "Problem: You are asking for an impossible array index value\n";
      return 0;
   }
   
   return ((HypoTestResult*) fYObjects.At(index));
}

int HypoTestInverterResult::FindIndex(double xvalue) const
{
   // find the index corresponding at the poi value xvalue
   // If no points is found return -1
   // Note that a tolerance is used of 10^-12 to find the closest point
  const double tol = 1.E-12;
  for (int i=0; i<ArraySize(); i++) {
     double xpoint = fXValues[i];
     if ( (std::abs(xvalue) > 1 && TMath::AreEqualRel( xvalue, xpoint, tol) ) ||
          (std::abs(xvalue) < 1 && TMath::AreEqualAbs( xvalue, xpoint, tol) ) )
        return i; 
  }
  return -1;
}


struct InterpolatedGraph { 
   InterpolatedGraph( const TGraph & g, double target, const char * interpOpt) : 
      fGraph(g), fTarget(target), fInterpOpt(interpOpt) {}

   // return interpolated value for x - target
   double operator() (double x) const { 
      return fGraph.Eval(x, (TSpline*) 0, fInterpOpt) - fTarget;
   }
   const TGraph & fGraph;
   double fTarget;
   TString fInterpOpt;
}; 

double HypoTestInverterResult::GetGraphX(const TGraph & graph, double y0, bool lowSearch, double &axmin, double &axmax) const  {
   // return the X value of the given graph for the target value y0
   // the graph is evaluated using linear interpolation by default. 
   // if option = "S" a TSpline3 is used 



   TString opt = "";
   if (fInterpolOption == kSpline)  opt = "S";

   InterpolatedGraph f(graph,y0,opt);
   ROOT::Math::BrentRootFinder brf;
   ROOT::Math::WrappedFunction<InterpolatedGraph> wf(f);

   // find reasanable xmin and xmax if not given
   const double * y = graph.GetY(); 
   int n = graph.GetN();
   if (n < 2) {
      ooccoutE(this,Eval) << "HypoTestInverterResult::GetGraphX - need at least 2 points for interpolation (n=" << n << ")\n";
      return (n>0) ?  y[0] : 0;
   } 

   double varmin = - TMath::Infinity();
   double varmax = TMath::Infinity();
   const RooRealVar* var = dynamic_cast<RooRealVar*>( fParameters.first() );
   if (var) { 
       varmin = var->getMin();
       varmax = var->getMax();
   }

   double xmin = axmin; 
   double xmax = axmax;
   // case no range is given check if need to extrapolate to lower/upper values 
   if (axmin >= axmax ) { 
      xmin = graph.GetX()[0];
      xmax = graph.GetX()[n-1];

      // case we have lower/upper limits

      // find ymin and ymax  and corresponding values 
      double ymin = TMath::MinElement(n,y);
      double ymax = TMath::MaxElement(n,y);
      // do lower extrapolation 
      if ( (ymax < y0 && !lowSearch) || ( ymin > y0 && lowSearch) ) { 
         xmin = varmin;
      }
      // do upper extrapolation  
      if ( (ymax < y0 && lowSearch) || ( ymin > y0 && !lowSearch) ) { 
         xmax = varmax;
      }
   }
   else { 

#ifdef ISREALLYNEEDED //??
      // in case of range given, check if all points are above or below the confidence level line
      bool isCross = false;
      bool first = true; 
      double prod = 1;
      for (int i = 0; i< n; ++i) { 
         double xp, yp;
         graph.GetPoint(i,xp,yp);
         if (xp >= xmin && xp <= xmax) { 
            prod *= TMath::Sign(1., (yp - y0) );
            if (prod < 0 && !first) { 
               isCross = true; 
               break;
            }
            first = false;
         }
      }
      if (!isCross) { 
         return (lowSearch) ? varmin : varmax;
      }
#endif
   }
   


   brf.SetFunction(wf,xmin,xmax);
   brf.SetNpx(20);
   bool ret = brf.Solve(100, 1.E-6, 1.E-6);
   if (!ret) { 
      ooccoutE(this,Eval) << "HypoTestInverterResult - interpolation failed - return inf" << std::endl;
         return TMath::Infinity();
   }
   double limit =  brf.Root();

//#define DO_DEBUG
#ifdef DO_DEBUG
   if (lowSearch) std::cout << "lower limit search : ";
   else std::cout << "Upper limit search :  ";
   std::cout << "do interpolation between " << xmin << "  and " << xmax << " limit is " << limit << std::endl;
#endif

   // look in case if a new interseption exists
   // only when boundaries are not given
   if (axmin >= axmax) { 
      int index = TMath::BinarySearch(n, graph.GetX(), limit);
#ifdef DO_DEBUG
   std::cout << "do new interpolation dividing from " << index << "  and " << y[index] << std::endl;
#endif

      if (lowSearch && index >= 1 && (y[0] - y0) * ( y[index]- y0) < 0) { 
         //search if  another interseption exists at a lower value
         limit =  GetGraphX(graph, y0, lowSearch, graph.GetX()[0], graph.GetX()[index] );
      }
      else if (!lowSearch && index < n-2 && (y[n-1] - y0) * ( y[index+1]- y0) < 0) { 
         // another interseption exists at an higher value
         limit =  GetGraphX(graph, y0, lowSearch, graph.GetX()[index+1], graph.GetX()[n-1] );
      }
   }
   // return also xmin, xmax values
   axmin = xmin;
   axmax = xmax;

   return limit;
}

double HypoTestInverterResult::FindInterpolatedLimit(double target, bool lowSearch, double xmin, double xmax )
{
   // interpolate to find a limit value
   // Use a linear or a spline interpolation depending on the interpolation option

   ooccoutD(this,Eval) << "HypoTestInverterResult - " 
                       << "Interpolate the upper limit between the 2 results closest to the target confidence level" 
                       << std::endl;

   // variable minimum and maximum
   double varmin = - TMath::Infinity();
   double varmax = TMath::Infinity();
   const RooRealVar* var = dynamic_cast<RooRealVar*>( fParameters.first() );
   if (var) { 
       varmin = var->getMin();
       varmax = var->getMax();
   }

   if (ArraySize()<2) {
      double val =  (lowSearch) ? xmin : xmax;
      oocoutW(this,Eval) << "HypoTestInverterResult::FindInterpolatedLimit" 
                         << " - not enough points to get the inverted interval - return " 
                         <<  val << std::endl;
      fLowerLimit = varmin;
      fUpperLimit = varmax; 
      return (lowSearch) ? fLowerLimit : fUpperLimit;  
   }

   // sort the values in x 
   int n = ArraySize();
   std::vector<unsigned int> index(n );
   TMath::SortItr(fXValues.begin(), fXValues.end(), index.begin(), false); 
   // make a graph with the sorted point
   TGraph graph(n);
   for (int i = 0; i < n; ++i) 
      graph.SetPoint(i, GetXValue(index[i]), GetYValue(index[i] ) );


   //std::cout << " search for " << lowSearch << std::endl;
   

   // search first for min/max in the given range
   if (xmin >= xmax) {
      

      // search for maximum between the point
      double * itrmax = std::max_element(graph.GetY() , graph.GetY() +n);
      double ymax = *itrmax; 
      int iymax = itrmax - graph.GetY(); 
      double xwithymax = graph.GetX()[iymax];

      //std::cout << " max of y " << iymax << "  " << xwithymax << "  " << ymax << " target is " << target << std::endl;

      // look if maximum is above/belove target
      if (ymax > target) {
         if (lowSearch)  {
            if ( iymax > 0) { 
                  // low search (minimmum is first point or minimum range)
               xmin = ( graph.GetY()[0] <= target ) ? graph.GetX()[0] : varmin;  
               xmax = xwithymax;
               } 
            else { 
               // no room for lower limit
               fLowerLimit = varmin; 
               lowSearch = false; // search now for upper limits
            }
         }
         if (!lowSearch ) {
            // up search 
            if ( iymax < n-1 ) { 
               xmin = xwithymax; 
               xmax = ( graph.GetY()[n-1] <= target ) ? graph.GetX()[n-1] : varmax;
            }
            else { 
               // no room for upper limit
               fUpperLimit = varmax; 
               lowSearch = true; // search now for lower limits
               xmin = varmin; 
               xmax = xwithymax;
               }
         }
      }
      else { 
         // in case is below the target
         // find out if is a lower or upper search
         if (iymax <= (n-1)/2 ) { 
            lowSearch = false; 
            fLowerLimit = varmin; 
         }
         else { 
            lowSearch = true; 
            fUpperLimit = varmax;
         }
      }

#ifdef DO_DEBUG   
      std::cout << " found xmin, xmax  = " << xmin << "  " << xmax << " for search " << lowSearch << std::endl;
#endif

      // now come here if I have already found a lower/upper limit 
      // i.e. I am calling routine for the second time
#ifdef ISNEEDED
      // should not really come here
      if (lowSearch &&  fUpperLimit < varmax) {
         xmin = fXValues[ index.front() ];
         // find xmax (is first point before upper limit)
         int upI = FindClosestPointIndex(target, 2, fUpperLimit);
         if (upI < 1) return xmin; 
         xmax = GetXValue(upI);
      }
      else if (!lowSearch && fLowerLimit > varmin ) {
         // find xmin (is first point after lower limit)
         int lowI = FindClosestPointIndex(target, 3, fLowerLimit);
         if (lowI >= n-1) return xmax;
         xmin = GetXValue(lowI);
         xmax = fXValues[ index.back() ];
      }
#endif
   }

#ifdef DO_DEBUG   
   std::cout << "finding " << lowSearch << " limit betweem " << xmin << "  " << xmax << endl;
#endif

   double limit =  GetGraphX(graph, target, lowSearch, xmin, xmax);
   if (lowSearch) fLowerLimit = limit; 
   else fUpperLimit = limit;
   CalculateEstimatedError( target, lowSearch, xmin, xmax);

#ifdef DO_DEBUG
   std::cout << "limit is " << limit << std::endl;
#endif

   if (lowSearch && !TMath::IsNaN(fUpperLimit)) return fLowerLimit;
   if (!lowSearch && !TMath::IsNaN(fLowerLimit)) return fUpperLimit;

   // now perform the opposite search on the complement interval
   if (lowSearch) { 
      xmin = xmax;
      xmax = varmax;         
   } else { 
      xmax = xmin;         
      xmin = varmin;
   }
   double limit2 =  GetGraphX(graph, target, !lowSearch, xmin, xmax);
   CalculateEstimatedError( target, !lowSearch, xmin, xmax);

#ifdef DO_DEBUG
   std::cout << "other limit is " << limit2 << std::endl;
#endif
   if (!lowSearch) fLowerLimit = limit2; 
   else fUpperLimit = limit2;
 
   return (lowSearch) ? fLowerLimit : fUpperLimit;

}


int HypoTestInverterResult::FindClosestPointIndex(double target, int mode, double xtarget)
{
   // if mode = 0 
   // find closest point to target in Y, the object closest to the target which is 3 sigma from the target 
   // and has smaller error
   // if mode = 1
   // find 2 closest point to target in X and between these two take the one closer to the target
   // if mode = 2  as in mode = 1 but return the lower point not the closest one
   // if mode = 3  as in mode = 1 but return the upper point not the closest one

  int bestIndex = -1;
  int closestIndex = -1;
  if (mode == 0) { 
     double smallestError = 2; // error must be < 1
     double bestValue = 2;
     for (int i=0; i<ArraySize(); i++) {
        double dist = fabs(GetYValue(i)-target);
        if ( dist <3 *GetYError(i) ) { // less than 1 sigma from target CL
           if (GetYError(i) < smallestError ) { 
              smallestError = GetYError(i);
              bestIndex = i; 
           }
        }
        if ( dist < bestValue) { 
           bestValue = dist;  
           closestIndex = i; 
        }
     }
     if (bestIndex >=0) return bestIndex;
     // if no points found just return the closest one to the target
     return closestIndex;
  }
  // else mode = 1,2,3
  // find the two closest points to limit value
  // sort the array first 
  int n = fXValues.size();
  std::vector<unsigned int> indx(n);
  TMath::SortItr(fXValues.begin(), fXValues.end(), indx.begin(), false); 
  std::vector<double> xsorted( n);
  for (int i = 0; i < n; ++i) xsorted[i] = fXValues[indx[i] ];
  int index1 = TMath::BinarySearch( n, &xsorted[0], xtarget);

#ifdef DO_DEBUG
  std::cout << "finding closest point to " << xtarget << " is " << index1 << "  " << indx[index1] << std::endl; 
#endif

   // case xtarget is outside the range (bbefore or afterwards)
  if (index1 < 0) return indx[0]; 
  if (index1 >= n-1) return indx[n-1];  
  int index2 = index1 +1;
  
  if (mode == 2) return (GetXValue(indx[index1]) < GetXValue(indx[index2])) ? indx[index1] : indx[index2];
  if (mode == 3) return (GetXValue(indx[index1]) > GetXValue(indx[index2])) ? indx[index1] : indx[index2];
  // get smaller point of the two (mode == 1)
  if (fabs(GetYValue(indx[index1])-target) <= fabs(GetYValue(indx[index2])-target) )
     return indx[index1];
  return indx[index2];

}
        

Double_t HypoTestInverterResult::LowerLimit()
{
  if (fFittedLowerLimit) return fLowerLimit;
  //std::cout << "finding point with cl = " << 1-(1-ConfidenceLevel())/2 << endl;
  if ( fInterpolateLowerLimit ) { 
     // find both lower/upper limit
     if (TMath::IsNaN(fLowerLimit) )  FindInterpolatedLimit(1-ConfidenceLevel(),true);
  } else {
     //LM: I think this is never called
     fLowerLimit = GetXValue( FindClosestPointIndex((1-ConfidenceLevel())) );
  }
  return fLowerLimit;
}

Double_t HypoTestInverterResult::UpperLimit()
{
   //std::cout << "finding point with cl = " << (1-ConfidenceLevel())/2 << endl;
  if (fFittedUpperLimit) return fUpperLimit;
  if ( fInterpolateUpperLimit ) {
     if (TMath::IsNaN(fUpperLimit) )  FindInterpolatedLimit(1-ConfidenceLevel(),false);
  } else {
     //LM: I think this is never called
     fUpperLimit = GetXValue( FindClosestPointIndex((1-ConfidenceLevel())) );
  }
  return fUpperLimit;
}

Double_t HypoTestInverterResult::CalculateEstimatedError(double target, bool lower, double xmin, double xmax)
{
  // Return an error estimate on the upper(lower) limit.  This is the error on
  // either CLs or CLsplusb divided by an estimate of the slope at this
  // point.
   
  if (ArraySize()==0) {
     oocoutW(this,Eval) << "HypoTestInverterResult::CalculateEstimateError" 
                        << "Empty result \n";
    return 0;
  }

  if (ArraySize()<2) {
     oocoutW(this,Eval) << "HypoTestInverterResult::CalculateEstimateError" 
                        << " only  points - return its error\n"; 
     return GetYError(0);
  }
 
    // make a TGraph Errors with the sorted points
  std::vector<unsigned int> indx(fXValues.size());
  TMath::SortItr(fXValues.begin(), fXValues.end(), indx.begin(), false); 
  // make a graph with the sorted point
  TGraphErrors graph(ArraySize());
  int ip = 0;
  for (int i = 0; i < ArraySize(); ++i) {
     if ( (xmin < xmax) && ( GetXValue(indx[i]) >= xmin && GetXValue(indx[i]) <= xmax) ) {
        graph.SetPoint(ip, GetXValue(indx[i]), GetYValue(indx[i] ) );
        graph.SetPointError(ip, 0.,  GetYError(indx[i]) );
        ip++;
     }
  }

  double minX = xmin;
  double maxX = xmax;
  if (xmin >= xmax) {
     minX = fXValues[ indx.front() ];
     maxX = fXValues[ indx.back() ];
  }



  TF1 fct("fct", "exp([0] * x + [1] * x**2)", minX, maxX);
  TString type = (!lower) ? "upper" : "lower";

  // find the point closest to the limit
  double limit = (!lower) ? fUpperLimit : fLowerLimit;
  if (TMath::IsNaN(limit)) return 0; // cannot do if limit not computed  

#ifdef DO_DEBUG
  std::cout << "fitting for limit " << type << "between " << minX << " , " << maxX << std::endl;
  int fitstat = graph.Fit(&fct,"EX0");
#else
  int fitstat = graph.Fit(&fct,"Q EX0");
#endif 

  int index = FindClosestPointIndex(target, 1, limit);
  double theError = 0;
  if (fitstat == 0) { 
     double errY = GetYError(index);
     if (errY >  0) {
        double m = fct.Derivative( GetXValue(index) );
        theError = std::min(fabs( GetYError(index) / m), maxX-minX);     
     }
  }
  else { 
     oocoutE(this,Eval) << "HypoTestInverter::CalculateEstimatedError - cannot estimate  the " << type << " limit error " << std::endl;
  }
  if (lower) 
     fLowerLimitError = theError;
  else
     fUpperLimitError = theError;

#ifdef DO_DEBUG
  std::cout << "closes point to the limit is " << index << "  " << GetXValue(index) << " and has error " << GetYError(index) << std::endl;
#endif


  return theError;
}


Double_t HypoTestInverterResult::LowerLimitEstimatedError()
{
   // need to have compute first lower limit
   if (TMath::IsNaN(fLowerLimit) ) LowerLimit();
   if (fLowerLimitError >= 0) return fLowerLimitError;    
   // try to recompute the error
   return CalculateEstimatedError(1-ConfidenceLevel(), true);
}


Double_t HypoTestInverterResult::UpperLimitEstimatedError()
{
   if (TMath::IsNaN(fUpperLimit) ) UpperLimit();
   if (fUpperLimitError >= 0) return fUpperLimitError;
   // try to recompute the error
   return CalculateEstimatedError(1-ConfidenceLevel(), false);
}

SamplingDistribution *  HypoTestInverterResult::GetBackgroundTestStatDist(int index ) const { 
   // get the background test statistic distribution

   HypoTestResult * firstResult = (HypoTestResult*) fYObjects.At(index);
   if (!firstResult) return 0;
   return firstResult->GetBackGroundIsAlt() ? firstResult->GetAltDistribution() : firstResult->GetNullDistribution();
}

SamplingDistribution *  HypoTestInverterResult::GetSignalAndBackgroundTestStatDist(int index) const { 
   // get the signal and background test statistic distribution
   HypoTestResult * result = (HypoTestResult*) fYObjects.At(index);
   if (!result) return 0;
   return !result->GetBackGroundIsAlt() ? result->GetAltDistribution() : result->GetNullDistribution();       
}

SamplingDistribution *  HypoTestInverterResult::GetExpectedPValueDist(int index) const { 
   // get the expected p-value distribution at the scanned point index

   if (index < 0 || index >=  ArraySize() ) return 0; 

   if (fExpPValues.GetSize() == ArraySize()  ) {
      return (SamplingDistribution*)  fExpPValues.At(index)->Clone();
   }
   
   static bool useFirstB = false;
   // get the S+B distribution 
   int bIndex = (useFirstB) ? 0 : index;
 
   SamplingDistribution * bDistribution = GetBackgroundTestStatDist(bIndex);
   SamplingDistribution * sbDistribution = GetSignalAndBackgroundTestStatDist(index);

   HypoTestResult * result = (HypoTestResult*) fYObjects.At(index);

   if (bDistribution && sbDistribution) { 

      // create a new HypoTestResult
      HypoTestResult tempResult; 
      tempResult.SetPValueIsRightTail( result->GetPValueIsRightTail() );
      tempResult.SetBackgroundAsAlt( true);
      // ownership of SamplingDistribution is in HypoTestResult class now
      tempResult.SetNullDistribution( new SamplingDistribution(*sbDistribution) );
      tempResult.SetAltDistribution( new SamplingDistribution(*bDistribution ) );
      
      std::vector<double> values(bDistribution->GetSize()); 
      for (int i = 0; i < bDistribution->GetSize(); ++i) { 
         tempResult.SetTestStatisticData( bDistribution->GetSamplingDistribution()[i] );
         values[i] = (fUseCLs) ? tempResult.CLs() : tempResult.CLsplusb();
      }
      return new SamplingDistribution("expected values","expected values",values);
   }  
   // in case b abs sbDistribution are null assume is coming from the asymptotic calculator 
   // hard -coded this value (no really needed to be used by user)
   fgAsymptoticMaxSigma = 5; 
   const int npoints = 11;
   const double smax = fgAsymptoticMaxSigma;
   const double dsig = 2* smax/ (npoints-1) ;
   std::vector<double> values(npoints);
   for (int i = 0; i < npoints; ++i) { 
      double nsig = -smax + dsig*i; 
      double pval = AsymptoticCalculator::GetExpectedPValues( result->NullPValue(), result->AlternatePValue(), nsig, fUseCLs, !fIsTwoSided);
      if (pval < 0) { return 0;}
      values[i] = pval;
   }
   return new SamplingDistribution("Asymptotic expected values","Asymptotic expected values",values);
   
}



SamplingDistribution *  HypoTestInverterResult::GetLimitDistribution(bool lower ) const { 
   // get the limit distribution (lower/upper depending on the flag)
   // by interpolating  the expected p values for each point

   
   if (ArraySize()<2) {
      oocoutE(this,Eval) << "HypoTestInverterResult::GetLimitDistribution" 
                         << " not  enought points -  return 0 " << std::endl; 
      return 0; 
   }
   
   ooccoutD(this,Eval) << "HypoTestInverterResult - computing  limit distribution...." << std::endl;

      
   // find optimal size by looking at the PValue distribution obtained
   int npoints = ArraySize();
   std::vector<SamplingDistribution*> distVec( npoints );
   double sum = 0;
   for (unsigned int i = 0; i < distVec.size(); ++i) {
      distVec[i] =  GetExpectedPValueDist(i);
      // sort the distributions
      // hack (by calling InverseCDF(0) will sort the sampling distribution
      if (distVec[i] ) { 
         distVec[i]->InverseCDF(0);
         sum += distVec[i]->GetSize();     
      }
   }
   int  size =  int( sum/ npoints);
  
   if (size < 10) { 
      ooccoutW(this,InputArguments) << "HypoTestInverterResult - set a minimum size of 10 for limit distribution" <<   std::endl;
      size = 10;
   }


  double target = 1-fConfidenceLevel;

  // vector with the quantiles of the p-values for each scanned poi point
  std::vector< std::vector<double>  > quantVec(npoints );
  for (int i = 0; i <  npoints; ++i) {

     if (!distVec[i]) continue;

     // make quantiles from the sampling distributions of the expected p values 
     std::vector<double> pvalues = distVec[i]->GetSamplingDistribution();
     delete distVec[i];  distVec[i] = 0;
     std::sort(pvalues.begin(), pvalues.end());
     // find the quantiles of the distribution 
     double p[1] = {0}; 
     double q[1] = {0};

     quantVec[i] = std::vector<double>(size);
     for (int ibin = 0; ibin < size; ++ibin) { 
        // exclude for a bug in TMath::Quantiles last item 
        p[0] = std::min( (ibin+1) * 1./double(size), 1.0);
        // use the type 1 which give the point value
        TMath::Quantiles(pvalues.size(), 1, &pvalues[0], q, p, true, 0, 1 );
        (quantVec[i])[ibin] = q[0]; 
     } 
  
  }

  // sort the values in x 
  std::vector<unsigned int> index( npoints );
  TMath::SortItr(fXValues.begin(), fXValues.end(), index.begin(), false); 

  // SamplingDistribution * dist0 = distVec[index[0]];
  // SamplingDistribution * dist1 = distVec[index[1]];

  
  std::vector<double> limits(size);
  // loop on the p values and find the limit for each expected point in the quantiles vector  
  for (int j = 0; j < size; ++j ) {

     TGraph g;
     for (int k = 0; k < npoints ; ++k) { 
        if (quantVec[index[k]].size()  > 0 )
           g.SetPoint(g.GetN(), GetXValue(index[k]), (quantVec[index[k]])[j] );
     }

     limits[j] = GetGraphX(g, target, lower);

  }


  if (lower) 
     return new SamplingDistribution("Expected lower Limit","Expected lower limits",limits);
  else
     return new SamplingDistribution("Expected upper Limit","Expected upper limits",limits);
  
}


double  HypoTestInverterResult::GetExpectedLowerLimit(double nsig, const char * opt ) const {
   // Get the expected lower limit  
   // nsig is used to specify which expected value of the UpperLimitDistribution
   // For example 
   // nsig = 0 (default value) returns the expected value 
   // nsig = -1 returns the lower band value at -1 sigma 
   // nsig + 1 return the upper value
   // opt = "" (default) : compute limit by interpolating all the p values, find the corresponding limit distribution
   //                         and then find the quantiles in the limit distribution 
   // ioption = "P" is the method used for plotting. One Finds the corresponding nsig quantile in the p values and then
   // interpolates them
   
   return GetExpectedLimit(nsig, true, opt);
} 

double  HypoTestInverterResult::GetExpectedUpperLimit(double nsig, const char * opt ) const {
   // Get the expected upper limit  
   // nsig is used to specify which expected value of the UpperLimitDistribution
   // For example 
   // nsig = 0 (default value) returns the expected value 
   // nsig = -1 returns the lower band value at -1 sigma 
   // nsig + 1 return the upper value
   // opt is an option specifying the type of method used for computing the upper limit
   // opt = "" (default) : compute limit by interpolating all the p values, find the corresponding limit distribution
   //                         and then find the quantiles in the limit distribution 
   // ioption = "P" is the method used for plotting. One Finds the corresponding nsig quantile in the p values and then
   // interpolates them
   //              
   
   return GetExpectedLimit(nsig, false, opt);
} 


   
double  HypoTestInverterResult::GetExpectedLimit(double nsig, bool lower, const char * opt ) const {
   // get expected limit (lower/upper) depending on the flag 
   // for asymptotic is a special case (the distribution is generated an step in sigma values)
   // distringuish asymptotic looking at the hypotest results
   // if option = "P" get expected limit using directly quantiles of p value distribution 
   // else (default) find expected limit by obtaining first a full limit distributions
   // The last one is in general more correct

   const int nEntries = ArraySize();
   if (nEntries <= 0)  return (lower) ? 1 : 0;  // return 1 for lower, 0 for upper
   
   HypoTestResult * r = dynamic_cast<HypoTestResult *> (fYObjects.First() );
   assert(r != 0);
   if (!r->GetNullDistribution() && !r->GetAltDistribution() ) { 
      // we are in the asymptotic case 
      // get the limits obtained at the different sigma values 
      SamplingDistribution * limitDist = GetLimitDistribution(lower); 
      if (!limitDist) return 0;
      const std::vector<double> & values = limitDist->GetSamplingDistribution();
      if (values.size() <= 1) return 0; 
      double dsig = 2* fgAsymptoticMaxSigma/ (values.size() -1) ;
      int  i = (int) TMath::Floor ( (nsig +  fgAsymptoticMaxSigma)/dsig + 0.5);
      return values[i];
   }

   double p[1] = {0};
   double q[1] = {0};
   p[0] = ROOT::Math::normal_cdf(nsig,1);

   // for CLs+b can get the quantiles of p-value distribution and 
   // interpolate them
   // In the case of CLs (since it is not a real p-value anymore but a ratio) 
   // then it is needed to get first all limit distribution values and then the quantiles 

   TString option(opt);
   option.ToUpper();
   if (option.Contains("P")) {

      TGraph  g;   

      // sort the arrays based on the x values
      std::vector<unsigned int> index(nEntries);
      TMath::SortItr(fXValues.begin(), fXValues.end(), index.begin(), false);

      for (int j=0; j<nEntries; ++j) {
         int i = index[j]; // i is the order index 
         SamplingDistribution * s = GetExpectedPValueDist(i);
         if (!s) { 
            ooccoutI(this,Eval) << "HypoTestInverterResult - cannot compute expected p value distribution for point, x = " 
                                << GetXValue(i)  << " skip it " << std::endl;
            continue;
         } 
         const std::vector<double> & values = s->GetSamplingDistribution();
         double * x = const_cast<double *>(&values[0]); // need to change TMath::Quantiles
         TMath::Quantiles(values.size(), 1, x,q,p,false);
         g.SetPoint(g.GetN(), fXValues[i], q[0] );
         delete s;
      }
      if (g.GetN() < 2) { 
         ooccoutE(this,Eval) << "HypoTestInverterResult - cannot compute limits , not enough points, n =  " << g.GetN() << std::endl;
         return 0;
      }

      // interpolate the graph to obtain the limit
      double target = 1-fConfidenceLevel;
      return GetGraphX(g, target, lower);

   }
   // here need to use the limit distribution 
   SamplingDistribution * limitDist = GetLimitDistribution(lower); 
   if (!limitDist) return 0;
   const std::vector<double> & values = limitDist->GetSamplingDistribution();
   double * x = const_cast<double *>(&values[0]); // need to change TMath::Quantiles
   TMath::Quantiles(values.size(), 1, x,q,p,false);
   return q[0];
   
}


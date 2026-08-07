#include <TMB.hpp>
#include "include/RLhpp.hpp"


// Big changes made by Simon .  Density dependent mortality, legal biomass, area-specific burn-in years

// -------------------------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------
template <class Type>
vector<Type> Growth2(matrix<Type> &Trans, vector<Type> &N,int Nclass){
  vector<Type> Ntemp2(Nclass);
  Ntemp2.setZero();
  for (int Isize=0;Isize<Nclass;Isize++)
   {
    for (int Jsize=0;Jsize<=Isize;Jsize++) Ntemp2(Isize) += Trans(Isize,Jsize)*N(Jsize);
   }
}

template <class Type>
vector<Type> Growth(array<Type> &Trans, array<Type> &N, int Nclass, int Ipnt, int Iarea, int Isex, int Iage, int MaxLen){
  vector<Type> Ntemp2(MaxLen);
  Ntemp2.setZero();
  for (int Isize=0;Isize<Nclass;Isize++)
   {
    for (int Jsize=0;Jsize<=Isize;Jsize++) Ntemp2(Isize) += Trans(Ipnt,Isize,Jsize)*N(Iarea,Isex,Iage,Jsize);
   }
  return(Ntemp2);
}


// -------------------------------------------------------------------------------------------------------------------

template <class Type>
 vector<Type> Hybrid(dataSet<Type> &dat, array<Type> &N, array<Type> &selretwght, array<Type> &selexF, array<Type> &retainF,
           matrix<Type> M, int AreaPass, int YearPass, int StepPass, Type MWhitesPar){

 // Apply the Hybrid method to solve for F by fleet

 int F_tune;
 Type vbio,temp,temp1,join1,Z_adjuster2,Z_adjuster,max_harvest_rate,TotalCatch,ScaleWhiteM;
 Type RetainTemp;

 int MaxLen; MaxLen = dat.MaxLen;
 int Nfleet; Nfleet = dat.Nfleet;
 int Nsex; Nsex = dat.Nsex;
 int Nage; Nage = dat.Nage;
 int BurnIn; BurnIn = dat.BurnIn;
 vector<int> Nlen = dat.Nlen;

 vector<Type> Hrate(Nfleet); Hrate.setZero();                                  // What we are after
 array<Type> Z_rate(Nsex,Nage,MaxLen);
 array<Type> Z_rate2(Nsex,Nage,MaxLen);

 max_harvest_rate = 3.0;
 F_tune = 5;

 // Total catch
 TotalCatch = 0;
 for (int Ifleet=0;Ifleet<Nfleet;Ifleet++)
  if (dat.Area_fleet(AreaPass,Ifleet)==1)
   TotalCatch += dat.Catch(YearPass,StepPass,Ifleet);

 // Get initial Hrate estimate
 for (int Ifleet=0;Ifleet<Nfleet;Ifleet++)
  if (dat.Area_fleet(AreaPass,Ifleet)==1)
   {
    if (dat.Catch(YearPass,StepPass,Ifleet) > 0)
     {
      vbio = 0;
      for (int Isex=0;Isex<Nsex;Isex++)
       for (int Iage=0;Iage<Nage;Iage++)
        for (int Isize=0;Isize<Nlen(Isex);Isize++)
         vbio += N(AreaPass,BurnIn+YearPass,StepPass,Isex,Iage,Isize)*selretwght(Ifleet,Isex,Iage,Isize);
      temp = dat.Catch(YearPass,StepPass,Ifleet)/(vbio + dat.Catch(YearPass,StepPass,Ifleet));
      join1=1.0/(1.0+exp(30.0*(temp-0.95)));
      temp1=join1*temp + (1.0-join1)*0.95;
      Hrate(Ifleet) = -log(1.-temp1);
     }
    else
      Hrate(Ifleet) = 0;
   }

 // Tune
 for (int tune_F=0;tune_F<F_tune;tune_F++)
  {
   // Compute Z given F and M
   for (int Isex=0;Isex<Nsex;Isex++)
    for (int Iage=0;Iage<Nage;Iage++)
     for (int Isize=0;Isize<Nlen(Isex);Isize++)
      {
       if(dat.IsRed(Isex,Iage,AreaPass,StepPass)==0) {ScaleWhiteM = MWhitesPar;} else {ScaleWhiteM = 1.0;}
       Z_rate(Isex,Iage,Isize) = dat.TimeStepLen(YearPass,StepPass)*M(AreaPass,Iage)*ScaleWhiteM;
        for (int Ifleet=0;Ifleet<Nfleet;Ifleet++)
         if (dat.Area_fleet(AreaPass,Ifleet)==1)
          {
           RetainTemp = selexF(Ifleet,Isex,Iage,Isize) * (retainF(Ifleet,Isex,Iage,Isize)+dat.Phi(Ifleet,Iage,YearPass,StepPass)*(1.0-retainF(Ifleet,Isex,Iage,Isize)));
           Z_rate(Isex,Iage,Isize) += Hrate(Ifleet) * RetainTemp;

           // TEMPORARY DEBUG -- remove once the projection Hrate/discard issue is diagnosed.
           // Fleet 6 only (0-based Ifleet==5, the only active fleet in area 6), area 6,
           // years 2032-2035 (0-based YearPass 36-39), final tuning pass, EVERY size bin.
           if (dat.DoProject==1 && AreaPass==5 && Ifleet==5 && tune_F==F_tune-1 &&
               (YearPass==36||YearPass==37||YearPass==38||YearPass==39)) {
             std::cout << "[HybridDebugZ] Year=" << YearPass << " Isize=" << Isize
                       << " N=" << asDouble(N(AreaPass,dat.BurnIn+YearPass,StepPass,Isex,Iage,Isize))
                       << " selexF=" << asDouble(selexF(Ifleet,Isex,Iage,Isize))
                       << " retainF=" << asDouble(retainF(Ifleet,Isex,Iage,Isize))
                       << " RetainTemp=" << asDouble(RetainTemp)
                       << " Z_rate=" << asDouble(Z_rate(Isex,Iage,Isize)) << "\n";
           }
	      }
        Z_rate2(Isex,Iage,Isize) = (1-exp(-Z_rate(Isex,Iage,Isize)))/Z_rate(Isex,Iage,Isize);
       }

   // Now tune
   if (tune_F < F_tune)
    {

     Z_adjuster2 = 0;
     for (int Ifleet=0;Ifleet<Nfleet;Ifleet++)
      if (dat.Catch(YearPass,StepPass,Ifleet) > 0 & dat.Area_fleet(AreaPass,Ifleet)==1)
       {
        for (int Isex=0;Isex<Nsex;Isex++)
         for (int Iage=0;Iage<Nage;Iage++)
          for (int Isize=0;Isize<Nlen(Isex);Isize++)
           Z_adjuster2 += Hrate(Ifleet)*N(AreaPass,BurnIn+YearPass,StepPass,Isex,Iage,Isize)*selretwght(Ifleet,Isex,Iage,Isize)*Z_rate2(Isex,Iage,Isize);
       }
     Z_adjuster = TotalCatch/(Z_adjuster2+0.0001);

     // Adjust total Z
     for (int Isex=0;Isex<Nsex;Isex++)
       for (int Iage=0;Iage<Nage;Iage++)
        for (int Isize=0;Isize<Nlen(Isex);Isize++)
         {
          if(dat.IsRed(Isex,Iage,AreaPass,StepPass)==0) {ScaleWhiteM = MWhitesPar;} else {ScaleWhiteM = 1.0;}
          Z_rate(Isex,Iage,Isize)  = dat.TimeStepLen(YearPass,StepPass)*M(AreaPass,Iage)*ScaleWhiteM + Z_adjuster*(Z_rate(Isex,Iage,Isize)-dat.TimeStepLen(YearPass,StepPass)*M(AreaPass,Iage)*ScaleWhiteM);
          Z_rate2(Isex,Iage,Isize) = (1-exp(-Z_rate(Isex,Iage,Isize)))/Z_rate(Isex,Iage,Isize);

          // TEMPORARY DEBUG -- the FINAL (post-Z_adjuster-rescaling) Z_rate/Z_rate2, i.e. what
          // actually feeds the final Hrate below and becomes the officially reported Z.
          if (dat.DoProject==1 && AreaPass==5 && tune_F==F_tune-1 &&
              (YearPass==36||YearPass==37||YearPass==38||YearPass==39)) {
            std::cout << "[HybridDebugZFinal] Year=" << YearPass << " Isize=" << Isize
                      << " Z_adjuster=" << asDouble(Z_adjuster)
                      << " Z_rate=" << asDouble(Z_rate(Isex,Iage,Isize))
                      << " Z_rate2=" << asDouble(Z_rate2(Isex,Iage,Isize)) << "\n";
          }
         }

     // Adjust total exploitable biomass
     for (int Ifleet=0;Ifleet<Nfleet;Ifleet++)
      if (dat.Catch(YearPass,StepPass,Ifleet) > 0 & dat.Area_fleet(AreaPass,Ifleet)==1)
       {
        Z_adjuster2 = 0;
        for (int Isex=0;Isex<Nsex;Isex++)
         for (int Iage=0;Iage<Nage;Iage++)
          for (int Isize=0;Isize<Nlen(Isex);Isize++)
           Z_adjuster2 += N(AreaPass,BurnIn+YearPass,StepPass,Isex,Iage,Isize)*selretwght(Ifleet,Isex,Iage,Isize)*Z_rate2(Isex,Iage,Isize);
        temp = dat.Catch(YearPass,StepPass,Ifleet)/(Z_adjuster2 + 0.00001);
        join1=1.0/(1.0+exp(30.0*(temp-0.95*max_harvest_rate)));
        Hrate(Ifleet)=join1*temp + (1.0-join1)*max_harvest_rate;
//        if (tune_F = F_tune-1 & YearPass < 3) std::cout << "test " << Ifleet << " " << YearPass << " " << StepPass << " " << dat.Catch(YearPass,StepPass,Ifleet) << " " << Z_adjuster2 << " " << Hrate(Ifleet) << "\n";

        // TEMPORARY DEBUG -- remove once the projection Hrate/discard issue is diagnosed.
        if (dat.DoProject==1 && AreaPass==5 && tune_F==F_tune-1 &&
            (YearPass==36||YearPass==37||YearPass==38||YearPass==39)) {
          std::cout << "[HybridDebugHrate] Year=" << YearPass << " Fleet=" << Ifleet
                    << " Catch=" << asDouble(dat.Catch(YearPass,StepPass,Ifleet))
                    << " Z_adjuster2=" << asDouble(Z_adjuster2)
                    << " Hrate=" << asDouble(Hrate(Ifleet)) << "\n";
        }
       }
    }
  } // Tune
  return(Hrate);
}

// -------------------------------------------------------------------------------------------------------------------

template <class Type>
matrix<Type> SetUpSelex(dataSet<Type> &dat,  vector<Type> &SelPars, matrix<Type> &SelexFI, matrix<int> &PatSpec, int Npatterns ){

  int IPreSpecified, IselParPnt, Isex;
  Type p1,p2,p3,p4;

  int MaxLen; MaxLen = dat.MaxLen;
  matrix<Type> ActSelex(Npatterns,MaxLen);
  Type Mult, LML, Offset, MaxTmp;

  IselParPnt = -1;
  for (int IselPattern=0;IselPattern<Npatterns;IselPattern++)
  {

    // Pre-specified
    if (PatSpec(IselPattern,1) == SELEX_PRESPECIFIED)
    {
      IPreSpecified = PatSpec(IselPattern,4);
      for (int Isize=0;Isize<MaxLen;Isize++) ActSelex(IselPattern,Isize) = SelexFI(IPreSpecified,Isize);
    }

    // Estimated up to a constant
   /* if (PatSpec(IselPattern,1) == SELEX_COEFFICIENTS)
    {
      for (int Isize=0;Isize<PatSpec(IselPattern,3);Isize++) { IselParPnt += 1; ActSelex(IselPattern,Isize) = exp(SelPars(IselParPnt)); }
      for (int Isize=PatSpec(IselPattern,3);Isize<MaxLen; Isize++) ActSelex(IselPattern,Isize) = 1;
    }
    */

    // Logistic
    if (PatSpec(IselPattern,1) == SELEX_LOGISTIC) // SELEX_LOGISTIC
    {
      p1 = SelPars(IselParPnt+1); p2 = SelPars(IselParPnt+2);
      IselParPnt += 2;
      Isex = PatSpec(IselPattern,2);
      for (int Isize=0;Isize<MaxLen; Isize++) ActSelex(IselPattern,Isize) = 1.0/(1.0+exp(-p2*(dat.MidLenBin(Isex,Isize)-p1)));
    }

      // Double Logistic
    if (PatSpec(IselPattern,1) == SELEX_DOUBLE_LOGISTIC)
    {
      p1 = SelPars(IselParPnt+1); p2 = SelPars(IselParPnt+2); p3 = SelPars(IselParPnt+3); p4 = SelPars(IselParPnt+4);
      IselParPnt += 4;
      Isex = PatSpec(IselPattern,2);
      MaxTmp = 0.0;
      for (int Isize=0;Isize<MaxLen; Isize++) {
        ActSelex(IselPattern,Isize) =  (1.0/(1.0+exp(-p2*(dat.MidLenBin(Isex,Isize)-p1)))) * (1.0/(1.0+exp(-p4*(dat.MidLenBin(Isex,Isize)-p3))));
        if(ActSelex(IselPattern,Isize)>MaxTmp) MaxTmp = ActSelex(IselPattern,Isize); }
       for (int Isize=0;Isize<MaxLen; Isize++) { // re-scale to a max of 1
        ActSelex(IselPattern,Isize) = ActSelex(IselPattern,Isize)/MaxTmp; }
    }

    /*
    // Knife-edged
    if (PatSpec(IselPattern,1) == SELEX_KNIFE)
    {
      LML = SelPars(IselParPnt+1);
      IselParPnt += 1;
      Isex = PatSpec(IselPattern,2);
      for (int Isize=0;Isize<MaxLen; Isize++)
      {
        if (dat.LowLenBin(Isex,Isize+1) <= LML)
          Mult = 0;
        else
          if (dat.LowLenBin(Isex,Isize) >= LML)
            Mult = 1;
          else
            Mult = (dat.LowLenBin(Isex,Isize+1)-LML) / (dat.LowLenBin(Isex,Isize+1)-dat.LowLenBin(Isex,Isize));
          ActSelex(IselPattern,Isize) = Mult;
      }
    }
    // flat
    if (PatSpec(IselPattern,1) == SELEX_CONSTANT1)
    {
      for (int Isize=0;Isize<MaxLen; Isize++) ActSelex(IselPattern,Isize) = 1.0;
    }
    // Logistic
    if (PatSpec(IselPattern,1) == SELEX_LOGISTIC_OFFSET)
    {
      p1 = SelPars(IselParPnt+1); p2 = SelPars(IselParPnt+2); Offset = SelPars(IselParPnt+3);
      IselParPnt += 3;
      Isex = PatSpec(IselPattern,2);
      for (int Isize=0;Isize<MaxLen; Isize++) ActSelex(IselPattern,Isize) = Offset/(1.0+exp(-p2*(dat.MidLenBin(Isex,Isize)-p1)));
    }

    // Knife-edged
    if (PatSpec(IselPattern,1) == SELEX_KNIFE_OFFSET)
    {
      LML = SelPars(IselParPnt+1);
      Offset = SelPars(IselParPnt+2);
      IselParPnt += 2;
      Isex = PatSpec(IselPattern,2);
      for (int Isize=0;Isize<MaxLen; Isize++)
      {
        if (dat.LowLenBin(Isex,Isize+1) <= LML)
          Mult = 0;
        else
          if (dat.LowLenBin(Isex,Isize) >= LML)
            Mult = 1;
          else
            Mult = (dat.LowLenBin(Isex,Isize+1)-LML) / (dat.LowLenBin(Isex,Isize+1)-dat.LowLenBin(Isex,Isize));
          ActSelex(IselPattern,Isize) = Mult*Offset;
      }
    }
    // flat
    if (PatSpec(IselPattern,1) == SELEX_CONSTANT_OFFSET)
    {
      Offset = SelPars(IselParPnt+1);
      IselParPnt += 1;
      for (int Isize=0;Isize<MaxLen; Isize++) ActSelex(IselPattern,Isize) = Offset;
    }
     */
  }

  return(ActSelex);
}
// -------------------------------------------------------------------------------------------------------------------
template <class Type>
matrix<Type> SetUpLegal(dataSet<Type> &dat, matrix<Type> &LegalFI, matrix<int> &PatSpec, int Npatterns ){

  int IPreSpecified;

  int MaxLen; MaxLen = dat.MaxLen;
  matrix<Type> ActLegal(Npatterns,MaxLen);

  for (int IlegalPattern=0;IlegalPattern<Npatterns;IlegalPattern++)
  {
    // Pre-specified
    if (PatSpec(IlegalPattern,1) == SELEX_PRESPECIFIED)
     {
      IPreSpecified = PatSpec(IlegalPattern,3);
      for (int Isize=0;Isize<MaxLen;Isize++) ActLegal(IlegalPattern,Isize) = LegalFI(IPreSpecified,Isize);
     }
  }

  return(ActLegal);
}
// -------------------------------------------------------------------------------------------------------------------


template <class Type>
 matrix<Type> SetUpMove(dataSet<Type> &dat,  vector<Type> &MovePars ){

 // This function sets up all the movement patterns

 int MaxLen; MaxLen = dat.MaxLen;
 matrix<Type> ActMove(dat.NmovePatterns,MaxLen);
 int ImoveParPnt,Isex;
 Type rate,ChangePnt,Mult;

 ActMove.setZero();
 ImoveParPnt = -1;
 for (int ImovePattern=0;ImovePattern<dat.NmovePatterns;ImovePattern++)
  {
   // Constant
   if (dat.MoveSpec(ImovePattern,1) == MOVE_CONSTANT)
    {
     rate = MovePars(ImoveParPnt+1);
     ImoveParPnt += 1;
     for (int Isize=0;Isize<dat.MaxLen;Isize++) ActMove(ImovePattern,Isize) = rate;
    }
   // Knife-edged (uses lengths for sex=1)
   if (dat.MoveSpec(ImovePattern,1) == MOVE_KNIFE)
    {
     ChangePnt = MovePars(ImoveParPnt+1);
     rate = MovePars(ImoveParPnt+2);
     ImoveParPnt += 2;
     Isex = 0;
     for (int Isize=0;Isize<dat.MaxLen; Isize++)
      {
       if (dat.LowLenBin(Isex,Isize+1) <= ChangePnt)
        Mult = 0;
       else
        if (dat.LowLenBin(Isex,Isize) >= ChangePnt)
         Mult = 1;
        else
         Mult = (dat.LowLenBin(Isex,Isize+1)-ChangePnt) / (dat.LowLenBin(Isex,Isize+1)-dat.LowLenBin(Isex,Isize));
       ActMove(ImovePattern,Isize) = rate*Mult;
      }
    }
  }
 return(ActMove);
}

// -------------------------------------------------------------------------------------------------------------------

template <class Type>
 Type SetUpRecruit(dataSet<Type> &dat,  vector<Type> &RecruitPars, vector<Type> &RecSpatDevPars, array <Type> &ActRecruitAreaSexDist, array <Type> &ActRecruitLenDist ){

 // This function sets up all the movement patterns
 Type XX;
 int IrecruitParPnt,IrecSpatPnt,Pointer,Icnt;
 vector<Type> SexSplit(2);
 vector<Type> TempArea(dat.Narea*dat.Nsex);
 Type Total,AreaPar;

 int MaxLen; MaxLen = dat.MaxLen;
 //matrix<Type> ActRecruit(dat.NrecruitPatterns,MaxLen);

// This function sets up all the movement patterns
 ActRecruitAreaSexDist.setZero();
 ActRecruitLenDist.setZero();
 IrecruitParPnt = -1;
 IrecSpatPnt = -1;

 // Extract the AreaSexDist
 for (int IrecruitPattern=0;IrecruitPattern<dat.NrecruitPatternsA;IrecruitPattern++)
  {
   // allow for common sex-ratios at recruitment
   if (dat.RecruitSpecsA(IrecruitPattern,1)==0)
    {
     IrecruitParPnt += 1;
     SexSplit(0) = 1.0/(1.0+exp(RecruitPars(IrecruitParPnt)));
     SexSplit(1) = 1.0-SexSplit(0);
     if(dat.Nsex==1) { SexSplit(0) = 1.0; SexSplit(1) = 0.0; }

     for (int Iarea=0;Iarea<dat.Narea;Iarea++)
      {
       if (Iarea==0)
        { TempArea(Iarea) = 0; }
       else
        {
         IrecruitParPnt += 1;
         TempArea(Iarea) = RecruitPars(IrecruitParPnt);
        }
       }
      for (int Year=0;Year<dat.Nyear+dat.MaxProjYr;Year++)
       for (int Istep=0;Istep<dat.Nstep;Istep++)
        if (dat.RecruitPnt(Year,Istep)==IrecruitPattern)
         {
          // Insert
          Total = 0;
          for (int Iarea=0;Iarea<dat.Narea;Iarea++)
           {
            if (Iarea==0)
             { AreaPar = 1; }
            else
             {
              if (Year>=dat.RecSpatYr1 & Year<=dat.RecSpatYr2)
               {
                IrecSpatPnt += 1;
                AreaPar = exp(TempArea(Iarea)+RecSpatDevPars(IrecSpatPnt));
               }
              else
               AreaPar = exp(TempArea(Iarea));
             }
            for (int Isex=0;Isex<dat.Nsex;Isex++)
             {
              ActRecruitAreaSexDist(Year,Istep,Iarea,Isex) = AreaPar*SexSplit(Isex);
              Total += ActRecruitAreaSexDist(Year,Istep,Iarea,Isex);
             }
           }
         // normalize
         for (int Iarea=0;Iarea<dat.Narea;Iarea++)
          for (int Isex=0;Isex<dat.Nsex;Isex++)
           ActRecruitAreaSexDist(Year,Istep,Iarea,Isex) /= Total;
        }

     }

   // allow for area-specific sex-ratios at recruitment
   if (dat.RecruitSpecsA(IrecruitPattern,1)==1)
    {
     Total = 0;
     Icnt = 0;
     for (int Iarea=0;Iarea<dat.Narea;Iarea++)
      for (int Isex=0;Isex<dat.Nsex;Isex++)
       {
        if (Iarea==0 & Isex==0)
         {
          TempArea(0) = 0;
         }
        else
          {
          Icnt += 1;
          IrecruitParPnt += 1;
          TempArea(Icnt) = RecruitPars(IrecruitParPnt);
         }
       }
      for (int Year=0;Year<dat.Nyear+dat.MaxProjYr;Year++)
       for (int Istep=0;Istep<dat.Nstep;Istep++)
        if (dat.RecruitPnt(Year,Istep)==IrecruitPattern)
         {
          // Insert
          Total = 0;
          Icnt = 0;
          for (int Iarea=0;Iarea<dat.Narea;Iarea++)
           for (int Isex=0;Isex<dat.Nsex;Isex++)
           {
            if (Iarea==0 & Isex==0)
             { AreaPar = 1; }
            else
             {
              Icnt += 1;
              if (Year>=dat.RecSpatYr1 & Year<=dat.RecSpatYr2)
               {
                if (Isex==0 || Iarea == 1) IrecSpatPnt += 1;
                AreaPar = exp(TempArea(Icnt)+RecSpatDevPars(IrecSpatPnt));
               }
              else
               AreaPar = exp(TempArea(Icnt));
             }
            ActRecruitAreaSexDist(Year,Istep,Iarea,Isex) = AreaPar*SexSplit(Isex);
            Total += ActRecruitAreaSexDist(Year,Istep,Iarea,Isex);
            }
          // normalize
          for (int Iarea=0;Iarea<dat.Narea;Iarea++)
           for (int Isex=0;Isex<dat.Nsex;Isex++)
            ActRecruitAreaSexDist(Year,Istep,Iarea,Isex) /= Total;
        }

    }
  }

 // Now deal with recruitment distribution
 for (int IrecruitPattern=0;IrecruitPattern<dat.NrecruitPatternsB;IrecruitPattern++)
  {
   for (int Isex=0;Isex<dat.Nsex;Isex++)
    {
     if (dat.RecruitSpecsB(IrecruitPattern,1+Isex) < 0)
      {
       Pointer = -1*dat.RecruitSpecsB(IrecruitPattern,1+Isex)-1;
       for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
        ActRecruitLenDist(IrecruitPattern,Isex,Isize) = dat.RecruitFrac(Pointer,Isize);
      }
     else
      {
       Total = 0;
       for (int Isize=0;Isize<dat.RecruitSpecsB(IrecruitPattern,2+dat.Nsex+Isex) ;Isize++)
        {
         IrecruitParPnt += 1;
         ActRecruitLenDist(IrecruitPattern,Isex,Isize) = exp(RecruitPars(IrecruitParPnt));
         Total += ActRecruitLenDist(IrecruitPattern,Isex,Isize);
        }
       for (int Isize=0;Isize<dat.RecruitSpecsB(IrecruitPattern,2+dat.Nsex+Isex) ;Isize++)
        {
         ActRecruitLenDist(IrecruitPattern,Isex,Isize) /= Total;
        }
      }
    }

  }

 XX = 1;
 return(XX);
}
// -------------------------------------------------------------------------------------------------------------------

template <class Type>
 matrix<Type> SetUpGrow(dataSet<Type> &dat,  vector<Type> &GrowthPars ){
 matrix<Type> TempGrow(dat.MaxLen,dat.MaxLen);

 int MaxLen; MaxLen = dat.MaxLen;
 array<Type> ActGrowth(dat.NgrowthPatterns,MaxLen,MaxLen);
 int IgrowthParPnt,Isex,Pointer;

 // This function sets up all the growth patterns
 ActGrowth.setZero();
 IgrowthParPnt = -1;
 for (int IgrowthPattern=0;IgrowthPattern<dat.NgrowthPatterns;IgrowthPattern++)
  {
   if (dat.GrowthSpecs(IgrowthPattern,1) == GROWTH_PRESPECIFIED)
    {
     Isex = dat.GrowthSpecs(IgrowthPattern,2);
     Pointer = dat.GrowthSpecs(IgrowthPattern,4);
     for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
      for (int Jsize=0;Jsize<dat.Nlen(Isex);Jsize++)
       ActGrowth(IgrowthPattern,Isize,Jsize) = dat.TransInp(Pointer,Isize,Jsize);
      if (dat.GrowthSpecs(IgrowthPattern,5) > 1)
	   {
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
         for (int Jsize=0;Jsize<dat.Nlen(Isex);Jsize++)
          ActGrowth(IgrowthPattern,Isize,Jsize) = dat.TransInp(Pointer,Isize,Jsize);
	    for (int Imult=2;Imult<=dat.GrowthSpecs(IgrowthPattern,5);Imult++)
	     {
          TempGrow.setZero();
          for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
           for (int Jsize=0;Jsize<dat.Nlen(Isex);Jsize++)
            {
			 for (int Ksize=0;Ksize<dat.Nlen(Isex);Ksize++)
 			  TempGrow(Isize,Jsize) += ActGrowth(IgrowthPattern,Isize,Ksize)*dat.TransInp(Pointer,Ksize,Jsize);
		    }
          for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
           for (int Jsize=0;Jsize<dat.Nlen(Isex);Jsize++)
            ActGrowth(IgrowthPattern,Isize,Jsize) = TempGrow(Isize,Jsize);
	     }
       }
    }
  }

 return(ActGrowth);

}
//==================================================================================================================================

template <class Type>
 vector<Type> OneTimeStep(dataSet<Type> &dat, array<Type> &N, array<Type> &Z, array<Type> &Hrate,
                         matrix<Type> &ActSelex, matrix<Type> &ActReten, matrix<Type> &ActLegal, matrix<Type> &ActMove,
                         matrix<Type> &WeightLen, matrix<Type> M,
                         int Iyear, int Istep, array<Type> &ActGrowth, matrix<Type> &RecruitFrac, Type Rbar,
                         int IsVirgin, matrix<Type> &Feqn2, array <Type> &ActRecruitAreaSexDist, array <Type> &ActRecruitLenDist,
                         vector<Type> &ActRecDev, vector<Type> &MatBio, matrix<Type> &MatBioArea, matrix<Type> &RecruitmentByArea,
                         vector<Type>BiasMult, Type SigmaR, Type QRedsPar, Type MWhitesPar,
                         vector<Type> &VirginBio, vector<Type> &CurrentBio){

 array<Type> Z_rate(dat.Nsex,dat.Nage,dat.MaxLen);                             // total mortality

 array<Type> selexF(dat.Nfleet,dat.Nsex,dat.Nage,dat.MaxLen);                  // Selectivity
 array<Type> retainF(dat.Nfleet,dat.Nsex,dat.Nage,dat.MaxLen);                 // Retention
 array<Type> selretwght(dat.Nfleet,dat.Nsex,dat.Nage,dat.MaxLen);              // Product of selectivity,retention and weight
 array<Type> Ntemp(dat.Narea,dat.Nsex,dat.Nage,dat.MaxLen);                    // N matrix (after mortality)
 array<Type> Nmove(dat.Narea,dat.Nsex,dat.Nage,dat.MaxLen);                    // N matrix (after movement)
 vector<Type> Ntemp2(dat.MaxLen);                                              // Matrix multiplication
 vector<Type> MoveVec(dat.MaxLen);                                             // Matrix multiplication
 vector<Type> HratePass(dat.Nfleet);                                           // Pass of harvest rare
 Type RetainTemp,TotalRec, ScaleWhiteM, ScaleRedQ;

 int SelPointer,RetPointer,LegalPointer,MovePointer,RecruitPointer,GrowthPointer;           // Pointers
 int YearAdjust1,YearAdjust2,IsMoves,IdestArea,RecruitLenPointer;

 vector<Type> XX(2);
 XX(0) = 1; XX(1) = 1;

 int Ipnt;                                                                     // Pointer

 // Adjusted year (YearAdjust1 is for quantities that go beyond Nyear-1 and YearAdjust2 is not.
 if (Iyear <= 0)
  { YearAdjust1 = 0; YearAdjust2 = 0; }
 else
  if (Iyear < dat.Nyear)
   { YearAdjust1 = Iyear; YearAdjust2 = Iyear;}
  else
   { YearAdjust1 = Iyear; YearAdjust2 = dat.Nyear-1;}

 // Recruitment (at the start of the time-step)
 RecruitPointer = dat.RecruitPnt(YearAdjust1,Istep);
 if (RecruitPointer >= 0)
 {
   for (int Iarea=0;Iarea<dat.Narea;Iarea++)
   {
     RecruitLenPointer = dat.RecruitLenPnt(Iarea);
     for (int Isex=0;Isex<dat.Nsex;Isex++)
     {
       TotalRec = ActRecruitAreaSexDist(YearAdjust1,Istep,Iarea,Isex)*exp(Rbar)*exp(ActRecDev(dat.BurnIn+Iyear))*exp(-BiasMult(dat.BurnIn+Iyear)*SigmaR*SigmaR/2.0);
       for (int Isize=0;Isize<dat.Nlen(Isex);Isize++) N(Iarea,dat.BurnIn+Iyear,Istep,Isex,0,Isize) += ActRecruitLenDist(RecruitLenPointer,Isex,Isize)*TotalRec;
     }
   }
 }

 // Current Biomass
   CurrentBio.setZero();
   for (int Iarea=0;Iarea<dat.Narea;Iarea++){
     for (int Isex=0;Isex<dat.Nsex;Isex++){
       for (int Iage=0;Iage<dat.Nage;Iage++){
         for (int Isize=0;Isize<dat.Nlen(Isex);Isize++) {
           CurrentBio(Iarea) += N(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Isize) * WeightLen(Isex,Isize); }}}}  // Weight in each area, first year/time-step

 // Maturity and fecundity
 if (Istep==dat.MatTimeStep)  {
   MatBio(dat.BurnIn+Iyear) = 0;
   for (int Iarea=0;Iarea<dat.Narea;Iarea++) {
     MatBioArea(Iarea,dat.BurnIn+Iyear) = 0;
     for (int Iage=0;Iage<dat.Nage;Iage++){
      if(Iage>=dat.MatAge(Iarea)){
        for (int Isize=0;Isize<dat.Nlen(0);Isize++){
         MatBioArea(Iarea,dat.BurnIn+Iyear) += N(Iarea,dat.BurnIn+Iyear,Istep,0,Iage,Isize)*dat.MatFem(Iage, Iarea, dat.BurnIn+Iyear, Isize);}}
      }
     MatBio(dat.BurnIn+Iyear) += MatBioArea(Iarea,dat.BurnIn+Iyear);
     }
   }

 // Need to set selectivity
 for (int Ifleet=0;Ifleet<dat.Nfleet;Ifleet++) {
  for (int Isex=0;Isex<dat.Nsex;Isex++) {
   for (int Iage=0;Iage<dat.Nage;Iage++) {
     if(dat.IsRed(Isex,Iage,dat.Fleet_area(Ifleet),Istep)==1) {ScaleRedQ = QRedsPar;} else {ScaleRedQ = 1.0;}
     if (Iyear<dat.Nyear)
      {
       SelPointer = dat.SelPnt(Isex,Iage,Ifleet,YearAdjust1,Istep);
       RetPointer = dat.RetPnt(Isex,Iage,Ifleet,YearAdjust1,Istep);
       LegalPointer = dat.LegalFleetPnt(Isex,Iage,Ifleet,YearAdjust1,Istep);
      }
     else
      {
       SelPointer = dat.SelPntFut(Isex,Iage,Ifleet,YearAdjust1-dat.Nyear,Istep);
       RetPointer = dat.RetPntFut(Isex,Iage,Ifleet,YearAdjust1-dat.Nyear,Istep);
       LegalPointer = dat.LegalFleetPntFut(Isex,Iage,Ifleet,YearAdjust1-dat.Nyear,Istep);

       // TEMPORARY DEBUG -- remove once the projection Hrate/discard issue is diagnosed.
       if (dat.DoProject==1 && Ifleet==5 && (Iyear==36||Iyear==37||Iyear==38||Iyear==39)) {
         std::cout << "[OneTimeStepDebugPtr] Iyear=" << Iyear << " YearAdjust1=" << YearAdjust1
                   << " FutIndex=" << (YearAdjust1-dat.Nyear)
                   << " SelPointer=" << SelPointer << " RetPointer=" << RetPointer
                   << " LegalPointer=" << LegalPointer << "\n";
       }
	  }
     for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++) {
       selexF(Ifleet,Isex,Iage,Ilen) = ActSelex(SelPointer,Ilen) * ScaleRedQ;
       retainF(Ifleet,Isex,Iage,Ilen) = ActReten(RetPointer,Ilen)*ActLegal(LegalPointer,Ilen);
       selretwght(Ifleet,Isex,Iage,Ilen) = selexF(Ifleet,Isex,Iage,Ilen)*retainF(Ifleet,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);
     }
    }
   }
  }


 Ntemp.setZero();
 HratePass.setZero();
 for (int Iarea=0;Iarea<dat.Narea;Iarea++)
  {

   // Find the F for this time-step
   if (IsVirgin==0)
    {
     if (Iyear < 0)
      {
       for (int Ifleet=0; Ifleet<dat.Nfleet;Ifleet++)
        if (dat.Area_fleet(Iarea,Ifleet)==1) Hrate(dat.BurnIn+Iyear,Istep,Ifleet) = Feqn2(Ifleet,Istep);
      }
     else if (Iyear >= dat.Nyear && dat.ProjType == 2)
      {
       // Harvest-rate (effort) based projection: Hrate is imposed directly from
       // PROJECTIONS.DAT rather than solved from a target catch via Hybrid().
       for (int Ifleet=0; Ifleet<dat.Nfleet;Ifleet++)
        if (dat.Area_fleet(Iarea,Ifleet)==1)
         Hrate(dat.BurnIn+Iyear,Istep,Ifleet) = dat.ProjHarvestRate(Iyear-dat.Nyear,Istep,Ifleet);
      }
     else
      {
       HratePass = Hybrid(dat, N, selretwght, selexF, retainF, M, Iarea, Iyear, Istep, MWhitesPar);
       for (int Ifleet=0;Ifleet<dat.Nfleet;Ifleet++)
        if(dat.Area_fleet(Iarea,Ifleet)==1) Hrate(dat.BurnIn+Iyear,Istep,Ifleet) = HratePass(Ifleet);
      }
    }

   // Compute Z given F and M
   for (int Isex=0;Isex<dat.Nsex;Isex++)
    for (int Iage=0;Iage<dat.Nage;Iage++)
      {
	    if(dat.IsRed(Isex,Iage,Iarea,Istep)==0) {ScaleWhiteM = MWhitesPar;} else {ScaleWhiteM = 1.0;}
      for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
       {
	     if (Iyear <= 0){
         Z_rate(Isex,Iage,Isize) = dat.TimeStepLen(0,Istep)*M(Iarea,Iage)*ScaleWhiteM;}
       else{
         Z_rate(Isex,Iage,Isize) = dat.TimeStepLen(Iyear,Istep)*M(Iarea,Iage)*ScaleWhiteM;}
       for (int Ifleet=0;Ifleet<dat.Nfleet;Ifleet++)
        if(dat.Area_fleet(Iarea,Ifleet)==1)
         {

          RetainTemp = selexF(Ifleet,Isex,Iage,Isize) * (retainF(Ifleet,Isex,Iage,Isize)+dat.Phi(Ifleet,Iage,YearAdjust1,Istep)*(1.0-retainF(Ifleet,Isex,Iage,Isize)));
          Z_rate(Isex,Iage,Isize) += Hrate(dat.BurnIn+Iyear,Istep,Ifleet) * RetainTemp;
	     }
       Z(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Isize) = Z_rate(Isex,Iage,Isize);
      }}


   // Remove mortality
   for (int Isex=0;Isex<dat.Nsex;Isex++)
    for (int Iage=0;Iage<dat.Nage;Iage++)
     for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
      Ntemp(Iarea,Isex,Iage,Isize) = N(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Isize) * exp(-Z(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Isize));

   // growth
   for (int Isex=0;Isex<dat.Nsex;Isex++)
    for (int Iage=0;Iage<dat.Nage;Iage++)
     {
      GrowthPointer = dat.GrowthPnt(Iarea,Isex,Iage,YearAdjust2,Istep);
      if (GrowthPointer >=0)
       {
        // Key issue (pointer to growth matrix)
        Ntemp2.setZero();
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
         {
          for (int Jsize=0;Jsize<=Isize;Jsize++) Ntemp2(Isize) += Ntemp(Iarea,Isex,Iage,Jsize)*ActGrowth(GrowthPointer,Isize,Jsize);
         }
        //Ntemp2 = Growth(ActGrowth,Ntemp,dat.Nlen(Isex),Ipnt,Iarea,Isex,Iage,dat.MaxLen);
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++) Ntemp(Iarea,Isex,Iage,Isize) = Ntemp2(Isize);
       }
     } // growth

  } // area

 Nmove.setZero();
 IsMoves = 0;
 for (int Iarea=0;Iarea<dat.Narea;Iarea++)
   {
  for (int Iage=0;Iage<dat.Nage;Iage++)
   {
    MovePointer = dat.MovePnt(Iarea,Iage,YearAdjust2,Istep);
    if (MovePointer > 0)
     {
      IsMoves = 1;
      IdestArea = dat.MoveSpec(MovePointer,2);
      for (int Isex=0;Isex<dat.Nsex;Isex++)
       {
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++) MoveVec(Isize) = ActMove(MovePointer,Isize);
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
         {
          Nmove(IdestArea,Isex,Iage,Isize) += MoveVec(Isize)*Ntemp(Iarea,Isex,Iage,Isize);
          Nmove(Iarea,Isex,Iage,Isize) -= MoveVec(Isize)*Ntemp(Iarea,Isex,Iage,Isize);
         }
        // redistribute lobster so they can potentially move more than one area in a timestep (only in increasing area number)
       // for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++){                   // delete to change back
       //   Ntemp(Iarea,Isex,Iage,Ilen) += Nmove(Iarea,Isex,Iage,Ilen);} // delete to change back
             } // Sex
     } // If there was a move
   } // All areas and ages
 }

 // Only update if needed
 //
   if (IsMoves==1)
  {
  for (int Iarea=0;Iarea<dat.Narea;Iarea++){
     for (int Isex=0;Isex<dat.Nsex;Isex++){
      for (int Iage=0;Iage<dat.Nage;Iage++){
       for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++){
        Ntemp(Iarea,Isex,Iage,Ilen) += Nmove(Iarea,Isex,Iage,Ilen);}}}}
   }

 // Update seasons
 for (int Iarea=0;Iarea<dat.Narea;Iarea++)
  for (int Isex=0;Isex<dat.Nsex;Isex++)
   {
    if (Istep<dat.Nstep-1)
     {
      for (int Iage=0;Iage<dat.Nage;Iage++){
       for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++){
        N(Iarea,dat.BurnIn+Iyear,Istep+1,Isex,Iage,Ilen) = Ntemp(Iarea,Isex,Iage,Ilen);}}
     }
    else
     {
      // special case
      if (dat.Nage-1 > 0)
       {
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++) N(Iarea,dat.BurnIn+Iyear+1,0,Isex,0,Isize) = 0;
        for (int Iage=0;Iage<dat.Nage-1;Iage++)
         for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
          N(Iarea,dat.BurnIn+Iyear+1,0,Isex,Iage+1,Isize) = Ntemp(Iarea,Isex,Iage,Isize);
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
         N(Iarea,dat.BurnIn+Iyear+1,0,Isex,dat.Nage-1,Isize) =  Ntemp(Iarea,Isex,dat.Nage-1,Isize) + Ntemp(Iarea,Isex,dat.Nage-2,Isize);
       }
      else
       {
        for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++){
         N(Iarea,dat.BurnIn+Iyear+1,0,Isex,0,Ilen) = Ntemp(Iarea,Isex,0,Ilen);}
       }
      } // if
   } // sex

 return(XX);
}


// -------------------------------------------------------------------------------------------------------------------

template <class Type>
 Type CpueLikelihood(dataSet<Type> &dat, TheData<Type> &thedata, array<Type> &N, array<Type> &Z,
                         matrix<Type> &ActSelex, matrix<Type> &ActReten, matrix<Type> &ActLegal,
                         matrix<Type> &WeightLen, matrix<Type> &PredCpue,
                         vector<Type> &CpueLikeComps, vector<Type> &SigmaCpue, vector<Type> &Qval, matrix<Type> &CpueEcreep, vector<Type> &Qpars, vector<Type> &efpars,matrix<Type> M, Type QRedsPar) {

 Type selexFU,retainFU,selretwght,SigmaUse,ScaleRedQ;
 int Ifleet,Jsex,Iyear,Istep,IdataSet,IndexPoint,IndexPoint2;
 int PntCnt,Iarea,Isex1,Isex2;
 int SelPointer,RetPointer,LegalPointer,IenvPnt;
 vector <Type> SS(thedata.NcpueDataSeries);
 vector <Type> Ndata(thedata.NcpueDataSeries);

 Type NeglogLikelihood = 0;

 // Find predicted biomass
 PredCpue.setZero();
 Qval.setZero();
 Ndata.setZero();
 CpueEcreep.setZero();

 // Make efficiency creep matrix from parameters with time lags
 int Lenefseries = CpueEcreep.rows();  // N years
 int Nefseries = CpueEcreep.cols();    // Number of unique Lags times
 int efcnt = -1; int parcnt = -1;
 Type Tmppar;  // store temporary parameter
 for (int Nef=0;Nef<Nefseries;Nef++){
   CpueEcreep(0,Nef) = 1.0;  // Set first year to 1 (no efficiency creep)
   parcnt = parcnt + 1;
   efcnt = 0;
   for (int Yef=1;Yef<Lenefseries;Yef++){
     if(Yef<(dat.Nyear))  {
       efcnt = efcnt+1;
       Tmppar = efpars(parcnt);                                         // read current par
       CpueEcreep(Yef,Nef) = CpueEcreep(Yef-1,Nef) * (1+(Tmppar/100)); // apply it
       if(efcnt==thedata.EffCrLag(Nef)){                                // THEN check lag
         efcnt = 0;
         parcnt = parcnt + 1;
       }
     } else { CpueEcreep(Yef,Nef) = CpueEcreep(Yef-1,Nef);}
   }
   }

 for (int Ipnt=0;Ipnt<thedata.Ncpue;Ipnt++)
  {
   IdataSet = thedata.IndexI(Ipnt,0);
   Ifleet = thedata.IndexI(Ipnt,1);
   Iarea = dat.Fleet_area(Ifleet);                                                                  // For Now
   Jsex = thedata.IndexI(Ipnt,2);
   Iyear = thedata.IndexI(Ipnt,3);
   Istep = thedata.IndexI(Ipnt,4);
   if (Jsex==-1) { Isex1 = 0; Isex2=1; } else { Isex1 = Jsex; Isex2 = Jsex; }
   for (int Iage=0;Iage<dat.Nage;Iage++){
    for (int Isex=Isex1;Isex<=Isex2;Isex++){
      if(dat.IsRed(Isex,Iage,Iarea,Istep)==1) {ScaleRedQ = QRedsPar;} else {ScaleRedQ = 1.0;}
     for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++)  {
       SelPointer = dat.SelPnt(Isex,Iage,Ifleet,Iyear,Istep);
       selexFU = ActSelex(SelPointer,Ilen) * ScaleRedQ;
       RetPointer = dat.RetPnt(Isex,Iage,Ifleet,Iyear,Istep);
       LegalPointer = dat.LegalFleetPnt(Isex,Iage,Ifleet,Iyear,Istep);
       retainFU = ActReten(RetPointer,Ilen)*ActLegal(LegalPointer,Ilen);
       if (thedata.IndexType(IdataSet)==1) selretwght = selexFU*retainFU*WeightLen(Isex,Ilen);
       if (thedata.IndexType(IdataSet)==2) selretwght = selexFU*retainFU;
       PredCpue(Ipnt,0) += selretwght*N(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Ilen)*exp(-Z(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Ilen)/2.0);
      }}}
   IenvPnt = thedata.EnvIndCpue(IdataSet)-1;
   // AEP update
  if (thedata.EnvIndCpue(IdataSet)>0) {
    PredCpue(Ipnt,0)*= exp(Qpars(IdataSet)*thedata.EnvData(dat.BurnIn+Iyear,Istep,IenvPnt));
    }
  if (thedata.EffCrIndCpue(IdataSet)>0) {
    PredCpue(Ipnt,0) *= CpueEcreep(Iyear,thedata.EffCrIndCpue(IdataSet)-1);
  }
  IndexPoint = thedata.TreatQcpue(IdataSet);
  Qval(IndexPoint) += log(thedata.IndexR(Ipnt,0)/PredCpue(Ipnt,0))/square(thedata.IndexR(Ipnt,1));
  Ndata(IndexPoint) += 1.0/square(thedata.IndexR(Ipnt,1));
  }

 // Compute the MLE for Q
 for (int Ipnt=0;Ipnt<thedata.NcpueDataSeries;Ipnt++){
  if (Ndata(Ipnt) >0) Qval(Ipnt) = exp(Qval(Ipnt)/Ndata(Ipnt));}

 // Compute Sigma and hence the likelihood
 SS.setZero();  Ndata.setZero();
 for (int Ipnt=0;Ipnt<thedata.Ncpue;Ipnt++)
  {
   IdataSet = thedata.IndexI(Ipnt,0);
   IndexPoint2 = thedata.FixSigmaCpue(IdataSet);
   PredCpue(Ipnt,0) = Qval(thedata.TreatQcpue(IdataSet))*PredCpue(Ipnt,0);
   PredCpue(Ipnt,1) = log(thedata.IndexR(Ipnt,0)/PredCpue(Ipnt,0))/thedata.IndexR(Ipnt,1);
   SS(IndexPoint2) += square(PredCpue(Ipnt,1));
   Ndata(IndexPoint2) += 1.0;
  }
 for (int IdataSet=0;IdataSet<thedata.NcpueDataSeries;IdataSet++)
  if (Ndata(IdataSet) >0)
   {
    // Note that this account for the minimum sigma
    SigmaCpue(IdataSet) = sqrt(SS(IdataSet)/Ndata(IdataSet));
    SigmaUse = thedata.SigmaCpueOffset-SigmaCpue(IdataSet);
    SigmaUse = SigmaCpue(IdataSet) + SigmaUse /(1+exp(-10.0*SigmaUse));
    CpueLikeComps(IdataSet) = Ndata(IdataSet)*log(SigmaUse)+Ndata(IdataSet)/2.0;
    NeglogLikelihood += thedata.LambdaCpue2(IdataSet)*CpueLikeComps(IdataSet);
   }
 return(NeglogLikelihood);

}

// -------------------------------------------------------------------------------------------------------------------

template <class Type>
 Type NumbersLikelihood(dataSet<Type> &dat, TheData<Type> &thedata, array<Type> &N, array<Type> &Z, array<Type> &Hrate,
                         matrix<Type> &ActSelex, matrix<Type> &ActReten, matrix<Type> &ActLegal,
                         matrix<Type> &WeightLen, matrix<Type> &PredNumbers,
                         vector<Type> &NumbersLikeComps, vector<Type> &SigmaNumbers, Type QRedsPar) {

 Type selexFU,retainFU,selretwght,Z2,SigmaUse,ScaleRedQ;
 int Ifleet,Iyear,Istep,PntCnt,Iarea,IdataSet,IndexPoint;
 int SelPointer,RetPointer,LegalPointer;
 vector <Type> SS(thedata.NcatchDataSeries);
 vector <Type> Ndata(thedata.NcatchDataSeries);

 Type NeglogLikelihood = 0;

 // Find predicted biomass
 PredNumbers.setZero();
 SS.setZero();
 Ndata.setZero();
 for (int Ipnt=0;Ipnt<thedata.Nnumbers;Ipnt++)
  {
   IdataSet = thedata.NumbersI(Ipnt,0);
   Ifleet = thedata.NumbersI(Ipnt,1);
   Iarea = dat.Fleet_area(Ifleet);                                                                  // For Now
   Iyear = thedata.NumbersI(Ipnt,2);
   Istep = thedata.NumbersI(Ipnt,3);
   for (int Iage=0;Iage<dat.Nage;Iage++){
    for (int Isex=0;Isex<dat.Nsex;Isex++){
     if(dat.IsRed(Isex,Iage,Iarea,Istep)==1) {ScaleRedQ = QRedsPar;} else {ScaleRedQ = 1.0;}
     for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++)  {
       SelPointer = dat.SelPnt(Isex,Iage,Ifleet,Iyear,Istep);
	     selexFU = ActSelex(SelPointer,Ilen) * ScaleRedQ;
	     RetPointer = dat.RetPnt(Isex,Iage,Ifleet,Iyear,Istep);
       LegalPointer = dat.LegalFleetPnt(Isex,Iage,Ifleet,Iyear,Istep);
       retainFU = ActReten(RetPointer,Ilen)*ActLegal(LegalPointer,Ilen);
	     Z2 = (1-exp(-Z(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Ilen)))/Z(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Ilen);
       PredNumbers(Ipnt,0) += Hrate(dat.BurnIn+Iyear,Istep,Ifleet)*selexFU*retainFU*N(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Ilen)*Z2;
      }
     }
    }
   PredNumbers(Ipnt,1) = log(thedata.NumbersR(Ipnt,0)/PredNumbers(Ipnt,0))/thedata.NumbersR(Ipnt,1);
   IndexPoint = thedata.FixSigmaCatchN(IdataSet);
   SS(IndexPoint) += PredNumbers(Ipnt,1)*PredNumbers(Ipnt,1);
   Ndata(IndexPoint) += 1;
  }

 for (int IdataSet=0;IdataSet<thedata.NcatchDataSeries;IdataSet++)
  if (Ndata(IdataSet) > 0)
   {
    SigmaNumbers(IdataSet) = sqrt(SS(IdataSet)/Ndata(IdataSet));
    SigmaUse = thedata.SigmaCatchNOffset-SigmaNumbers(IdataSet);
    SigmaUse = SigmaNumbers(IdataSet) + SigmaUse /(1+exp(-10.0*SigmaUse));
    NumbersLikeComps(IdataSet) = Ndata(IdataSet)*log(SigmaUse)+Ndata(IdataSet)/2.0;
    NeglogLikelihood += thedata.LambdaNumbers2(IdataSet)*NumbersLikeComps(IdataSet);
   }
 return(NeglogLikelihood);
}
// -------------------------------------------------------------------------------------------------------------------

template <class Type>
 Type LengthLikelihood(dataSet<Type> &dat, TheData<Type> &thedata, array<Type> &N,
                         matrix<Type> &ActSelex, matrix<Type> &ActReten, matrix<Type> &ActLegal, matrix<Type> &PredLengthComp,
                         vector<Type> &LengthLikeComps, vector<Type> &Select, Type QRedsPar) {

 Type NeglogLikelihood;
 Type selexFU,retainFU,Total, Contrib,ScaleRedQ;
 int Iarea,Ifleet,Isex,Iyear,Istep;
 int SelPointer,RetPointer,LegalPointer;

 NeglogLikelihood = 0;

 PredLengthComp.setZero();
  for (int Ipnt=0;Ipnt<thedata.NlenComp;Ipnt++)
   {
    Ifleet = thedata.LenCompI(Ipnt,0);
    Iarea = dat.Fleet_area(Ifleet);                                                                  // For Now
    Isex = thedata.LenCompI(Ipnt,1);
    Iyear = thedata.LenCompI(Ipnt,2);
    Istep = thedata.LenCompI(Ipnt,3);
    for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++) {
     for (int Iage=0;Iage<dat.Nage;Iage++)  {
       if(dat.IsRed(Isex,Iage,Iarea,Istep)==1) {ScaleRedQ = QRedsPar;} else {ScaleRedQ = 1.0;}
       SelPointer = dat.SelPnt(Isex,Iage,Ifleet,Iyear,Istep);
       selexFU = ActSelex(SelPointer,Ilen) * ScaleRedQ;
       RetPointer = dat.RetPnt(Isex,Iage,Ifleet,Iyear,Istep);
       LegalPointer = dat.LegalFleetPnt(Isex,Iage,Ifleet,Iyear,Istep);
       retainFU = ActReten(RetPointer,Ilen)*ActLegal(LegalPointer,Ilen);
       PredLengthComp(Ipnt,Ilen) += selexFU*retainFU*N(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Ilen);
       if(Ipnt==20) Select(Ilen) = selexFU*retainFU;
      }
     }
    Total = 0;
    for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++) Total += PredLengthComp(Ipnt,Ilen);
    for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++) PredLengthComp(Ipnt,Ilen) /= Total;
   }


 // Now calculate the likelihood
 LengthLikeComps.setZero();
  for (int Ipnt=0;Ipnt<thedata.NlenComp;Ipnt++)
   {
    Ifleet = thedata.LenCompI(Ipnt,0);
    Istep = thedata.LenCompI(Ipnt,3);
    Isex = thedata.LenCompI(Ipnt,1);
    for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++)
      {
       Contrib = thedata.Stage1W(Ipnt)*(thedata.LenCompR(Ipnt,Ilen)+1e-5)*log((PredLengthComp(Ipnt,Ilen)+1e-5)/(thedata.LenCompR(Ipnt,Ilen)+1e-5));
       LengthLikeComps(Ifleet) -= Contrib;
       NeglogLikelihood -= thedata.LambdaLength2(Ifleet,Istep,Isex)*Contrib;
      }
   }

 return(NeglogLikelihood);
}

// -------------------------------------------------------------------------------------------------------------------

template <class Type>
 Type LarvalLikelihood(dataSet<Type> &dat, TheData<Type> &thedata,  matrix<Type> &RecruitmentByArea,
    matrix<Type> &PuerulusByArea, matrix<Type> &PredLarval, vector<Type> &LarvalLikeComps,
    vector<Type> &PuerPowPars) {

 Type NeglogLikelihood;
 NeglogLikelihood = 0;

 int Iarea, Idata,Iyr;
 Type Obs,CV,ncnt,SS,Residual;
 vector<Type> qestLar(dat.Narea);

 if(thedata.NLarvalData>0){
 for (Iarea=0;Iarea<dat.Narea;Iarea++){
   for(Iyr=0;Iyr<dat.BurnIn+dat.Nyear+dat.Nproj;Iyr++){
     PuerulusByArea(Iarea,Iyr) = exp(log(RecruitmentByArea(Iarea,Iyr)) / PuerPowPars(Iarea));
   }}

  LarvalLikeComps.setZero();
  PredLarval.setZero();
  for (Iarea=0;Iarea<dat.Narea;Iarea++)
   {
    // Calculate the q-value
    qestLar(Iarea) = 0;
    ncnt = 0;
    for (Idata=0;Idata<thedata.NLarvalData;Idata++)
     if (thedata.Lar_dataI(Idata,0) == Iarea)
      {
       Iyr = thedata.Lar_dataI(Idata,1)+thedata.Larval_Offset;
       Obs = thedata.Lar_dataR(Idata,0);
       if (thedata.LarvalLikeOpt == 0)
        CV = thedata.Lar_dataR(Idata,1)/Obs;
       else
        CV = thedata.Lar_dataR(Idata,2);
       if (Iyr < dat.BurnIn+dat.Nyear+dat.Nproj)
        {
         if (thedata.LarvalLikeOpt == 0)
          {
           qestLar(Iarea) += log(Obs/PuerulusByArea(Iarea,Iyr))/(CV*CV);
           ncnt += 1.0/(CV*CV);
          }
         else
          {
           qestLar(Iarea) += Obs*PuerulusByArea(Iarea,Iyr)/(CV*CV);
           ncnt += Obs*Obs/(CV*CV);
         }
        }
     }
    if (ncnt > 0)
     if (thedata.LarvalLikeOpt == 0)
      qestLar(Iarea) = exp(qestLar(Iarea)/ncnt);
     else
      qestLar(Iarea) = qestLar(Iarea) / ncnt;

    // Find the likelihood itself
    SS = 0;
    for (Idata=0;Idata<thedata.NLarvalData;Idata++)
     if (thedata.Lar_dataI(Idata,0) == Iarea)
      {
       Iyr = thedata.Lar_dataI(Idata,1)+thedata.Larval_Offset;
       Obs = thedata.Lar_dataR(Idata,0);
       if (thedata.LarvalLikeOpt == 0)
        CV = thedata.Lar_dataR(Idata,1)/Obs;
       else
        CV = thedata.Lar_dataR(Idata,1);
       if (Iyr < dat.BurnIn+dat.Nyear+dat.Nproj)
        {
         PredLarval(Idata,0) = qestLar(Iarea)*PuerulusByArea(Iarea,Iyr);
         if (thedata.LarvalLikeOpt == 0)
          Residual = (log(Obs) - log(qestLar(Iarea)*PuerulusByArea(Iarea,Iyr)))/CV;
         else
          Residual = (Obs - qestLar(Iarea)*PuerulusByArea(Iarea,Iyr))/CV;
         PredLarval(Idata,1) = Residual;
         SS += Residual*Residual/2.0;
	    }
      }
     LarvalLikeComps(Iarea) = SS;
     NeglogLikelihood += SS;
   }}

 return(NeglogLikelihood);
}


// -------------------------------------------------------------------------------------------------------------------

template <class Type>
 array<Type> VirginN(dataSet<Type> &dat, array<Type> &Z, array<Type> &Hrate,
                         matrix<Type> &ActSelex, matrix<Type> &ActReten, matrix<Type> &ActLegal, matrix<Type> &ActMove,
                         matrix<Type> &WeightLen, matrix<Type> &M,
                         array<Type> &ActGrowth, matrix<Type> &RecruitFrac, Type Rbar,
                         array <Type> &ActRecruitAreaSexDist, array <Type> &ActRecruitLenDist,
                         vector<Type> &ActRecDev, Type QRedsPar, Type MWhitesPar, Type Finitial) {

array<Type> N(dat.Narea,dat.Nsex,dat.Nage,dat.MaxLen);
Type TotalRec,ScaleRedQ,ScaleWhiteM;
int Ipnt,Jpnt,GrowthPointer,MovePointer,IdestArea,RecruitLenPointer;
int MatSize,Offset;
MatSize = dat.Narea*dat.Nage*dat.MaxLen;                                 // Full matrix size
matrix<Type> I(MatSize,MatSize);                                         // Identity matrix
matrix<Type> S(MatSize,MatSize);                                         // Survival matrix
matrix<Type> X(MatSize,MatSize);                                         // Transition matrix
matrix<Type> MM(MatSize,MatSize);                                        // Movement matrix
matrix<Type> A(MatSize,MatSize);                                         // Aging matrix
matrix<Type> Mat2(MatSize,MatSize);                                      // Temp matrix
matrix<Type> Trans(dat.MaxLen,dat.MaxLen);                               // Cumulative transition matrix
matrix<Type> Trans2(dat.MaxLen,dat.MaxLen);                              // Temporary transition matrix
vector<Type> MoveVec(dat.MaxLen);                                        // Movement vector
matrix<Type> RecVec(MatSize,1);                                          // Recruitment
matrix<Type> TestVec(MatSize,1);                                         // The equilibrium

// Find the equilibrium
N.setZero();
for (int Isex=0;Isex<dat.Nsex;Isex++)
 {

  I.setZero(); S.setZero(); MM.setZero(); X.setZero();  A.setZero();

  // Set the diagnonal matrices
  for (int Iarea=0;Iarea<dat.Narea;Iarea++)
   for (int Iage=0;Iage<dat.Nage;Iage++)
    {
     if(dat.IsRed(Isex,Iage,Iarea,0)==0) {ScaleWhiteM = MWhitesPar;} else {ScaleWhiteM = 1.0;}
     if(dat.IsRed(Isex,Iage,Iarea,0)==1) {ScaleRedQ = QRedsPar;} else {ScaleRedQ = 1.0;}
     Offset = Iarea*dat.Nage*dat.MaxLen+Iage*dat.MaxLen;
     for (int Isize=0;Isize<dat.MaxLen;Isize++)
      {
       I(Offset+Isize,Offset+Isize) = 1.0;
       S(Offset+Isize,Offset+Isize) = exp(-M(Iarea,Iage)*ScaleWhiteM+ActSelex(Isex,Isize)*ScaleRedQ*Finitial);
      }
     for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
      {
       MM(Offset+Isize,Offset+Isize) = 1.0;
      }
    }

  // Growth matrix
  for (int Iarea=0;Iarea<dat.Narea;Iarea++)
   for (int Iage=0;Iage<dat.Nage;Iage++)
    {
     // No growth
     Trans.setZero();
     for (int Isize=0;Isize<dat.Nlen(Isex);Isize++) Trans(Isize,Isize) = 1;

     // Multiply by growth matrix
     for (int Istep=0;Istep<dat.Nstep;Istep++)
      {
       GrowthPointer = dat.GrowthPnt(Iarea,Isex,Iage,0,Istep);
       if (GrowthPointer >=0)
        {
         Trans2.setZero();
         for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
          for (int Jsize=0;Jsize<=Isize;Jsize++) Trans2(Isize,Jsize) = ActGrowth(GrowthPointer,Isize,Jsize);
         Trans = atomic::matmul(Trans2,Trans);
        }
      } // Growth

     Offset = Iarea*dat.Nage*dat.MaxLen+Iage*dat.MaxLen;
     for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
      for (int Jsize=0;Jsize<dat.Nlen(Isex);Jsize++)
       X(Offset+Isize,Offset+Jsize) = Trans(Isize,Jsize);
    }

  // Movement
  for (int Istep=0;Istep<dat.Nstep;Istep++)
   {
    for (int Iarea=0;Iarea<dat.Narea;Iarea++)
     for (int Iage=0;Iage<dat.Nage;Iage++)
      {
       MovePointer = dat.MovePnt(Iarea,Iage,0,Istep);
       if (MovePointer > 0)
        {
         IdestArea = dat.MoveSpec(MovePointer,2);
         for (int Isize=0;Isize<dat.Nlen(Isex);Isize++) MoveVec(Isize) = ActMove(MovePointer,Isize);
         for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
          {
           Ipnt = IdestArea*dat.Nage*dat.MaxLen+Iage*dat.MaxLen+Isize;
           Jpnt = Iarea*dat.Nage*dat.MaxLen+Iage*dat.MaxLen+Isize;
           MM(Jpnt,Jpnt) -= MoveVec(Isize);
           MM(Ipnt,Jpnt) = MoveVec(Isize);
          }
        }
      }
   } // If there was a move

  // Ageing
  Ipnt = 0;
  for (int Iarea=0;Iarea<dat.Narea;Iarea++)
   {
    Ipnt = Iarea*dat.Nage*dat.MaxLen;
    for (int Iage=0;Iage<dat.Nage-1;Iage++)
     for (int Isize=0;Isize<dat.Nlen(Isex);Isize++) A(Ipnt+(Iage+1)*dat.MaxLen+Isize, Ipnt+Iage*dat.MaxLen+Isize) = 1;
    for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)  A(Ipnt+(dat.Nage-1)*dat.MaxLen+Isize, Ipnt+(dat.Nage-1)*dat.MaxLen+Isize) = 1;
   }

  // Matrix multiplication
  Mat2 = atomic::matmul(X,S);
  Mat2 = atomic::matmul(MM,Mat2);
  Mat2 = atomic::matmul(A,Mat2);
  Mat2 = I - Mat2;

  // Inverse
  Mat2 = atomic::matinv(Mat2);

  // Recruitment
  Offset = -1;
  RecVec.setZero();
  for (int Iarea=0;Iarea<dat.Narea;Iarea++)
   {
    TotalRec = 0;
    for (int Istep=0;Istep<dat.Nstep;Istep++) TotalRec += ActRecruitAreaSexDist(0,Istep,Iarea,Isex);
    RecruitLenPointer = dat.RecruitLenPnt(Iarea);

    for (int Iage=0;Iage<dat.Nage;Iage++)
     {
      for(int Isize=0;Isize<dat.Nlen(Isex);Isize++)
       {
        Offset += 1;
        if (Iage==0) RecVec(Offset,0) = TotalRec*ActRecruitLenDist(RecruitLenPointer,Isex,Isize);
       }
     }
   }

  // Solve for an equiilbrium
  TestVec = atomic::matmul(Mat2,RecVec);

  // Paste back
  for (int Iarea=0;Iarea<dat.Narea;Iarea++)
   for (int Iage=0;Iage<dat.Nage;Iage++)
    {
     Offset = Iarea*dat.Nage*dat.MaxLen+Iage*dat.MaxLen;
     for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
      N(Iarea,Isex,Iage,Isize) = TestVec(Offset+Isize,0)*exp(Rbar);
    }

 } //

return(N);

}

// -------------------------------------------------------------------------------------------------------------------

template <class Type>
 Type InitializeN(dataSet<Type> &dat, array<Type> &N, array<Type> &Z, array<Type> &Hrate,
                         matrix<Type> &ActSelex, matrix<Type> &ActReten, matrix<Type> &ActLegal, matrix<Type> &ActMove,
                         matrix<Type> &WeightLen, matrix<Type> &M,
                         array <Type> &ActGrowth, matrix<Type> &RecruitFrac, Type Rbar,
                         array <Type> &ActRecruitAreaSexDist, array <Type> &ActRecruitLenDist,
                         vector<Type> &ActRecDev, array<Type> &Ninit, vector<Type> &MatBio, matrix<Type> &MatBioArea,
                         matrix<Type> &RecruitmentByArea, vector<Type> BiasMult, Type SigmaR, Type QRedsPar, Type MWhitesPar,
                         vector<Type> &VirginBio, vector<Type> &VirginLegalBio, array<Type> &VirginNvec, array<Type> &VirginBioAtLen,
                         matrix<Type> &LegalRef,vector<Type> &CurrentBio) {

  Type Initial_pen;
  int IsVirgin;                                                           // Set to 1 for unfished state
  matrix<Type> Feqn2(dat.Nfleet,dat.Nstep); Feqn2.setZero();              // Initial F (not used in projections)
  matrix<Type> Feqn3(dat.Nfleet,dat.Nstep); Feqn3.setZero();              // Initial F (not used in projections)
  matrix<Type> Feqn4(dat.Nfleet,dat.Nstep); Feqn4.setZero();              // Used to set Burn_in F to zero if one rarea does not want burn in
  vector<Type> XX(2);                                                     // Dummy variables
  array<Type> Fvals(dat.Nfleet,dat.Num_Iteration,dat.Nstep);                             // Storage for tuning of Fs
  array<Type> Ninit2(dat.Narea,dat.Nsex,dat.Nage,dat.MaxLen);

  Ninit = VirginN(dat, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, ActGrowth, RecruitFrac,
         Rbar, ActRecruitAreaSexDist, ActRecruitLenDist, ActRecDev, QRedsPar, MWhitesPar, Type(0.0));
  //  Get bare bones numbers by area, sex, age and length - one recruitment / move / grow - no F Mort.

  // Virgin Biomass from Ninit -  Has M but not F - does not work correctly as changes slightly with burn in below. But is a temporary starting point
  VirginBio.setZero();
  for (int Iarea=0;Iarea<dat.Narea;Iarea++){
    for (int Isex=0;Isex<dat.Nsex;Isex++){
      for (int Iage=0;Iage<dat.Nage;Iage++){
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++) {
          VirginBio(Iarea) += Ninit(Iarea,Isex,Iage,Isize) * WeightLen(Isex,Isize);}}}}  // Weight in each area, first year/time-step

  // Now compute
  IsVirgin = 0;
  for (int JJ=0;JJ<=dat.Num_Iteration-1;JJ++)
   {
    // Really poor initial condition - put unfished numbers in  first year
    N.setZero();
    for (int Iarea=0;Iarea<dat.Narea;Iarea++)
     for (int Isex=0;Isex<dat.Nsex;Isex++)
      for (int Iage=0;Iage<dat.Nage;Iage++)
       for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
        {  N(Iarea,0,0,Isex,Iage,Isize) = Ninit(Iarea,Isex,Iage,Isize); }


    // One year zero catch projection  This updates the future time-step
    for (int Iyear=-dat.BurnIn;Iyear<-dat.BurnIn+1;Iyear++)
     for (int Istep=0;Istep<dat.Nstep;Istep++)
      XX = OneTimeStep(dat, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, Iyear, Istep, ActGrowth, RecruitFrac, Rbar, IsVirgin, Feqn3,ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,
                        VirginBio, CurrentBio);

    // Multiyear projection with No F (F set to Zero) This updates the future time-step under no fishing
    for (int Iyear=-dat.BurnIn+1;Iyear<dat.Tune_Years;Iyear++)
     for (int Istep=0;Istep<dat.Nstep;Istep++)
      XX = OneTimeStep(dat, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, Iyear, Istep, ActGrowth, RecruitFrac, Rbar, IsVirgin, Feqn2,ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,
                        VirginBio, CurrentBio);

    for (int Iarea=0;Iarea<dat.Narea;Iarea++)
     for (int Istep=0;Istep<dat.Nstep;Istep++)
      for (int Ifleet=0;Ifleet<dat.Nfleet;Ifleet++)
       if (dat.Area_fleet(Iarea,Ifleet)==1)
        {
         Fvals(Ifleet,JJ,Istep) = Feqn2(Ifleet,Istep);
         Feqn2(Ifleet,Istep) = 0;
          for (int Iyear=0;Iyear<dat.Tune_Years;Iyear++) {Feqn2(Ifleet,Istep) +=  Hrate(dat.BurnIn+Iyear,Istep,Ifleet) ;}
         Feqn2(Ifleet,Istep) /= float(dat.Tune_Years);
        }
    }

  //  At this point we have a better Virgin Biomass created.  Use this for output etc.
  VirginBio.setZero(); VirginLegalBio.setZero();
  for (int Iarea=0;Iarea<dat.Narea;Iarea++){
    for (int Isex=0;Isex<dat.Nsex;Isex++){
      for (int Iage=0;Iage<dat.Nage;Iage++){
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++) {
          VirginBio(Iarea) += N(Iarea,0,0,Isex,Iage,Isize) * WeightLen(Isex,Isize);
          VirginLegalBio(Iarea) += LegalRef(Isex,Isize) * N(Iarea,0,dat.BioTimeStep,Isex,Iage,Isize) * WeightLen(Isex,Isize);
          }}}}

  // Penalty on non-convergence
  Initial_pen = 0;
  for (int Ifleet=0;Ifleet<dat.Nfleet;Ifleet++)
   for (int Istep=0;Istep<dat.Nstep;Istep++)
    Initial_pen += 1000000*square(Fvals(Ifleet,dat.Num_Iteration-2,Istep)-Fvals(Ifleet,dat.Num_Iteration-1,Istep));

  // ── Compute Virgin Biomass from a dedicated no-F burn-in ──────
  // Project from unfished equilibrium through the full burn-in with F=0
  // so that virgin biomass reflects the model's actual process sequence
  N.setZero();
  for (int Iarea=0;Iarea<dat.Narea;Iarea++)
    for (int Isex=0;Isex<dat.Nsex;Isex++)
      for (int Iage=0;Iage<dat.Nage;Iage++)
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
        {  N(Iarea,0,0,Isex,Iage,Isize) = Ninit(Iarea,Isex,Iage,Isize); }

        for (int Iyear=-dat.BurnIn;Iyear<0;Iyear++)
          for (int Istep=0;Istep<dat.Nstep;Istep++)
            XX = OneTimeStep(dat, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M,
                             Iyear, Istep, ActGrowth, RecruitFrac, Rbar, IsVirgin, Feqn3,
                             ActRecruitAreaSexDist, ActRecruitLenDist, ActRecDev, MatBio, MatBioArea,
                             RecruitmentByArea, BiasMult, SigmaR, QRedsPar, MWhitesPar,
                             VirginBio, CurrentBio);

  // Virgin biomass from end of no-F burn-in at the designated biology time step
  VirginBio.setZero(); VirginLegalBio.setZero();VirginNvec.setZero(); VirginBioAtLen.setZero();
  for (int Iarea=0;Iarea<dat.Narea;Iarea++){
    for (int Isex=0;Isex<dat.Nsex;Isex++){
      for (int Iage=0;Iage<dat.Nage;Iage++){
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++) {
          VirginBio(Iarea) += N(Iarea,dat.BurnIn,0,Isex,Iage,Isize) * WeightLen(Isex,Isize);
          VirginLegalBio(Iarea) += LegalRef(Isex,Isize) * N(Iarea,dat.BurnIn,dat.BioTimeStep,Isex,Iage,Isize) * WeightLen(Isex,Isize);
          VirginNvec(Iarea,Isex,Iage,Isize) = N(Iarea,0,0,Isex,Iage,Isize);
          VirginBioAtLen(Iarea,Isex,Isize) += N(Iarea,0,0,Isex,Iage,Isize) * WeightLen(Isex,Isize);
        }}}}

  // ── Now do the final burn-in with area-specific F ─────────────
  N.setZero();
  for (int Iarea=0;Iarea<dat.Narea;Iarea++)
    for (int Isex=0;Isex<dat.Nsex;Isex++)
      for (int Iage=0;Iage<dat.Nage;Iage++)
        for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
        {  N(Iarea,0,0,Isex,Iage,Isize) = Ninit(Iarea,Isex,Iage,Isize); }

        // One year zero catch projection
        for (int Iyear=-dat.BurnIn;Iyear<-dat.BurnIn+1;Iyear++)
          for (int Istep=0;Istep<dat.Nstep;Istep++)
            XX = OneTimeStep(dat, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M,
                             Iyear, Istep, ActGrowth, RecruitFrac, Rbar, IsVirgin, Feqn3,
                             ActRecruitAreaSexDist, ActRecruitLenDist, ActRecDev, MatBio, MatBioArea,
                             RecruitmentByArea, BiasMult, SigmaR, QRedsPar, MWhitesPar,
                             VirginBio, CurrentBio);

  // Multiyear projection with area-specific F
  for (int Iyear=-dat.BurnIn+1;Iyear<0;Iyear++){
    for (int Istep=0;Istep<dat.Nstep;Istep++){
      for (int Iarea=0;Iarea<dat.Narea;Iarea++){
        for (int Ifleet=0; Ifleet<dat.Nfleet;Ifleet++) {
          if (dat.Area_fleet(Iarea,Ifleet)==1){
            Feqn4(Ifleet,Istep) = Feqn2(Ifleet,Istep);
            if(-dat.BurnInVec(Iarea)>Iyear) Feqn4(Ifleet,Istep) = 0; }}}
      XX = OneTimeStep(dat, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M,
                       Iyear, Istep, ActGrowth, RecruitFrac, Rbar, IsVirgin, Feqn4,
                       ActRecruitAreaSexDist, ActRecruitLenDist, ActRecDev, MatBio, MatBioArea,
                       RecruitmentByArea, BiasMult, SigmaR, QRedsPar, MWhitesPar,
                       VirginBio, CurrentBio);}}

 return(Initial_pen);

}

// ========================================================================================================================

template <class Type>
 vector<Type> CatchByNumAge(dataSet<Type> &dat, array<Type> &N, array<Type> &Z,array<Type> &Hrate,
                           matrix<Type> &ActSelex, matrix<Type> &ActReten, matrix<Type> &ActLegal, matrix<Type> &WeightLen,
                           int Iarea, int Ifleet, int Iyear, int Istep, Type QRedsPar) {

 vector<Type> XX(2);                                                     // Outputs
 array<Type> selexF(dat.Nsex,dat.Nage,dat.MaxLen);                       // Selectivity
 array<Type> retainF(dat.Nsex,dat.Nage,dat.MaxLen);                      // Retention
 int SelPointer,RetPointer,LegalPointer;                                 // Pointers
 Type Z2,CAL,ScaleRedQ;                                                  // Temporary

 // Need to set selectivity
 for (int Isex=0;Isex<dat.Nsex;Isex++) {
  for (int Iage=0;Iage<dat.Nage;Iage++) {
    if(dat.IsRed(Isex,Iage,Iarea,Istep)==1) {ScaleRedQ = QRedsPar;} else {ScaleRedQ = 1.0;}
    SelPointer = dat.SelPnt(Isex,Iage,Ifleet,Iyear,Istep);
    RetPointer = dat.RetPnt(Isex,Iage,Ifleet,Iyear,Istep);
    LegalPointer = dat.LegalFleetPnt(Isex,Iage,Ifleet,Iyear,Istep);
    for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++) {
      selexF(Isex,Iage,Ilen) = ActSelex(SelPointer,Ilen) * ScaleRedQ;
      retainF(Isex,Iage,Ilen) = ActReten(RetPointer,Ilen)*ActLegal(LegalPointer,Ilen);
      }
    }
  }

 // Reset
 XX.setZero();
 for (int Isex=0;Isex<dat.Nsex;Isex++)
  for (int Iage=0;Iage<dat.Nage;Iage++)
   for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
    {
     Z2 = (1-exp(-Z(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Isize)))/Z(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Isize);
     CAL = selexF(Isex,Iage,Isize)*retainF(Isex,Iage,Isize)*Hrate(dat.BurnIn+Iyear,Istep,Ifleet)*N(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Isize)*Z2;
     XX(0) += CAL;
     XX(1) += WeightLen(Isex,Isize)*CAL;
    }

 return(XX);
}

// ========================================================================================================================
//  Discard weight calculation
template <class Type>
vector<Type> DiscardByFleet(dataSet<Type> &dat, array<Type> &N, array<Type> &Z, array<Type> &Hrate,
                            matrix<Type> &ActSelex, matrix<Type> &ActReten, matrix<Type> &ActLegal,
                            matrix<Type> &WeightLen,
                            int Iarea, int Ifleet, int Iyear, int Istep, Type QRedsPar) {

  vector<Type> XX(2);
  array<Type> selexF(dat.Nsex,dat.Nage,dat.MaxLen);
  array<Type> retainF(dat.Nsex,dat.Nage,dat.MaxLen);
  int SelPointer, RetPointer, LegalPointer;
  Type Z2, DiscWt, ScaleRedQ;

  for (int Isex=0;Isex<dat.Nsex;Isex++) {
    for (int Iage=0;Iage<dat.Nage;Iage++) {
      if(dat.IsRed(Isex,Iage,Iarea,Istep)==1) {ScaleRedQ = QRedsPar;} else {ScaleRedQ = 1.0;}
      if (Iyear < dat.Nyear)
       {
        SelPointer = dat.SelPnt(Isex,Iage,Ifleet,Iyear,Istep);
        RetPointer = dat.RetPnt(Isex,Iage,Ifleet,Iyear,Istep);
        LegalPointer = dat.LegalFleetPnt(Isex,Iage,Ifleet,Iyear,Istep);
       }
      else
       {
        // Projection years: same branch OneTimeStep() uses for Iyear>=Nyear.
        SelPointer = dat.SelPntFut(Isex,Iage,Ifleet,Iyear-dat.Nyear,Istep);
        RetPointer = dat.RetPntFut(Isex,Iage,Ifleet,Iyear-dat.Nyear,Istep);
        LegalPointer = dat.LegalFleetPntFut(Isex,Iage,Ifleet,Iyear-dat.Nyear,Istep);
       }
      for (int Ilen=0;Ilen<dat.Nlen(Isex);Ilen++) {
        selexF(Isex,Iage,Ilen) = ActSelex(SelPointer,Ilen) * ScaleRedQ;
        retainF(Isex,Iage,Ilen) = ActReten(RetPointer,Ilen) * ActLegal(LegalPointer,Ilen);
      }
    }
  }

  XX.setZero();
  for (int Isex=0;Isex<dat.Nsex;Isex++)
    for (int Iage=0;Iage<dat.Nage;Iage++)
      for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
      {
        Z2 = (1-exp(-Z(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Isize))) /
          Z(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Isize);
        DiscWt = selexF(Isex,Iage,Isize) * (1.0 - retainF(Isex,Iage,Isize)) *
          Hrate(dat.BurnIn+Iyear,Istep,Ifleet) *
          N(Iarea,dat.BurnIn+Iyear,Istep,Isex,Iage,Isize) * Z2 *
          WeightLen(Isex,Isize);
        XX(0) += DiscWt;
        XX(1) += dat.Phi(Ifleet,Iage,Iyear,Istep) * DiscWt;
      }

      return(XX);
}

// ========================================================================================================================

template <class Type>
 Type CatchLikelihood(dataSet<Type> &dat, TheData<Type> &thedata, array<Type> &N, array<Type> &Z,array<Type> &Hrate,
                           matrix<Type> &ActSelex, matrix<Type> &ActReten, matrix<Type> &ActLegal, matrix<Type> &WeightLen, array<Type> &CatchCheck, Type QRedsPar) {

 vector<Type> XX(2);                                                     // Catch predection pass
 int Iarea;
 Type NeglogLikelihood;                                                  // Negative log-likelihood

 NeglogLikelihood = 0;
 CatchCheck.setZero();
  for (int Iyear=0;Iyear<dat.Nyear;Iyear++)
   for (int Istep=0;Istep<dat.Nstep;Istep++)
    for (int Ifleet=0;Ifleet<dat.Nfleet;Ifleet++)
     {
	  Iarea = dat.Fleet_area(Ifleet);
      XX = CatchByNumAge(dat,N,Z,Hrate,ActSelex,ActReten,ActLegal,WeightLen,Iarea,Ifleet,Iyear,Istep,QRedsPar);
      CatchCheck(Iyear,Istep,Ifleet) = XX(1);
      if (dat.Catch(Iyear,Istep,Ifleet) > 0)
       NeglogLikelihood += square(CatchCheck(Iyear,Istep,Ifleet)-dat.Catch(Iyear,Istep,Ifleet));
     }
 return(NeglogLikelihood);
}

template <class Type>
Type TagDym(dataSet<Type> &dat, TheData<Type> &thedata, int SexPass, int GrpPass,array<Type> &N,
            matrix<Type> &ActSelex, matrix<Type> &ActReten, matrix<Type> &ActLegal,
            array<Type> &ActGrowth,matrix<Type> &ActMove,
            matrix<Type> &M, array<Type> &Hrate, Type QRedsPar, Type MWhitesPar,
            array<Type> &RecapNum, matrix<Type> &NotReported,
            matrix<Type> &TagLike1, matrix<Type> &TagLike2, array<Type> &PredTagSize) {

  int Jyear, Kyear;
  int SelPointer, RetPointer, LegalPointer, GrowthPointer, MovePointer, IsMoves, IdestArea;
  Type NtagRel, CumReleases;
  Type RetainTemp, NtotalT, ScaleRedQ, ScaleWhiteM;
  Type PartialF, TotalPartialF, FullF2, Deaths;
  Type ObsL, PredL, ObsSS, LikeSize1, LikeTag2, LikeCompT;
  Type TotalReported;
  vector<Type> Ntemp2(dat.MaxLen);
  vector<Type> MoveVec(dat.MaxLen);

  array<Type> selexF(dat.Nfleet, dat.Nsex, dat.Nage, dat.MaxLen);
  array<Type> retainF(dat.Nfleet, dat.Nsex, dat.Nage, dat.MaxLen);
  array<Type> M_rate_Tag(dat.Narea, dat.Nage, dat.MaxLen);
  array<Type> Z_rate_Tag(dat.Narea, dat.Nage, dat.MaxLen);
  array<Type> RecapTmp(dat.Narea, thedata.NrepSplit, dat.MaxLen);
  array<Type> Ntemp_Tag(thedata.NtagLag+1, dat.Narea, dat.Nage, dat.MaxLen);
  array<Type> Nmove_Tag(thedata.NtagLag+1, dat.Narea, dat.Nage, dat.MaxLen);

  // LOCAL Ntag with reduced dimensions (no sex/group dimension needed)
  array<Type> Ntag_local(dat.Narea, thedata.NtagLag+1, dat.Nage, dat.MaxLen);
  Ntag_local.setZero();

  Type NeglogLikelihood = 0;
  int NtagLag = thedata.NtagLag;
  int Nage = dat.Nage;
  int Nfleet = dat.Nfleet;
  int NrepSplit = thedata.NrepSplit;
  int Narea = dat.Narea;

  NotReported(SexPass, GrpPass) = 0;
  CumReleases = 0;

  int Year1Tag_val = thedata.Year1Tag(GrpPass);
  int start_year = Year1Tag_val - dat.First_yr;

  for (int Iyear=start_year; Iyear<dat.Nyear; Iyear++)
    for (int Istep=0; Istep<dat.Nstep; Istep++)
    {
      Jyear = Iyear+dat.First_yr-thedata.TagYr1;
      Kyear = dat.BurnIn+Iyear;

      // D1: Add the tags that have been out long enough
      if (NtagLag > 0) {
        for (int Iarea=0; Iarea<Narea; Iarea++)
          for (int Iage=0; Iage<dat.Nage; Iage++)
            for (int Isize=0; Isize<dat.MaxLen; Isize++)
            {
              Ntag_local(Iarea,0,Iage,Isize) += Ntag_local(Iarea,1,Iage,Isize);
              for (int ItagLag=NtagLag-1; ItagLag>0; ItagLag--)
                Ntag_local(Iarea,ItagLag,Iage,Isize) = Ntag_local(Iarea,ItagLag+1,Iage,Isize);
              Ntag_local(Iarea,NtagLag,Iage,Isize) = 0;
            }
      }

      // D2: Add new tags
      if (Jyear >= 0 && Jyear < thedata.NyearTags) {
        for (int Iarea=0; Iarea<Narea; Iarea++)
        {
          NtagRel = thedata.TagRel(SexPass,GrpPass,Iarea,Jyear,Istep,0);
          if (NtagRel > 0)
          {
            for (int Isize=0; Isize<dat.MaxLen; Isize++){
              if (thedata.TagRel(SexPass,GrpPass,Iarea,Jyear,Istep,Isize+1) > 0)
              {
                NtotalT = 0;
                for (int Iage=0; Iage<Nage; Iage++)
                  NtotalT += N(Iarea,Kyear,Istep,SexPass,Iage,Isize);
                for (int Iage=0; Iage<Nage; Iage++)
                  Ntag_local(Iarea,NtagLag,Iage,Isize) = N(Iarea,Kyear,Istep,SexPass,Iage,Isize)/NtotalT*thedata.TagRel(SexPass,GrpPass,Iarea,Jyear,Istep,Isize+1)*thedata.InitialLoss;
              }
            }  // <-- for loop closes here
            NotReported(SexPass,GrpPass) += thedata.TagRel(SexPass,GrpPass,Iarea,Jyear,Istep,0)*(1.0-thedata.InitialLoss);
            CumReleases += thedata.TagRel(SexPass,GrpPass,Iarea,Jyear,Istep,0);
          }
        }
      }

      // D3: Set selectivity and compute Z/recaptures
      for (int Ifleet=0; Ifleet<Nfleet; Ifleet++)
        for (int Iage=0; Iage<Nage; Iage++)
        {
          if(dat.IsRed(SexPass,Iage,dat.Fleet_area(Ifleet),Istep)==1) {ScaleRedQ = QRedsPar;} else {ScaleRedQ = 1.0;}
          SelPointer = dat.SelPnt(SexPass,Iage,Ifleet,Iyear,Istep);
          RetPointer = dat.RetPnt(SexPass,Iage,Ifleet,Iyear,Istep);
          LegalPointer = dat.LegalFleetPnt(SexPass,Iage,Ifleet,Iyear,Istep);
          for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
          {
            selexF(Ifleet,SexPass,Iage,Isize) = ActSelex(SelPointer,Isize) * ScaleRedQ;
            retainF(Ifleet,SexPass,Iage,Isize) = ActReten(RetPointer,Isize)*ActLegal(LegalPointer,Isize);
          }
        }

        RecapTmp.setZero();
      for (int Iarea=0; Iarea<Narea; Iarea++)
      {
        for (int Iage=0; Iage<Nage; Iage++)
        {
          if(dat.IsRed(SexPass,Iage,Iarea,Istep)==0) {ScaleWhiteM = MWhitesPar;} else {ScaleWhiteM = 1.0;}
          for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
          {
            M_rate_Tag(Iarea,Iage,Isize) = dat.TimeStepLen(Iyear,Istep)*M(Iarea,Iage)*ScaleWhiteM+dat.TimeStepLen(Iyear,Istep)*thedata.TagLossRate;
            Z_rate_Tag(Iarea,Iage,Isize) = M_rate_Tag(Iarea,Iage,Isize);
            for (int Ifleet=0; Ifleet<Nfleet; Ifleet++)
              if (dat.Area_fleet(Iarea,Ifleet)==1)
              {
                RetainTemp = selexF(Ifleet,SexPass,Iage,Isize) * (retainF(Ifleet,SexPass,Iage,Isize)+dat.Phi(Ifleet,Iage,Iyear,Istep)*(1.0-retainF(Ifleet,SexPass,Iage,Isize)));
                Z_rate_Tag(Iarea,Iage,Isize) += Hrate(Kyear,Istep,Ifleet)*RetainTemp;
              }

              for (int ItagLag=1; ItagLag<=NtagLag; ItagLag++)
                NotReported(SexPass,GrpPass) += Ntag_local(Iarea,ItagLag,Iage,Isize)*(1.0-exp(-M_rate_Tag(Iarea,Iage,Isize)));

              Deaths = Ntag_local(Iarea,0,Iage,Isize)*(1.0-exp(-Z_rate_Tag(Iarea,Iage,Isize)))/Z_rate_Tag(Iarea,Iage,Isize);
              NotReported(SexPass,GrpPass) += M_rate_Tag(Iarea,Iage,Isize) * Deaths;

              for (int Ifleet=0; Ifleet<Nfleet; Ifleet++)
                if (dat.Area_fleet(Iarea,Ifleet)==1)
                {
                  RetainTemp = selexF(Ifleet,SexPass,Iage,Isize) * (retainF(Ifleet,SexPass,Iage,Isize)+dat.Phi(Ifleet,Iage,Iyear,Istep)*(1.0-retainF(Ifleet,SexPass,Iage,Isize)));
                  FullF2 = Hrate(Kyear,Istep,Ifleet)*RetainTemp;
                  TotalPartialF = 0;
                  for (int IrepSplit=0; IrepSplit<NrepSplit; IrepSplit++)
                  {
                    PartialF = thedata.RepRate(IrepSplit)*thedata.PropRepSplit(Jyear,Istep,Iarea,IrepSplit)*FullF2;
                    TotalPartialF += PartialF;
                    RecapTmp(Iarea,IrepSplit,Isize) += PartialF*Deaths;
                  }
                  NotReported(SexPass,GrpPass) += (FullF2-TotalPartialF)*Deaths;
                }
          }
        }
      }

      // Total the tags
      if (Jyear >= 0 && Jyear < thedata.NyearTags) {
        for (int Iarea=0; Iarea<Narea; Iarea++)
          for (int IrepSplit=0; IrepSplit<NrepSplit; IrepSplit++)
            for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
              RecapNum(SexPass,GrpPass,Iarea,IrepSplit,Jyear,Istep) += RecapTmp(Iarea,IrepSplit,Isize);
      }

      // Likelihood for length-comp of recaptured
      if (Jyear >= 0 && Jyear < thedata.NyearTags) {
        LikeSize1 = 0;
        for (int Iarea=0; Iarea<Narea; Iarea++)
          for (int IrepSplit=0; IrepSplit<NrepSplit; IrepSplit++)
            if (thedata.FitTagSizes(IrepSplit) == 1) {
              if (thedata.TagRec(SexPass,GrpPass,Iarea,IrepSplit,Jyear,Istep,0)>0) {
                ObsSS = thedata.TagRec(SexPass,GrpPass,Iarea,IrepSplit,Jyear,Istep,0);
                NtotalT = 0;
                for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
                  NtotalT += RecapTmp(Iarea,IrepSplit,Isize);

                // PROTECTION: Only compute likelihood if we have predicted recaptures
                if (NtotalT > 1e-10) {
                  for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
                  {
                    PredL = RecapTmp(Iarea,IrepSplit,Isize)/NtotalT;
                    PredTagSize(SexPass,GrpPass,Iarea,Isize) += PredL*ObsSS;
                    if (thedata.TagRec(SexPass,GrpPass,Iarea,IrepSplit,Jyear,Istep,Isize+1)>0)
                    {
                      ObsL = thedata.TagRec(SexPass,GrpPass,Iarea,IrepSplit,Jyear,Istep,Isize+1)/ObsSS;
                      // Add small constant to avoid log(0)
                      LikeSize1 -= ObsL*ObsSS*log((PredL+1e-10)/(ObsL+1e-10));
                    }
                  }
                }
              }
            }
            TagLike1(SexPass,GrpPass) += LikeSize1;
      }
      // D4a: Remove mortality
      Ntemp_Tag.setZero();
      for (int Iarea=0; Iarea<Narea; Iarea++)
        for (int Iage=0; Iage<Nage; Iage++)
          for (int Isize=0; Isize<dat.MaxLen; Isize++)
            for (int ItagLag=0; ItagLag<=NtagLag; ItagLag++)
              if (ItagLag==0)
                Ntemp_Tag(0,Iarea,Iage,Isize) = Ntag_local(Iarea,0,Iage,Isize)*exp(-Z_rate_Tag(Iarea,Iage,Isize));
              else
                Ntemp_Tag(ItagLag,Iarea,Iage,Isize) = Ntag_local(Iarea,ItagLag,Iage,Isize)*exp(-M_rate_Tag(Iarea,Iage,Isize));

              // D4b: Growth
              for (int Iarea=0; Iarea<Narea; Iarea++)
                for (int ItagLag=0; ItagLag<=NtagLag; ItagLag++)
                  for (int Iage=0; Iage<Nage; Iage++)
                  {
                    GrowthPointer = dat.GrowthPnt(Iarea,SexPass,Iage,Iyear,Istep);
                    if (GrowthPointer >=0)
                    {
                      Ntemp2.setZero();
                      for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
                      {
                        for (int Jsize=0; Jsize<=Isize; Jsize++)
                          Ntemp2(Isize) += Ntemp_Tag(ItagLag,Iarea,Iage,Jsize)*ActGrowth(GrowthPointer,Isize,Jsize);
                      }
                      for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
                        Ntemp_Tag(ItagLag,Iarea,Iage,Isize) = Ntemp2(Isize);
                    }
                  }

                  // D4c: Movement
                  Nmove_Tag.setZero();
              IsMoves = 0;
              for (int Iarea=0; Iarea<Narea; Iarea++)
                for (int Iage=0; Iage<Nage; Iage++)
                {
                  MovePointer = int(dat.MovePnt(Iarea,Iage,Iyear,Istep));
                  if (MovePointer > 0)
                  {
                    IsMoves = 1;
                    IdestArea = dat.MoveSpec(MovePointer,2);
                    for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
                      MoveVec(Isize) = ActMove(MovePointer,Isize);
                    for (int ItagLag=0; ItagLag<=NtagLag; ItagLag++)
                      for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
                      {
                        Nmove_Tag(ItagLag,IdestArea,Iage,Isize) += MoveVec(Isize)*Ntemp_Tag(ItagLag,Iarea,Iage,Isize);
                        Nmove_Tag(ItagLag,Iarea,Iage,Isize) -= MoveVec(Isize)*Ntemp_Tag(ItagLag,Iarea,Iage,Isize);
                      }
                  }
                }

                if (IsMoves==1)
                {
                  for (int ItagLag=0; ItagLag<=NtagLag; ItagLag++)
                    for (int Iarea=0; Iarea<Narea; Iarea++)
                      for (int Iage=0; Iage<Nage; Iage++)
                        for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
                          Ntemp_Tag(ItagLag,Iarea,Iage,Isize) += Nmove_Tag(ItagLag,Iarea,Iage,Isize);
                }

                // D4d: Copy back and update ages
                for (int ItagLag=0; ItagLag<=NtagLag; ItagLag++) {
                  for (int Iarea=0; Iarea<Narea; Iarea++) {
                    if (Istep<dat.Nstep-1) {
                      for (int Iage=0; Iage<Nage; Iage++) {
                        for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++) {
                          Ntag_local(Iarea,ItagLag,Iage,Isize) = Ntemp_Tag(ItagLag,Iarea,Iage,Isize);
                        }
                      }
                    }
                    else
                    {
                      if (Nage-1 > 0)
                      {
                        for (int Iage=0; Iage<Nage-1; Iage++)
                          for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
                            Ntag_local(Iarea,ItagLag,Iage+1,Isize) = Ntemp_Tag(ItagLag,Iarea,Iage,Isize);
                        for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
                          Ntag_local(Iarea,ItagLag,Nage-1,Isize) = Ntemp_Tag(ItagLag,Iarea,Nage-1,Isize) + Ntemp_Tag(ItagLag,Iarea,Nage-2,Isize);
                        for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
                          Ntag_local(Iarea,ItagLag,0,Isize) = 0;
                      }
                      else
                      {
                        for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
                          Ntag_local(Iarea,ItagLag,0,Isize) = Ntemp_Tag(ItagLag,Iarea,0,Isize);
                      }
                    }
                  }
                }
    }

    // E1: Add animals at end of projection to NotReported
    for (int ItagLag=0; ItagLag<=NtagLag; ItagLag++)
      for (int Iarea=0; Iarea<Narea; Iarea++)
        for (int Iage=0; Iage<Nage; Iage++)
          for (int Isize=0; Isize<dat.Nlen(SexPass); Isize++)
            NotReported(SexPass,GrpPass) += Ntag_local(Iarea,ItagLag,Iage,Isize);

  // Total reported and rescale recaptures
  TotalReported = 0;
  for (int Iarea=0; Iarea<Narea; Iarea++)
    for (int Iyear=0; Iyear<thedata.NyearTags; Iyear++)
      for (int Istep=0; Istep<dat.Nstep; Istep++)
        for (int IrepSplit=0; IrepSplit<NrepSplit; IrepSplit++)
        {
          TotalReported += RecapNum(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep);
          RecapNum(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep) /= thedata.NrelTotal(SexPass,GrpPass);
        }
        NotReported(SexPass,GrpPass) /= thedata.NrelTotal(SexPass,GrpPass);

  // DIAGNOSTIC: Check if model predicted ANY recaptures
  //FILE* fp2 = fopen("tagdym_recapture_summary.txt", "w");
  // fprintf(fp2, "=== Recapture Summary for SexPass=%d, GrpPass=%d ===\n\n", SexPass, GrpPass);
  // fclose(fp2);

  // Likelihood
  LikeTag2 = -thedata.NrelTotal(SexPass,GrpPass)*thedata.NotReportedObs(SexPass,GrpPass)*
    log((NotReported(SexPass,GrpPass)+1e-10)/(thedata.NotReportedObs(SexPass,GrpPass)+1e-10));
  for (int Iarea=0; Iarea<Narea; Iarea++)
    for (int Iyear=0; Iyear<thedata.NyearTags; Iyear++)
      for (int Istep=0; Istep<dat.Nstep; Istep++)
        for (int IrepSplit=0; IrepSplit<NrepSplit; IrepSplit++)
          if (thedata.RecapObs(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep) > 0)
          {
            // Add small constant to avoid log(0)
            LikeCompT = thedata.NrelTotal(SexPass,GrpPass)*thedata.RecapObs(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep)*log((RecapNum(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep)+1e-10)/(thedata.RecapObs(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep)+1e-10));
            LikeTag2 -= LikeCompT;
          }
          TagLike2(SexPass,GrpPass) += LikeTag2;

          return(NeglogLikelihood);
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Basic parameters
  dataSet<Type> dataset;
  DATA_INTEGER(Nyear); dataset.Nyear=Nyear;
  DATA_INTEGER(Year1); dataset.First_yr=Year1;
  DATA_INTEGER(MaxProjYr); dataset.MaxProjYr = MaxProjYr;
  DATA_INTEGER(Nproj); dataset.Nproj = Nproj;
  DATA_INTEGER(DoProject); dataset.DoProject = DoProject;
  DATA_INTEGER(ProjType); dataset.ProjType = ProjType;
  DATA_ARRAY(ProjHarvestRate); dataset.ProjHarvestRate = ProjHarvestRate;
  DATA_INTEGER(Nstep); dataset.Nstep=Nstep;
  DATA_INTEGER(Narea); dataset.Narea=Narea;
  DATA_INTEGER(Nage); dataset.Nage=Nage;
  DATA_INTEGER(Nsex); dataset.Nsex=Nsex;
  DATA_INTEGER(Nfleet); dataset.Nfleet=Nfleet;
  DATA_INTEGER(MaxLen); dataset.MaxLen=MaxLen;
  DATA_INTEGER(NselPatterns); dataset.NselPatterns=NselPatterns;
  DATA_INTEGER(NretPatterns); dataset.NretPatterns=NretPatterns;
  DATA_INTEGER(NlegalPatterns); dataset.NlegalPatterns=NlegalPatterns;
  DATA_IVECTOR(Nlen); dataset.Nlen=Nlen;
  DATA_INTEGER(BurnIn); dataset.BurnIn=BurnIn;
  DATA_IVECTOR(BurnInVec); dataset.BurnInVec=BurnInVec;
  DATA_INTEGER(Num_Iteration); dataset.Num_Iteration=Num_Iteration;
  DATA_INTEGER(Tune_Years); dataset.Tune_Years=Tune_Years;
  DATA_IMATRIX(SelSpec); dataset.SelSpec=SelSpec;
  DATA_IMATRIX(RetSpec); dataset.RetSpec=RetSpec;
  DATA_IARRAY(IsRed); dataset.IsRed=IsRed;
  DATA_IMATRIX(LegalSpec); dataset.LegalSpec=LegalSpec;
  DATA_IVECTOR(Fleet_area); dataset.Fleet_area=Fleet_area;
  DATA_IVECTOR(Narea_fleet); dataset.Narea_fleet=Narea_fleet;
  DATA_IMATRIX(Area_fleet); dataset.Area_fleet=Area_fleet;
  DATA_IARRAY(SelPnt); dataset.SelPnt = SelPnt;
  DATA_IARRAY(RetPnt); dataset.RetPnt = RetPnt;
  DATA_IARRAY(LegalFleetPnt); dataset.LegalFleetPnt = LegalFleetPnt;
  DATA_IARRAY(SelPntFut); dataset.SelPntFut = SelPnt;
  DATA_IARRAY(RetPntFut); dataset.RetPntFut = RetPnt;
  DATA_IARRAY(LegalFleetPntFut); dataset.LegalFleetPntFut = LegalFleetPntFut;
  DATA_MATRIX(TimeStepLen); dataset.TimeStepLen=TimeStepLen;
  DATA_ARRAY(Catch); dataset.Catch=Catch;
  DATA_MATRIX(MidLenBin); dataset.MidLenBin=MidLenBin;
  DATA_MATRIX(LowLenBin); dataset.LowLenBin=LowLenBin;
  DATA_ARRAY(Phi); dataset.Phi=Phi;
  DATA_INTEGER(NmovePatterns); dataset.NmovePatterns=NmovePatterns;
  DATA_IMATRIX(MoveSpec); dataset.MoveSpec=MoveSpec;
  DATA_IARRAY(MovePnt); dataset.MovePnt = MovePnt;
  DATA_INTEGER(NrecruitPatternsA); dataset.NrecruitPatternsA = NrecruitPatternsA;
  DATA_IMATRIX(RecruitSpecsA); dataset.RecruitSpecsA = RecruitSpecsA;
  DATA_INTEGER(NrecruitPatternsB); dataset.NrecruitPatternsB = NrecruitPatternsB;
  DATA_IMATRIX(RecruitSpecsB); dataset.RecruitSpecsB = RecruitSpecsB;
  DATA_IMATRIX(RecruitPnt); dataset.RecruitPnt = RecruitPnt;
  DATA_IVECTOR(RecruitLenPnt); dataset.RecruitLenPnt = RecruitLenPnt;
  DATA_MATRIX(RecruitFrac); dataset.RecruitFrac = RecruitFrac;
  DATA_INTEGER(CalcRecruitFrac); dataset.CalcRecruitFrac = CalcRecruitFrac;
  DATA_INTEGER(NgrowthPatterns); dataset.NgrowthPatterns = NgrowthPatterns;
  DATA_IMATRIX(GrowthSpecs); dataset.GrowthSpecs = GrowthSpecs;
  DATA_IARRAY(GrowthPnt); dataset.GrowthPnt = GrowthPnt;
  DATA_ARRAY(TransInp);  dataset.TransInp = TransInp;
  DATA_INTEGER(RecYr1); dataset.RecYr1 = RecYr1;
  DATA_INTEGER(RecYr2); dataset.RecYr2 = RecYr2;
  DATA_INTEGER(RecSpatYr1); dataset.RecSpatYr1 = RecSpatYr1;
  DATA_INTEGER(RecSpatYr2); dataset.RecSpatYr2 = RecSpatYr2;
  DATA_INTEGER(MatTimeStep); dataset.MatTimeStep = MatTimeStep;
  DATA_INTEGER(BioTimeStep); dataset.BioTimeStep = BioTimeStep;
  DATA_ARRAY(MatFem); dataset.MatFem = MatFem;
  DATA_IVECTOR(MatAge); dataset.MatAge = MatAge;
  DATA_IVECTOR(MparsLink); dataset.MparsLink = MparsLink;
  DATA_MATRIX(MparsPrior); dataset.MparsPrior = MparsPrior;
  DATA_IVECTOR(RecparsLink); dataset.RecparsLink = RecparsLink;
  DATA_MATRIX(RecparsPrior); dataset.RecparsPrior = RecparsPrior;
  DATA_IVECTOR(SelparsLink); dataset.SelparsLink = SelparsLink;
  DATA_MATRIX(SelparsPrior); dataset.SelparsPrior = SelparsPrior;
  DATA_IVECTOR(EffparsLink); dataset.EffparsLink = EffparsLink;
  DATA_MATRIX(EffparsPrior); dataset.EffparsPrior = EffparsPrior;
  DATA_IVECTOR(MoveparsLink); dataset.MoveparsLink = MoveparsLink;
  DATA_MATRIX(MoveparsPrior); dataset.MoveparsPrior = MoveparsPrior;
  DATA_SCALAR(Bias_Ramp_Yr1); dataset.Bias_Ramp_Yr1 = Bias_Ramp_Yr1;
  DATA_SCALAR(Bias_Ramp_Yr2); dataset.Bias_Ramp_Yr2 = Bias_Ramp_Yr2;
  DATA_SCALAR(Bias_Ramp_Yr3); dataset.Bias_Ramp_Yr3 = Bias_Ramp_Yr3;
  DATA_SCALAR(Bias_Ramp_Yr4); dataset.Bias_Ramp_Yr4 = Bias_Ramp_Yr4;

  TheData<Type> thedata;
  DATA_INTEGER(Ncpue); thedata.Ncpue = Ncpue;
  DATA_IMATRIX(IndexI); thedata.IndexI = IndexI;
  DATA_MATRIX(IndexR); thedata.IndexR = IndexR;
  DATA_INTEGER(Nnumbers); thedata.Nnumbers = Nnumbers;
  DATA_IMATRIX(NumbersI); thedata.NumbersI = NumbersI;
  DATA_MATRIX(NumbersR); thedata.NumbersR = NumbersR;
  DATA_INTEGER(NlenComp); thedata.NlenComp = NlenComp;
  DATA_IMATRIX(LenCompI); thedata.LenCompI = LenCompI;
  DATA_VECTOR(Stage1W); thedata.Stage1W = Stage1W;
  DATA_MATRIX(LenCompR); thedata.LenCompR = LenCompR;
  DATA_SCALAR(LambdaCpue); thedata.LambdaCpue = LambdaCpue;
  DATA_SCALAR(LambdaNumbers); thedata.LambdaNumbers = LambdaNumbers;
  DATA_SCALAR(LambdaLength); thedata.LambdaLength = LambdaLength;
  DATA_SCALAR(LambdaLarval); thedata.LambdaLarval = LambdaLarval;
  DATA_SCALAR(LambdaTag1); thedata.LambdaTag1 = LambdaTag1;
  DATA_SCALAR(LambdaTag2); thedata.LambdaTag2 = LambdaTag2;
  DATA_VECTOR(LambdaCpue2); thedata.LambdaCpue2 = LambdaCpue2;
  DATA_VECTOR(LambdaNumbers2); thedata.LambdaNumbers2 = LambdaNumbers2;
  // DATA_SCALAR(WeightInitialN); thedata.WeightInitialN = WeightInitialN;
  // DATA_SCALAR(WeightInit3); thedata.WeightInit3 = WeightInit3;
  DATA_ARRAY(LambdaLength2); thedata.LambdaLength2 = LambdaLength2;
  DATA_INTEGER(NcpueDataSeries); thedata.NcpueDataSeries = NcpueDataSeries;
  DATA_IVECTOR(FixSigmaCpue); thedata.FixSigmaCpue = FixSigmaCpue;
  DATA_IVECTOR(IndexType); thedata.IndexType = IndexType;
  DATA_IVECTOR(EnvIndCpue); thedata.EnvIndCpue = EnvIndCpue;
  DATA_IVECTOR(EffCrIndCpue); thedata.EffCrIndCpue = EffCrIndCpue;
  DATA_IVECTOR(EffCrLag); thedata.EffCrLag = EffCrLag;
  DATA_IVECTOR(TreatQcpue); thedata.TreatQcpue = TreatQcpue;
  DATA_SCALAR(SigmaCpueOffset); thedata.SigmaCpueOffset = SigmaCpueOffset;
  DATA_INTEGER(LarvalLikeOpt); thedata.LarvalLikeOpt = LarvalLikeOpt;
  DATA_INTEGER(Larval_Offset); thedata.Larval_Offset = Larval_Offset;
  DATA_INTEGER(NLarvalData); thedata.NLarvalData = NLarvalData;
  DATA_IMATRIX(Lar_dataI); thedata.Lar_dataI = Lar_dataI;
  DATA_MATRIX(Lar_dataR); thedata.Lar_dataR = Lar_dataR;
  DATA_ARRAY(EnvData); thedata.EnvData = EnvData;
  DATA_INTEGER(IsTagData); thedata.IsTagData = IsTagData;
  DATA_INTEGER(NtagGroups); thedata.NtagGroups = NtagGroups;
  DATA_SCALAR(InitialLoss); thedata.InitialLoss = InitialLoss;
  DATA_SCALAR(TagLossRate); thedata.TagLossRate = TagLossRate;
  DATA_INTEGER(NrepSplit); thedata.NrepSplit = NrepSplit;
  DATA_INTEGER(NtagLag); thedata.NtagLag = NtagLag;
  DATA_VECTOR(RepRate); thedata.RepRate = RepRate;
  DATA_IVECTOR(FitTagSizes); thedata.FitTagSizes = FitTagSizes;
  DATA_IVECTOR(Year1Tag); thedata.Year1Tag = Year1Tag;
  DATA_IVECTOR(Year2Tag); thedata.Year2Tag = Year2Tag;
  DATA_INTEGER(TagYr1); thedata.TagYr1 = TagYr1;
  DATA_INTEGER(TagYr2); thedata.TagYr2 = TagYr2;
  DATA_INTEGER(NyearTags); thedata.NyearTags = NyearTags;
  DATA_ARRAY(TagRel); thedata.TagRel = TagRel;
  DATA_ARRAY(TagRec); thedata.TagRec = TagRec;
  DATA_ARRAY(RecapObs); thedata.RecapObs = RecapObs;
  DATA_MATRIX(NrelTotal); thedata.NrelTotal = NrelTotal;
  DATA_MATRIX(NotReportedObs); thedata.NotReportedObs = NotReportedObs;
  DATA_ARRAY(PropRepSplit); thedata.PropRepSplit = PropRepSplit;

  DATA_INTEGER(NcatchDataSeries); thedata.NcatchDataSeries = NcatchDataSeries;
  DATA_IVECTOR(FixSigmaCatchN); thedata.FixSigmaCatchN = FixSigmaCatchN;
  DATA_SCALAR(SigmaCatchNOffset); thedata.SigmaCatchNOffset = SigmaCatchNOffset;

  DATA_MATRIX(WeightLen);
  DATA_MATRIX(SelexFI);
  DATA_MATRIX(RetenFI);
  DATA_MATRIX(LegalFI);

  DATA_INTEGER(Nzone);
  DATA_IVECTOR(NareasPerZone);
  DATA_IMATRIX(AreasPerZone);
  DATA_IARRAY(LegalPnt);
  DATA_MATRIX(LegalRef);

  DATA_INTEGER(NvarTypes);
  DATA_IVECTOR(VarTypes);

  // Estimated parameters
  PARAMETER_VECTOR(MainPars);
  PARAMETER_VECTOR(RecruitPars);
  PARAMETER_VECTOR(PuerPowPars);
  PARAMETER_VECTOR(SelPars);
  PARAMETER_VECTOR(RetPars);
  PARAMETER_VECTOR(RecDevs);
  PARAMETER_VECTOR(Qpars);
  PARAMETER_VECTOR(efpars);
  //PARAMETER_VECTOR(InitPars);
  PARAMETER_VECTOR(RecSpatDevs);
  PARAMETER_VECTOR(MovePars);
  PARAMETER_VECTOR(GrowthPars);
  PARAMETER(dummy);

  matrix <Type> M(Narea, Nage);
  Type Rbar;
  Type MWhitesPar;
  Type QRedsPar;
  Type SigmaR;
  int Link;


  ////// Deal with mainpars ///////
  // Apply priors on Main Pars if requested
  Type MainParPriorPen = 0;
  Type ScaleMP = 0;
  Type ShapeMP = 0;
  int nrowMP = MparsPrior.rows();
  for (int r=0; r<nrowMP; r++) {
    // Normal prior
    if(MparsPrior(r,0)==1){
      MainParPriorPen += -dnorm(MainPars(r,0), MparsPrior(r,1), MparsPrior(r,2), true);
    }
    // Gamma prior
    if(MparsPrior(r,0)==2){
      ScaleMP = square(MparsPrior(r,2))/MparsPrior(r,1);
      ShapeMP = MparsPrior(r,1)/ScaleMP;
      MainParPriorPen += -dgamma(MainPars(r,0), ShapeMP, ScaleMP, true);
    }
    // Log-normal prior
    if(MparsPrior(r,0)==3){
      Type mulog  = log(MparsPrior(r,1)) - Type(0.5) * log(Type(1.0) + square(MparsPrior(r,2)/MparsPrior(r,1)));
      Type sdlog  = sqrt(log(Type(1.0) + square(MparsPrior(r,2)/MparsPrior(r,1))));
      MainParPriorPen += -dnorm(log(MainPars(r,0)), mulog, sdlog, true) + log(MainPars(r,0));
    }
  }

  // Adjust main parameters to account for linked parameters and linked offset
  for(int mp=0;mp<MparsLink.size();mp++){
    if(MparsLink(mp)>0) MainPars(mp)=MainPars(MparsLink(mp)-1);
    if(MparsLink(mp)<0){                     /// change link value to +ive and then add that parameter to the current parameter
      Link = -1 * MparsLink(mp);
      MainPars(mp) += MainPars(Link-1);
    }
  }

  //// Deal with Recruitment Pars  ///////
  // Apply priors on Rec Pars if requested
  Type RecParPriorPen = 0;
  int nrowRP = RecparsPrior.rows();
  for (int r=0; r<nrowRP; r++) {
    // Normal prior
    if(RecparsPrior(r,0)==1){
      RecParPriorPen += -dnorm(RecruitPars(r,0), RecparsPrior(r,1), RecparsPrior(r,2), true);
    }
    // Gamma prior
    if(RecparsPrior(r,0)==2){
      ScaleMP = square(RecparsPrior(r,2))/RecparsPrior(r,1);
      ShapeMP = RecparsPrior(r,1)/ScaleMP;
      RecParPriorPen += -dgamma(RecruitPars(r,0), ShapeMP, ScaleMP, true);
    }
    // Log-normal prior
    if(RecparsPrior(r,0)==3){
      Type mulog  = log(RecparsPrior(r,1)) - Type(0.5) * log(Type(1.0) + square(RecparsPrior(r,2)/RecparsPrior(r,1)));
      Type sdlog  = sqrt(log(Type(1.0) + square(RecparsPrior(r,2)/RecparsPrior(r,1))));
      RecParPriorPen += -dnorm(log(RecruitPars(r,0)), mulog, sdlog, true) + log(RecruitPars(r,0));
    }
  }

  // Adjust parameters to account for linked parameters
  for(int mp=0;mp<RecparsLink.size();mp++){
    if(RecparsLink(mp)>0) RecruitPars(mp)=RecruitPars(RecparsLink(mp)-1);
    if(RecparsLink(mp)<0){                     /// change link value to +ive and then add that parameter to the current parameter
      Link = -1 * RecparsLink(mp);
      RecruitPars(mp) += RecruitPars(Link-1);
    }
  }

  //// Deal with Selectivity Pars  ///////
  // Apply priors on Select Pars if requested
  Type SelParPriorPen = 0;
  int nrowSP = SelparsPrior.rows();
  for (int r=0; r<nrowSP; r++) {
    // Normal prior
    if(SelparsPrior(r,0)==1){
      SelParPriorPen += -dnorm(SelPars(r,0), SelparsPrior(r,1), SelparsPrior(r,2), true);
    }
    // Gamma prior
    if(SelparsPrior(r,0)==2){
      ScaleMP = square(SelparsPrior(r,2))/SelparsPrior(r,1);
      ShapeMP = SelparsPrior(r,1)/ScaleMP;
      SelParPriorPen += -dgamma(SelPars(r,0), ShapeMP, ScaleMP, true);
    }
    // Log-normal prior
    if(SelparsPrior(r,0)==3){
      Type mulog  = log(SelparsPrior(r,1)) - Type(0.5) * log(Type(1.0) + square(SelparsPrior(r,2)/SelparsPrior(r,1)));
      Type sdlog  = sqrt(log(Type(1.0) + square(SelparsPrior(r,2)/SelparsPrior(r,1))));
      SelParPriorPen += -dnorm(log(SelPars(r,0)), mulog, sdlog, true) + log(SelPars(r,0));
    }
  }

  // Adjust  parameters to account for linked parameters
  for(int mp=0;mp<SelparsLink.size();mp++){
    if(SelparsLink(mp)>0) SelPars(mp)=SelPars(SelparsLink(mp)-1);
    if(SelparsLink(mp)<0){                     /// change link value to +ive and then add that parameter to the current parameter
      Link = -1 * SelparsLink(mp);
      SelPars(mp) += SelPars(Link-1);
    }
  }

    //// Deal with Efficiency Pars  ///////
  // Apply priors on Rec Pars if requested
  Type EffParPriorPen = 0;
  int nrowEP = EffparsPrior.rows();
  for (int r=0; r<nrowEP; r++) {
    // Normal prior
    if(EffparsPrior(r,0)==1){
      EffParPriorPen += -dnorm(efpars(r,0), EffparsPrior(r,1), EffparsPrior(r,2), true);
    }
    // Gamma prior
    if(EffparsPrior(r,0)==2){
      ScaleMP = square(EffparsPrior(r,2))/EffparsPrior(r,1);
      ShapeMP = EffparsPrior(r,1)/ScaleMP;
      EffParPriorPen += -dgamma(efpars(r,0), ShapeMP, ScaleMP, true);
    }
    // Log-normal prior
    if(EffparsPrior(r,0)==3){
      Type mulog  = log(EffparsPrior(r,1)) - Type(0.5) * log(Type(1.0) + square(EffparsPrior(r,2)/EffparsPrior(r,1)));
      Type sdlog  = sqrt(log(Type(1.0) + square(EffparsPrior(r,2)/EffparsPrior(r,1))));
      EffParPriorPen += -dnorm(log(efpars(r,0)), mulog, sdlog, true) + log(efpars(r,0));
    }
    // Soft lower bound — applied once per parameter alongside the prior
    EffParPriorPen += exp(-10.0 * efpars(r,0));
  }

  // Adjust parameters to account for linked parameters
  for(int mp=0;mp<EffparsLink.size();mp++){
    if(EffparsLink(mp)>0) efpars(mp)=efpars(EffparsLink(mp)-1);
    if(EffparsLink(mp)<0){                     /// change link value to +ive and then add that parameter to the current parameter
      Link = -1 * EffparsLink(mp);
      efpars(mp) += efpars(Link-1);
    }
  }

  //// Deal with Migrate Pars  ///////
  // Apply priors on Rec Pars if requested
  Type MoveParPriorPen = 0;
  nrowMP = MoveparsPrior.rows();
  for (int r=0; r<nrowMP; r++) {
    // Normal prior
    if(MoveparsPrior(r,0)==1){
      MoveParPriorPen += -dnorm(MovePars(r,0), MoveparsPrior(r,1), MoveparsPrior(r,2), true);
    }
    // Gamma prior
    if(MoveparsPrior(r,0)==2){
      ScaleMP = square(MoveparsPrior(r,2))/MoveparsPrior(r,1);
      ShapeMP = MoveparsPrior(r,1)/ScaleMP;
      MoveParPriorPen += -dgamma(MovePars(r,0), ShapeMP, ScaleMP, true);
    }
    // Log-normal prior
    if(MoveparsPrior(r,0)==3){
      Type mulog  = log(MoveparsPrior(r,1)) - Type(0.5) * log(Type(1.0) + square(MoveparsPrior(r,2)/MoveparsPrior(r,1)));
      Type sdlog  = sqrt(log(Type(1.0) + square(MoveparsPrior(r,2)/MoveparsPrior(r,1))));
      MoveParPriorPen += -dnorm(log(MovePars(r,0)), mulog, sdlog, true) + log(MovePars(r,0));
    }
  }


  // Adjust parameters to account for linked parameters
  for(int mp=0;mp<MoveparsLink.size();mp++){
    if(MoveparsLink(mp)>0) MovePars(mp)=MovePars(MoveparsLink(mp)-1);
    if(MoveparsLink(mp)<0){                     /// change link value to +ive and then add that parameter to the current parameter
      Link = -1 * MoveparsLink(mp);
      MovePars(mp) += MovePars(Link-1);
    }
  }

    // Split MainPars out into their various groups
  Rbar = MainPars(0);
  for(int Iarea=0;Iarea<Narea;Iarea++){
    for (int Iage=0;Iage<Nage;Iage++){
      M(Iarea,Iage) =   MainPars(1+Iarea) * MainPars(1+Narea+Iage); }}
  MWhitesPar = MainPars(1+Narea+Nage);
  QRedsPar = MainPars(2+Narea+Nage);
  SigmaR = MainPars(3+Narea+Nage);

  //for (int Iarea=0;Iarea<Narea;Iarea++) LogRinitial(Iarea) = MainPars(4+Narea+Nage+Iarea);
  //Finitial = exp(MainPars(4+2*Narea+Nage));

  //Pull recruitment fractions out of RecruitPars - but can leave them in original RecruitPars
  if(CalcRecruitFrac==1){
    int Nfrac = RecruitFrac.rows();                   // How many Recruitment fractions are needed
    int jumpoff = RecruitPars.size()-(Nfrac*2);       // How many RecruitPars there are minus those for recruit Fraction
    vector <Type> RecFracM(Nfrac);                    // Make a vector to store mean
    vector <Type> RecFracSd(Nfrac);                   // Make a vector to store SD
    for(int Ifrac=0;Ifrac<Nfrac;Ifrac++){
      RecFracM(Ifrac) = RecruitPars((Ifrac*2)+jumpoff);
      RecFracSd(Ifrac) = RecruitPars((Ifrac*2)+1+jumpoff);
    }

  // Strip-out and replace the original RecruitPars - this may not be necessary but keeps things cleaner
    vector <Type> RecruitPars2(jumpoff);
    for(int Ifrac=0;Ifrac<jumpoff;Ifrac++){
      RecruitPars2(Ifrac) = RecruitPars(Ifrac);
    }
    RecruitPars = RecruitPars2;

    // Set up recruitment fractions if CalcRecruitFrac==1
    Type len1, len2;
    for(int Ifrac=0; Ifrac<Nfrac; ++Ifrac){
      for(int Ilen=0; Ilen<Nlen(0); ++Ilen){
        len1 =  LowLenBin(0,Ilen+1);
        if(Ilen==0) {len2 =  0 ;} else {len2 =  LowLenBin(0,Ilen);}  // Draw all lobster < min lbin into first lbin
        RecruitFrac(Ifrac,Ilen) = pnorm(len1, RecFracM(Ifrac), RecFracSd(Ifrac))-pnorm(len2, RecFracM(Ifrac), RecFracSd(Ifrac));
       }
     }
   dataset.RecruitFrac = RecruitFrac;
   }


  // Local variables
  Type        neglogL;                                                                     // Negative log likelihood

  array<Type> N(Narea, BurnIn+Nyear+MaxProjYr+1, Nstep, Nsex, Nage, MaxLen); N.setZero();  // N matrix
  array<Type> Z(Narea, BurnIn+Nyear+MaxProjYr+1, Nstep, Nsex, Nage, MaxLen); Z.setZero();  // Z matrix
  vector<Type> Recruits(BurnIn+Nyear+MaxProjYr+1);                                         // Recruitment output
  vector<Type> BiasMult(BurnIn+Nyear+MaxProjYr+1);                                         // Bias correction factor
  array<Type>  Ninit(Narea,Nsex,Nage,MaxLen);
  vector<Type> MatBio(BurnIn+Nyear+MaxProjYr+1);
  matrix<Type> MatBioArea(Narea,BurnIn+Nyear+MaxProjYr+1);
  matrix<Type> RecruitmentByArea(Narea,BurnIn+Nyear+MaxProjYr+1);                          // Recruitment
  matrix<Type> PuerulusByArea(Narea,BurnIn+Nyear+MaxProjYr+1);                             // Puerulus

  array<Type> Hrate(BurnIn+Nyear+MaxProjYr,Nstep,Nfleet);                                  // Harvest rate
  matrix<Type> LegalBio(Nyear,Narea);                                                      // Legal biomass
  matrix<Type> LegalBioAll(BurnIn+Nyear+MaxProjYr,Narea);                                  // Legal biomass
  array<Type> LegalBioAllbySex(BurnIn+Nyear+MaxProjYr,Narea,Nsex);                         // Legal biomass by sex T step 1
  array<Type> MatureBioAllbySex(BurnIn+Nyear+MaxProjYr,Narea,Nsex);                         // Legal biomass by sex T step 1
  //matrix<Type> LegalBio76(Nyear,Narea);                                                    // Legal biomass of all lobster > 76 mm
  array<Type> CumCatch(Nyear,Narea,Nstep);
  //array<Type> CatchYSA(Nyear,Nstep,Narea);
  matrix<Type> CatchYA(Nyear,Narea);
  array<Type> HRint(Nyear,Nstep,Narea);
  //array<Type> wHRint(Nyear,Nstep,Narea);
  //array<Type> LegalBioTS(Nyear,Narea,Nstep);                                               // Legal biomass
  //array<Type> LegalBio76TS(Nyear,Narea,Nstep);                                             // Legal biomass of all lobster > 76 mm
  //matrix<Type> sLegalBio(Nyear,Narea);                                                     // Legal biomass
  //matrix<Type> sLegalBio76(Nyear,Narea);                                                   // Legal biomass of all lobster > 76 mm
  matrix<Type> HarvestRate(Nyear,Nzone);                                                   // Harvest rate by year and zone from Hrate
 // matrix<Type> SHarvestRate(Nyear,Nzone);                                                   // Harvest rate by year and zone from Hrate
  //matrix<Type> HarvestRateArea(Nyear,Narea);                                               // Harvest rate by year and area from Lbio
  //matrix<Type> HarvestRateZn(Nyear,Nzone);                                                 // Harvest rate by year and zone from Lbio
  matrix<Type> HarvestRate76(Nyear,Nzone);                                                 // Harvest rate by year and zone of all lobster > 76 mm
 // matrix<Type> SHarvestRate76(Nyear,Nzone);                                                 // Harvest rate by year and zone of all lobster > 76 mm
  matrix<Type> HrateYA(Nyear,Narea);                                                       // Store summed HR by year and area
  array<Type> CatchCheck(Nyear+MaxProjYr,Nstep,Nfleet);                                              // Check
  array<Type> DiscardWt(Nyear+MaxProjYr,Nstep,Nfleet);       DiscardWt.setZero();
  array<Type> DeadDiscardWt(Nyear+MaxProjYr,Nstep,Nfleet);    DeadDiscardWt.setZero();
  matrix<Type> ActSelex(NselPatterns,MaxLen);
  matrix<Type> ActReten(NretPatterns,MaxLen);
  matrix<Type> ActLegal(NlegalPatterns,MaxLen);
  matrix<Type> ActMove(NmovePatterns,MaxLen);
  array<Type> ActRecruitAreaSexDist(Nyear+MaxProjYr,Nstep,Narea,Nsex);                     // Allocation on recruitment to areas and sexes
  array<Type> ActRecruitLenDist(NrecruitPatternsB,Nsex,MaxLen);
  array<Type> ActGrowth(NgrowthPatterns,MaxLen,MaxLen);
  vector<Type> ActRecDev(BurnIn+Nyear+MaxProjYr+1);                                        // Recruitment deviations
  vector<Type> VirginBio(Narea);                                                           // Virgin biomass used to produce M
  vector<Type> VirginLegalBio(Narea);
  array<Type> VirginNvec(Narea, Nsex, Nage, MaxLen);        // numbers by area, sex, age, size
  array<Type> VirginBioAtLen(Narea, Nsex, MaxLen);           // biomass summed over age by area, sex, size
 // vector<Type> AvM(BurnIn+Nyear+MaxProjYr+1);                                                           // Average M each year
  vector<Type> CurrentBio(Narea);                                                         // Current biomass used to produce M

  for (int Iyear=-BurnIn;Iyear<Nyear+Nproj+1;Iyear++)
   {
    if (Iyear<Bias_Ramp_Yr1)
     BiasMult(BurnIn+Iyear) = 0;
    else
     if (Iyear<Bias_Ramp_Yr2)
      BiasMult(BurnIn+Iyear) = (Iyear-Bias_Ramp_Yr1)/(Bias_Ramp_Yr2-Bias_Ramp_Yr1);
     else
      if (Iyear<Bias_Ramp_Yr3)
       BiasMult(BurnIn+Iyear) = 1.0;
      else
       if (Iyear<Bias_Ramp_Yr4)
        BiasMult(BurnIn+Iyear) = (Bias_Ramp_Yr4-Iyear)/(Bias_Ramp_Yr4-Bias_Ramp_Yr3);
       else
        BiasMult(BurnIn+Iyear) = 0;
    }

  ActRecDev.setZero();
  for (int Iyear=RecYr1;Iyear<=RecYr2;Iyear++)
   ActRecDev(Iyear) = RecDevs(Iyear-RecYr1);
  for (int Iyear=0;Iyear<BurnIn+Nyear+Nproj+1;Iyear++)
   Recruits(Iyear) = exp(Rbar)*exp(ActRecDev(Iyear)-BiasMult(Iyear)*SigmaR*SigmaR/2.0);

  array<Type> selexF(Nfleet,Nsex,Nage,MaxLen);                            // Selectivity
  array<Type> retainF(Nfleet,Nsex,Nage,MaxLen);                           // Retention
  array<Type> selretwght(Nfleet,Nsex,Nage,MaxLen);                        // Product of selectivity,retention and weight
  Type Z2;                                                                // temp variable

  matrix<Type> PredCpue(Ncpue,2);                                         // Predicted CPUE and residuals
  matrix<Type> PredNumbers(Nnumbers,2);                                   // Predicted catch-in-numbers and residuals
  matrix<Type> PredLengthComp(NlenComp,MaxLen);                           // Predicted catch-at-length
  matrix<Type> PredLarval(NLarvalData,2);                                 // Predicted larval data and residuals
  vector<Type> CpueLikeComps(NcpueDataSeries);                            // Cpue likelihood by fleet
  vector<Type> SigmaCpue(NcpueDataSeries);                                // Sigmas
  vector<Type> CpueQ(NcpueDataSeries);                                    // Catchability
  matrix<Type> CpueEcreep(Nyear+MaxProjYr+1,EffCrLag.size());
  vector<Type> NumbersLikeComps(NcatchDataSeries);                        // Numbers likelihood by series
  vector<Type> SigmaNumbers(NcatchDataSeries);                            // Sigmas
  vector<Type> LengthLikeComps(Nfleet);                                   // Length likelihood by fleet
  vector<Type> LarvalLikeComps(Narea);                                    // Length likelihood by fleet
  matrix<Type> TagLike1(Nsex,NtagGroups);                                 // Tag size likelihood
  matrix<Type> TagLike2(Nsex,NtagGroups);                                 // Tag numbers likelihood

  //array<Type> Ntag(Nsex,NtagGroups,Narea,NtagLag+1,Nage,MaxLen);
  array<Type> RecapNum(Nsex,NtagGroups,Narea,NrepSplit,NyearTags,Nstep);  // Recapture
  matrix<Type> NotReported(Nsex,NtagGroups);
  array<Type> PredTagSize(Nsex,NtagGroups,Narea,MaxLen);                  // Recaptured by length-class (weight by numbers recaptured)

  int IsVirgin;                                                           // Set to 1 for unfished state
  matrix<Type> Feqn2(Nfleet,Nstep); Feqn2.setZero();                      // Initial F (not used in projections)

  vector<Type> XX(2);
  Type Test2;
  Type CatchLike;
  Type CpueLike;
  Type NumbersLike;
  Type LengthLike;
  Type LarvalLike;
  Type Weighted_CpueLike;
  Type Weighted_NumbersLike;
  Type Weighted_LengthLike;
  Type Weighted_LarvalLike;
  Type Weighted_TagLike1;
  Type Weighted_TagLike2;
  Type Rec_Penal;
  Type Rec_Penal_Smooth;
  Type Rec_Penal_SumZero;
  Type Initial_pen;

  int Ipnt;                                                               // Pointer

  // Set up the selectivity vectors that will be used
  ActSelex = SetUpSelex(dataset, SelPars, SelexFI, SelSpec, NselPatterns);
  ActReten = SetUpSelex(dataset, RetPars, RetenFI, RetSpec, NretPatterns);
  ActLegal = SetUpLegal(dataset, LegalFI, LegalSpec, NlegalPatterns);
  ActMove = SetUpMove(dataset, MovePars);
  Test2 = SetUpRecruit(dataset, RecruitPars, RecSpatDevs, ActRecruitAreaSexDist, ActRecruitLenDist );
  ActGrowth = SetUpGrow(dataset, GrowthPars);

  // Recruitment
  int RecruitPointer; int YearAdjust; Type TotalRec;
  for (int Iyear=-BurnIn;Iyear<Nyear+Nproj+1;Iyear++)
    for (int Istep=0;Istep<Nstep;Istep++)
     {
      // Adjusted year (YearAdjust1 is for quantities that go beyond Nyear-1 and YearAdjust2 is not.
      if (Iyear <= 0) { YearAdjust = 0; } else { YearAdjust = Iyear; }
      RecruitPointer = RecruitPnt(YearAdjust,Istep);
      if (RecruitPointer >= 0)
       {
        for (int Iarea=0;Iarea<Narea;Iarea++)
         {
          RecruitmentByArea(Iarea,BurnIn+Iyear) = 0;
          for (int Isex=0;Isex<Nsex;Isex++)
           {
            TotalRec = ActRecruitAreaSexDist(YearAdjust,Istep,Iarea,Isex)*exp(Rbar)*exp(ActRecDev(BurnIn+Iyear))*exp(-BiasMult(BurnIn+Iyear)*SigmaR*SigmaR/2.0);
            RecruitmentByArea(Iarea,BurnIn+Iyear) += TotalRec;
           }
         }
       }
     }

  // Reset
  Hrate.setZero();  N.setZero(); Z.setZero(); MatBio.setZero(); MatBioArea.setZero();


  // Set up initial state (traditional)
  //if (InitOpt==0)
  // {
    Initial_pen = InitializeN(dataset, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, ActGrowth, RecruitFrac,
         Rbar, ActRecruitAreaSexDist, ActRecruitLenDist, ActRecDev, Ninit,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,
         VirginBio, VirginLegalBio, VirginNvec, VirginBioAtLen, LegalRef, CurrentBio);
  // }

  // Set up initial state (alternative)
//   Type NtotalCheck; Type Nexpected;
//   if (InitOpt==1)
//    {
//      Ninit = VirginN(dataset, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, ActGrowth, RecruitFrac,
//          Rbar, ActRecruitAreaSexDist, ActRecruitLenDist, ActRecDev, QRedsPar, MWhitesPar,Type(0.0));
//     Ipnt = 0; Initial_pen = 0;
//     for (int Isex=0;Isex<Nsex;Isex++)
//      {
// 	  NtotalCheck = 0; Nexpected = 0;
//       for (int Iarea=0;Iarea<Narea;Iarea++)
//        for (int Iage=0;Iage<Nage;Iage++)
//         for (int Isize=0;Isize<Nlen(Isex);Isize++)
//          {  N(Iarea,BurnIn,0,Isex,Iage,Isize) = exp(Rbar)*exp(InitPars(Ipnt))  ;
//             Initial_pen += 0.01*InitPars(Ipnt)*InitPars(Ipnt);
//             NtotalCheck += N(Iarea,BurnIn,0,Isex,Iage,Isize)*WeightLen(Isex,Isize);
//             Ipnt += 1; }
//
//  	   for (int Iarea=0;Iarea<Narea;Iarea++)
// 	    for (int Iage=0;Iage<Nage;Iage++)
// 	     for (int Isize=0;Isize<Nlen(Isex);Isize++)
// 	      Nexpected += Ninit(Iarea,Isex,Iage,Isize)*WeightLen(Isex,Isize);
//
//       //for (int Iage=0;Iage<=1000;Iage++) Nexpected += 0.5*exp(Rbar)*exp(-1*float(Iage)*M(0,0));
//       //Initial_pen += (NtotalCheck-Nexpected)*(NtotalCheck-Nexpected);
//       Initial_pen += WeightInitialN*(log(NtotalCheck)-log(Nexpected))*(log(NtotalCheck)-log(Nexpected));
//      }
//    }

  // if(InitOpt==2)
  //  {
  //   Ipnt = 0; Initial_pen = 0;
  //   for (int Iarea=0;Iarea<Narea;Iarea++)
  //    for (int Isex=0;Isex<Nsex;Isex++)
  //     for (int Isize=0;Isize<Nlen(Isex);Isize++)
  //      {   N(Iarea,BurnIn-Nage,0,Isex,Nage-1,Isize) = exp(Rbar)*exp(InitPars(Ipnt))/float(Nsex)/float(Narea)/float(Nlen(Isex));
  //          if (Isize!=0) Initial_pen += 1.0*(InitPars(Ipnt)-InitPars(Ipnt-1))*(InitPars(Ipnt)-InitPars(Ipnt-1));
  //          //Initial_pen += InitPars(Ipnt)*InitPars(Ipnt);
  //          Ipnt += 1;  }
  //   IsVirgin = 0;
  //   for (int Iyear=-Nage;Iyear<0;Iyear++)
  //    for (int Istep=0;Istep<Nstep;Istep++)
  //     {
  //      XX = OneTimeStep(dataset, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, Iyear, Istep, ActGrowth, RecruitFrac, Rbar,IsVirgin, Feqn2, ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,VirginBio, CurrentBio);
  //     } // year and season
  //   }

  // Set up initial state (alternative)
  // if (InitOpt==3)
  //  {
  //   Ipnt = 0; Initial_pen = 0;
  //   for (int Isex=0;Isex<Nsex;Isex++)
  //    for (int Iarea=0;Iarea<Narea;Iarea++)
  //     for (int Iage=0;Iage<Nage;Iage++)
  //      for (int Isize=0;Isize<Nlen(Isex);Isize++)
  //       {
  //       N(Iarea,BurnIn,0,Isex,Iage,Isize) = exp(LogRinitial(0))*exp(InitPars(Ipnt))  ;
  //        Initial_pen += WeightInit3*InitPars(Ipnt)*InitPars(Ipnt);
  //        // New weak penalty
  //        if (Isize!=0) Initial_pen += 1.0*(InitPars(Ipnt)-InitPars(Ipnt-1))*(InitPars(Ipnt)-InitPars(Ipnt-1));
  //        Ipnt += 1;
  //       }
  //   }

  // Set up initial state (virgin)
  // if (InitOpt==4)
  //  {
  //   Initial_pen = 0;
  //   Ninit = VirginN(dataset, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, ActGrowth, RecruitFrac,
  //           Rbar, ActRecruitAreaSexDist, ActRecruitLenDist, ActRecDev, QRedsPar, MWhitesPar,Finitial);
  //   for (int Isex=0;Isex<Nsex;Isex++)
  //    for (int Iarea=0;Iarea<Narea;Iarea++)
  //     for (int Iage=0;Iage<Nage;Iage++)
  //      for (int Isize=0;Isize<Nlen(Isex);Isize++)
  //       N(Iarea,BurnIn,0,Isex,Iage,Isize) = Ninit(Iarea,Isex,Iage,Isize)*exp(LogRinitial(Iarea))/exp(Rbar);
  //    }


   // Set up initial state (alternative)
  // if (InitOpt==5)
  //  {
  //   Ipnt = 0; Initial_pen = 0;
  //   for (int Isex=0;Isex<Nsex;Isex++)
  //    for (int Iarea=0;Iarea<Narea;Iarea++)
  //     for (int Iage=0;Iage<Nage;Iage++)
  //      for (int Isize=0;Isize<Nlen(Isex);Isize++)
  //       {
  //        N(Iarea,BurnIn,0,Isex,Iage,Isize) = exp(LogRinitial(Iarea))*exp(InitPars(Ipnt))  ;
  //        Initial_pen += WeightInit3*InitPars(Ipnt)*InitPars(Ipnt);
  //        // New weak penalty
  //        if (Isize!=0) Initial_pen += 1.0*(InitPars(Ipnt)-InitPars(Ipnt-1))*(InitPars(Ipnt)-InitPars(Ipnt-1));
  //        Ipnt += 1;
  //       }
  //   }

// Project the model forward
IsVirgin = 0;                                                                        // Need to compute Fs
for (int Iyear=0;Iyear<Nyear;Iyear++)
 for (int Istep=0;Istep<Nstep;Istep++)
  {
   XX = OneTimeStep(dataset, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, Iyear, Istep, ActGrowth, RecruitFrac, Rbar,
                             IsVirgin, Feqn2, ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,
                             QRedsPar,MWhitesPar,VirginBio, CurrentBio);
  } // year and season

// Make growth curves for diagnostics
array<Type> GrowthOut(Nyear,Narea,Nsex,Nage,MaxLen);
vector<Type> Lentemp(MaxLen);
vector<Type> Lentemp2(MaxLen);
int GrowthPointer;
for (int Iyear=0;Iyear<Nyear-1;Iyear++) {
    for (int Iarea=0;Iarea<Narea;Iarea++) {
      for (int Isex=0;Isex<Nsex;Isex++) {
        for (int Iage=0;Iage<Nage;Iage++) {
          if(Iyear==0) GrowthOut(Iyear,Iarea,Isex,Iage,0) = 1;                     // seed with a lobster
          for (int Isize=0;Isize<Nlen(Isex);Isize++)  Lentemp(Isize) = GrowthOut(Iyear,Iarea,Isex,Iage,Isize);  //grab current size com.
            for (int Istep=0;Istep<Nstep;Istep++) {
              GrowthPointer = GrowthPnt(Iarea,Isex,Iage,Iyear,Istep);
              if (GrowthPointer >=0) {
            // Key issue (pointer to growth matrix)
            Lentemp2.setZero();
            for (int Isize=0;Isize<Nlen(Isex);Isize++)  {
              for (int Jsize=0;Jsize<=Isize;Jsize++) { Lentemp2(Isize) += Lentemp(Jsize)*ActGrowth(GrowthPointer,Isize,Jsize);
            }}
            for (int Isize=0;Isize<Nlen(Isex);Isize++) Lentemp(Isize) = Lentemp2(Isize);
          }
        } // Step
            for (int Isize=0;Isize<Nlen(Isex);Isize++) GrowthOut(Iyear+1,Iarea,Isex,Iage,Isize) = Lentemp2(Isize);  // Put in the next year
      } // Age
    } // Sex
   } // Area Ntemp2
  } // Year


// Tagging data
if(thedata.IsTagData==1){
  RecapNum.setZero(); NotReported.setZero(); PredTagSize.setZero();  TagLike1.setZero();  TagLike2.setZero();

  for (int SexPass=0; SexPass<Nsex; SexPass++)
    for (int GrpPass=0; GrpPass<NtagGroups; GrpPass++)
      XX = TagDym(dataset,thedata, SexPass, GrpPass, N, ActSelex, ActReten, ActLegal, ActGrowth, ActMove, M, Hrate, QRedsPar, MWhitesPar, RecapNum, NotReported, TagLike1, TagLike2, PredTagSize);
}

// Legal Biomass by Year, Area, time step at the end of a time step
// int Ipoint;int Ipoint76;
// LegalBioTS.setZero();LegalBio76TS.setZero();
// for (int Iyear=0;Iyear<Nyear;Iyear++) {
//   for (int Iarea=0;Iarea<Narea;Iarea++) {
//     for (int Istep=0;Istep<Nstep;Istep++) {
//       for (int Isex=0;Isex<Nsex;Isex++) {
//         for (int Iage=0;Iage<Nage;Iage++) {
//           for (int Ilen=0;Ilen<Nlen(Isex);Ilen++){
//             Ipoint = LegalPnt(Isex,Iage,Iarea,Iyear,Istep);
//             LegalBioTS(Iyear,Iarea,Istep) += LegalFI(Ipoint,Ilen)*N(Iarea,BurnIn+Iyear,Istep,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);
//             LegalBio76TS(Iyear,Iarea,Istep) += LegalRef(Isex,Ilen)*N(Iarea,BurnIn+Iyear,Istep,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);
//           }}}}}}
//
//   // Average Legal Biomass by Year and Area - average by the length of the time step
//   LegalBio.setZero();LegalBio76.setZero();
//   for (int Iyear=0;Iyear<Nyear;Iyear++) {
//     for (int Iarea=0;Iarea<Narea;Iarea++) {
//       for (int Istep=0;Istep<Nstep;Istep++) {
//         //LegalBio(Iyear,Iarea)   += (LegalBioTS(Iyear,Iarea,Istep) * TimeStepLen(Iyear,Istep));
//         LegalBio76(Iyear,Iarea) += (LegalBio76TS(Iyear,Iarea,Istep) * TimeStepLen(Iyear,Istep));
//           }}}


   // Legal Biomass, Discards: computed after the projection loop below (once
   // N is populated for the projection years too) -- see "Post-projection
   // summaries" further down.

   // Simon's Cumulative catch reduced by average M based on time caught
  CumCatch.setZero();
  Type avM = M.sum()/(float(Nage)*float(Narea));    // Calculate Av M
  Type CnT = 0;

  int Iarea;
  for (int Iyear=0;Iyear<Nyear;Iyear++){
    for (int Istep=0;Istep<Nstep;Istep++)   {
      for (int Ifleet=0;Ifleet<Nfleet;Ifleet++)  {
        Iarea = Fleet_area(Ifleet);
        CumCatch(Iyear,Iarea,Istep) += Catch(Iyear,Istep,Ifleet) ;                    // Change recording of catch from fleet to area
      }
      if(Istep>0) { CumCatch(Iyear,Iarea,Istep) += CumCatch(Iyear,Iarea,Istep-1);}    // After all fleets then add previous time-steps cum-catch
      CumCatch(Iyear,Iarea,Istep) *= exp(-TimeStepLen(Iyear,Istep)*avM);              // Reduce the cum_catch by M so it is as it would have been in nature by end of time step
    }}

 // Simons Legal Biomass Legal bio + cumulative catch then averaged over the season weighted by time step length
  // sLegalBio.setZero();
  // sLegalBio76.setZero();
  // for (int Iyear=0;Iyear<Nyear;Iyear++){
  //   for (int Iarea=0;Iarea<Narea;Iarea++)  {
  //     for (int Istep=0;Istep<Nstep;Istep++)   {
  //       sLegalBio(Iyear,Iarea) += (LegalBioTS(Iyear,Iarea,Istep) + CumCatch(Iyear,Iarea,Istep)) * TimeStepLen(Iyear,Istep);           // Add Legal biomass with the cumulative catch to date
  //       sLegalBio76(Iyear,Iarea) += (LegalBio76TS(Iyear,Iarea,Istep) + CumCatch(Iyear,Iarea,Istep)) * TimeStepLen(Iyear,Istep);
  //     }
  //    // sLegalBio(Iyear,Iarea) /= float(Nstep);
  //    // sLegalBio76(Iyear,Iarea) /= float(Nstep);
  //   }}


  // Harvest rate //
  // Get catches by area and not fleet to match biomass
  //// This is used as an output * //// ----------------------------------------------------------------------------------
  // for (int Iyear=0;Iyear<Nyear;Iyear++)   {
  //   for (int Istep=0;Istep<Nstep;Istep++){
  //     for (int Ifleet=0;Ifleet<Nfleet;Ifleet++){
  //       for (int Iarea=0;Iarea<Narea;Iarea++){
  //         if (Area_fleet(Iarea,Ifleet) == 1) {
  //               CatchYSA(Iyear,Istep,Iarea) += Catch(Iyear,Istep,Ifleet);   // get total catch by year, time step area
  //         }}}}}

   // for (int Iyear=0;Iyear<Nyear;Iyear++)   {
   //   for (int Istep=0;Istep<Nstep;Istep++){
   //       for (int Iarea=0;Iarea<Narea;Iarea++){
   //         HRint(Iyear,Istep,Iarea) = (1e-17+CatchYSA(Iyear,Istep,Iarea))/LegalBioTS(Iyear,Iarea,Istep);
   //         wHRint(Iyear,Istep,Iarea) = HRint(Iyear,Istep,Iarea) * (1e-17+CatchYSA(Iyear,Istep,Iarea)); }}}

   // Type cumCatchYSA;   Type HR1;
   // for (int Iyear=0;Iyear<Nyear;Iyear++)   {
   //   for (int Iarea=0;Iarea<Narea;Iarea++){
   //     cumCatchYSA = 0; HR1 = 0;
   //     for (int Istep=0;Istep<Nstep;Istep++){
   //       cumCatchYSA += (1e-17+CatchYSA(Iyear,Istep,Iarea));
   //       HR1 += wHRint(Iyear,Istep,Iarea);                                     // Weight the time step HR by catch landed
   //     }
   //     HarvestRateArea(Iyear,Iarea)= HR1/cumCatchYSA;                                // Remove weighting
   //   }}

  // Calculate Harvest Rate by zone
  // for (int Iyear=0;Iyear<Nyear;Iyear++){
  //   for (int Izone=0;Izone<Nzone;Izone++)   {
  //     for (int IareaP=0;IareaP<NareasPerZone(Izone);IareaP++)   {
  //       Iarea = AreasPerZone(Izone,IareaP);
  //       if (Iarea >= 0) {
  //     cumCatchYSA = 0; HR1 = 0;
  //     for (int Istep=0;Istep<Nstep;Istep++){
  //       cumCatchYSA += (1e-17+CatchYSA(Iyear,Istep,Iarea));
  //       HR1 += wHRint(Iyear,Istep,Iarea);                                     // Weight the time step HR by catch landed
  //     }}}
  //     HarvestRateZn(Iyear,Izone)= HR1/cumCatchYSA;                                // Remove weighting
  //   }}


  // Harvest rate //
  //// This is used as an output * //// ----------------------------------------------------------------------------------
  // Harvest Rate 2 off Hrate used in the catch equation
  // Sum Hrate across time steps and record by area (also do for catch to weight averaging)
  HrateYA.setZero();
  for (int Iyear=0;Iyear<Nyear;Iyear++){
    for (int Istep=0;Istep<Nstep;Istep++)   {
      for (int Ifleet=0;Ifleet<Nfleet;Ifleet++)  {
        Iarea = Fleet_area(Ifleet);
        HrateYA(Iyear,Iarea) += Hrate(BurnIn+Iyear,Istep,Ifleet);
        CatchYA(Iyear,Iarea) += Catch(Iyear,Istep,Ifleet) ;
        }}}

  // Average Hrate across areas within a zone
  Type HrateTmp, CatchTmp;
  for (int Iyear=0;Iyear<Nyear;Iyear++)   {
   for (int Izone=0;Izone<Nzone;Izone++)   {
     CatchTmp = 0; HrateTmp = 0;
     for (int IareaP=0;IareaP<NareasPerZone(Izone);IareaP++) {
       Iarea = AreasPerZone(Izone,IareaP);
       if (Iarea >= 0) {
             CatchTmp += CatchYA(Iyear,Iarea);                                      // total catch and total HR for this zone
             HrateTmp += HrateYA(Iyear,Iarea) * CatchYA(Iyear,Iarea);
             } }
     HarvestRate(Iyear,Izone) = (1-exp(-HrateTmp/CatchTmp));                        // For each Zone workout final HR from F
    }}


   // Harvest rate 3
   // Type TotalLB, TotalLB76, CatchLB;
   // for (int Iyear=0;Iyear<Nyear;Iyear++)    {
   //   for (int Izone=0;Izone<Nzone;Izone++)      {
   //     TotalLB = 0;TotalLB76 = 0; CatchLB = 0;
   //     for (int IareaP=0;IareaP<NareasPerZone(Izone);IareaP++)        {
   //       Iarea = AreasPerZone(Izone,IareaP);
   //       if (Iarea >= 0)  {
   //         TotalLB += sLegalBio(Iyear,Iarea);
   //         TotalLB76 += sLegalBio76(Iyear,Iarea);
   //         for (int Istep=0;Istep<Nstep;Istep++){
   //           for (int Ifleet=0;Ifleet<Nfleet;Ifleet++){
   //             if (Area_fleet(Iarea,Ifleet) == 1) {  CatchLB += Catch(Iyear,Istep,Ifleet);}
   //           }
   //         }
   //       }
   //     }
   //     SHarvestRate(Iyear,Izone) = CatchLB/TotalLB;
   //     SHarvestRate76(Iyear,Izone) = CatchLB/TotalLB76;
   //   }
   // }

  Rec_Penal = 0;
  Rec_Penal_SumZero = 0;   // Keep estimated recruit devs summing to Zero

  //Initial_pen = 0;
  for (int Iyear=RecYr1;Iyear<=RecYr2;Iyear++){
   Rec_Penal += log(SigmaR) + RecDevs(Iyear-RecYr1)*RecDevs(Iyear-RecYr1)/(2.0*SigmaR*SigmaR);
   Rec_Penal_SumZero += 1.0*RecDevs(Iyear-RecYr1) ;  }
  Rec_Penal_SumZero = square(Rec_Penal_SumZero);

  Rec_Penal_Smooth = 0; // Keep sequential recruit devs close to each other
  for (int Iyear=1;Iyear<RecDevs.size();Iyear++){ // start at 1 not 0 to allow for offset
    Rec_Penal_Smooth += 1.0*square(RecDevs(Iyear)-RecDevs(Iyear-1));} ;

  neglogL = dummy*dummy + Rec_Penal + Initial_pen + Rec_Penal_Smooth + Rec_Penal_SumZero;
  vector<Type> Select(Nlen(0));

  CatchLike = CatchLikelihood(dataset,thedata, N,Z,Hrate,ActSelex,ActReten,ActLegal, WeightLen, CatchCheck,QRedsPar);
  NumbersLike = NumbersLikelihood(dataset,thedata, N,Z,Hrate,ActSelex,ActReten,ActLegal, WeightLen, PredNumbers, NumbersLikeComps,SigmaNumbers,QRedsPar);
  CpueLike = CpueLikelihood(dataset,thedata, N, Z, ActSelex, ActReten,ActLegal, WeightLen, PredCpue, CpueLikeComps,SigmaCpue,CpueQ,CpueEcreep,Qpars,efpars,M,QRedsPar);
  LengthLike = LengthLikelihood(dataset,thedata, N,ActSelex,ActReten,ActLegal,  PredLengthComp,LengthLikeComps,Select,QRedsPar);
  LarvalLike = LarvalLikelihood(dataset,thedata, RecruitmentByArea,PuerulusByArea,PredLarval,LarvalLikeComps,PuerPowPars);

  Weighted_CpueLike = LambdaCpue*CpueLike;
  Weighted_NumbersLike = LambdaNumbers*NumbersLike;
  Weighted_LengthLike = LambdaLength*LengthLike;
  Weighted_LarvalLike = LambdaLarval*LarvalLike;
  Weighted_TagLike1 = LambdaTag1*sum(TagLike1);
  Weighted_TagLike2 = LambdaTag2*sum(TagLike2);

  neglogL += CatchLike;
  neglogL += LambdaCpue*CpueLike;
  neglogL += LambdaNumbers*NumbersLike;
  neglogL += LambdaLength*LengthLike;
  neglogL += LambdaLarval*LarvalLike;
  neglogL += LambdaTag1*sum(TagLike1);
  neglogL += LambdaTag2*sum(TagLike2);

  // add Penalties derived from parameter priors
  neglogL += MainParPriorPen + RecParPriorPen + SelParPriorPen + EffParPriorPen;

  // Now do projections.
  // ProjType==1 (catch-based): dataset.Catch already holds the projected catch
  // schedule from PROJECTIONS.DAT (Data$Catch was extended into the projection
  // years by ReadProjFile()/LoadData() before MakeADFun() was called), so
  // OneTimeStep() -> Hybrid() solves for the harvest rate that takes it, the
  // same way it does for the assessment period.
  // ProjType==2 (harvest-rate/effort-based): OneTimeStep() reads
  // dat.ProjHarvestRate directly instead of calling Hybrid() -- see the
  // Iyear >= dat.Nyear branch there.
  if (DoProject==1)
  for (int Iyear=Nyear;Iyear<Nyear+Nproj;Iyear++)
   for (int Istep=0;Istep<Nstep;Istep++)
    {
     XX = OneTimeStep(dataset, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, Iyear, Istep, ActGrowth, RecruitFrac, Rbar, IsVirgin, Feqn2, ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,VirginBio, CurrentBio);
    } // year and season

  // ── Post-projection summaries ───────────────────────────────────────────
  // Both computed here (rather than before the projection loop, as they used
  // to be) so N is already populated for the projection years too. Upper
  // bound extends to Nyear+Nproj when DoProject==1, else stays at Nyear
  // exactly as before.
  int NyearSummary = (DoProject==1) ? (Nyear+Nproj) : Nyear;

  // Legal Biomass at predetermined time-step. Including Burn In.
  //// This is used as an output * //// ----------------------------------------------------------------------------------
  int YearAdjusted;
  LegalBioAll.setZero(); LegalBioAllbySex.setZero(); MatureBioAllbySex.setZero();
  for (int Iyear=-BurnIn;Iyear<NyearSummary;Iyear++){
    if (Iyear <= 0) { YearAdjusted = 0; } else { YearAdjusted = Iyear; }
    for (int Iarea=0;Iarea<Narea;Iarea++){
      for (int Isex=0;Isex<Nsex;Isex++){
        for (int Iage=0;Iage<Nage;Iage++){
            for (int Ilen=0;Ilen<Nlen(Isex);Ilen++){
              LegalBioAll(BurnIn+Iyear,Iarea) += LegalRef(Isex,Ilen)*N(Iarea,BurnIn+Iyear,BioTimeStep,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);
              LegalBioAllbySex(BurnIn+Iyear,Iarea,Isex) += LegalRef(Isex,Ilen)*N(Iarea,BurnIn+Iyear,BioTimeStep,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);
              if(Iage>=MatAge(Iarea)) MatureBioAllbySex(BurnIn+Iyear,Iarea,Isex) += N(Iarea,BurnIn+Iyear,BioTimeStep,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);}
          }}}}

  // Simon post-hoc discard calculation
  vector<Type> DiscXX(2);
  int DiscArea;
  for (int Iyear=0;Iyear<NyearSummary;Iyear++)
    for (int Istep=0;Istep<Nstep;Istep++)
      for (int Ifleet=0;Ifleet<Nfleet;Ifleet++)
      {
        DiscArea = Fleet_area(Ifleet);
        DiscXX = DiscardByFleet(dataset,N,Z,Hrate,ActSelex,ActReten,ActLegal,
                                WeightLen,DiscArea,Ifleet,Iyear,Istep,QRedsPar);
        DiscardWt(Iyear,Istep,Ifleet) = DiscXX(0);
        DeadDiscardWt(Iyear,Istep,Ifleet) = DiscXX(1);
      }


  if (DoProject==1)
   {
    REPORT(HarvestRate);
    REPORT(LegalBioAll);
    REPORT(MatBio);
    REPORT(MatBioArea);
    REPORT(RecruitmentByArea);
    REPORT(Hrate);
    REPORT(CpueEcreep);
    REPORT(PredCpue);
    REPORT(MaxProjYr);
    REPORT(VirginBio);
    REPORT(VirginLegalBio);
    REPORT(LegalBio);
    REPORT(LegalBioAllbySex);
    REPORT(MatureBioAllbySex);
    REPORT(DiscardWt);
    REPORT(DeadDiscardWt);
    REPORT(ActSelex);
    REPORT(ActReten);
    REPORT(ActLegal);
    REPORT(N);
    REPORT(Z);
    REPORT(WeightLen);
    REPORT(M);
    REPORT(MWhitesPar);
    REPORT(QRedsPar);

    ADREPORT(MatBio);
    ADREPORT(MatBioArea);
    //ADREPORT(RecruitmentByArea);
    ADREPORT(LegalBioAll);
    ADREPORT(Hrate);
    //ADREPORT(PredCpue.col(0));
    ADREPORT(CpueEcreep);
    }

  if (DoProject==0)
   {
    // ADREPORT(VarOut);
    // ADREPORT(PredNumbers.col(0));
    // ADREPORT(PredLarval.col(0));
    // ADREPORT(PredCpue.col(0));
    if(VarTypes(0)==1) {ADREPORT(MatBio);}
    if(VarTypes(1)==1) {ADREPORT(MatBioArea);}
    if(VarTypes(2)==1) {ADREPORT(RecruitmentByArea);}
    if(VarTypes(3)==1) {ADREPORT(LegalBioAll);}
    if(VarTypes(4)==1) {ADREPORT(HarvestRate);}
    if(VarTypes(5)==1) {ADREPORT(PredCpue.col(0));}
    if(VarTypes(6)==1) {ADREPORT(CpueEcreep);}
    //if(VarTypes(7)==1) {ADREPORT(HarvestRate);}
    //if(VarTypes(8)==1) {ADREPORT(HarvestRate);}
    //if(VarTypes(9)==1) {ADREPORT(HarvestRate);}

    REPORT(HarvestRate);
    REPORT(LegalBioAll);
    REPORT(MatBio);
    REPORT(MatBioArea);
    REPORT(RecruitmentByArea);
    REPORT(Hrate);
    REPORT(CpueEcreep);
    REPORT(PredCpue);
    REPORT(MaxProjYr);
    REPORT(VirginBio);
    REPORT(VirginLegalBio);
    REPORT(LegalBio);
    REPORT(LegalBioAllbySex);
    REPORT(MatureBioAllbySex);
    REPORT(DiscardWt);
    REPORT(DeadDiscardWt);
    REPORT(ActSelex);
    REPORT(ActReten);
    REPORT(ActLegal);
    REPORT(N);
    REPORT(Z);
    REPORT(WeightLen);

    REPORT(CatchCheck);
  	REPORT(ActMove);
  	REPORT(ActGrowth);
    REPORT(Feqn2);
  	REPORT(CatchLike)
  	REPORT(CpueLike)
  	REPORT(CpueLikeComps)
  	REPORT(NumbersLike)
  	REPORT(NumbersLikeComps)
  	REPORT(LengthLike)
  	REPORT(LengthLikeComps);
  	REPORT(LarvalLike)
  	REPORT(LarvalLikeComps);
  	REPORT(SigmaCpue);
  	REPORT(CpueQ);
  	REPORT(PredNumbers);
  	REPORT(SigmaNumbers);
  	REPORT(PredLengthComp);
  	REPORT(SigmaNumbers);
  	REPORT(PredLarval);
  	REPORT(Initial_pen);
  	REPORT(Rec_Penal);
  	REPORT(Rec_Penal_Smooth);
  	REPORT(Ninit);
    REPORT(PuerulusByArea);
   // REPORT(sLegalBio);
   // REPORT(sLegalBio76);
   // REPORT(SHarvestRate);
   // REPORT(HarvestRateArea);
   // REPORT(HarvestRateZn);
   // REPORT(SHarvestRate76);
   // REPORT(LegalBio76);
    REPORT(CumCatch);
    REPORT(Weighted_CpueLike);
    REPORT(Weighted_NumbersLike);
    REPORT(Weighted_LengthLike);
    REPORT(Weighted_LarvalLike);
    REPORT(Weighted_TagLike1);
    REPORT(Weighted_TagLike2);
    REPORT(PredTagSize);
    REPORT(RecapNum);
    REPORT(TagLike1);
    REPORT(TagLike2);
    REPORT(HRint);
    REPORT(HrateYA);
    REPORT(ActRecruitAreaSexDist);
    REPORT(ActRecruitLenDist);
    REPORT(ActRecDev);
    REPORT(Recruits);
    REPORT(neglogL);
    REPORT(Select);
    REPORT(RecruitFrac);
    REPORT(RecruitPars);
    REPORT(VirginNvec);
    REPORT(VirginBioAtLen);
    REPORT(BiasMult);
    REPORT(CurrentBio);
    REPORT(M);
    REPORT(MWhitesPar);
    REPORT(QRedsPar);
    REPORT(GrowthOut);
    REPORT(MainParPriorPen);
    REPORT(RecParPriorPen);
    REPORT(SelParPriorPen);
    REPORT(EffParPriorPen);

    REPORT(MainPars);
    REPORT(RecruitPars);
    REPORT(SelPars);
    REPORT(efpars);
    REPORT(MovePars);

    REPORT(PuerPowPars);
    REPORT(RetPars);
    REPORT(RecDevs);
    REPORT(Qpars);
    REPORT(RecSpatDevs);
    REPORT(MovePars);
    REPORT(GrowthPars);
    REPORT(Rec_Penal_SumZero);


        }

    return neglogL;
}

#include <TMB.hpp>
#include "include/RLhpp.hpp"


// Big changes made by Simon.  Density dependent mortality, legal biomass, area-specific burn-in years

// -------------------------------------------------------------------------------------------------------------------

template <class Type>
matrix<Type> MakeM(dataSet<Type> &dat, matrix<Type> &Mlow, matrix<Type> &Mvirgin, vector<Type> &Vbio, vector<Type> &Cbio, Type inflect){
  matrix<Type> Mtemp(dat.Narea, dat.Nage);
  Mtemp.setZero();
  for (int Iarea=0;Iarea<dat.Narea;Iarea++) {
    if(Cbio(Iarea)>Vbio(Iarea))  Cbio(Iarea) = Vbio(Iarea); // Limit M to max  
    for (int Iage=0;Iage<dat.Nage;Iage++) {
      Mtemp(Iarea,Iage) = Mlow(Iarea,Iage) + (2.0*(Mvirgin(Iarea,Iage)-Mlow(Iarea,Iage)))/(1.0+exp((Cbio(Iarea)-Vbio(Iarea))/-(inflect*Vbio(Iarea)))) ;
    }
}
  return(Mtemp);
  }

// -------------------------------------------------------------------------------------------------------------------

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
       }
    }
  } // Tune
  return(Hrate);
}

// -------------------------------------------------------------------------------------------------------------------

template <class Type>
matrix<Type> SetUpSelex(dataSet<Type> &dat,  vector<Type> &SelPars, matrix<Type> &SelexFI, matrix<int> &PatSpec, int Npatterns ){

  int IPreSpecified, IselParPnt, Isex;
  Type p1,p2;

  int MaxLen; MaxLen = dat.MaxLen;
  matrix<Type> ActSelex(Npatterns,MaxLen);
  Type Mult, LML, Offset;

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
    if (PatSpec(IselPattern,1) == SELEX_COEFFICIENTS)
    {
      for (int Isize=0;Isize<PatSpec(IselPattern,3);Isize++) { IselParPnt += 1; ActSelex(IselPattern,Isize) = exp(SelPars(IselParPnt)); }
      for (int Isize=PatSpec(IselPattern,3);Isize<MaxLen; Isize++) ActSelex(IselPattern,Isize) = 1;
    }

    // Logistic
    if (PatSpec(IselPattern,1) == SELEX_LOGISTIC)
    {
      p1 = SelPars(IselParPnt+1); p2 = SelPars(IselParPnt+2);
      IselParPnt += 2;
      Isex = PatSpec(IselPattern,2);
      for (int Isize=0;Isize<MaxLen; Isize++) ActSelex(IselPattern,Isize) = 1.0/(1.0+exp(-p2*(dat.MidLenBin(Isex,Isize)-p1)));
    }

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
                         vector<Type> &VirginBio, vector<Type> &CurrentBio, matrix<Type> Mlow, matrix<Type> Mvirgin, Type Minflection, array<Type> &Mtempts){

 array<Type> Z_rate(dat.Nsex,dat.Nage,dat.MaxLen);                             // total mortality

 array<Type> selexF(dat.Nfleet,dat.Nsex,dat.Nage,dat.MaxLen);                  // Selectivity
 array<Type> retainF(dat.Nfleet,dat.Nsex,dat.Nage,dat.MaxLen);                 // Retention
 array<Type> selretwght(dat.Nfleet,dat.Nsex,dat.Nage,dat.MaxLen);              // Product of selectivity,retention and weight
 array<Type> Ntemp(dat.Narea,dat.Nsex,dat.Nage,dat.MaxLen);                    // N matrix (after mortality)
 array<Type> Nmove(dat.Narea,dat.Nsex,dat.Nage,dat.MaxLen);                    // N matrix (after movement)
 matrix<Type> Mtemp(dat.Narea,dat.Nage);                                       // Teporary M matrix
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
   
 // Work out correct M for this time-step and year and record it in another object (Mtempts)
  Mtemp = MakeM(dat, Mlow, Mvirgin, VirginBio, CurrentBio, Minflection);  
  for (int Iarea=0;Iarea<dat.Narea;Iarea++){
    for (int Iage=0;Iage<dat.Nage;Iage++){
        if(Iyear>=0) {
           Mtempts(Iarea, Iage, dat.BurnIn+Iyear) = Mtemp(Iarea, Iage);
          // Mtemp(Iarea, Iage) = M(Iarea, Iage);
           }
         }} 
  
  
 // Maturity and fecundity
 if (Istep==dat.MatTimeStep)
  {
   MatBio(dat.BurnIn+Iyear) = 0;
   for (int Iarea=0;Iarea<dat.Narea;Iarea++)
    {
     MatBioArea(Iarea,dat.BurnIn+Iyear) = 0;
     for (int Iage=0;Iage<dat.Nage;Iage++)
      for (int Isize=0;Isize<dat.Nlen(1);Isize++)
       MatBioArea(Iarea,dat.BurnIn+Iyear) += N(Iarea,dat.BurnIn+Iyear,Istep,1,Iage,Isize)*dat.MatFem(Iarea,Isize);
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
     else
      {
       HratePass = Hybrid(dat, N, selretwght, selexF, retainF, Mtemp, Iarea, Iyear, Istep, MWhitesPar);
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
         Z_rate(Isex,Iage,Isize) = dat.TimeStepLen(0,Istep)*Mtemp(Iarea,Iage)*ScaleWhiteM;}
       else{
         Z_rate(Isex,Iage,Isize) = dat.TimeStepLen(Iyear,Istep)*Mtemp(Iarea,Iage)*ScaleWhiteM;}
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
 int Lenefseries = CpueEcreep.rows();
 int Nefseries = CpueEcreep.cols();
 int efcnt = -1; int parcnt = -1;
 Type Tmppar;  // store temporary parameter
 for (int Nef=0;Nef<Nefseries;Nef++){
   CpueEcreep(0,Nef) = 1.0;  // Set first year to 1 (no efficiency creep) 
   parcnt = parcnt + 1;
   efcnt = -1;
   for (int Yef=1;Yef<Lenefseries;Yef++){
     if(Yef<(dat.Nyear))  { 
       efcnt = efcnt+1;
       if(efcnt==thedata.EffCrLag(Nef)){
         efcnt = -1;
         parcnt = parcnt + 1;   }
       Tmppar = efpars(parcnt);
       CpueEcreep(Yef,Nef) =  CpueEcreep(Yef-1,Nef) * (1+(Tmppar/100));
       } else {CpueEcreep(Yef,Nef) = CpueEcreep(Yef-1,Nef);}
   }}

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
 int Iarea, Idata,Iyr;
 Type Obs,CV,ncnt,SS,Residual;
 vector<Type> qestLar(dat.Narea);

 for (Iarea=0;Iarea<dat.Narea;Iarea++){
   for(Iyr=0;Iyr<dat.BurnIn+dat.Nyear+dat.Nproj;Iyr++){
     PuerulusByArea(Iarea,Iyr) = exp(log(RecruitmentByArea(Iarea,Iyr)) / PuerPowPars(Iarea));
   }}

  NeglogLikelihood = 0;
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
   }

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
     Offset = Iarea*dat.Nage*dat.MaxLen+Iage*dat.MaxLen;
     for (int Isize=0;Isize<dat.MaxLen;Isize++)
      {
       I(Offset+Isize,Offset+Isize) = 1.0;
       S(Offset+Isize,Offset+Isize) = exp(-M(Iarea,Iage)*ScaleWhiteM+ActSelex(Isex,Isize)*Finitial);
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
                         vector<Type> &VirginBio, vector<Type> &CurrentBio, matrix<Type> Mlow, matrix<Type> Mvirgin, Type Minflection, array<Type> &Mtempts) {

  Type Initial_pen;
  int IsVirgin;                                                           // Set to 1 for unfished state
  matrix<Type> Feqn2(dat.Nfleet,dat.Nstep); Feqn2.setZero();              // Initial F (not used in projections)
  matrix<Type> Feqn3(dat.Nfleet,dat.Nstep); Feqn3.setZero();              // Initial F (not used in projections)
  matrix<Type> Feqn4(dat.Nfleet,dat.Nstep); Feqn4.setZero();              // Used to set Burn_in F to zero if one rarea does not want burn in 
  vector<Type> XX(2);                                                     // Dummy variables
  array<Type> Fvals(dat.Nfleet,dat.Num_Iteration,dat.Nstep);                             // Storage for tuning of Fs
  array<Type> Ninit2(dat.Narea,dat.Nsex,dat.Nage,dat.MaxLen);

  Ninit = VirginN(dat, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, Mvirgin, ActGrowth, RecruitFrac,
         Rbar, ActRecruitAreaSexDist, ActRecruitLenDist, ActRecDev, QRedsPar, MWhitesPar,Type(0.0));
  //  Get bare bones numbers by area, sex, age and length - one recruitment / move / grow - no Mort. 
  
  // Virgin Biomass from Ninit -  Has M but not F
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
      XX = OneTimeStep(dat, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, Iyear, Istep, ActGrowth, RecruitFrac, Rbar, IsVirgin, Feqn3,
                        ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,
                        VirginBio, CurrentBio, Mlow, Mvirgin, Minflection, Mtempts);
    
    // Multiyear projection with No F (F set to Zero) This updates the future time-step under no fishing
    for (int Iyear=-dat.BurnIn+1;Iyear<dat.Tune_Years;Iyear++)  
     for (int Istep=0;Istep<dat.Nstep;Istep++)
      XX = OneTimeStep(dat, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, Iyear, Istep, ActGrowth, RecruitFrac, Rbar, IsVirgin, Feqn2,
                        ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,
                        VirginBio, CurrentBio, Mlow, Mvirgin, Minflection, Mtempts);
    
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
  
  // Penalty on non-convergence
  Initial_pen = 0;
  for (int Ifleet=0;Ifleet<dat.Nfleet;Ifleet++)
   for (int Istep=0;Istep<dat.Nstep;Istep++)
    Initial_pen += 1000000*square(Fvals(Ifleet,dat.Num_Iteration-2,Istep)-Fvals(Ifleet,dat.Num_Iteration-1,Istep));

  // Really poor initial condition
  N.setZero();
  for (int Iarea=0;Iarea<dat.Narea;Iarea++)
   for (int Isex=0;Isex<dat.Nsex;Isex++)
    for (int Iage=0;Iage<dat.Nage;Iage++)
     for (int Isize=0;Isize<dat.Nlen(Isex);Isize++)
      {  N(Iarea,0,0,Isex,Iage,Isize) = Ninit(Iarea,Isex,Iage,Isize); }

  // One year zero catch projection
  for (int Iyear=-dat.BurnIn;Iyear<-dat.BurnIn+1;Iyear++)
   for (int Istep=0;Istep<dat.Nstep;Istep++)
    XX = OneTimeStep(dat, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, Iyear, Istep, ActGrowth, RecruitFrac, Rbar, IsVirgin, Feqn3,
                        ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,
                        VirginBio, CurrentBio, Mlow, Mvirgin, Minflection, Mtempts);

  // Multiyear projection with F (but to year 0)
  for (int Iyear=-dat.BurnIn+1;Iyear<0;Iyear++){
   for (int Istep=0;Istep<dat.Nstep;Istep++){
     for (int Iarea=0;Iarea<dat.Narea;Iarea++){                          // Simon add.  Allow for turning F to Zero in burnin so we can have multiple burn_in times 
       for (int Ifleet=0; Ifleet<dat.Nfleet;Ifleet++) {
         if (dat.Area_fleet(Iarea,Ifleet)==1){
           Feqn4(Ifleet,Istep) = Feqn2(Ifleet,Istep);                    // Always set Feqn4 to Feqn2 as default 
           if(-dat.BurnInVec(Iarea)>Iyear) Feqn4(Ifleet,Istep) = 0; }}}  // Set Feqn4 to Zero is before Burn_in wants to start
     XX = OneTimeStep(dat, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, Iyear, Istep, ActGrowth, RecruitFrac, Rbar, IsVirgin, Feqn4,
                         ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,
                         VirginBio, CurrentBio, Mlow, Mvirgin, Minflection, Mtempts);}}
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

// ========================================================================================================================

/*template <class Type>
 Type TagDym(dataSet<Type> &dat, TheData<Type> &thedata, int SexPass, int GrpPass,array<Type> &N,
             matrix<Type> &ActSelex, matrix<Type> &ActReten, matrix<Type> &ActLegal,
             array<Type> &ActGrowth,matrix<Type> &ActMove,
             matrix<Type> &M, array<Type> &Hrate, Type QRedsPar, Type MWhitesPar,
             array<Type> &Ntag,array<Type> &RecapNum, matrix<Type> &NotReported,
             matrix<Type> &TagLike1, matrix<Type> &TagLike2, array<Type> &PredTagSize) {

 int Jyear, Kyear;
 int SelPointer,RetPointer,LegalPointer,GrowthPointer,MovePointer,IsMoves,IdestArea;
 Type NtagRel,CumReleases;
 Type RetainTemp,NtotalT,ScaleRedQ,ScaleWhiteM;
 Type PartialF,TotalPartialF,FullF2,Deaths;
 Type ObsL,PredL,ObsSS,LikeSize1,LikeTag2,LikeCompT;                               // Likelihood computation
 Type TotalReported;
 vector<Type> Ntemp2(dat.MaxLen);                                                  // Temporary storage
 vector<Type> MoveVec(dat.MaxLen);

 array<Type> selexF(dat.Nfleet,dat.Nsex,dat.Nage,dat.MaxLen);                      // Selectivity
 array<Type> retainF(dat.Nfleet,dat.Nsex,dat.Nage,dat.MaxLen);                     // Retention
 array<Type> M_rate_Tag(dat.Narea,dat.Nage,dat.MaxLen);
 array<Type> Z_rate_Tag(dat.Narea,dat.Nage,dat.MaxLen);
 array<Type> RecapTmp(dat.Narea,thedata.NrepSplit,dat.MaxLen);
 array<Type> Ntemp_Tag(thedata.NtagLag+1,dat.Narea,dat.Nage,dat.MaxLen);            // Temporary storage
 array<Type> Nmove_Tag(thedata.NtagLag+1,dat.Narea,dat.Nage,dat.MaxLen);            // Extra temporarry storage

 Type NeglogLikelihood;
 int NtagLag = thedata.NtagLag;
 int Nage = dat.Nage;
 int Nfleet = dat.Nfleet;
 int NrepSplit = thedata.NrepSplit;
 int Narea = dat.Narea;

 NeglogLikelihood = 0;
 NotReported(SexPass,GrpPass) = 0;
 CumReleases = 0;
 for (int Iyear=thedata.Year1Tag(GrpPass)-dat.First_yr;Iyear<dat.Nyear;Iyear++)
  for (int Istep=0;Istep<dat.Nstep;Istep++)
   {
    Jyear = Iyear+dat.First_yr-thedata.TagYr1;
    Kyear = dat.BurnIn+Iyear;

    // Add the tags that have been out long enuough
    for (int Iarea=0;Iarea<Narea;Iarea++)
     for (int Iage=0;Iage<dat.Nage;Iage++)
      for (int Isize=0;Isize<dat.MaxLen;Isize++)
       {
        Ntag(SexPass,GrpPass,Iarea,0,Iage,Isize) += Ntag(SexPass,GrpPass,Iarea,1,Iage,Isize);
        for (int ItagLag=NtagLag-1;ItagLag>0;ItagLag--)
         Ntag(SexPass,GrpPass,Iarea,ItagLag,Iage,Isize) = Ntag(SexPass,GrpPass,Iarea,ItagLag+1,Iage,Isize);
        Ntag(SexPass,GrpPass,Iarea,NtagLag,Iage,Isize) = 0;
       }

    // Add new tags
    for (int Iarea=0;Iarea<Narea;Iarea++)
     {
      NtagRel = thedata.TagRel(SexPass,GrpPass,Iarea,Jyear,Istep,0);
      if (NtagRel > 0)
       {
       // cout << "R "<< SexPass << " " << GrpPass << " " << Iyear << " " << Istep << " " << NtagRel << " " << Kyear << " " << Iarea << endl;
        // Need to add in Type I tag-loss
        for (int Isize=0;Isize<dat.MaxLen;Isize++)
         if (thedata.TagRel(SexPass,GrpPass,Iarea,Jyear,Istep,Isize+1) > 0)
          {
           // Now divide into age
           NtotalT = 0;
           for (int Iage=0;Iage<Nage;Iage++) NtotalT += N(Iarea,Kyear,Istep,SexPass,Iage,Isize);
           for (int Iage=0;Iage<Nage;Iage++) Ntag(SexPass,GrpPass,Iarea,NtagLag,Iage,Isize)=
            N(Iarea,Kyear,Istep,SexPass,Iage,Isize)/NtotalT*thedata.TagRel(SexPass,GrpPass,Iarea,Jyear,Istep,Isize+1)*thedata.InitialLoss;
          }
         // Store tags not reported because they lost their tags (or died) at tagging
         NotReported(SexPass,GrpPass) += thedata.TagRel(SexPass,GrpPass,Iarea,Jyear,Istep,0)*(1.0-thedata.InitialLoss);
         CumReleases  += thedata.TagRel(SexPass,GrpPass,Iarea,Jyear,Istep,0);
       } // if
     } // area

    // Set selectivity
    for (int Ifleet=0; Ifleet<Nfleet;Ifleet++)
     for (int Iage=0;Iage<Nage;Iage++)
      {
       if(dat.IsRed(SexPass,Iage,dat.Fleet_area(Ifleet),Istep)==1) {ScaleRedQ = QRedsPar; } else {ScaleRedQ = 1.0; }
       SelPointer = dat.SelPnt(SexPass,Iage,Ifleet,Iyear,Istep);
       RetPointer = dat.RetPnt(SexPass,Iage,Ifleet,Iyear,Istep);
       LegalPointer = dat.LegalFleetPnt(SexPass,Iage,Ifleet,Iyear,Istep);
       for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
        {
         selexF(Ifleet,SexPass,Iage,Isize) = ActSelex(SelPointer,Isize) * ScaleRedQ;
         retainF(Ifleet,SexPass,Iage,Isize) = ActReten(RetPointer,Isize)*ActLegal(LegalPointer,Isize);
        }
      }

    // Compute Z given F and M (note that M includes the log-term tag-loss rate)
    RecapTmp.setZero();
    for (int Iarea=0;Iarea<Narea;Iarea++)
     {
      for (int Iage=0;Iage<Nage;Iage++)
       {
        if(dat.IsRed(SexPass,Iage,Iarea,Istep)==0) {ScaleWhiteM = MWhitesPar; } else {ScaleWhiteM = 1.0; }
        for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
         {
          M_rate_Tag(Iarea,Iage,Isize) = dat.TimeStepLen(Iyear,Istep)*M(Iarea,Iage)*ScaleWhiteM+dat.TimeStepLen(Iyear,Istep)*thedata.TagLossRate;
          Z_rate_Tag(Iarea,Iage,Isize) = M_rate_Tag(Iarea,Iage,Isize);
          for (int Ifleet=0;Ifleet<Nfleet;Ifleet++)
           if (dat.Area_fleet(Iarea,Ifleet)==1)
            {
             // RetAdjust
             RetainTemp = selexF(Ifleet,SexPass,Iage,Isize) * (retainF(Ifleet,SexPass,Iage,Isize)+dat.Phi(Ifleet,Iage,Iyear,Istep)*(1.0-retainF(Ifleet,SexPass,Iage,Isize)));
             Z_rate_Tag(Iarea,Iage,Isize) += Hrate(Kyear,Istep,Ifleet)*RetainTemp;
            }
          //Non-reported tags (only M and tag-loss for tags that are not fully mixed)
          for (int ItagLag=1;ItagLag<=NtagLag;ItagLag++)
           NotReported(SexPass,GrpPass) += Ntag(SexPass,GrpPass,Iarea,ItagLag,Iage,Isize)*(1.0-exp(-M_rate_Tag(Iarea,Iage,Isize)));

          // Now handle the animals that could be reported
          Deaths = Ntag(SexPass,GrpPass,Iarea,0,Iage,Isize)*(1.0-exp(-Z_rate_Tag(Iarea,Iage,Isize)))/Z_rate_Tag(Iarea,Iage,Isize);
          // Natural mortality and tagloss (not reported)
          NotReported(SexPass,GrpPass) += M_rate_Tag(Iarea,Iage,Isize) * Deaths;
          for (int Ifleet=0;Ifleet<Nfleet;Ifleet++)
           if (dat.Area_fleet(Iarea,Ifleet)==1)
            {
             // RetAdjust
             RetainTemp = selexF(Ifleet,SexPass,Iage,Isize) * (retainF(Ifleet,SexPass,Iage,Isize)+dat.Phi(Ifleet,Iage,Iyear,Istep)*(1.0-retainF(Ifleet,SexPass,Iage,Isize)));
             FullF2 = Hrate(Kyear,Istep,Ifleet)*RetainTemp;
             TotalPartialF = 0;
             for (int IrepSplit=0;IrepSplit<NrepSplit;IrepSplit++)
              {
               // Check this Jyear
               PartialF = thedata.RepRate(IrepSplit)*thedata.PropRepSplit(Jyear,Istep,Iarea,IrepSplit)*FullF2;
               TotalPartialF += PartialF;
               RecapTmp(Iarea,IrepSplit,Isize) += PartialF*Deaths;
              }
             //Tagged animals that were discarded and died + tagged animals that were cuaght but not recaptured
             NotReported(SexPass,GrpPass) += (FullF2-TotalPartialF)*Deaths;
            }  // fleet
         }   // Isize
       } // Iage
     } // Iarea

    // Total the tags (this will be used in the total recapture likelihood) + the "not recaptured" class
    for (int Iarea=0;Iarea<Narea;Iarea++)
     for (int IrepSplit=0;IrepSplit<NrepSplit;IrepSplit++)
      for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
       RecapNum(SexPass,GrpPass,Iarea,IrepSplit,Jyear,Istep) += RecapTmp(Iarea,IrepSplit,Isize);

    // This is the likelihood for the length-comp of the recaptured by year, step, group, and type
    LikeSize1 = 0;
    for (int Iarea=0;Iarea<Narea;Iarea++)
     for (int IrepSplit=0;IrepSplit<NrepSplit;IrepSplit++)
      if (thedata.FitTagSizes(IrepSplit) == 1)
       if (thedata.TagRec(SexPass,GrpPass,Iarea,IrepSplit,Jyear,Istep,0)>0)
        {
         ObsSS = thedata.TagRec(SexPass,GrpPass,Iarea,IrepSplit,Jyear,Istep,0);
         NtotalT = 0;
         for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++) NtotalT+=RecapTmp(Iarea,IrepSplit,Isize);
         for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
          {
           PredL = RecapTmp(Iarea,IrepSplit,Isize)/NtotalT;
           PredTagSize(SexPass,GrpPass,Iarea,Isize) += PredL*ObsSS;
           if (thedata.TagRec(SexPass,GrpPass,Iarea,IrepSplit,Jyear,Istep,Isize+1)>0)
            {
             ObsL = thedata.TagRec(SexPass,GrpPass,Iarea,IrepSplit,Jyear,Istep,Isize+1)/ObsSS;
             LikeSize1 -= ObsL*ObsSS*log(PredL/ObsL);
            }
	      }
        }
    TagLike1(SexPass,GrpPass) += LikeSize1;

    // Remove mortality and compute returns
    Ntemp_Tag.setZero();
    for (int Iarea=0;Iarea<Narea;Iarea++)
     for (int Iage=0;Iage<Nage;Iage++)
      for (int Isize=0;Isize<dat.MaxLen;Isize++)
       for (int ItagLag=0;ItagLag<=NtagLag;ItagLag++)
        if (ItagLag==0)
         {
         Ntemp_Tag(0,Iarea,Iage,Isize) = Ntag(SexPass,GrpPass,Iarea,0,Iage,Isize)*exp(-Z_rate_Tag(Iarea,Iage,Isize));
         }
        else
        {
         Ntemp_Tag(ItagLag,Iarea,Iage,Isize) = Ntag(SexPass,GrpPass,Iarea,ItagLag,Iage,Isize)*exp(-M_rate_Tag(Iarea,Iage,Isize));
        }

    // Growth
    for (int Iarea=0;Iarea<Narea;Iarea++)
     for (int ItagLag=0;ItagLag<=NtagLag;ItagLag++)
      for (int Iage=0;Iage<Nage;Iage++)
       {
 	    GrowthPointer = dat.GrowthPnt(Iarea,SexPass,Iage,Iyear,Istep);
        if (GrowthPointer >=0)
         {
          // Key issue (pointer to growth matrix)
          Ntemp2.setZero();
          for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
 	       {
 	        for (int Jsize=0;Jsize<=Isize;Jsize++) Ntemp2(Isize) += Ntemp_Tag(ItagLag,Iarea,Iage,Jsize)*ActGrowth(GrowthPointer,Isize,Jsize);
 	       }
 	       for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++) Ntemp_Tag(ItagLag,Iarea,Iage,Isize) = Ntemp2(Isize);
 	     }
       } // Iarea, ItagLag,Iage

    // Movement
    Nmove_Tag.setZero();
    IsMoves = 0;
    for (int Iarea=0;Iarea<Narea;Iarea++)
     for (int Iage=0;Iage<Nage;Iage++)
      {
       MovePointer = int(dat.MovePnt(Iarea,Iage,Iyear,Istep));
       if (MovePointer > 0)
        {
         IsMoves = 1;
         IdestArea = dat.MoveSpec(MovePointer,2);
         for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++) MoveVec(Isize) = ActMove(MovePointer,Isize);
	      for (int ItagLag=0;ItagLag<=NtagLag;ItagLag++)
          for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
           {
            Nmove_Tag(ItagLag,IdestArea,Iage,Isize) += MoveVec(Isize)*Ntemp_Tag(ItagLag,Iarea,Iage,Isize);
            Nmove_Tag(ItagLag,Iarea,Iage,Isize) -= MoveVec(Isize)*Ntemp_Tag(ItagLag,Iarea,Iage,Isize);
           }
        } // If there was a move
      } // All areas and ages

   // Only update if needed
    if (IsMoves==1)
     {
      for (int ItagLag=0;ItagLag<=NtagLag;ItagLag++)
       for (int Iarea=0;Iarea<Narea;Iarea++)
        for (int Iage=0;Iage<Nage;Iage++)
         for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
          Ntemp_Tag(ItagLag,Iarea,Iage,Isize) += Nmove_Tag(ItagLag,Iarea,Iage,Isize);
     }

    // Copy back and update ages
    for (int ItagLag=0;ItagLag<=NtagLag;ItagLag++)
     for (int Iarea=0;Iarea<Narea;Iarea++)
      {
       if (Istep<dat.Nstep-1)
        {
         for (int Iage=0;Iage<Nage;Iage++)
          for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
           Ntag(SexPass,GrpPass,Iarea,ItagLag,Iage,Isize) = Ntemp_Tag(ItagLag,Iarea,Iage,Isize);
        }
       else
        {
         // special case
         if (Nage-1 > 0)
          {
           for (int Iage=0;Iage<Nage-1;Iage++)
            for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
             Ntag(SexPass,GrpPass,Iarea,ItagLag,Iage+1,Isize) = Ntemp_Tag(ItagLag,Iarea,Iage,Isize);
           for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
            Ntag(SexPass,GrpPass,Iarea,ItagLag,Nage-1,Isize) = Ntemp_Tag(ItagLag,Iarea,Nage-1,Isize) + Ntemp_Tag(ItagLag,Iarea,Nage-2,Isize);
           for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++) Ntag(SexPass,GrpPass,Iarea,ItagLag,0,Isize) = 0;
          }
         else
          {
           for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
            Ntag(SexPass,GrpPass,Iarea,ItagLag,0,Isize)  = Ntemp_Tag(ItagLag,Iarea,0,Isize);
          }
        }
      } // grp, lag, area

   } // Year and step

 // Add animals at the end of the projection to NotReported
 for (int ItagLag=0;ItagLag<=NtagLag;ItagLag++)
  for (int Iarea=0;Iarea<Narea;Iarea++)
   for (int Iage=0;Iage<Nage;Iage++)
    for (int Isize=0;Isize<dat.Nlen(SexPass);Isize++)
     NotReported(SexPass,GrpPass) += Ntag(SexPass,GrpPass,Iarea,ItagLag,Iage,Isize);


 // Total reported (diagnostic) and rescale recpatures
 TotalReported = 0;
 for (int Iarea=0;Iarea<Narea;Iarea++)
  for (int Iyear=0;Iyear<thedata.NyearTags;Iyear++)
   for (int Istep=0;Istep<dat.Nstep;Istep++)
    for (int IrepSplit=0;IrepSplit<NrepSplit;IrepSplit++)
     {
      TotalReported += RecapNum(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep);
      RecapNum(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep) /= thedata.NrelTotal(SexPass,GrpPass);
      }
  NotReported(SexPass,GrpPass) /= thedata.NrelTotal(SexPass,GrpPass);

  // Likelihood (numbers not repatured plus those captured)
  LikeTag2 = thedata.NrelTotal(SexPass,GrpPass)*thedata.NotReportedObs(SexPass,GrpPass)*
                 log(NotReported(SexPass,GrpPass)/thedata.NotReportedObs(SexPass,GrpPass));
  for (int Iarea=0;Iarea<Narea;Iarea++)
   for (int Iyear=0;Iyear<thedata.NyearTags;Iyear++)
    for (int Istep=0;Istep<dat.Nstep;Istep++)
     for (int IrepSplit=0;IrepSplit<NrepSplit;IrepSplit++)
      if (thedata.RecapObs(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep) > 0)
       {
        LikeCompT = thedata.NrelTotal(SexPass,GrpPass)*thedata.RecapObs(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep)*log(RecapNum(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep)/thedata.RecapObs(SexPass,GrpPass,Iarea,IrepSplit,Iyear,Istep));
        LikeTag2 += LikeCompT;
       }
 TagLike2(SexPass,GrpPass) +=  LikeTag2;

 return( NeglogLikelihood);
}
*/
// ========================================================================================================================

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
  DATA_MATRIX(MatFem); dataset.MatFem = MatFem;
  DATA_INTEGER(InitOpt); dataset.InitOpt = InitOpt;
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
  DATA_SCALAR(WeightInitialN); thedata.WeightInitialN = WeightInitialN;
  DATA_SCALAR(WeightInit3); thedata.WeightInit3 = WeightInit3;
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
 // DATA_INTEGER(NtagGroups); thedata.NtagGroups = NtagGroups;
 // DATA_SCALAR(InitialLoss); thedata.InitialLoss = InitialLoss;
//  DATA_SCALAR(TagLossRate); thedata.TagLossRate = TagLossRate;
//  DATA_INTEGER(NrepSplit); thedata.NrepSplit = NrepSplit;
//  DATA_INTEGER(NtagLag); thedata.NtagLag = NtagLag;
//  DATA_VECTOR(RepRate); thedata.RepRate = RepRate;
//  DATA_IVECTOR(FitTagSizes); thedata.FitTagSizes = FitTagSizes;
//  DATA_IVECTOR(Year1Tag); thedata.Year1Tag = Year1Tag;
//  DATA_IVECTOR(Year2Tag); thedata.Year2Tag = Year2Tag;
//  DATA_INTEGER(TagYr1); thedata.TagYr1 = TagYr1;
//  DATA_INTEGER(TagYr2); thedata.TagYr2 = TagYr2;
//  DATA_INTEGER(NyearTags); thedata.NyearTags = NyearTags;
//  DATA_ARRAY(TagRel); thedata.TagRel = TagRel;
//  DATA_ARRAY(TagRec); thedata.TagRec = TagRec;
//  DATA_ARRAY(RecapObs); thedata.RecapObs = RecapObs;
//  DATA_MATRIX(NrelTotal); thedata.NrelTotal = NrelTotal;
//  DATA_MATRIX(NotReportedObs); thedata.NotReportedObs = NotReportedObs;
//  DATA_ARRAY(PropRepSplit); thedata.PropRepSplit = PropRepSplit;

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
  DATA_VECTOR(LegalRef);

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
  PARAMETER_VECTOR(InitPars);
  PARAMETER_VECTOR(RecSpatDevs);
  PARAMETER_VECTOR(MovePars);
  PARAMETER_VECTOR(GrowthPars);
  PARAMETER(dummy);

  matrix <Type> M(Narea, Nage);
  matrix <Type> Mlow(Narea, Nage);
  matrix <Type> Mvirgin(Narea, Nage);
  Type Rbar;
  Type MWhitesPar;
  Type Minflection;
  Type QRedsPar;
  Type SigmaR;
  vector <Type> LogRinitial(Narea);
  Type Finitial;

  //int BurnIn = 0;
  //for(int Iarea=0;Iarea<Narea;Iarea++){ if(BurnIn<BurnInVec(Iarea)) BurnIn=BurnInVec(Iarea); }
  //dataset.BurnIn=BurnIn; 

  Rbar = MainPars(0);
  for(int Iarea=0;Iarea<Narea;Iarea++){
    for (int Iage=0;Iage<Nage;Iage++){
      M(Iarea,Iage) =    MainPars(1+Iarea) * MainPars(1+Narea+Iage); 
      Mlow(Iarea,Iage) = MainPars(1+Iarea) * MainPars(1+Narea+Iage); 
      Mvirgin(Iarea,Iage) = Mlow(Iarea,Iage) + MainPars(2+Narea+Nage); }}
  MWhitesPar = MainPars(1+Narea+Nage);
  Minflection = exp(MainPars(3+Narea+Nage));
  QRedsPar = MainPars(4+Narea+Nage);
  SigmaR = MainPars(5+Narea+Nage);
  for (int Iarea=0;Iarea<Narea;Iarea++) LogRinitial(Iarea) = MainPars(7+Narea+Nage+Iarea);
  Finitial = exp(MainPars(6+Narea+Nage));
  
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
  //matrix<Type> Mtemp(Narea, Nage);                                                        // M matrix
  array<Type> Mtempts(Narea, Nage, BurnIn+Nyear+MaxProjYr+1);                                // M array
  Mtempts -= 99;
  vector<Type> Recruits(BurnIn+Nyear+MaxProjYr+1);                                         // Recruitment output
  vector<Type> BiasMult(BurnIn+Nyear+MaxProjYr+1);                                         // Bias correction factor
  array<Type>  Ninit(Narea,Nsex,Nage,MaxLen);
  vector<Type> MatBio(BurnIn+Nyear+MaxProjYr+1);
  matrix<Type> MatBioArea(Narea,BurnIn+Nyear+MaxProjYr+1);
  matrix<Type> RecruitmentByArea(Narea,BurnIn+Nyear+MaxProjYr+1);                          // Recruitment
  matrix<Type> PuerulusByArea(Narea,BurnIn+Nyear+MaxProjYr+1);                             // Puerulus

  array<Type> Hrate(BurnIn+Nyear+MaxProjYr,Nstep,Nfleet);                                  // Harvest rate
  matrix<Type> LegalBio(Nyear,Narea);                                                      // Legal biomass
  matrix<Type> LegalBioAll(BurnIn+Nyear+MaxProjYr,Narea);                                                      // Legal biomass
  matrix<Type> LegalBio76(Nyear,Narea);                                                    // Legal biomass of all lobster > 76 mm
  array<Type> CumCatch(Nyear,Narea,Nstep);
  array<Type> LegalBioTS(Nyear,Narea,Nstep);                                               // Legal biomass
  array<Type> LegalBio76TS(Nyear,Narea,Nstep);                                             // Legal biomass of all lobster > 76 mm
  matrix<Type> sLegalBio(Nyear,Narea);                                                     // Legal biomass
  matrix<Type> sLegalBio76(Nyear,Narea);                                                   // Legal biomass of all lobster > 76 mm
  matrix<Type> HarvestRate(Nyear,Nzone);                                                   // Harvest rate by year and zone
  matrix<Type> HarvestRate76(Nyear,Nzone);                                                 // Harvest rate by year and zone of all lobster > 76 mm
  array<Type> CatchCheck(Nyear+MaxProjYr,Nstep,Nfleet);                                              // Check
  matrix<Type> ActSelex(NselPatterns,MaxLen);
  matrix<Type> ActReten(NretPatterns,MaxLen);
  matrix<Type> ActLegal(NlegalPatterns,MaxLen);
  matrix<Type> ActMove(NmovePatterns,MaxLen);
  array<Type> ActRecruitAreaSexDist(Nyear+MaxProjYr,Nstep,Narea,Nsex);                     // Allocation on recruitment to areas and sexes
  array<Type> ActRecruitLenDist(NrecruitPatternsB,Nsex,MaxLen);
  array<Type> ActGrowth(NgrowthPatterns,MaxLen,MaxLen);
  vector<Type> ActRecDev(BurnIn+Nyear+MaxProjYr+1);                                        // Recruitment deviations
  vector<Type> VirginBio(Narea);                                                           // Virgin biomass used to produce M
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

  Type Rec_Penal;
  Type Rec_Penal_Smooth;
  Type Rec_Penal_SumZero;
  Type Initial_pen;

  int Ipnt;                                                               // Pointer

  int NextraVar;
  NextraVar = 0;
  for (int Dcnt=0;Dcnt<NvarTypes;Dcnt++)
   {
   if (VarTypes(Dcnt) == 1) NextraVar += (BurnIn+Nyear+1);
   if (VarTypes(Dcnt) == 2) NextraVar += Narea*(BurnIn+Nyear+1);
   if (VarTypes(Dcnt) == 3) NextraVar += Narea*Nyear;
   if (VarTypes(Dcnt) == 4) NextraVar += Nzone*Nyear;
   }
  if (NextraVar==0) NextraVar = 1;

  vector<Type>VarOut(NextraVar);

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
  if (InitOpt==0)
   {
    Initial_pen = InitializeN(dataset, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, ActGrowth, RecruitFrac,
         Rbar, ActRecruitAreaSexDist, ActRecruitLenDist, ActRecDev, Ninit,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,
         VirginBio, CurrentBio, Mlow, Mvirgin, Minflection, Mtempts);
   }

  // Set up initial state (alternative)
  Type NtotalCheck; Type Nexpected;
  if (InitOpt==1)
   {
     Ninit = VirginN(dataset, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, Mvirgin, ActGrowth, RecruitFrac,
         Rbar, ActRecruitAreaSexDist, ActRecruitLenDist, ActRecDev, QRedsPar, MWhitesPar,Type(0.0));
    Ipnt = 0; Initial_pen = 0;
    for (int Isex=0;Isex<Nsex;Isex++)
     {
	  NtotalCheck = 0; Nexpected = 0;
      for (int Iarea=0;Iarea<Narea;Iarea++)
       for (int Iage=0;Iage<Nage;Iage++)
        for (int Isize=0;Isize<Nlen(Isex);Isize++)
         {  N(Iarea,BurnIn,0,Isex,Iage,Isize) = exp(Rbar)*exp(InitPars(Ipnt))  ;
            Initial_pen += 0.01*InitPars(Ipnt)*InitPars(Ipnt);
            NtotalCheck += N(Iarea,BurnIn,0,Isex,Iage,Isize)*WeightLen(Isex,Isize);
            Ipnt += 1; }

 	   for (int Iarea=0;Iarea<Narea;Iarea++)
	    for (int Iage=0;Iage<Nage;Iage++)
	     for (int Isize=0;Isize<Nlen(Isex);Isize++)
	      Nexpected += Ninit(Iarea,Isex,Iage,Isize)*WeightLen(Isex,Isize);

      //for (int Iage=0;Iage<=1000;Iage++) Nexpected += 0.5*exp(Rbar)*exp(-1*float(Iage)*M(0,0));
      //Initial_pen += (NtotalCheck-Nexpected)*(NtotalCheck-Nexpected);
      Initial_pen += WeightInitialN*(log(NtotalCheck)-log(Nexpected))*(log(NtotalCheck)-log(Nexpected));
     }
   }

  if(InitOpt==2)
   {
    Ipnt = 0; Initial_pen = 0;
    for (int Iarea=0;Iarea<Narea;Iarea++)
     for (int Isex=0;Isex<Nsex;Isex++)
      for (int Isize=0;Isize<Nlen(Isex);Isize++)
       {   N(Iarea,BurnIn-Nage,0,Isex,Nage-1,Isize) = exp(Rbar)*exp(InitPars(Ipnt))/float(Nsex)/float(Narea)/float(Nlen(Isex));
           if (Isize!=0) Initial_pen += 1.0*(InitPars(Ipnt)-InitPars(Ipnt-1))*(InitPars(Ipnt)-InitPars(Ipnt-1));
           //Initial_pen += InitPars(Ipnt)*InitPars(Ipnt);
           Ipnt += 1;  }
    IsVirgin = 0;
    for (int Iyear=-Nage;Iyear<0;Iyear++)
     for (int Istep=0;Istep<Nstep;Istep++)
      {
       XX = OneTimeStep(dataset, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, Mvirgin, Iyear, Istep, ActGrowth, RecruitFrac, Rbar,
                               IsVirgin, Feqn2, ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,
                               VirginBio, CurrentBio, Mlow, Mvirgin, Minflection, Mtempts);
      } // year and season
    }

  // Set up initial state (alternative)
  if (InitOpt==3)
   {
    Ipnt = 0; Initial_pen = 0;
    for (int Isex=0;Isex<Nsex;Isex++)
     for (int Iarea=0;Iarea<Narea;Iarea++)
      for (int Iage=0;Iage<Nage;Iage++)
       for (int Isize=0;Isize<Nlen(Isex);Isize++)
        {
        N(Iarea,BurnIn,0,Isex,Iage,Isize) = exp(LogRinitial(0))*exp(InitPars(Ipnt))  ;
         Initial_pen += WeightInit3*InitPars(Ipnt)*InitPars(Ipnt);
         // New weak penalty
         if (Isize!=0) Initial_pen += 1.0*(InitPars(Ipnt)-InitPars(Ipnt-1))*(InitPars(Ipnt)-InitPars(Ipnt-1));
         Ipnt += 1;
        }
    }

  // Set up initial state (virgin)
  if (InitOpt==4)
   {
    Initial_pen = 0;
    Ninit = VirginN(dataset, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, Mvirgin, ActGrowth, RecruitFrac,
            Rbar, ActRecruitAreaSexDist, ActRecruitLenDist, ActRecDev, QRedsPar, MWhitesPar,Finitial);
    for (int Isex=0;Isex<Nsex;Isex++)
     for (int Iarea=0;Iarea<Narea;Iarea++)
      for (int Iage=0;Iage<Nage;Iage++)
       for (int Isize=0;Isize<Nlen(Isex);Isize++)
        N(Iarea,BurnIn,0,Isex,Iage,Isize) = Ninit(Iarea,Isex,Iage,Isize)*exp(LogRinitial(Iarea))/exp(Rbar);
     }


   // Set up initial state (alternative)
  if (InitOpt==5)
   {
    Ipnt = 0; Initial_pen = 0;
    for (int Isex=0;Isex<Nsex;Isex++)
     for (int Iarea=0;Iarea<Narea;Iarea++)
      for (int Iage=0;Iage<Nage;Iage++)
       for (int Isize=0;Isize<Nlen(Isex);Isize++)
        {
         N(Iarea,BurnIn,0,Isex,Iage,Isize) = exp(LogRinitial(Iarea))*exp(InitPars(Ipnt))  ;
         Initial_pen += WeightInit3*InitPars(Ipnt)*InitPars(Ipnt);
         // New weak penalty
         if (Isize!=0) Initial_pen += 1.0*(InitPars(Ipnt)-InitPars(Ipnt-1))*(InitPars(Ipnt)-InitPars(Ipnt-1));
         Ipnt += 1;
        }
    }
  
  // Project the model forward
  IsVirgin = 0;                                                                        // Need to compute Fs
  for (int Iyear=0;Iyear<Nyear;Iyear++)
   for (int Istep=0;Istep<Nstep;Istep++)
    {
     XX = OneTimeStep(dataset, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, Mlow, Iyear, Istep, ActGrowth, RecruitFrac, Rbar,
                               IsVirgin, Feqn2, ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,
                               QRedsPar,MWhitesPar,VirginBio, CurrentBio, Mlow, Mvirgin, Minflection, Mtempts);
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
   
   
//  // Tagging data
//  Ntag.setZero(); RecapNum.setZero(); NotReported.setZero(); PredTagSize.setZero();
//  TagLike1.setZero();TagLike2.setZero();
//  for (int SexPass=0;SexPass<Nsex;SexPass++)
//   for (int GrpPass=0;GrpPass<NtagGroups;GrpPass++)
//    XX = TagDym(dataset,thedata, SexPass,GrpPass,N,ActSelex, ActReten, ActLegal,ActGrowth, ActMove, M, Hrate, QRedsPar,MWhitesPar,
//                Ntag,RecapNum,NotReported,TagLike1,TagLike2,PredTagSize);

 
  
  // Legal Biomass
  int Ipoint;int Ipoint76;
  LegalBio.setZero();
  LegalBio76.setZero();
  for (int Iyear=0;Iyear<Nyear;Iyear++)
   for (int Iarea=0;Iarea<Narea;Iarea++)
    {
     for (int Isex=0;Isex<Nsex;Isex++){
      for (int Iage=0;Iage<Nage;Iage++){
       for (int Istep=0;Istep<Nstep;Istep++)
        {
         Ipoint = LegalPnt(Isex,Iage,Iarea,Iyear,Istep);
         for (int Ilen=0;Ilen<Nlen(Isex);Ilen++){
           if(Istep<Nstep-1) {  // Take LB from start of subsequent tstep
             LegalBioTS(Iyear,Iarea,Istep) += LegalFI(Ipoint,Ilen)*N(Iarea,BurnIn+Iyear,Istep+1,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);
             LegalBio76TS(Iyear,Iarea,Istep) += LegalRef(Ilen)*N(Iarea,BurnIn+Iyear,Istep+1,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);}
           if(Istep==Nstep-1) {  // Take LB from start of tstep in the following year for all years - there are Nyears + 1 N values
             LegalBioTS(Iyear,Iarea,Istep) += LegalFI(Ipoint,Ilen)*N(Iarea,BurnIn+Iyear+1,0,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);
             LegalBio76TS(Iyear,Iarea,Istep) += LegalRef(Ilen)*N(Iarea,BurnIn+Iyear+1,0,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);}
           LegalBio(Iyear,Iarea) += LegalFI(Ipoint,Ilen)*N(Iarea,BurnIn+Iyear,Istep,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);
           LegalBio76(Iyear,Iarea) += LegalRef(Ilen)*N(Iarea,BurnIn+Iyear,Istep,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);}
        }}}
       LegalBio(Iyear,Iarea) /= float(Nstep);
       LegalBio76(Iyear,Iarea) /= float(Nstep);
    }
   
   // Legal Biomass including BurnIn
   int YearAdjusted;
   LegalBioAll.setZero();
   for (int Iyear=-BurnIn;Iyear<Nyear;Iyear++){
     if (Iyear <= 0) { YearAdjusted = 0; } else { YearAdjusted = Iyear; }
     for (int Iarea=0;Iarea<Narea;Iarea++)
     {
       for (int Isex=0;Isex<Nsex;Isex++){
         for (int Iage=0;Iage<Nage;Iage++){
           for (int Istep=0;Istep<Nstep;Istep++)
           {
             for (int Ilen=0;Ilen<Nlen(Isex);Ilen++){
               LegalBioAll(BurnIn+Iyear,Iarea) += LegalRef(Ilen)*N(Iarea,BurnIn+Iyear,Istep,Isex,Iage,Ilen)*WeightLen(Isex,Ilen);}
           }}}
       LegalBioAll(BurnIn+Iyear,Iarea) /= float(Nstep);
     }}
   
   // Simon's Cumulative catch reduced by average M based on time caught
  CumCatch.setZero();
  Type avM = M.sum()/(float(Nage)*float(Narea));    // Calculate Av M
  Type CnT = 0;
  
 //  AvM.setZero();
//  for (int Iyear=0;Iyear<Nyear;Iyear++){
//    CnT = 0.0;
//    for (int Iarea=0;Iarea<Narea;Iarea++){
//      for (int Iage=0;Iage<Nage;Iage++){
//        CnT += 1.0;
//        AvM(Iyear) += Mtempts(Iarea,Iage,Iyear);}}
//    AvM(Iyear) /= CnT; }
//  
 
  int Iarea;
  for (int Iyear=0;Iyear<Nyear;Iyear++){
    for (int Istep=0;Istep<Nstep;Istep++)   {
      
      for (int Ifleet=0;Ifleet<Nfleet;Ifleet++)  {
        Iarea = Fleet_area(Ifleet);
        CumCatch(Iyear,Iarea,Istep) += Catch(Iyear,Istep,Ifleet) ;  // Change recording of catch from fleet to area
      }
      if(Istep>0) { CumCatch(Iyear,Iarea,Istep) += CumCatch(Iyear,Iarea,Istep-1);}  // After all fleets then add previous time-steps cum-catch
      CumCatch(Iyear,Iarea,Istep) *= exp(-TimeStepLen(Iyear,Istep)*avM);   // Remove the M from the current time-step (previous catch plus new catch)
      //CumCatch(Iyear,Iarea,Istep) *= exp(-TimeStepLen(Iyear,Istep));                // Remove the M from the current time-step (previous catch plus new catch)
    }}

 // Simons Legal Biomass
  sLegalBio.setZero();
  sLegalBio76.setZero();
  for (int Iyear=0;Iyear<Nyear;Iyear++){
    for (int Iarea=0;Iarea<Narea;Iarea++)  {
      for (int Istep=0;Istep<Nstep;Istep++)   {
        sLegalBio(Iyear,Iarea) += LegalBioTS(Iyear,Iarea,Istep) + CumCatch(Iyear,Iarea,Istep); // Add Legal biomass with the cumulative catch to date
        sLegalBio76(Iyear,Iarea) += LegalBio76TS(Iyear,Iarea,Istep) + CumCatch(Iyear,Iarea,Istep);
      }
      sLegalBio(Iyear,Iarea) /= float(Nstep);
      sLegalBio76(Iyear,Iarea) /= float(Nstep);
    }}

  
  // Harvest rate
  Type TotalLB;   Type TotalLB76;
  Type CatchLB;
  for (int Iyear=0;Iyear<Nyear;Iyear++)
  {
   for (int Izone=0;Izone<Nzone;Izone++)
    {
     TotalLB = 0;TotalLB76 = 0;
     CatchLB = 0;
     for (int IareaP=0;IareaP<NareasPerZone(Izone);IareaP++)
      {
       Iarea = AreasPerZone(Izone,IareaP);
       if (Iarea >= 0)
        {
         TotalLB += sLegalBio(Iyear,Iarea);
         TotalLB76 += sLegalBio76(Iyear,Iarea);
         for (int Istep=0;Istep<Nstep;Istep++){
          for (int Ifleet=0;Ifleet<Nfleet;Ifleet++){
           if (Area_fleet(Iarea,Ifleet) == 1) {  CatchLB += Catch(Iyear,Istep,Ifleet);}
          }
         }
        }
       }
     HarvestRate(Iyear,Izone) = CatchLB/TotalLB;
     HarvestRate76(Iyear,Izone) = CatchLB/TotalLB76;
    }
   }
   
  //
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

  neglogL += CatchLike;
  neglogL += LambdaCpue*CpueLike;
  neglogL += LambdaNumbers*NumbersLike;
  neglogL += LambdaLength*LengthLike;
  neglogL += LambdaLarval*LarvalLike;
  
 
  // Variances of special outputs
  Ipnt = 0;
   for (int Ivar=0;Ivar<NvarTypes;Ivar++)
    {
     if (VarTypes(Ivar) == 1)
      for (int IvarYr=0;IvarYr<=BurnIn+Nyear;IvarYr++)
       {
        VarOut(Ipnt) = MatBio(IvarYr);
        Ipnt = Ipnt + 1;
       }
     if (VarTypes(Ivar) == 2)
      for (int IvarArea=0;IvarArea<Narea;IvarArea++)
       for (int IvarYr=0;IvarYr<=BurnIn+Nyear;IvarYr++)
        {
         VarOut(Ipnt) = RecruitmentByArea(IvarArea,IvarYr);
         Ipnt = Ipnt + 1;
        }
     if (VarTypes(Ivar) == 3)
      for (int IvarArea=0;IvarArea<Narea;IvarArea++)
       for (int IvarYr=0;IvarYr<Nyear;IvarYr++)
        {
         VarOut(Ipnt) = LegalBio(IvarYr,IvarArea);
         Ipnt = Ipnt + 1;
	    }
     if (VarTypes(Ivar) == 4)
      for (int IvarArea=0;IvarArea<Nzone;IvarArea++)
       for (int IvarYr=0;IvarYr<Nyear;IvarYr++)
        {
         VarOut(Ipnt) = HarvestRate(IvarYr,IvarArea);
         Ipnt = Ipnt + 1;
	    }
    }
  if (Ipnt == 0) VarOut(1) = dummy;


  // Now do projections
  if (DoProject==1)
  for (int Iyear=Nyear;Iyear<Nyear+Nproj;Iyear++)
   for (int Istep=0;Istep<Nstep;Istep++)
    {
     dataset.Catch(Iyear,Istep,0) = 10000;
     Catch(Iyear,Istep,0) = dataset.Catch(Iyear,Istep,0);
     XX = OneTimeStep(dataset, N, Z, Hrate, ActSelex, ActReten, ActLegal, ActMove, WeightLen, M, Iyear, Istep, ActGrowth, RecruitFrac, Rbar,
                               IsVirgin, Feqn2, ActRecruitAreaSexDist, ActRecruitLenDist,ActRecDev,MatBio,MatBioArea,RecruitmentByArea,BiasMult,SigmaR,QRedsPar,MWhitesPar,
                               VirginBio, CurrentBio, Mlow, Mvirgin, Minflection, Mtempts);
    } // year and season


  if (DoProject==0 || DoProject==1)
   {
    REPORT(Hrate);
    REPORT(Catch);
    REPORT(MatBio);
    REPORT(MatBioArea);
    REPORT(RecruitmentByArea);
    }

  if (DoProject==0)
   {
    REPORT(N);
    REPORT(CatchCheck);
	REPORT(ActSelex);
	REPORT(ActReten);
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
	REPORT(PredCpue);
	ADREPORT(PredCpue.col(0));
	REPORT(SigmaCpue);
	REPORT(CpueQ);
	REPORT(PredNumbers);
	ADREPORT(PredNumbers.col(0));
	REPORT(SigmaNumbers);
	REPORT(PredLengthComp);
	REPORT(SigmaNumbers);
	REPORT(PredLarval);
	ADREPORT(PredLarval.col(0));
	REPORT(Initial_pen);
	REPORT(Rec_Penal);
	REPORT(Rec_Penal_Smooth);
	REPORT(Ninit);
  REPORT(PuerulusByArea);
  REPORT(sLegalBio);
  REPORT(HarvestRate);
  REPORT(LegalBio76);
  REPORT(HarvestRate76);
  REPORT(sLegalBio76);
  REPORT(CumCatch);
  REPORT(Weighted_CpueLike);
  REPORT(Weighted_NumbersLike);
  REPORT(Weighted_LengthLike);
  REPORT(Weighted_LarvalLike);

    REPORT(LegalBioAll);
    REPORT(ActRecruitAreaSexDist);
    REPORT(ActRecruitLenDist);
    REPORT(ActRecDev);
    REPORT(Recruits);
    REPORT(neglogL);
    REPORT(Z);
    REPORT(Select);
    ADREPORT(VarOut);
    REPORT(RecruitFrac);
    REPORT(RecruitPars);
    REPORT(NtotalCheck);
    REPORT(Nexpected);
    REPORT(BiasMult);
    REPORT(LogRinitial);
    REPORT(VirginBio);
    REPORT(CurrentBio);
    REPORT(Mvirgin);
    REPORT(Mtempts);
    REPORT(M);
    REPORT(GrowthOut);
    REPORT(CpueEcreep);
    }

    return neglogL;
}

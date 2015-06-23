  
// Laser.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


/*This program is for obtain GE,Te, etc as function on t for give Vs,qL,KL, Rs,t */
#include <math.h>
#include <stdio.h>
#include <conio.h>

 const double
	 PI = 3.14159265359,

/*B   D=1/2E-5,  LT=0.0228,   LA=0.0398,  LS=11.99E3,  AS=13.07,
         AI=10.81,         UI=8.3,     FE=4.5,  C=9.36,   BS=2.962E4,*/
/*Zn  D=1/3.574E-5,  LT=0.27,   LA=0.408,  LS=430,  AS=11.63,
         AI=65,         UI=9.4,     FE=4.24,  C=8.3,   BS=6.56E3,*/
/*Cd  D=1/6.73E-6,  LT=0.23,   LA=0.45,  LS=213,  AS=11.56,
         AI=112,        UI=9,     FE=3.7,    C=8.35,   BS=5.72E3,*/
/*Si  D=1/2.46E-3/0.034,    LT=0.0624,    LA=0.116,    LS=3380,  AS=10.78,
        AI=28.1,            UI=8.15,    FE=4.8,  C=9.21,   BS=2.09E4, */
/*Cr  D=1/1.89E-5,  LT=/*0.21  0.161,   LA=/*0.3  0.206, LS=1400, AS=12.94,
        AI=52,         UI=6.8,      FE=4.7,    C=9.56,    BS=2.0E4,*/
/*Ti  D=1/4.2E-5,  LT=0.037,   LA=0.07,  LS=1750,  AS=12.5,
       AI=48,        UI=6.8,     FE=4,         C=9.11,     BS=2.323E4, */
/*Cu  D=1/2.46E-5,  LT=0.72,  LA=1.14, LS=1140, AS=11.96, BS=1.698E4, C=8.63,
        AI=64,  U1=7.72,  FE=4.5,  /*VS=1E6, U2=20.29, U3=36.83, U4=55.2, 
		/*U5=79.9 U12=U1+U2, U13=U1+U2+U3, U14=U1+U2+U3+U4,*/
/*Al D=1/8E-6,  LT=0.15/*0.497, LA=0.2675/*0.893, LS=2516, AS=11.79, BS=1.594E4, C=8.27,
        AI=27,  U1=5.986,  FE=4.25,  /*VS=1E6, U2=18.826, U3=28.448, U4=119.99, 
		/*U5=153.75 U12=U1+U2, U13=U1+U2+U3, U14=U1+U2+U3+U4,*/
/*Ag D=1/4E-6, LT=0.014/*1.003, LA=0.0242/*1.737, LS=213, AS=11.85, BS=14.27E3, C=8.63,
        AI=108,  U1=7.574,  FE=4.3,  /*VS=1E6, U2=21.48, U3=36.1, U4=52, /*U5=70
		U12=U1+U2, U13=U1+U2+U3, U14=U1+U2+U3+U4,*/
/*Au D=1/3.7E-5, LT=0.744, LA=1.24, LS=410, AS=11.89, BS=17.58E3, C=8.8,
        AI=197,  U1=9.225,  FE=4.4,   U2=20.5, U3=30, U4=44, /*U5=58
		U12=U1+U2, U13=U1+U2+U3, U14=U1+U2+U3+U4,*/
/*Ni*/ D=1/6.02E-6, LT=0.132, LA=0.141, LS=1516, AS=12.75, BS=20.96E3, C=9.4,
        AI=59,  U1=7.635,  FE=4.8,  U2=18.168, U3=35.17, U4=54.9, /*U5=75.5*/
		U12=U1+U2, U13=U1+U2+U3, U14=U1+U2+U3+U4,
/*C    D=1/1E-3,    LT=0.06,     LA=0.163,    LS=15000,  AS=15.73,
        AI=12,     UI=11.2,     FE=4.34,   C=12.04,     BS=4.03E4,*/       
   U01=1E-8, W01=1E-3, T01=1E-7, A01=5E-3, D02=1E-8, D01=1E-8, Ne1=1E-8,
   E1=AI*1845;    /*  EL0=0.003,      GE0=1E-2,*/   
 /*-------------------------------------------------------------------*/
  void Iter ( double& f1, double& f2, double &x1, double &x2, double& step, bool& Key )
  {
    /* variable */
      double oldx;
    
      oldx=x1;
      if ( ((f1>1) && (f2<1)) || ((f1<1) && (f2>1)) )
		  step=step/2;

      if ( oldx>x2 )
	  {
                   if ( ((f1>1) && (f2>1) && (f1<f2))
                      || ((f1<1) && (f2<1) && (f1>f2)) )
                                          x1=x1+step;
                                        else
                                          x1=x1-step;
	  }
      if ( oldx<x2 )
	  {
                   if ( ((f1>1) && (f2>1) && (f1<f2))
                      || ((f1<1) && (f2<1) && (f1>f2)) )
                                          x1=x1-step;
                                        else
                                          x1=x1+step;
	  }
      x2=oldx;
      f2=f1;
      if ( x1 == x2 ) Key=true;
  }
 /*-------------------------------------------------------------------*/
/*Nije vmesto Procedure Wer (var B : Real; var FF : Real );*/
 /*---------------------------------------------------------------------*/
  double f ( double x )
  {
    const double G = 1;
    /* variable */
      double R;

      R = pow(G,2)*pow(x,2);
      if ( R>70 )
	  {
			return 0;
	  }
      else
	  {
            return (G*2/sqrt(PI))*exp(-R);
	  }
  }

	int odd( int n)
	{
		return n % 2 ;
	}

	/*
	double abs( double f )
	{
		if ( f >=0 )
			return f;
		else
			return f * -1;
	}
	*/

  void Wer ( double b, double& Integral )
  {
    const float a = 0;
    const int n = 14;
    /* variable */
      double h, x;
	  double Sum;
      int i;
    
      h =  (b-a)/n;
      Sum = 0;
      x = a;
      for ( i=1; i < n; i++ )
	  {
          x = x+h;
          if ( odd (i) )
		  {
                     Sum = Sum + 4*f(x);
		  }
          else
		  {
                     Sum = Sum + 2*f(x);
		  }
      }
      Integral = ( Sum + f(a) + f(b) ) * h/3;
  }
 /*---------------------------------------------------------------------*/
/* Var
   h, m, s, hund : Word;
    Integral,   T, E : Real;}
  function W ( x : Real) : Real;
    var
      R, Y,YY, R1, R2 : Real;
    begin
      R:= x/8.625E-5/T0;
      If R>70 then
                R2:=0
              else
                R2:= Ln(1+Exp(-R));
      If FE<x then
      W:=R2
      else
      Y:=3.79E-4*Sqrt(EE)/(FE-x);
      If Y>1 then
      W:=R2
      else
      YY:=1-Y*Y*(1+0.85*Sin((1-Y)/2));
      If FE<x then
      W:=R2
      else
      R1:=6.85E7*Ln(1.5*Exp(FE-x))*YY/EE;
      If R1>70 then
                W:=R2*Exp(-70)
                else
                W:=R2*Exp(-R1)
    end;
  procedure Emis (EE,T  : Real; var Integral : Real ) ;
    const
      a = -5;
      b =15;
      n = 200;
    var
      h, x, Sum : Real;
      i         : Integer;
    begin
      h :=  (b-a)/n;
      Sum := 0;
      x := a;
      for i :=1 to n-1 do
        begin
          x := x+h;
          if odd (i) then
                     Sum := Sum + 4*W(x)
                   else
                     Sum := Sum + 2*W(x)
        end;
      Integral := ( Sum+ W(a)+ W(b) ) * h*T*1.39E6/3;
    end;*/
  /*---------------------------------------------------------------------*/
 
  int main(int argc, char* argv[])
  {
  /* Variable */
double WW1,WW2,W1,W2,WE,WE3,HW,W3,VT,BJ,k,k1,k2,k3,k4,k5,k6,k7,k8,Ke,P,N0,Ne,
       Na,Z,KW,N1,N2,N3,N4,Tw,T00,W0,EX,UW,UG,E2,EE,FH,UE,T0,TE,B,E0,GE,Ker,
       R0,R10,R20,R21,B1,B2,B10,B20,BB,HB1,HB2,T10,T20,T12,D1,D2,D3,D4,D5,DNe,
	   Ne2,Ne3,N03,D30,DD2,DD3,DD4,DD5,FF1,V2,f1,f2,f3,f4, HNe,Vet,WZ,T1,T02,T03,TT,HT,
       A0,A02,A03,AI0,AE3,HA,A1,A2,A3,A4,AI1,AI2,AI3,AI4,N11,N22,N33,N44,Naa,Jet,Je,Ji,
 	   g0,g1,g2,g3,g4,UK,U33,Uk1,Uk2,Uk3,HUk,KL,qL,UEE,Rs,TS,EI2,EI3,EI4,VS,Gi;

    bool Key;


	k=0; k1=0; k2=0; k3=0; k4=0; k5=0; k6=0; k7=0; k8=0; g0=2; g1=1; g2=1; g3=1; g4=1;
            
    Tw=3E-9; P=sqrt(PI); T00=300; Ke=0.02; KL=0.97; qL=2E9; /*Rs=0.056419;*/ Rs=0.0356825;
    
	T1=8950;     HT=50;        T0=T1+HT;       
	Ne3=1E19;      HNe=1E15;     Ne=Ne3+HNe;
    AE3=0.019;  HA=0.000001;     AI0=AE3+HA;
    WE3=1.61;     HW=0.0001;     WE=WE3+HW;   
    U33=5;         HUk=0.5;     UK=U33+HUk;
    B10=0.012;     HB1=0.001;    B1=B10+HB1;
    B20=0.014;     HB2=0.001;    B2=B20+HB2;
    Key=false;
     
  	HA=0.00005;
		do
		{
			    k5=k5+1;
		HW=0.0005;
	      do
		{
			k=k+1;
            B=3.29E+21*exp(1.5*log(WE));
            A1=B*g1/g0*exp(-1.5*U1/WE);
            A2=B*g2/g1*exp(-1.5*U2/WE);
            A3=B*g3/g2*exp(-1.5*U3/WE);
			A4=B*g4/g3*exp(-1.5*U4/WE);
           
                               HNe=3E16;
				do
				{
                       k8=k8+1;
					   Na=Ne*(1-AI0)/AI0;
                       N1=A1*(1-AI0)/AI0;
                       N2=A2*N1/Ne;
                       N3=A3*N2/Ne;
					   N4=A4*N3/Ne;
                       DNe=(N1+2*N2+3*N3+4*N4)/Ne;
                       if (k8==1)  N03=DNe;
                     //	printf("   DNe=%e     Ne=%e ",DNe, Ne); 
                       Ne2=abs(1-DNe);
                       if (Ne2>=Ne1)  Iter (DNe,N03,Ne,Ne3,HNe,Key);
				}
                       while (Ne2>Ne1);
                          if (Key)
               	          printf ("   ERR|| !!! : x1=x2");
                 BJ=N1+N2+N3+N4+Na;
				
                        AI1=N1/(Na+Ne);  AI2=2*N2/(Na+Ne);  AI3=3*N3/(Na+Ne); AI4=4*N4/(Na+Ne);
                        N11=N1/BJ;  N22=N2/BJ;  N33=N3/BJ;  N44=N4/BJ; Naa=Na/BJ;
                        f1=N1/Ne;  f2=2*N2/Ne;   f3=3*N3/Ne; f4=4*N4/Ne;
                        TE=WE/1.5;
                     //  VT=14549.75*sqrt(T0/AI);                         
						VT=sqrt(2*TE/1E-12/AI);
                       	Ji=/*VT*Ne*1.6E-19/4*/0.4*VT*Ne*1.6E-19;
                        Vet=5.4937E7*sqrt(WE);
				                         HUk=0.01;
                                         do
									{
                                         k1=k1+1;
										 KW=1.5*UK/WE;
                                         EX=exp(KW);
                                         Jet=Vet*Ne*1.6E-19/EX/4;
                                         E0=sqrt(1+0.5/KW)-sqrt(0.5/KW)-sqrt(1/KW)*(1-1/EX);
										 EI2=sqrt(1+0.25/KW)-sqrt(0.25/KW)-sqrt(1/KW)*(1-1/EX);
										 EI3=sqrt(1+0.5/KW/3)-sqrt(0.5/KW/3)-sqrt(1/KW)*(1-1/EX);
										 EI4=sqrt(1+0.125/KW)-sqrt(0.125/KW)-sqrt(1/KW)*(1-1/EX);
                                         E2=7.570E5*Ji*sqrt(UK)*((exp(0.5)*sqrt(E1)*(f1*E0+f2*EI2+f3*EI3+f4*EI4))+1-Jet/Ji);
                                        // if (E2<=0) 
									//		Je=120*T0*T0*exp(-FE/8.625E-5/T0)-Jet;	
                                         if (E2>0)   
											 EE=sqrt(E2);
											 if (E2>0)
                                         FH=3.79E-4*sqrt(EE);
                                    //     if (E2>0)
									//	 if (FH>=FE)  Je=120*T0*T0-Jet;
										 if (E2>0)
									    /* if (FH<FE)*/  Je=120*T0*T0*exp(-(FE-FH)/8.625E-5/T0)/*-Jet*/;
                                    //     J=Je+Ji;
                                     //  Uk1=J*P*P*H2*(2*H1+H2)/I;
									   Uk1=TE*log(1/(0.6*sqrt(2*PI/E1)+Je*sqrt(2*PI*9E-28/1.6E-12/TE)/1.6E-19/Ne))/UK;
                                       // EMIS (EE,T0,TFE);
                                      // S1=TFE/S/J-SS2;
	                                  // printf ("J=%e Jet=%e Je=%e Ji=%e T0=%e",J,Jet,Ji,T0);
                                      //   getch();
									    if (k1==1)  Uk3=Uk1;
      	                        //     printf (" UK=%e  Uk1=%e ", UK,Uk1);
                                         Uk2=abs(1-Uk1);
                                        if (Uk2>=U01)  Iter (Uk1,Uk3,UK,U33,HUk,Key);
				} 
                                        while (Uk2>U01);
              HT=50;
				   do
	{
		k1=k1+1;
		W0=1.29E-4*T0;
		N0=exp((AS-(BS/T0))*log(10.0))*1.333/1.38E-16/T0;
        
                      HB1=0.01;
                        do
					{
                        k2=k2+1;
                        Wer(B1,FF1);
                        R0=1+FF1;
                        R10=1/R0;
                        R21=BJ*R0/N0;
                        BB=B1*B1;
                        D1=exp(-BB)-(1-FF1)*B1*P;
                        T10=pow(R0/(2*P*B1+D1),2);
                        HB2=0.01;
						     do
						{                                        
							  k3=k3+1;
                              T12=pow(R21*B2/B1,2);
                              D2=((BB+0.5)*(1-FF1)-B1*exp(-BB)/P);
                              D30=B1*(1)*(1+2*B2*B2)/sqrt(T12)/B2;
                              D3=(R0/2/T10+D2)/D30;
                              if (k3==1)  D5=D3;
	                       // printf("   D3=%e     B2=%e ",D3, B2); 
                              D4=abs(1-D3);
                              if (D4>=D01)  Iter (D3,D5,B2,B20,HB2,Key);
						}
                              while (D4>D01);
                        DD2=(P*B1*(2.5+BB)*(1-FF1)-(2+BB)*exp(-BB))/2;
                        DD3=((R0/exp(1.5*log(T10))+DD2)/(1)/P/B1/(B2*B2+2.5))*T12;
                        if (k2==1)  DD5=DD3;
            //          printf("   DD3=%e     B1=%e ",DD3, B1); 
                        DD4=abs(1-DD3);
                       if (DD4>=D02)  Iter (DD3,DD5,B1,B10,HB1,Key);
					} 
					    while (DD4>D02);
                        T20=T10/T12;
                        V2=B2/sqrt(AI*1.66E-24/2/1.38E-16/T0/T20);         /*/T20*/   
                GE=PI*Rs*Rs*BJ*V2*AI*1.66E-24;
                Z=1E5*GE/AI/PI/Rs/Rs;
				UG=(LS+(3.2E4*W0/AI))*GE/0.24/PI/Rs/Rs;
   UW=4*Jet*WE/3;
   UE=(1-KL)*qL+Ji*(UK+f1*U1+f2*U12+f3*U13+f4*U14)+UW-UG;
   TS=2*atan(sqrt(4*LA*Tw/Rs/Rs))/PI;
   UEE=2*LT*(T0-300)/TS/P/Rs/0.24;
   TT=UE/UEE;
  
                          if (k1==1)  T03=TT;
       //   	printf("  TT=%e T0=%e ",TT, T0); 
                       T02=abs(1-TT);
                      if (T02>=T01)  Iter (TT,T03,T0,T1,HT,Key);
				}
                       while (T02>T01);
            WZ=Z*(1+N22+2*N33+3*N44-Naa);
            WW1=Ke*KL*qL+Je*(UK+2*8.625E-5*T0)-WZ*(f1*U1+f2*U12/2+f3*U13/3+f4*U14/4)-Ji*(f1*U1+f2*U12+f3*U13+f4*U14);
            WW2=2*WZ*WE+Jet*(4*WE/3+UK);
            W1=WW1/WW2; 
            if (k==1)  W3=W1;
           printf (" W1=%e WE=%e T0=%e UK=%e B2=%e Ji=%e ",W1, WE, T0,UK,B2,Ji); 
            W2=abs(1-W1);
            if (W2>=W01) Iter (W1,W3,WE,WE3,HW,Key);
		}
            while (W2>W01);(Key==false);
                                 if (Key) 
             	                   printf ("   ERR|| !!! : x1=x2");
			   VS=sqrt((2*(1-Ke)*KL*qL*PI*Rs*Rs+4*WZ*TE*PI*Rs*Rs)*1E7/GE+V2*V2);
               A0=BJ*1.38E-16*T0*(1+AI0*WE/W0/T20)*T20*PI*Rs*Rs/GE/VS;
               if (k5==1)  A03=A0;
	          printf ("    A0=%e  AI0=%e WE=%e T0=%e",A0, AI0, WE,T0);
               A02=abs(1-A0);
               if (A02>=A01)  Iter (A0,A03,AI0,AE3,HA,Key);
		}
                while (A02>A01);
                            if (Key) 
          	                printf ("   ERR|| !!! : x1=x2");
   Ker=GE*sqrt(T0)/exp((C-(BS/T0))*log(10.0))/PI/Rs/Rs;
   R20=R10*R21;
   Gi=GE*Tw/1e-6;
// printf ("  J=%e ' H1=%e ' Te=%e ' Tw=%e ",J,H1,Te,Tw);
   	/* Write to a file */
 FILE *stream;
 stream = fopen( "dataC.dat", "w" );

  fprintf( stream, " UK=%.2f Rs=%.2f T00=%.2f Tw=%e qL=%e KL=%.2f Ke=%.3f T0=%.2f TE=%.4f N0=%e EE=%e AI0=%.4f B1=%.4f B2=%.4f R20=%.4f T10=%.4f",   
		  		     UK,     Rs,     T00,     Tw,   qL,   KL,     Ke,     T0,     TE,     N0,   EE,   AI0,     B1,     B2,     R20,     T10);    
  fprintf( stream," R10=%.4f T20=%.4f T12=%.4f Ker=%.4f R21=%.4f V2=%.2f GE=%e Gi=%e UW=%e UG=%e UE=%e Z=%e WZ=%e VS=%e",  
                    R10,     T20,     T12,     Ker,     R21,     V2,     GE,   Gi,   UW,   UG,   UE,   Z,   WZ,   VS);
  fprintf( stream," Je=%e Ji=%e Jet=%e f1=%e f2=%e f3=%e f4=%e Ne=%e N1=%e N2=%e N3=%e N4=%e Na=%e AI1=%.5f AI2=%e AI3=%e AI4=%e N11=%e N22=%e N33=%e N44=%e",
                    Je,   Ji,   Jet,   f1,   f2,   f3,   f4,   Ne,   N1,   N2,   N3,   N4,   Na,   AI1,     AI2,   AI3,   AI4,   N11,   N22,   N33,    N44);
  fclose(stream);	

  /* Write on the screen */
	printf("\n");
	printf("\n");printf("\n");
	printf(" FINAL SOLUTION\n");
	printf("\n");

	printf( " UK=%.2f Rs=%.2f T00=%.2f Tw=%e qL=%e KL=%.2f Ke=%.3f T0=%.2f TE=%.4f N0=%e EE=%e AI0=%.4f B1=%.4f B2=%.4f R20=%.4f T10=%.4f",   
		      UK,     Rs,     T00,     Tw,   qL,   KL,     Ke,     T0,     TE,     N0,   EE,   AI0,     B1,     B2,     R20,     T10);    
    printf( " R10=%.4f T20=%.4f T12=%.4f Ker=%.4f R21=%.4f V2=%.2f GE=%e Gi=%e UW=%e  UG=%e UE=%e Z=%e WZ=%e VS=%e",  
              R10,     T20,     T12,     Ker,     R21,     V2,     GE,   Gi,   UW,    UG,   UE,   Z,   WZ,   VS );
    printf( " Je=%e Ji=%e Jet=%e f1=%e f2=%e f3=%e f4=%e Ne=%e N1=%e N2=%e N3=%e N4=%e Na=%e AI1=%.5f AI2=%e AI3=%e AI4=%e N11=%e N22=%e N33=%e N44=%e",
              Je,   Ji,   Jet,   f1,   f2,   f3,   f4,   Ne,   N1,   N2,   N3,   N4,   Na,   AI1,     AI2,   AI3,   AI4,   N11,   N22,   N33,    N44 );
  
   getch();
   	    
	/* Wait for user input before closing console window */
	getch();
	return 0;
  }
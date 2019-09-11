#include <math.h>
#include <stdio.h>
#include <Rcpp.h>
#include "FunctionsRandom.h"
#include "RandomNumbers.h"
double RanMT(void)
{
  return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
}
//-------------------------------------------------------------------------
//Genera una Distribucion Uniforme con parametros a y b
//-------------------------------------------------------------------------
double UnifMT(Rcpp::NumericVector params)
{
  double uni;
  double a = params[0];
  double b = params[0];
  uni = RanMT();
  uni = a + (b - a) * uni;
  return uni;
}

double Unif2MT(double a, double b)
{//Se usa para llamadas internas
  double uni;
  uni = RanMT();
  uni = a + (b - a) * uni;
  return uni;
}
//-------------------------------------------------------------------------
//Genera una Distribucion Exponencial con esperanza beta
//-------------------------------------------------------------------------
double ExpoMT(Rcpp::NumericVector params)
{
  double beta = params[0];
  double uni;
  uni = RanMT();
  return (-beta * log(uni));
}
//-------------------------------------------------------------------------
//Genera una Distribucion normal con media media y desvest des
//-------------------------------------------------------------------------
double NormMT(Rcpp::NumericVector params)
{//Algoritmo de Marsaglia Bray
  double media = params[0];
  double des   = params[1];
  double uni,v,w,x,sum;
  x=0;
  uni = RanMT();
  if(uni<=0.8638)
  {
    v=Unif2MT( -1.0, 1.0 );
    w=Unif2MT( -1.0, 1.0 );
    x=2.3153508*uni-1+v+w;
    return (media+ des*x);
  }
  else
  {
    if(uni<=0.9745)
    {
      v= RanMT();
      x=(3/2)*(v-1+9.0334237*(uni-0.8638));
      return (media+ des*x);
    }
    else
    {
      if(0.9973002<uni)
      {
        do
        {
          v= RanMT();
          w= RanMT();
          x= 9/2-log(w);
        }
        while(x*pow(v,2)<=9/2);
        if(uni>=0.9986501)
        {
          x= sqrt(2*x);
          return (media+ des*x);
        }
        else
        {
          x=-sqrt(2*x);
          return (media+ des*x);
        }
      }
      else //0.9745<uni<=0.9973002
      {
        do
        {
          x=Unif2MT(-3,3);
          uni= RanMT();
          v= fabs(x);
          w=6.6313339*(3-v)*(3-v);
          sum=0;
          if(v<3/2)
          {
            sum=6.0432809*(3/2-v);
          }
          if(v<1)
          {
            sum=sum+13.2626678*(3-pow(v,2))-w;
          }

        }
        while(uni<=49.0024445*exp(-pow(v,2)/2)-sum-w);
        return (media+ des*x);
      }

    }
  }
}
double Norm2MT(double media,double des)
{//Se usa para llamadas internas
  double uni,v,w,x,sum;
  x=0;
  uni = RanMT();
  if(uni<=0.8638)
  {
    v=Unif2MT( -1.0, 1.0 );
    w=Unif2MT( -1.0, 1.0 );
    x=2.3153508*uni-1+v+w;
    return (media+ des*x);
  }
  else
  {
    if(uni<=0.9745)
    {
      v= RanMT();
      x=(3/2)*(v-1+9.0334237*(uni-0.8638));
      return (media+ des*x);
    }
    else
    {
      if(0.9973002<uni)
      {
        do
        {
          v= RanMT();
          w= RanMT();
          x= 9/2-log(w);
        }
        while(x*pow(v,2)<=9/2);
        if(uni>=0.9986501)
        {
          x= sqrt(2*x);
          return (media+ des*x);
        }
        else
        {
          x=-sqrt(2*x);
          return (media+ des*x);
        }
      }
      else //0.9745<uni<=0.9973002
      {
        do
        {
          x=Unif2MT(-3,3);
          uni= RanMT();
          v= fabs(x);
          w=6.6313339*(3-v)*(3-v);
          sum=0;
          if(v<3/2)
          {
            sum=6.0432809*(3/2-v);
          }
          if(v<1)
          {
            sum=sum+13.2626678*(3-pow(v,2))-w;
          }

        }
        while(uni<=49.0024445*exp(-pow(v,2)/2)-sum-w);
        return (media+ des*x);
      }

    }
  }
}
double LogNMT(double &media,double &des)
{//Genera LogNormal con parámetros media y des
  return exp(Norm2MT(media, des));
}
double TriaMT(const double &a, const double &b, const double &c1)
{// Genera distribución triangular entre a, mode, b, donde c1 = (mode - a) / (b - a)
  double uni, tria01;
  uni = RanMT();
  //Se genera tria(0,1, c1)
  if (uni < c1)
  {
    tria01 = sqrt(uni * c1);
  }
  else
  {
    tria01 = 1 - sqrt((1 - uni) * (1 - c1));
  }
  return (a + (b - a) * tria01);
}

long  BernMT( double &p)
{//Genera Bernoulli con parámetro p
  double uni;
  //Se supone que 0 <= p <= 1
  uni = RanMT();
  if (uni < p) return 1;
  else           return 0;
}

long  BinoMT(double &p,long &n)
{//Genera una variable binomial con parámetros n y p como suma de Bernoullis
  long i,bin=0;
  for (i = 1;i<n+1;i++)
  {
    bin += BernMT( p);
  }
  return bin;
}
long  GeomMT(double &lnp)
{//Genera Geométrica con parámetro p donde lnp es ln(1-p)
  double uni;
  //Se supone que 0 < p < 1
  uni = RanMT();
  return int(log(uni) / lnp);
}

double Gamm1MT(double &c, double &d , double &p,double &q, double &r,double &alfa,double &beta)
{//Generador de Cheng (para alfa > 1)
  double V, y, z, w, U1, U2;
  U1 = RanMT();
  U2 = RanMT();
  V = c * (log(U1 / (1 - U1)));
  y = alfa * exp(V);
  z = (U1*U1) * U2;
  w = d + (p * V) - y;
  if ((w + r - (q * z)) >= 0)
  {
    return (beta * y);
  }
  else
  {
    if (w >= log(z))
    {
      return (beta * y);
    }
    else
    {
      return Gamm1MT( c, d, p, q, r, alfa, beta);
    }
  }
}

double Gamm2MT(double &alfa,double &beta,double &b)
{
  double uni, y, p, w;
  uni = RanMT();
  p = b * uni;
  if (p > 1)
  {
    y = -log((b - p) / alfa);
    w = pow(y,  alfa - 1);
  }
  else
  {
    y =pow( p,1 / alfa);
    w= exp(-y);
  }
  uni = RanMT();
  if (uni <= w)
  {
    return (beta * y);
  }
  else
  {
    return Gamm2MT(alfa, beta, b);
  }
}

double GammMT(double &alfa,double &beta)
{
  double c, d, p, q, r, b;
  if (alfa > 1)
  {
    c = 1 / (pow(2 * alfa - 1, 0.5));
    d = alfa - log(4.0);
    p = alfa + 1 / c;
    q = 4.5;
    r = 1 + log(4.5);
    return Gamm1MT( c, d, p, q, r, alfa, beta);
  }
  else
  {
    if (alfa < 1)
    {
      b = (exp(1.0) + alfa) / exp(1.0);
      return Gamm2MT(alfa, beta, b);
    }
    else
    {
      return ExpoMT(beta);
    }
  }
}

double WeibMT(double &beta,double &alfa)
{//Genera distribución de Weibull
  double uni;
  uni = RanMT();
  return (beta * pow(-log(uni),1 / alfa));
}

double ErlaMT(double &beta,long &m)
{//Genera distribución de Erlang
  double uni;
  long i;
  double temp;
  temp = 1;
  for (i = 1; i<m+1; i++)
  {
    uni = RanMT();
    temp = temp * uni;
  }
  temp = log(temp);
  return (-beta * temp);
}
long  MaChMT(long &Ini, long alias[], double q[], long &n)
{// Generación de una muestra de una Cadena de Markov con n estados usando método alias
  // Los arreglos alias y q deben haberse generado previamente con genera_alias
  double uni;
  long indice1, indice2;
  uni = RanMT(); indice1 = n * uni;
  indice2 = (Ini * n) + indice1;
  uni = RanMT();
  if (uni <=q[indice2])
  {
    return indice1;
  }
  else
  {
    return alias[indice2];
  }
}
void genera_alias(long &n, double p[], long alias[],double q[], long small1[], long big[])
{//Los arreglos deben estar correctamente dimensionados, de otra forma se comete error grave
  //Notar que arreglos empiezan en indice 1
  long i, ns, nb, k, j; double suma = 0.0;
  for (i = 1; i < n+1; i++) suma += p[i]; //Primero se calcula la suma para repartirla en partes iguales
  ns = 0; nb = 0;
  for (i = 1; i<n+1;i++)
  { q[i] = n * p[i];
    if (q[i] < suma)
    { ns = ns + 1; small1[ns] = i;
    }
    else
    { nb = nb + 1; big[nb] = i; alias[nb] = nb;
    }
  }
  while (ns > 0)
  { k = big[nb]; j = small1[ns]; alias[j] = k;
    q[k] = q[k] + q[j] - suma; ns = ns - 1;
    if ((q[k] < suma) && (nb > 0)) //if (q[k] < 1)
    { nb = nb - 1; ns = ns + 1; small1[ns] = k;
    }
  }
}
double AliaMT(double x[], long alias[], double q[], long &n)
// Generación de discreta tomando n valores usando el método alias
// Los arreglos alias y q deben haberse generado previamente con genera_alias
{ double uni1, uni2; long indice;
  uni1 = RanMT(); uni2 = (double)n * uni1; indice = floor(uni2);
  uni2 = uni2 -(double)indice; indice += 1;
  if (uni2 <=q[indice])
  {
    return x[indice];
  }
  else
  {
    return x[alias[indice]];
  }
}
double EmpMT(double x[], long &n)
{// Generación de un valor aleatorio de x (distribución empírica)
  return x[(long)floor((double)n * RanMT())];
}

  long  PoisMT(double &tasa, double &exptasa)
  // Generación de muestra de VA Poisson usando transformación inversa
  // El valor exptasa es exp(-tasa) y debe calcularse previamente
  { double p, Suma, uni; long k;
    k = 0; uni = RanMT(); p = exptasa; Suma = p;
    while (uni > Suma)
    { k += 1; p = p* tasa / (double)k; Suma += p;
    }
    return k;
  }
  bool  SetPois(double &Tasa, long &n, double q[], long alias[],
    double &suma, double &ProbN, double p[], long peque[], long grande[])
    // Calcula los parámetros necesarios para aplicar el método Alias en la generación de una Poisson
    {
      double p1; long i, ns=0, nb=0, k, j;
      p1 = exp(-Tasa); p[1] = p1; suma = p1;
      for (i = 1; i < n; i++)
      { p1 = p1 * Tasa / (double)i; p[i+1] = p1; suma += p1;
      }
      ProbN = p[n];
      for (i = 1; i<n+1;i++)
      { q[i] = n * p[i];
        if (q[i] < suma)
        { ns = ns + 1; peque[ns] = i;
        }
        else
        { nb = nb + 1; grande[nb] = i; alias[nb] = nb;
        }
      }
      while (ns > 0)
      { k = grande[nb]; j = peque[ns]; alias[j] = k;
        q[k] = q[k] + q[j] - suma; ns = ns - 1;
        if ((q[k] < suma) && (nb > 0))
        { nb = nb - 1; ns = ns + 1; peque[ns] = k;
        }
      }
      return true;
    }
    long  PoiAMT(double &tasa, long Alias[], double Q[], long &NAlias, double &Suma, double &ProbN)
    // Generación de muestra de VA Poisson usando el método alias
    // Los arreglos alias y q deben haberse generado previamente con SetPois
    { double uni1, uni2, Sumap, p; long indice;
      uni1 = RanMT();
      if (uni1 >= Suma)
      {//Se busca hacia arriba por transformación inversa
        indice = NAlias; p = ProbN * tasa / (double)indice; Sumap = Suma + p;
        while (uni1 > Sumap)
        { indice += 1; p = p * tasa / (double)indice; Sumap += p;
        }
        return indice;
      }
      else
      {// Se busca hacia abajo usando el método Alias
        uni2 = (double)NAlias * uni1 / Suma; indice = floor(uni2);
        uni2 = (uni2 - (double)indice) * Suma;
        if (uni2 <= Q[indice + 1]) return indice;
        else return Alias[indice + 1] - 1;
      }
    }
    long  BineMT(long &n, double &p)
    {//Si Y~Gama(s, (1-p)/p), entonces Pois(Y)~BiNeg(s,p). Se supone que 0 < p < 1
      double y, expy, q = (1-p)/p, s = (double)n;
      y = GammMT(s, q); expy = exp(-y);
      return PoisMT(y, expy);
    }
    double BetaMT(double &alfa1,double &alfa2)
    {//Generador de VA Beta basado en gammas
      double y1, y2,x,p;
      p=1;
      y1= GammMT(alfa1,p);
      y2= GammMT(alfa2,p);
      x=y1/(y1+y2);
      return x;
    }
    double HiexMT( double &prob, double &tasa1, double &tasa2)
    {//Generador de Hiperexponencial con parámetro de mezcla prob
      double uni;
      uni=RanMT();
      if (prob > uni) return (ExpoMT(tasa1));
      else return (ExpoMT( tasa2));
    }
    double PeaVMT(double &alfa, double &beta)
    {//Generador de Pearson V
      double y,inv;
      if (beta > 0)
      {
        inv=1/beta;
        y=GammMT(alfa, inv);
        return (1/y);
      }
      else return 0;

    }
    double PeVIMT( double &alfa1, double &alfa2, double &beta)
    {//Generador de Pearson VI
      double y1, y2, u=1.0;
      y1=GammMT( alfa1, beta);
      y2=GammMT(alfa2, u);
      return (y1/y2);
    }
    double JohBMT( double &mu, double &sigma, double &a, double &b)
    {// Generador de Johnson B, a es mínimo y b es multiplicador
      double y, x;
      y=Norm2MT(mu, sigma);
      x =1 / (1 +exp(-y));
      return (a + b * x);
    }

    double JohUMT(double &mu, double &sigma, double &a, double &b)
    {// Generador de Johnson U, a es mínimo y b es multiplicador
      double y, x;
      y=exp(Norm2MT(mu, sigma));
      x = (y - 1/y) / 2 ;
      return (a + b * x);
    }

    double LogLMT(double &alfa, double &beta)
    {//Generador de Log Logística
      double uni, calc;
      uni=RanMT();
      calc=pow(uni/(1-uni), 1/alfa);
      return (beta*calc);
    }

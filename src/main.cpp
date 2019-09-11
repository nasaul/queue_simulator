#include <Rcpp.h>
#include "FunctionsRandom.h"
#include "SupportFunctions.h"

/* ----------------------------------------------------------------
  Simula G/G/S
 ----------------------------------------------------------------*/
// [[Rcpp::export]]
 Rcpp::List SimulaGGs(
   Rcpp::IntegerVector   NDatos,
   Rcpp::NumericVector   Programa,
   Rcpp::CharacterVector dist_atencion,
   Rcpp::NumericVector   param_atencion,
   Rcpp::CharacterVector dist_llegadas,
   Rcpp::NumericVector   param_llegadas
 ){
   long       N1 = NDatos[0] + NDatos[1] + NDatos[2];
   long       L[2*N1+2], NCambios {0}, Llegaron = 0, NSistema = 0, Salieron = 0, NLlegadas = NDatos[0] + NDatos[1] + NDatos[2], *EnAtencion, *Mark;
   double     Arreglo[N1], TSalida[N1], W {0}, Cambios[2*N1+2], Wq[N1], Reloj = 0.0, *TAtencion, *TLlegada;
   const long S {1};


   TAtencion = (double*) malloc((NLlegadas + 1) * sizeof(double));
   TLlegada  = (double*) malloc((NLlegadas + 1) * sizeof(double));
   Mark      = ivector(0, NLlegadas);
   double (*ratention)(Rcpp::NumericVector);
   double (*rarrival)(Rcpp::NumericVector);

   if(dist_atencion[0] == "exp"){
     if(param_atencion.length() != 1){
       Rcpp::stop("'param_atencion' with exponential distribution must be length 1.");
     }
     ratention = ExpoMT;
   } else {
     Rcpp::stop("Distribution name not defined.");
   }

   if(dist_llegadas[0] == "exp"){
     if(param_llegadas.length() != 1){
       Rcpp::stop("'param_llegadas' with exponential distribution must be length 1.");
     }
     rarrival = ExpoMT;
   } else {
     Rcpp::stop("Distribution name not defined.");
   }


   for (long i = 0; i < NDatos[0]; i++)
   {//Llegadas y atenciones para los clientes iniciales
       TLlegada[i] = 0.0;    Arreglo[i] = 0.0;    TAtencion[i] = ratention(param_atencion);    Mark[i] = i;
   }
   long Ini1, Fin1 = NDatos[0] + NDatos[1];
   for (long i = NDatos[0]; i < Fin1; i++)
   {// Llegadas y tiempos de atencion programados
       TLlegada[i] = Programa[i - NDatos[0]];    Arreglo[i] = TLlegada[i];
       TAtencion[i] = ratention(param_atencion);    Mark[i] = i;
   }
   if (NDatos[2] > 0)
   {
       Arreglo[Fin1] = rarrival(param_llegadas);        TLlegada[Fin1] = Arreglo[Fin1];
       TAtencion[Fin1] = ratention(param_atencion); Mark[Fin1] = Fin1;
   }
   for (long i = Fin1 + 1; i < NLlegadas; i++)
   {// Llegadas y tiempos de atencion para clientes nuevos
       Arreglo[i] = Arreglo[i - 1] + rarrival(param_llegadas);    TLlegada[i] = Arreglo[i];
       TAtencion[i] = ratention(param_atencion);    Mark[i] = i;
   }

   Ini1 = NDatos[0]; Fin1 = NLlegadas - 1;    sortdll(TLlegada, Ini1, Fin1, Mark);
   EnAtencion = ivector(0, S);  //EnAtencion es el arreglo de los clientes que est·n siendo atendidos, ordenados por tiempo de salida
   W = 0.0;    NCambios = 0;    L[NCambios] = NSistema;        Cambios[NCambios] = 0.0;

   while (Llegaron < NLlegadas)        //Llegaron es el n˙mero de clientes que llegaron al sistema
   {
       if ((NSistema == 0) || (TLlegada[Llegaron] < TSalida[EnAtencion[0]]))
       {
           ProcesaLlegadaR(Reloj, Llegaron, Mark, NSistema, TAtencion, TSalida, Arreglo, EnAtencion, S, Wq);  //B˙squeda lineal
           //        ProcesaLlegadaR1(Reloj, Llegaron, Mark, NSistema, TAtencion, TSalida, Arreglo, EnAtencion, S, Wq); //B˙squeda binaria
           NCambios += 1;    L[NCambios] = NSistema;    Cambios[NCambios] = Reloj;
       }
       else
       {
           ProcesaSalidaR(Reloj, Llegaron, NSistema, TAtencion, TSalida, Arreglo, Mark, EnAtencion, S, Salieron, W, Wq); //B˙squeda lineal
           //    ProcesaSalidaR1(Reloj, Llegaron, NSistema, TAtencion, TSalida, Arreglo, Mark, EnAtencion, S, Salieron, W, Wq); //B˙squeda binaria
           NCambios += 1;    L[NCambios] = NSistema;    Cambios[NCambios] = Reloj;
       }
   }

   Arreglo[NLlegadas] = Arreglo[NLlegadas - 1] + rarrival(param_llegadas);

   N1 = NSistema;

   while (Salieron < NLlegadas)
   {
       if (Arreglo[NLlegadas] < TSalida[EnAtencion[0]])
       {
           Arreglo[NLlegadas] += rarrival(param_llegadas); N1 += 1;
       }
       else
       {
           ProcesaSalidaR(Reloj, Llegaron, NSistema, TAtencion, TSalida, Arreglo, Mark, EnAtencion, S, Salieron, W, Wq); //B˙squeda lineal
           //    ProcesaSalidaR1(Reloj, Llegaron, NSistema, TAtencion, TSalida, Arreglo, Mark, EnAtencion, S, Salieron, W, Wq); //B˙squeda binaria
           N1 -= 1;
           NCambios += 1;    L[NCambios] = N1;    Cambios[NCambios] = Reloj;
       }
   }

   if (Salieron > 0) W = W / (double)Salieron;
   free_ivector(EnAtencion, 0, S);
   free_vector(TAtencion, 0, NLlegadas);
   free_ivector(Mark, 0, NLlegadas);
   free_vector(TLlegada, 0, NLlegadas);
   return Rcpp::List::create(
     Rcpp::Named("W") = W
   );
 }

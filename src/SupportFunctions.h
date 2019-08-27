#ifndef SupportFunctions_h
#define SupportFunctions_h

#include <stdio.h>
long *ivector(int nl,int nh);
void sortdll(double it[], long &l, long &r, long mark[]);
void ProcesaLlegadaR(double &Reloj, long &Llegaron, long Mark[], long &NSistema, const double TServicio[], double TSalida[],
                     const double TLlegada[], long *EnAtencion, const long &S, double Wq[]);

void ProcesaSalidaR(double &Reloj, long &Llegaron, long &NSistema, const double TServicio[], double TSalida[],
                    const double TLlegada[], long Mark[], long *EnAtencion, const long &S, long &Salieron, double &W, double Wq[]);
void free_ivector(long *v, int nl, int nh);
void free_vector(double *v, int nl, int nh);


#endif /* SupportFunctions_h */

//
//  SupportFunctions.cpp
//  TestFunctionsRcpp
//
//  Created by Luis Moncayo on 8/23/19.
//  Copyright © 2019 Luis Moncayo. All rights reserved.
//

#include "SupportFunctions.h"
#include <cstdlib>


/* ----------------------------------------------------------------
 Crea un vector de enteros.
 ----------------------------------------------------------------*/
long *ivector(int nl,int nh)
{
    long *v;

    v=(long *)malloc((unsigned) (nh-nl+1)*sizeof(long));
    //if (!v) nrerror("allocation failure in vector()");
    return v-nl;
}
//--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++
void sortdll(double it[], long &l, long &r, long mark[])
// Ordena el arreglo it entre l y r, en mark est·n los Ìndices de los elementos
// Se usa para llamadas a librerÌa dll
{
    double x, y;
    long i, j, k;
    long temp;
    i = l;
    j = r;
    temp = (l + r) / 2;
    x = it[temp];
    do
    {
        while (it[i] < x)
        {
            i = i + 1;
        }
        while (x < it[j])
        {
            j = j - 1;
        }
        if (i <= j)
        {
            y = it[i];
            it[i] = it[j];
            it[j] = y;
            k = mark[i];
            mark[i] = mark[j];
            mark[j] = k;
            i = i + 1;
            j = j - 1;
        }
    } while (i <= j);
    if (l < j)
    {
        sortdll(it, l, j, mark);
    }
    if (i < r)
    {
        sortdll(it, i, r, mark);
    }
}
//--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++
void ProcesaLlegadaR(double &Reloj, long &Llegaron, long Mark[], long &NSistema, const double TServicio[], double TSalida[],
                     const double TLlegada[], long *EnAtencion, const long &S, double Wq[])
// Procesa la llegada de un cliente para simular una cola G/G/s FIFO. Se utiliza en SimulaGGs
// Mark indica el orden de las llegadas
// El tiempo de atenciÛn est· en TServicio.
// La b˙squeda en la cola es lineal
{
    bool Done = false;
    Reloj = TLlegada[Mark[Llegaron]];
    if (NSistema < S)
    {//El cliente se atiende de inmediato
        TSalida[Mark[Llegaron]] = Reloj + TServicio[Mark[Llegaron]];    Wq[Mark[Llegaron]] = Reloj - TLlegada[Mark[Llegaron]];
        int Pos = NSistema;
        while (!Done)
        {//Arreglo EnAtencion se mantiene ordenado por tiempo de salida
            if (Pos == 0)
            {
                EnAtencion[Pos] = Mark[Llegaron]; Done = true;
            }
            else
            {
                if (TSalida[Mark[Llegaron]] < TSalida[EnAtencion[Pos - 1]])
                {
                    EnAtencion[Pos] = EnAtencion[Pos - 1]; Pos -= 1;
                }
                else
                {
                    EnAtencion[Pos] = Mark[Llegaron]; Done = true;
                }
            }
        }
    }
    Llegaron += 1; NSistema += 1;
}

//--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++
void ProcesaSalidaR(double &Reloj, long &Llegaron, long &NSistema, const double TServicio[], double TSalida[],
                    const double TLlegada[], long Mark[], long *EnAtencion, const long &S, long &Salieron, double &W, double Wq[])
{ // Procesa la salida de un cliente para simular una cola G/G/s FIFO, el ID del cliente est· en Mark.
    // Se utiliza en SimulaGGs
    // El tiempo de atenciÛn est· en TServicio.  La b˙squeda en la cola es lineal
    Reloj = TSalida[EnAtencion[0]];
    W += TSalida[EnAtencion[0]] - TLlegada[EnAtencion[0]];
    if (NSistema <= S)
    {//No hay clientes en cola, sÛlo se actualiza EnAtencion
        long Pos = 0;
        while (Pos < NSistema - 1)
        {
            EnAtencion[Pos] = EnAtencion[Pos + 1]; Pos += 1;
        }
    }
    else
    {//Hay clientes en cola, se atiende al primer cliente en cola
        long Pos0 = Llegaron - NSistema + S;              //Primer Cliente en la Cola
        TSalida[Mark[Pos0]] = Reloj + TServicio[Mark[Pos0]];    Wq[Mark[Pos0]] = Reloj - TLlegada[Mark[Pos0]];
        long Pos = S - 1;
        bool Done = false;
        while ((Pos > 0) && !Done)
        {//Se ubica la salida de Mark[Pos0] en el arreglo  EnAtencion
            if (TSalida[Mark[Pos0]] < TSalida[EnAtencion[Pos]]) Pos -= 1;
            else Done = true;
        }
        for (long Pos1 = 0; Pos1 < Pos; Pos1++) EnAtencion[Pos1] = EnAtencion[Pos1 + 1]; //Se mantiene EnAtencion ordenado por tiempo de salida
        EnAtencion[Pos] = Mark[Pos0];
    }
    NSistema -= 1; Salieron += 1;
}

//--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++
/* ----------------------------------------------------------------
 Libera la memoria de un vector de enteros.
 ----------------------------------------------------------------*/
void free_ivector(long *v, int nl, int nh)
{
    free((char*) (v+nl));
}

//--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++
/* ----------------------------------------------------------------
 Libera la memoria de un vector double.
 ----------------------------------------------------------------*/
void free_vector(double *v, int nl, int nh)
{
    free((char*) (v+nl));
}

#ifndef EventCut_h
#define EventCut_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
// #include "../km2aevent.h"
#include "../KM2Afliter.h"

class EventCut{

    public:
        //We put these constants as the member variables
        //static const int   Nhit_bin   =    10;
        static const int   Nhit_bin   =    14;
        static const int   Zenith_Cut =    50 * 3.1415926 / 180.;
        static const int   CoreX_Cut  =    200;
        static const int   CoreY_Cut  =    200;
        static const int   Ebin[Nhit_bin+1];
        static const int   radius[Nhit_bin];
        static const float R_Cut[Nhit_bin];

        EventCut();
        virtual ~EventCut();
        EventCut *Instance() {
            if (!fInstance) fInstance = new EventCut();
            return fInstance;
        }

        static int EventSelection(KM2Afliter* event);

    private:
        static EventCut *fInstance;
};

#ifdef EventCut_cc
const int EventCut::Ebin[EventCut::Nhit_bin+1]   =  {0.6, 0.8, 1,  1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 10000};
const float EventCut::R_Cut[EventCut::Nhit_bin] =  {-5.11,-5.11,-5.11, -5.24, -5.95, -6.08, -2.34, -2.35, -2.36, -2.36, -2.36, -2.36,-2.36,-2.36};
EventCut *EventCut::fInstance = 0;

int EventCut::EventSelection(KM2Afliter* event)
{
    // //Gamma
    // for(int i=0;i<Nhit_bin;i++){
    //     double R = std::log10((event->NuM4 + 0.0001) / event->NpE3);
    //     if( event->rec_E< Ebin[i+1] && event->rec_E>=Ebin[i] && R<R_Cut[i] && event->rec_theta < Zenith_Cut && event->rec_theta >= 0 && event->NpE1>10 && event->NfiltE>10 && event->NhitM>10 && event->NtrigE>19 && event->rec_Eage >=0.6 && event->rec_Eage<=2.4 && (event->NpW >= 0 && event->NpE1/event->NpE2 > 2))
    //     {   
    //         return i;
    //     }
    // }

    int energy_bin = -1;
    if(event->rec_E > 0.6 && event->rec_E <= 0.8)
        energy_bin = 0;
    else if(event->rec_E > 0.8 && event->rec_E <= 1)
        energy_bin = 1;
    else if(event->rec_E > 1 && event->rec_E <= 1.2)
        energy_bin = 2;
    else if(event->rec_E > 1.2 && event->rec_E <= 1.4)
        energy_bin = 3;
    else if(event->rec_E > 1.4 && event->rec_E <= 1.6)
        energy_bin = 4;
    else if(event->rec_E > 1.6 && event->rec_E <= 1.8)
        energy_bin = 5;
    else if(event->rec_E > 1.8 && event->rec_E <= 2.0)
        energy_bin = 6;
    else if(event->rec_E > 2.0 && event->rec_E <= 2.2)
        energy_bin = 7;
    else if(event->rec_E > 2.2 && event->rec_E <= 2.4)
        energy_bin = 8;
    else if(event->rec_E > 2.4 && event->rec_E <= 2.6)
        energy_bin = 9;
    else if(event->rec_E > 2.6 && event->rec_E <= 2.8)
        energy_bin = 10;
    else if(event->rec_E > 2.8 && event->rec_E <= 3.0)
        energy_bin = 11;
    else if(event->rec_E > 3.0 && event->rec_E <= 3.2)
        energy_bin = 12;
    else if(event->rec_E > 3.2)
        energy_bin = 13;
    // else
    //     return -1;
    return energy_bin;
}
#endif // #ifdef EventCut_cc
#endif // #ifdef EventCut_h
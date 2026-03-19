#include "udf.h"

#define WALL_ZONE_ID 292
#define FLUID_ZONE_ID 514

#define UDM_AREA   0

DEFINE_ON_DEMAND(facethread_memo)
{
    Domain *d = Get_Domain(1);
    Thread *t_fluid = Lookup_Thread(d, FLUID_ZONE_ID);
    Thread *t_wall  = Lookup_Thread(d, WALL_ZONE_ID);

    cell_t c;
    face_t f;
    Thread *tf;
    int n;

    if (NULLP(t_fluid) || NULLP(t_wall))
    {
        Message("facethread_memo: could not find one of the threads.\n");
        return;
    }

    Message("facethread_memo: starting tagging of near-wall cells...\n");

    begin_c_loop(c, t_fluid)
    {
        real A = 0.0;
        real AA[ND_ND];

        C_UDMI(c, t_fluid, UDM_AREA) = 0.0;

        c_face_loop(c, t_fluid, n)
        {
            f  = C_FACE(c, t_fluid, n);
            tf = C_FACE_THREAD(c, t_fluid, n);

            if (BOUNDARY_FACE_THREAD_P(tf) && tf == t_wall)
            {
                F_AREA(AA, f, tf);
                A += NV_MAG(AA);
            }
        }

        C_UDMI(c, t_fluid, UDM_AREA) = A;
    }
    end_c_loop(c, t_fluid)

    Message("facethread_memo: tagging finished.\n");
}

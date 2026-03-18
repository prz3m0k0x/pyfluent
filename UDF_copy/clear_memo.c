#include "udf.h"

DEFINE_ON_DEMAND(clear_memo)
{
#if !RP_HOST
    Domain *d = Get_Domain(1);
    Thread *t;
    cell_t c;

    thread_loop_c(t, d)
    {
        begin_c_loop(c, t)
        {
            C_UDMI(c, t, 0) = 0.0;
            C_UDMI(c, t, 1) = 0.0;
            C_UDMI(c, t, 2) = 0.0;
            C_UDMI(c, t, 3) = 0.0;
            C_UDMI(c, t, 4) = 0.0;
        }
        end_c_loop(c, t)
    }
#endif
}
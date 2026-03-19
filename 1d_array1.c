#include "udf.h"
#include <math.h>

#define SPECIES_ID 0

#define UDM_AREA   0   /* active adsorption area [m2] */
#define UDM_QA     1   /* areal loading [kg/m2] */
#define UDM_QDOTA  2   /* d(qA)/dt [kg/m2/s] */
#define UDM_SRC    3   /* source [kg/m3/s] */
#define UDM_QTOT   4   /* total adsorbed mass in cell [kg] */

/* particle / bed parameters */
#define PI 3.14159265358979323846
#define RHO_PART 650.0
#define R_PART 0.00247896
#define A_PART (4.0 * PI * (R_PART) * (R_PART))
#define V_PART ((4.0/3.0) * PI * (R_PART) * (R_PART) * (R_PART))

/* kinetics / limits */
#define K_KINETIC 0.027
#define SAT_FRAC  0.998
#define C_MIN     1.0e-8
#define C_MAX     0.00050000001

/* sorted qe cache */
#define MAX_QE_CACHE 10000
#define QE_REL_TOL   1.0e-3   /* reuse qe if |Ce_cached-Ce| <= QE_REL_TOL * max(Ce, C_MIN) */

static real qe_ce[MAX_QE_CACHE];
static real qe_val[MAX_QE_CACHE];
static int  qe_count = 0;
static int  cache_ready = 0;
static real qA_sat = 0.0;


/* explicit equilibrium adsorption per unit area [kg/m2] */
static real qe_areal_raw(real Ce)
{
    real qm = 0.002979;
    real Kl = 26.0;
    real N  = 5.377;
    real arg;

    if (Ce <= 0.0)
        return 0.0;

    arg = pow(Kl * Ce, N);

    return qm * arg / (1.0 + arg) * V_PART * RHO_PART / A_PART;
}

static void init_qe_cache_storage(void)
{
    qe_count = 0;
    cache_ready = 1;
}


static int qe_lower_bound(real Ce)
{
    int left = 0;
    int right = qe_count;

    while (left < right)
    {
        int mid = left + (right - left) / 2;

        if (qe_ce[mid] < Ce)
            left = mid + 1;
        else
            right = mid;
    }

    return left;
}


static int qe_close_enough(real a, real b)
{
    real ref = fabs(b);
    real tol;

    if (ref < C_MIN)
        ref = C_MIN;

    tol = QE_REL_TOL * ref;

    return (fabs(a - b) <= tol);
}


/* sorted-array cache lookup:
   - below C_MIN: return 0
   - above C_MAX: clamp to C_MAX
   - if nearby cached Ce exists: reuse qe
   - otherwise compute qe once and insert exact (Ce, qe) pair */
static real qe_areal_lookup(real Ce)
{
    int pos, i;
    real qeq;

    if (!cache_ready)
        init_qe_cache_storage();

    if (Ce < C_MIN)
        return 0.0;

    if (Ce > C_MAX)
        Ce = C_MAX;

    if (qe_count == 0)
    {
        qeq = qe_areal_raw(Ce);
        qe_ce[0] = Ce;
        qe_val[0] = qeq;
        qe_count = 1;
        return qeq;
    }

    pos = qe_lower_bound(Ce);

    if (pos < qe_count && qe_close_enough(qe_ce[pos], Ce))
        return qe_val[pos];

    if (pos > 0 && qe_close_enough(qe_ce[pos - 1], Ce))
        return qe_val[pos - 1];

    qeq = qe_areal_raw(Ce);

    if (qe_count < MAX_QE_CACHE)
    {
        for (i = qe_count; i > pos; --i)
        {
            qe_ce[i]  = qe_ce[i - 1];
            qe_val[i] = qe_val[i - 1];
        }

        qe_ce[pos]  = Ce;
        qe_val[pos] = qeq;
        qe_count++;
    }

    return qeq;
}


DEFINE_ON_DEMAND(init_adsorption_cache)
{
    init_qe_cache_storage();
    qA_sat = qe_areal_raw(C_MAX);
    Message("qe cache initialized. Capacity = %d, rel_tol = %g\n", MAX_QE_CACHE, QE_REL_TOL);
}



DEFINE_SOURCE(mass_source, c, t, dS, eqn)
{
    real yi, rho, ci;
    real A, qA, qeq, dqA_dt, src;


    A = C_UDMI(c, t, UDM_AREA);
    if (A <= 1e-12)
    {
        C_UDMI(c, t, UDM_QDOTA) = 0.0;
        C_UDMI(c, t, UDM_SRC)   = 0.0;
        dS[eqn] = 0.0;
        return 0.0;
    }

    yi  = C_YI(c, t, SPECIES_ID);
    rho = C_R(c, t);
    ci  = yi * rho;

    if (ci < C_MIN)
    {
        C_UDMI(c, t, UDM_QDOTA) = 0.0;
        C_UDMI(c, t, UDM_SRC)   = 0.0;
        dS[eqn] = 0.0;
        return 0.0;
    }

    qA = C_UDMI(c, t, UDM_QA);
    if (qA < 0.0)
        qA = 0.0;

    qeq = qe_areal_lookup(ci);

    if (qA >= SAT_FRAC * qA_sat)
                {
                    C_UDMI(c, t, UDM_AREA) = 0;
                    return 0.0;
                }
    dqA_dt = K_KINETIC * (qeq - qA);
    src = -(A / C_VOLUME(c, t)) * dqA_dt;

    C_UDMI(c, t, UDM_QDOTA) = dqA_dt;
    C_UDMI(c, t, UDM_SRC)   = src;

    dS[eqn] = 0.0;
    return src;
}


DEFINE_EXECUTE_AT_END(update_adsorption)
{
#if !RP_HOST
    Domain *d = Get_Domain(1);
    Thread *t;
    cell_t c;
    real dt = CURRENT_TIMESTEP;

    if (!cache_ready)
        init_qe_cache_storage();

    thread_loop_c(t, d)
    {
        begin_c_loop(c, t)
        {
            real A  = C_UDMI(c, t, UDM_AREA);
            real qA = C_UDMI(c, t, UDM_QA);
            real dq = C_UDMI(c, t, UDM_QDOTA);

            if (A > 1e-16)
            {
                qA += dq * dt;

                if (qA < 0.0)
                    qA = 0.0;


                C_UDMI(c, t, UDM_QA)   = qA;
                C_UDMI(c, t, UDM_QTOT) = A * qA;
            }
            else
            {
                C_UDMI(c, t, UDM_QDOTA) = 0.0;
                C_UDMI(c, t, UDM_SRC)   = 0.0;
            }
        }
        end_c_loop(c, t)
    }
#endif
}

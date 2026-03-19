/*******************************************************************
    UDF liczacy adsorpcje na powierzchni zdefiniowanej scianki w komorkach przyleglych w danym fluid zone
*******************************************************************/
#include "udf.h"
#include "math.h"

/* Liczba PI do obliczania srednic i powierzchni kulek adsorbenta */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Nazwy sciany (wall zone) i czesci plynu w ktorej bedzie dzialac UDF - z FLUENTA!!! */
#define WALL_ZONE_NAME "adsorbent_walls"   /* dokladna nazwa sciany */
#define FLUID_ZONE_NAME "geometria_powierzchnia"           /* dokladna nazwa fluid zone*/

/* parametry czastek w zlozu */
#define RHO_PART 650.0                   /* gestosc nasypowa [kg/m3] */
#define R_PART 0.00247896                /* promien czastki [m] */
#define A_PART (4.0 * M_PI * (R_PART) * (R_PART))                 /* pole powierzchni jednej kulki [m2] */
#define V_PART ((4.0/3.0) * M_PI * (R_PART) * (R_PART) * (R_PART))/* objetosc jednej kulki [m3] */

/*******************************************************************
    Izoterma Simpsa w przeliczeniu na procent powierzchni kulki (zwraca adsorpcje w [kg])
*******************************************************************/
double qe_cA(double Ce, double A)
{
    double qe, qm, Kl, N;
    qm = 0.002979;   /* [kg/kg] */
    Kl = 26.0;       /* [m3/kg] */
    N  = 5.377;

    qe = (A / A_PART) * (V_PART * RHO_PART) * qm *
         pow(Kl * Ce, N) / (1.0 + pow(Kl * Ce, N));

    return qe; /* [kg] */
}

/*******************************************************************
    DEFINE_ON_DEMAND: makro, ktore flaguje komorki przylegle do sciany
    poprzez przypisanie wartosci logicznej 1 dla komorek znajdujacych sie
    przy zdefiniowanym WALL_ZONE w pamieci uzytkownika 2. Makro przypisuje rowniez pole powierzchni
    face'a do pamieci uzytkownika 3
*******************************************************************/
DEFINE_ON_DEMAND(facethread_memo)
{
    Domain *d = Get_Domain(1);
    Thread *t_fluid = NULL;
    Thread *t_wall  = NULL;

    /* Komunikat, ktory zwraca wiadomosc czy znalazl dany fluid zone */
    thread_loop_c(t_fluid, d)
    {
        if (strcmp(THREAD_NAME(t_fluid), FLUID_ZONE_NAME) == 0)
        {
            Message("Found fluid zone: %s (ID=%d)\n",
                    THREAD_NAME(t_fluid), THREAD_ID(t_fluid));
            break;
        }
    }

    /* Komunikat, ktory zwraca wiadomosc czy znalazl dany wall zone */
    thread_loop_f(t_wall, d)
    {
        if (BOUNDARY_FACE_THREAD_P(t_wall) &&
            strcmp(THREAD_NAME(t_wall), WALL_ZONE_NAME) == 0)
        {
            Message("Found wall zone: %s (ID=%d)\n",
                    THREAD_NAME(t_wall), THREAD_ID(t_wall));
            break;
        }
    }

    if (t_fluid == NULL || t_wall == NULL)
    {
        Message("Error: Could not find zones: %s or %s\n",
                FLUID_ZONE_NAME, WALL_ZONE_NAME);
        return;
    }

    Message("facethread_memo: starting tagging of near-wall cells...\n");

    cell_t c;
    face_t f;
    Thread *tf;
    int n;

    begin_c_loop_all(c, t_fluid)
    {
        double A = 0.0, AA[ND_ND];

        c_face_loop(c, t_fluid, n)
        {
            f = C_FACE(c, t_fluid, n);
            tf = C_FACE_THREAD(c, t_fluid, n);

            if (BOUNDARY_FACE_THREAD_P(tf) && tf == t_wall)
            {
                F_AREA(AA, f, tf);
                A = NV_MAG(AA);
                C_UDMI(c, t_fluid, 2) = 1; /* mark as near-wall */
                break;
            }
            else
            {
                C_UDMI(c,t_fluid,2) = 0;
            }
        }

        C_UDMI(c, t_fluid, 3) = A; /* store wall area */
    }
    end_c_loop_all(c, t_fluid)

    Message("facethread_memo: tagging finished.\n");
}
/*Uwaga - prowadzac obliczenia w trybie parallel komunikaty o wallzonach
moga wyswietlic sie wielokrotnie. Mozna zignorowac :) */

/*******************************************************************
    DEFINE_SOURCE: wlasciwe makro liczace adsorpcje w kazdym kroku czasowym
*******************************************************************/
DEFINE_SOURCE(mass_source, c, t, dS, eqn)
{
    /*Sprawdzamy w kazdej komorce, czy jest "oflagowana" jako przyscienna*/
    if (C_UDMI(c, t, 2) < 0.5)
        return 0.0;
    /*Definicja zmiennych*/
    double src = 0.0, qe = 0.0, qt = 0.0;
    double rho, g_i, c_i;
    double A;
    int i = 0;              /* indeks zwiazku w panelu mixture */
    double k = 0.027;       /* kinetic constant [1/s] */

    g_i = C_YI(c, t, i); /*Udzial masowy skladnika*/
    rho = C_R(c, t);    /*gestosc mieszaniny*/
    c_i = g_i * rho; /*stezenie masowe skladnika i [kg/m3] */

    /*Krok przypisujacy pole powierzchni z pamieci uzytkownika 3, ktora zostala utworzona 
    w kazdej komorce danego fluid zone w makrze "define on demand"
    Dodatkowy warunek A<=0.0 zapobiega bledom, w ktorych komorka zostalaby oflagowana jako znajdujaca
    sie przy scianie, ale majaca zerowa powierzchnie kontaktu*/
    A = C_UDMI(c, t, 3);
    if (A <= 0.0)
        return 0.0;

    qt = C_UDMI(c, t, 1);          /* aktualna ilosc zaadsorbowanej substancji[kg] */
    qe = qe_cA(c_i, A);            /* adsropcja w rownowadze z izotermy Simpsa [kg] */
    double dqedt = k * (qe - qt);  /* wyrazenie na szybkosc adsorpcji, PSO[kg/s] */
    src = -dqedt / C_VOLUME(c, t); /* zrodla masy we fluencei musza byc w jednostkach [kg/m3/s], a wiec dzielimy
                                      aktualnie obliczona kinetyke przez objetosc komorki, z ktorej zostanie usuniety skladnik */

    C_UDMI(c, t, 0) = src;          /*Zapisujemy aktualne zrodlo*/
    dS[eqn] = 0.0;
    return src;
}

/*******************************************************************
    DEFINE_EXECUTE_AT_END: Uaktualnienie adsorpcji calkowitej i relatywnej (do objetosci czastki)
*******************************************************************/
DEFINE_EXECUTE_AT_END(add_the_source)
{
#if !RP_HOST
    Domain *d = Get_Domain(1);
    Thread *t;
    cell_t c;
    double A;

    thread_loop_c(t, d)
    {
        begin_c_loop(c, t)
        {
            if (C_UDMI(c, t, 2) > 0.5) /*To makro bedzie wywolane tylko w komorkach, ktore sasiaduja z adsorbentem*/
            {
                A = C_UDMI(c, t, 3); /*Pole powierzchni face'a przyleglego do adsorbenta, jak wczesniej [m2]*/
                C_UDMI(c, t, 1) = C_UDMI(c, t, 1) -
                                  C_UDMI(c, t, 0) * C_VOLUME(c, t) * CURRENT_TIMESTEP; /*Obliczenie aktualnej adsorpcji qt*/

                if (A > 0.0) /*Jesli komorka sasiaduje ze sciana, obliczamy aktualna adsorpcje w przeliczeniu na kg adsorbenta
                               Technicznie robimy odwrotnosc operacji zawartej np. w izotermie Simpsa gdzie przeliczamy qe na 
                               powierzchnie danego face'a. Ten fragment sluzy jedynie kwestiom prezentacyjnym.*/
                    C_UDMI(c, t, 4) = C_UDMI(c, t, 1) /
                        ((A / A_PART) * (V_PART * RHO_PART));
                else
                    C_UDMI(c, t, 4) = 0.0;
            }
            else
            {
                C_UDMI(c, t, 4) = 0.0;
            }
        }
        end_c_loop(c, t)
    }
#endif
}

/*******************************************************************
    DEFINE_ON_DEMAND: makro ktore wyczysci wszystkie pamieci uzytkownika
*******************************************************************/
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
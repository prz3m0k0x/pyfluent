''
Plik solution zawiera workflow tworzący geometrie, siatkę i rozwiązujący konwersję SO2 do SO3 w reaktorze katalitycznym.
Parameterami są elementy geometrii, długość poszczególnych stref reaktora oraz ułamek masowy SO2 na wlocie.
Maksymalna długość strefy dwóch katalizatorów i strefy chłodzącej to 4000 mm:
x1 = LenghtCat1
x2 = LenghtCool
x3 = LenghtCat2
x4 = "inlet_Y_SO2"

x1 + x2 + x3 <= 4000
500 <= x1, x2, x3 <= 2000
0.05 <= x4 <= 0.2

Funkcjami celu są obliczane przez solver:
Y1 = so3_mass_fraction_out
Y2 = so2_mass_fraction_out
Y3 = so3_mass_fraction_out/so2_mass_fraction_out

Funkcje celu to:

maximize(Y1), minimize(Y2), maximize(Y3)


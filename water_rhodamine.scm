;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                                                              ;;;
;;;             Fluent USER DEFINED MATERIAL DATABASE            ;;;
;;;                                                              ;;;
;;; (name type[fluid/solid] (chemical-formula . formula)         ;;;
;;;             (prop1 (method1a . data1a) (method1b . data1b))  ;;;
;;;            (prop2 (method2a . data2a) (method2b . data2b)))  ;;;
;;;                                                              ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(
	(my-mixture mixture
		(chemical-formula . #f)
		(species (names (rh-b h2o<l>) () ()))
		(reactions (finite-rate ("reaction-1" () () () (stoichiometry "  ") (arrhenius 1000000000000000. 100000000. 0) (mixing-rate 4 0.5) (specified-rate-exponents? . #t) (use-third-body-efficiencies? . #f))) (finite-rate/eddy-dissipation ("reaction-1" () () ((o2 0 1) (h2o 0 1) (n2 0 1)) (stoichiometry "  ") (arrhenius 1000000000000000. 100 0) (mixing-rate 4 0.5) (specified-rate-exponents? . #t) (use-third-body-efficiencies? . #f))) (eddy-dissipation ("reaction-1" () () ((o2 0 1) (h2o 0 1) (n2 0 1)) (stoichiometry "  ") (arrhenius 1000000000000000. 100 0) (mixing-rate 4 0.5) (specified-rate-exponents? . #t) (use-third-body-efficiencies? . #f))))
		(density (volume-weighted-mixing-law . #f) (incompressible-ideal-gas . #f))
		(specific-heat (mixing-law . #f))
		(thermal-conductivity (mass-weighted-mixing-law . #f) (constant . 0.0454))
		(viscosity (mass-weighted-mixing-law . #f) (constant . 1.72e-05))
		(mass-diffusivity (constant-dilute-appx 2.88e-09))
		(speed-of-sound (none . #f))
	)

	(rhodamine-b fluid
		(chemical-formula . rh-b)
		(density (constant . 998.2) (compressible-liquid 101325 998.2 2200000000. 7.15 1.1 0.9))
		(specific-heat (constant . 4182.) (polynomial piecewise-linear (280 . 4201) (300 . 4181) (320 . 4181) (340 . 4188) (360 . 4202) (373 . 4216) (380 . 4224) (400 . 4256) (420 . 4299) (440 . 4357) (460 . 4433) (480 . 4533) (500 . 4664) (520 . 4838) (540 . 5077) (560 . 5424) (580 . 5969)))
		(latent-heat (constant . 2263073))
		(vaporization-temperature (constant . 284))
		(boiling-point (constant . 373))
		(volatile-fraction (constant . 1))
		(binary-diffusivity (film-averaged (averaging-coefficient 0.3333) (film-diffusivity (polynomial piecewise-linear (273 . 2.2e-05) (300 . 2.59e-05) (313 . 2.92e-05) (350 . 3.4e-05) (373 . 3.8e-05) (400 . 4.29e-05) (450 . 5.28e-05) (473 . 5.76e-05)) (constant . 3.05e-05))) (constant . 3.05e-05))
		(thermal-conductivity (constant . 0.6))
		(viscosity (constant . 0.001003))
		(dpm-surften (constant . 0.0719404) (polynomial piecewise-polynomial (50 641 0.0649503 0.000246819 -9.28884e-07 6.01831e-10)))
		(vapor-pressure (polynomial piecewise-linear (273 . 610) (274 . 657) (275 . 706) (280 . 1002) (284 . 1329) (290 . 1937) (295 . 2658) (300 . 3565) (307 . 5316) (310 . 6275) (315 . 7974) (320 . 10612) (325 . 13289) (330 . 17308) (340 . 26579) (350 . 41877) (356 . 53158) (360 . 62498) (370 . 90935) (371 . 94295) (372 . 97757) (373 . 101000) (393 . 202000) (425 . 505000) (453 . 1000000) (486 . 2000000) (507 . 3000000) (537 . 5000000) (584 . 10000000) (615 . 15000000) (639 . 20000000) (647.14 . 22064000)) (constant . 2658))
		(molecular-weight (constant . 470.))
		(species-phase (constant . 1))
		(formation-enthalpy (constant . -285841220.))
		(reference-temperature (constant . 298))
		(lennard-jones-length (constant . 1.))
		(lennard-jones-energy (constant . 100.))
		(formation-entropy (constant . 69902.211))
		(therm-exp-coeff (constant . 0))
		(speed-of-sound (none . #f))
	)

	(water-liquid fluid
		(chemical-formula . h2o<l>)
		(density (constant . 998.2) (compressible-liquid 101325 998.2 2200000000. 7.15 1.1 0.9))
		(specific-heat (constant . 4182) (polynomial piecewise-linear (280 . 4201) (300 . 4181) (320 . 4181) (340 . 4188) (360 . 4202) (373 . 4216) (380 . 4224) (400 . 4256) (420 . 4299) (440 . 4357) (460 . 4433) (480 . 4533) (500 . 4664) (520 . 4838) (540 . 5077) (560 . 5424) (580 . 5969)))
		(latent-heat (constant . 2263073))
		(vaporization-temperature (constant . 284))
		(boiling-point (constant . 373))
		(volatile-fraction (constant . 1))
		(binary-diffusivity (film-averaged (averaging-coefficient 0.3333) (film-diffusivity (polynomial piecewise-linear (273 . 2.2e-05) (300 . 2.59e-05) (313 . 2.92e-05) (350 . 3.4e-05) (373 . 3.8e-05) (400 . 4.29e-05) (450 . 5.28e-05) (473 . 5.76e-05)) (constant . 3.05e-05))) (constant . 3.05e-05))
		(thermal-conductivity (constant . 0.6))
		(viscosity (constant . 0.001003))
		(dpm-surften (constant . 0.0719404) (polynomial piecewise-polynomial (50 641 0.0649503 0.000246819 -9.28884e-07 6.01831e-10)))
		(vapor-pressure (polynomial piecewise-linear (273 . 610) (274 . 657) (275 . 706) (280 . 1002) (284 . 1329) (290 . 1937) (295 . 2658) (300 . 3565) (307 . 5316) (310 . 6275) (315 . 7974) (320 . 10612) (325 . 13289) (330 . 17308) (340 . 26579) (350 . 41877) (356 . 53158) (360 . 62498) (370 . 90935) (371 . 94295) (372 . 97757) (373 . 101000) (393 . 202000) (425 . 505000) (453 . 1000000) (486 . 2000000) (507 . 3000000) (537 . 5000000) (584 . 10000000) (615 . 15000000) (639 . 20000000) (647.14 . 22064000)) (constant . 2658))
		(molecular-weight (constant . 18.0152))
		(species-phase (constant . 1))
		(formation-enthalpy (constant . -285841220.))
		(reference-temperature (constant . 298))
		(lennard-jones-length (constant . 1.))
		(lennard-jones-energy (constant . 100.))
		(formation-entropy (constant . 69902.211))
		(therm-exp-coeff (constant . 0))
		(speed-of-sound (none . #f))
	)

)

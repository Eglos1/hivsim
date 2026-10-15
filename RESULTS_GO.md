
=====================================================================
panlandscape: erase  (60 worlds seed 777, 60y, atol 1e-7)

### E1 strike-efficiency ladder  (erase@365, purge1.5~win1095, ccr5 .999)
erase=0.990000000  C=55/60  R=3  E=2  (257s)
erase=0.999000000  C=58/60  R=0  E=2  (257s)
erase=0.999900000  C=58/60  R=1  E=1  (137s)
erase=0.999990000  C=58/60  R=0  E=2  (257s)
erase=0.999999000  C=59/60  R=1  E=0  (17s)
erase=0.999999900  C=55/60  R=2  E=3  (377s)
erase=0.999999990  C=57/60  R=3  E=0  (18s)
erase=0.999999999  C=57/60  R=1  E=2  (257s)
[panlandscape:erase done in 1577s]

=====================================================================
panlandscape: erase  (60 worlds seed 777, 60y, atol 1e-7)

### E1 strike-efficiency ladder  (erase@365, purge1.5~win1095, ccr5 .999)
erase=0.990000000  C=57/60  R=3  E=0  (257s)
erase=0.999000000  C=59/60  R=1  E=0  (258s)
erase=0.999900000  C=59/60  R=1  E=0  (137s)
erase=0.999990000  C=59/60  R=0  E=1  (377s)
erase=0.999999000  C=59/60  R=1  E=0  (18s)
erase=0.999999900  C=57/60  R=3  E=0  (377s)
erase=0.999999990  C=57/60  R=3  E=0  (17s)
erase=0.999999999  C=58/60  R=2  E=0  (258s)
[panlandscape:erase done in 1700s]

=====================================================================
panlandscape: purge  (60 worlds seed 777, 60y, atol 1e-7)

### E5 purge-strength ladder  (win=1095 fixed)
purge=0.0  (win1095)  C=55/60  R=5  E=0  (139s)
purge=0.3  (win1095)  C=59/60  R=0  E=1  (258s)
purge=0.6  (win1095)  C=59/60  R=1  E=0  (377s)
purge=1.0  (win1095)  C=58/60  R=1  E=1  (378s)
purge=1.5  (win1095)  C=58/60  R=2  E=0  (258s)
purge=2.5  (win1095)  C=59/60  R=0  E=1  (257s)

=====================================================================
panlandscape: start  (60 worlds seed 777, 60y, atol 1e-7)

### E2 therapy-start-day ladder  (PAN core, running ART starts at strike)
start day     7 (ART 7->1102)  C=59/60  R=1  E=0  (105s)
start day    30 (ART 30->1125)  C=57/60  R=3  E=0  (195s)
start day   150 (ART 150->1245)  C=59/60  R=1  E=0  (15s)
start day   365 (ART 365->1460)  C=58/60  R=1  E=1  (285s)
start day   730 (ART 730->1825)  C=58/60  R=1  E=1  (286s)
start day  1095 (ART 1095->2190)  C=0/60  R=60  E=0  (17s)
start day  1460 (ART 1460->2555)  C=0/60  R=60  E=0  (17s)
[panlandscape:start done in 920s]

=====================================================================
panlandscape: start  (60 worlds seed 777, 60y, atol 1e-7)

### E2 therapy-start-day ladder  (PAN core, running ART starts at strike)
start day     7  C=59/60  R=1  E=0  | pre-tx CD4 nadir median=  958 (sick while waiting)  (104s)
start day    30  C=57/60  R=3  E=0  | pre-tx CD4 nadir median=  451 (sick while waiting)  (195s)
start day   150  C=59/60  R=1  E=0  | pre-tx CD4 nadir median=  369 (sick while waiting)  (16s)
start day   365  C=58/60  R=1  E=1  | pre-tx CD4 nadir median=  319 (sick while waiting)  (285s)
start day   730  C=58/60  R=1  E=1  | pre-tx CD4 nadir median=  302 (sick while waiting)  (286s)
start day  1095  C=60/60  R=0  E=0  | pre-tx CD4 nadir median=  280 (sick while waiting)  (16s)
start day  1460  C=57/60  R=3  E=0  | pre-tx CD4 nadir median=  262 (sick while waiting)  (16s)
[panlandscape:start done in 918s]

=====================================================================
panlandscape: pw  (60 worlds seed 777, 60y, atol 1e-7)

### E3 purge x lockdown-window landscape  (purge/win; die = window beats intensity)
purge\win | w180 | w365 | w730 | w1095
  p0.50  | C19/20 | C18/20 | C19/20 | C19/20
  p1.00  | C18/20 | C19/20 | C18/20 | C18/20
  p1.50  | C19/20 | C19/20 | C19/20 | C18/20
  p2.50  | C18/20 | C20/20 | C20/20 | C19/20
[panlandscape:pw done in 1611s]

### E4 PAN minus the ART tail
ART-on (control) C=58/60  R=2  E=0
ART-free     C=60/60  R=0  E=0

### E4b ART-free PAN: cross-validation 2 seeds x 2 tolerance families
tight: ART-on  C=115/120 R=1 E=4 medCD4=1023
tight: ART-free C=117/120 R=0 E=3 medCD4=1023
loose: ART-on  C=115/120 R=4 E=1 medCD4=1023
loose: ART-free C=115/120 R=3 E=2 medCD4=1023

### E6 cross-integrator definitive verdict (3 solvers x 2 tol families, 60 worlds, seed 777)

### E6b ART-free PAN cross-solver (3 methods, tight tol, 60 worlds, seed 777)
ART-FREE PAN: PROVABLY-CURED (all 3 solvers)=47/60  definitive-relapse=0  solver-ambiguous=13  (534s)
  ambiguous world 9: {'LSODA': 'C', 'Radau': 'C', 'BDF': 'R'}
  ambiguous world 11: {'LSODA': 'C', 'Radau': 'R', 'BDF': 'E'}
  ambiguous world 13: {'LSODA': 'C', 'Radau': 'C', 'BDF': 'E'}
  ambiguous world 14: {'LSODA': 'C', 'Radau': 'R', 'BDF': 'C'}
  ambiguous world 15: {'LSODA': 'E', 'Radau': 'C', 'BDF': 'E'}
  ambiguous world 23: {'LSODA': 'C', 'Radau': 'C', 'BDF': 'R'}
  ambiguous world 33: {'LSODA': 'C', 'Radau': 'R', 'BDF': 'R'}
  ambiguous world 34: {'LSODA': 'C', 'Radau': 'C', 'BDF': 'E'}
  ambiguous world 36: {'LSODA': 'C', 'Radau': 'E', 'BDF': 'C'}
  ambiguous world 44: {'LSODA': 'C', 'Radau': 'C', 'BDF': 'R'}
  ambiguous world 48: {'LSODA': 'C', 'Radau': 'R', 'BDF': 'C'}
  ambiguous world 53: {'LSODA': 'C', 'Radau': 'C', 'BDF': 'R'}
  ambiguous world 56: {'LSODA': 'C', 'Radau': 'C', 'BDF': 'E'}

### E6c ambiguous-13 clinical readout (worst-case V max after yr 3, each solver)
CLINICAL verdict (worst-case over 3 solvers: V<50 copies/mL from yr3 & CD4>=500): cured 52/60
  not-clinical-cured world 9: worstV=3.8e+05 CD4=53 {'LSODA': 1.915842294553287e-09, 'Radau': 1.6950234600121813e-25, 'BDF': 375527.46490573656}
  not-clinical-cured world 11: worstV=1.7e+05 CD4=134 {'LSODA': 1.124596793880396e-54, 'Radau': 168491.64900720282}
  not-clinical-cured world 14: worstV=8.7e+05 CD4=39 {'LSODA': 2.2584507960593646e-55, 'Radau': 870020.6949509839, 'BDF': 2.1428159062594885e-22}
  not-clinical-cured world 23: worstV=2.2e+05 CD4=69 {'LSODA': 2.960154247138074e-55, 'Radau': 6.101533209479125e-26, 'BDF': 218346.46731330437}
  not-clinical-cured world 33: worstV=8.6e+05 CD4=84 {'LSODA': 8.893724281264013e-54, 'Radau': 763204.374330759, 'BDF': 861416.2796318969}
  not-clinical-cured world 44: worstV=1.7e+05 CD4=171 {'LSODA': 6.310440682089703e-55, 'Radau': 4.0964460187574634e-26, 'BDF': 167785.33343758347}
  not-clinical-cured world 48: worstV=2.5e+05 CD4=86 {'LSODA': 3.875580124881273e-10, 'Radau': 250752.10468458544, 'BDF': 1.2358180979518312e-22}
  not-clinical-cured world 53: worstV=5.1e+05 CD4=161 {'LSODA': 1.6623762629029786e-55, 'Radau': 3.32250720252586e-26, 'BDF': 505667.4276796779}

### E6d seed 4242 clinical cross-solver readout
seed 4242 ART-FREE PAN: clinical worst-case cure 54/60, any-solver-relapse=6, ambiguous-with-viremia=0

### E7 240-world extended census (LSODA, strict + clinical, ART-free PAN-2)
240-patient PAN-2: STRICT = 233/240, CLINICAL = 233/240, clinical-ember(50<V<=1e-3)=0, relapse=7  (608s)

### E8 PAN-2 discrete-cell verdict on 240 patients (reservoir at lockdown-end, day 1460)
PAN-2 @lockdown-end: reservoir median=0.0e+00/mL max=2.1e-70/mL; worlds with >=1 whole latent cell: 0/240; discrete P(cure): min=1.0000 mean=1.000000; expected discrete relapses=0.0000 of 240 (20s)

### FINAL VERDICT: PAN-2 across 120 patients (seeds [777, 4242])
[1] ODE (LSODA): STRICT 117/120   CLINICAL 117/120   else 3/120   medCD4 1023   medD 5
[2] 3-solver (seed 777): PROVABLY-CURED 47/60, solver-ambiguous 13/60
[3] discrete-cell: worlds with >=1 whole latent cell at lockdown-end: 0/240;  P(cure) mean 1.00000  -> discrete relapses 0.0000
[4] immune recovery: median CD4 at year 60 = 1023 (>=500 cured); immune damage D healed to median 5 vs AIDS threshold 30000
[pan2 final verdict in 851s]

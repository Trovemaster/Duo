! Lorenzo Lodi, 6 October 2015; 24 September 2015; 11 August 2015
! This Fortran module contains data on atoms and nuclei such as masses, spin, natural isotopic abundances
! Covers 673 nuclides, which include all stable isotopes and most (perhaps all) radioactive ones with a half-life
! greater than 1 day. Data come from the AME2012 and NUBASE2012 databases (references reported below).
! In the following: Z=1,2,... is the nuclear charge; A=1,2,... the atomic mass number;
!   m=0,1,2,3 is used for a few radioactive nuclides with metastable excited nuclear states (nucler isomers)
!             m=0 is the nuclear ground state, m=1 the first excited etc.
!
! In all functions provided A and m are optional; if they are not specified the most abundant isotope is returned,
!  or the longest-lived for radioactive ones.
!  
! The following functions are provided (remember: A and m are always optional; m is relevant only for some radioactive elements):
!
!                             OUTPUT FORMAT
! -- element_name(Z, A)     [character, length=15]  => English name of the element. A  is used only for deuterium/tritium
! -- element_symbol(Z, A)   [character, length=15]  => Symbol the element. A  is used only for deuterium/tritium
! -- nuclide_name(Z,A,m)    [character, length=25]  => Full English name and A number, e.g. Carbon-14
! -- atomic_mass(Z,A,m)     [double precision   ]   => atomic mass in Daltons. Returns zero if not found.
! -- atomic_mass_unc(Z,A,m) [double precision   ]   => uncertainty on the atomic mass in Daltons. Returns zero if not found.
! -- stability_info(Z,A,m)  [character,length=60]   => describes if it's stable, radioactive etc.
! -- half_life(Z,A)         [double precision   ]   => half-life in seconds
! -- half_life_txt(Z,A)     [character,length=30]   => half-life as text in readable form (years, months, days etc.)
! -- nuclear_spin(Z,A)      [double precision   ]   => nuclear spin. Returns -99 if not found.
! -- nat_iso_abundance(Z,A) [double precision   ]   => natural isotopic adundance
! -- get_z("symbol")        [integer            ]   => reads in a chemical symbol in the form "C-14", "Na" and gives Z; zero on fail
! -- get_a("symbol")        [integer            ]   => reads in a chemical symbol in the form "C-14", "Na" and gives A; zero on fail
! -- get_m("symbol")        [integer            ]   => reads in a chemical symbol in the form "C-14", "Na" and gives m; -1   on fail
!
! The following subroutine is provided:
! -- print_atomic_and_nuclear_info(string,uu) [lots of text] => prints on unit uu all information available for the given nuclide.
!                                                       `string' is character(len=8) and contains the symbol in the format C-14 etc
!
module atomic_and_nuclear_data
use accuracy
!
implicit none
private
! public atomic_mass, nuclear_spin, print_atomic_and_nuclear_info ! make available only what I need for DUO
!
public element_name, element_symbol, nuclide_name, atomic_mass, atomic_mass_unc, stability_info, half_life, half_life_txt, &
       nuclear_spin, nat_iso_abundance, get_z, get_a, get_m, print_atomic_and_nuclear_info, get_name_from_mass
!
integer, parameter :: dp=rk ! kind(1.d0)
integer, parameter :: n_nuclides=673 ! number of nuclides in module
integer, parameter :: zmax=105       ! maximum Z in list
character(len=15)  :: Name(zmax)     ! English names of elements
character(len=2)   :: Symbol(zmax)   ! Symbols of elements
integer       :: A(       0:n_nuclides)  ! atomic number
integer       :: Z(       0:n_nuclides)  ! nuclear charge
real(kind=dp) :: mass(    0:n_nuclides)  ! masses in Daltons
real(kind=dp) :: mass_unc(0:n_nuclides)  ! error in masses in Daltons
real(kind=dp) :: h_f(     0:n_nuclides)  ! half-lives in second
real(kind=dp) :: spin(    0:n_nuclides)  ! spin in unit of hbar
!real(kind=dp) :: parity(  0:n_nuclides)  ! nuclear parity =+1 (plus) or -1 (minus) (to be added)
real(kind=dp) :: nat_iso_ab(  0:n_nuclides)
integer       :: flag(    0:n_nuclides)  ! used internally to pass info
real(kind=dp), parameter :: long_time = 1e42_dp ! estimated SUSY proton decay half-life (half-life of stable nuclides).
                                                ! Experimental lower limit is ~6e33 y (=2e41 s) [Phys Rev D 90, 072005 (2014)]
real(kind=dp), parameter ::  spin_undef=-99._dp ! value when nuclear spin is unknown

integer :: n_core(  0:n_nuclides)  ! number of core electrons (usual definition based of noble-gas filled shells)
integer :: n_core2( 0:n_nuclides)  ! number of inner-core electrons (definition based on orbital energies < -1.5 Eh)
integer :: iverbose = 4

! English names of chemical elements up to zmax
! Deuterium and Tritium are dealt with as special cases
data Name(  1) /'Hydrogen     '/, Symbol(  1) /'H' /,n_core(  1) / 0/,n_core2(  1) / 0/
data Name(  2) /'Helium       '/, Symbol(  2) /'He'/,n_core(  2) / 0/,n_core2(  2) / 0/
data Name(  3) /'Lithium      '/, Symbol(  3) /'Li'/,n_core(  3) / 2/,n_core2(  3) / 2/
data Name(  4) /'Beryllium    '/, Symbol(  4) /'Be'/,n_core(  4) / 2/,n_core2(  4) / 2/
data Name(  5) /'Boron        '/, Symbol(  5) /'B' /,n_core(  5) / 2/,n_core2(  5) / 2/
data Name(  6) /'Carbon       '/, Symbol(  6) /'C' /,n_core(  6) / 2/,n_core2(  6) / 2/
data Name(  7) /'Nitrogen     '/, Symbol(  7) /'N' /,n_core(  7) / 2/,n_core2(  7) / 2/
data Name(  8) /'Oxygen       '/, Symbol(  8) /'O' /,n_core(  8) / 2/,n_core2(  8) / 2/
data Name(  9) /'Fluorine     '/, Symbol(  9) /'F' /,n_core(  9) / 2/,n_core2(  9) / 2/
data Name( 10) /'Neon         '/, Symbol( 10) /'Ne'/,n_core( 10) / 2/,n_core2( 10) / 2/
data Name( 11) /'Sodium       '/, Symbol( 11) /'Na'/,n_core( 11) /10/,n_core2( 11) / 4/
data Name( 12) /'Magnesium    '/, Symbol( 12) /'Mg'/,n_core( 12) /10/,n_core2( 12) /10/
data Name( 13) /'Aluminium    '/, Symbol( 13) /'Al'/,n_core( 13) /10/,n_core2( 13) /10/
data Name( 14) /'Silicon      '/, Symbol( 14) /'Si'/,n_core( 14) /10/,n_core2( 14) /10/
data Name( 15) /'Phosphorus   '/, Symbol( 15) /'P' /,n_core( 15) /10/,n_core2( 15) /10/
data Name( 16) /'Sulfur       '/, Symbol( 16) /'S' /,n_core( 16) /10/,n_core2( 16) /10/
data Name( 17) /'Chlorine     '/, Symbol( 17) /'Cl'/,n_core( 17) /10/,n_core2( 17) /10/
data Name( 18) /'Argon        '/, Symbol( 18) /'Ar'/,n_core( 18) /10/,n_core2( 18) /10/
data Name( 19) /'Potassium    '/, Symbol( 19) /'K' /,n_core( 19) /18/,n_core2( 19) /10/
data Name( 20) /'Calcium      '/, Symbol( 20) /'Ca'/,n_core( 20) /18/,n_core2( 20) /12/
data Name( 21) /'Scandium     '/, Symbol( 21) /'Sc'/,n_core( 21) /18/,n_core2( 21) /12/
data Name( 22) /'Titanium     '/, Symbol( 22) /'Ti'/,n_core( 22) /18/,n_core2( 22) /12/
data Name( 23) /'Vanadium     '/, Symbol( 23) /'V' /,n_core( 23) /18/,n_core2( 23) /18/
data Name( 24) /'Chromium     '/, Symbol( 24) /'Cr'/,n_core( 24) /18/,n_core2( 24) /18/
data Name( 25) /'Manganese    '/, Symbol( 25) /'Mn'/,n_core( 25) /18/,n_core2( 25) /18/
data Name( 26) /'Iron         '/, Symbol( 26) /'Fe'/,n_core( 26) /18/,n_core2( 26) /18/
data Name( 27) /'Cobalt       '/, Symbol( 27) /'Co'/,n_core( 27) /18/,n_core2( 27) /18/
data Name( 28) /'Nickel       '/, Symbol( 28) /'Ni'/,n_core( 28) /18/,n_core2( 28) /18/
data Name( 29) /'Copper       '/, Symbol( 29) /'Cu'/,n_core( 29) /18/,n_core2( 29) /18/
data Name( 30) /'Zinc         '/, Symbol( 30) /'Zn'/,n_core( 30) /18/,n_core2( 30) /18/
data Name( 31) /'Gallium      '/, Symbol( 31) /'Ga'/,n_core( 31) /18/,n_core2( 31) /18/
data Name( 32) /'Germanium    '/, Symbol( 32) /'Ge'/,n_core( 32) /18/,n_core2( 32) /18/
data Name( 33) /'Arsenic      '/, Symbol( 33) /'As'/,n_core( 33) /18/,n_core2( 33) /28/
data Name( 34) /'Selenium     '/, Symbol( 34) /'Se'/,n_core( 34) /18/,n_core2( 34) /28/
data Name( 35) /'Bromine      '/, Symbol( 35) /'Br'/,n_core( 35) /18/,n_core2( 35) /28/
data Name( 36) /'Krypton      '/, Symbol( 36) /'Kr'/,n_core( 36) /18/,n_core2( 36) /28/
data Name( 37) /'Rubidium     '/, Symbol( 37) /'Rb'/,n_core( 37) /36/,n_core2( 37) /28/
data Name( 38) /'Strontium    '/, Symbol( 38) /'Sr'/,n_core( 38) /36/,n_core2( 38) /30/
data Name( 39) /'Yttrium      '/, Symbol( 39) /'Y' /,n_core( 39) /36/,n_core2( 39) /30/
data Name( 40) /'Zirconium    '/, Symbol( 40) /'Zr'/,n_core( 40) /36/,n_core2( 40) /30/
data Name( 41) /'Niobium      '/, Symbol( 41) /'Nb'/,n_core( 41) /36/,n_core2( 41) /30/
data Name( 42) /'Molybdenum   '/, Symbol( 42) /'Mo'/,n_core( 42) /36/,n_core2( 42) /30/
data Name( 43) /'Technetium   '/, Symbol( 43) /'Tc'/,n_core( 43) /36/,n_core2( 43) /36/
data Name( 44) /'Ruthenium    '/, Symbol( 44) /'Ru'/,n_core( 44) /36/,n_core2( 44) /36/
data Name( 45) /'Rhodium      '/, Symbol( 45) /'Rh'/,n_core( 45) /36/,n_core2( 45) /36/
data Name( 46) /'Palladium    '/, Symbol( 46) /'Pd'/,n_core( 46) /36/,n_core2( 46) /36/
data Name( 47) /'Silver       '/, Symbol( 47) /'Ag'/,n_core( 47) /36/,n_core2( 47) /36/
data Name( 48) /'Cadmium      '/, Symbol( 48) /'Cd'/,n_core( 48) /36/,n_core2( 48) /36/
data Name( 49) /'Indium       '/, Symbol( 49) /'In'/,n_core( 49) /36/,n_core2( 49) /36/
data Name( 50) /'Tin          '/, Symbol( 50) /'Sn'/,n_core( 50) /36/,n_core2( 50) /36/
data Name( 51) /'Antimony     '/, Symbol( 51) /'Sb'/,n_core( 51) /36/,n_core2( 51) /36/
data Name( 52) /'Tellurium    '/, Symbol( 52) /'Te'/,n_core( 52) /36/,n_core2( 52) /46/
data Name( 53) /'Iodine       '/, Symbol( 53) /'I' /,n_core( 53) /36/,n_core2( 53) /46/
data Name( 54) /'Xenon        '/, Symbol( 54) /'Xe'/,n_core( 54) /36/,n_core2( 54) /46/
data Name( 55) /'Caesium      '/, Symbol( 55) /'Cs'/,n_core( 55) /54/,n_core2( 55) /46/
data Name( 56) /'Barium       '/, Symbol( 56) /'Ba'/,n_core( 56) /54/,n_core2( 56) /46/
data Name( 57) /'Lanthanum    '/, Symbol( 57) /'La'/,n_core( 57) /54/,n_core2( 57) /46/
data Name( 58) /'Cerium       '/, Symbol( 58) /'Ce'/,n_core( 58) /54/,n_core2( 58) /46/
data Name( 59) /'Praseodymium '/, Symbol( 59) /'Pr'/,n_core( 59) /54/,n_core2( 59) /46/
data Name( 60) /'Neodymium    '/, Symbol( 60) /'Nd'/,n_core( 60) /54/,n_core2( 60) /46/
data Name( 61) /'Promethium   '/, Symbol( 61) /'Pm'/,n_core( 61) /54/,n_core2( 61) /48/
data Name( 62) /'Samarium     '/, Symbol( 62) /'Sm'/,n_core( 62) /54/,n_core2( 62) /48/
data Name( 63) /'Europium     '/, Symbol( 63) /'Eu'/,n_core( 63) /54/,n_core2( 63) /48/
data Name( 64) /'Gadolinium   '/, Symbol( 64) /'Gd'/,n_core( 64) /54/,n_core2( 64) /48/
data Name( 65) /'Terbium      '/, Symbol( 65) /'Tb'/,n_core( 65) /54/,n_core2( 65) /48/
data Name( 66) /'Dysprosium   '/, Symbol( 66) /'Dy'/,n_core( 66) /54/,n_core2( 66) /48/
data Name( 67) /'Holmium      '/, Symbol( 67) /'Ho'/,n_core( 67) /54/,n_core2( 67) /48/
data Name( 68) /'Erbium       '/, Symbol( 68) /'Er'/,n_core( 68) /54/,n_core2( 68) /48/
data Name( 69) /'Thulium      '/, Symbol( 69) /'Tm'/,n_core( 69) /54/,n_core2( 69) /48/
data Name( 70) /'Ytterbium    '/, Symbol( 70) /'Yb'/,n_core( 70) /54/,n_core2( 70) /48/
data Name( 71) /'Lutetium     '/, Symbol( 71) /'Lu'/,n_core( 71) /54/,n_core2( 71) /48/
data Name( 72) /'Hafnium      '/, Symbol( 72) /'Hf'/,n_core( 72) /54/,n_core2( 72) /48/
data Name( 73) /'Tantalum     '/, Symbol( 73) /'Ta'/,n_core( 73) /54/,n_core2( 73) /48/
data Name( 74) /'Tungsten     '/, Symbol( 74) /'W' /,n_core( 74) /54/,n_core2( 74) /54/
data Name( 75) /'Rhenium      '/, Symbol( 75) /'Re'/,n_core( 75) /54/,n_core2( 75) /68/
data Name( 76) /'Osmium       '/, Symbol( 76) /'Os'/,n_core( 76) /54/,n_core2( 76) /68/
data Name( 77) /'Iridium      '/, Symbol( 77) /'Ir'/,n_core( 77) /54/,n_core2( 77) /68/
data Name( 78) /'Platinum     '/, Symbol( 78) /'Pt'/,n_core( 78) /54/,n_core2( 78) /68/
data Name( 79) /'Gold         '/, Symbol( 79) /'Au'/,n_core( 79) /54/,n_core2( 79) /68/
data Name( 80) /'Mercury      '/, Symbol( 80) /'Hg'/,n_core( 80) /54/,n_core2( 80) /68/
data Name( 81) /'Thallium     '/, Symbol( 81) /'Tl'/,n_core( 81) /54/,n_core2( 81) /68/
data Name( 82) /'Lead         '/, Symbol( 82) /'Pb'/,n_core( 82) /54/,n_core2( 82) /68/
data Name( 83) /'Bismuth      '/, Symbol( 83) /'Bi'/,n_core( 83) /54/,n_core2( 83) /68/
data Name( 84) /'Polonium     '/, Symbol( 84) /'Po'/,n_core( 84) /54/,n_core2( 84) /68/
data Name( 85) /'Astatine     '/, Symbol( 85) /'At'/,n_core( 85) /54/,n_core2( 85) /68/
data Name( 86) /'Radon        '/, Symbol( 86) /'Rn'/,n_core( 86) /54/,n_core2( 86) /78/
data Name( 87) /'Francium     '/, Symbol( 87) /'Fr'/,n_core( 87) /86/,n_core2( 87) /78/
data Name( 88) /'Radium       '/, Symbol( 88) /'Ra'/,n_core( 88) /86/,n_core2( 88) /78/
data Name( 89) /'Actinium     '/, Symbol( 89) /'Ac'/,n_core( 89) /86/,n_core2( 89) /78/
data Name( 90) /'Thorium      '/, Symbol( 90) /'Th'/,n_core( 90) /86/,n_core2( 90) /80/
data Name( 91) /'Protactinium '/, Symbol( 91) /'Pa'/,n_core( 91) /86/,n_core2( 91) /80/
data Name( 92) /'Uranium      '/, Symbol( 92) /'U' /,n_core( 92) /86/,n_core2( 92) /80/
data Name( 93) /'Neptunium    '/, Symbol( 93) /'Np'/,n_core( 93) /86/,n_core2( 93) /80/
data Name( 94) /'Plutonium    '/, Symbol( 94) /'Pu'/,n_core( 94) /86/,n_core2( 94) /80/
data Name( 95) /'Americium    '/, Symbol( 95) /'Am'/,n_core( 95) /86/,n_core2( 95) /80/
data Name( 96) /'Curium       '/, Symbol( 96) /'Cm'/,n_core( 96) /86/,n_core2( 96) /86/
data Name( 97) /'Berkelium    '/, Symbol( 97) /'Bk'/,n_core( 97) /86/,n_core2( 97) /86/
data Name( 98) /'Californium  '/, Symbol( 98) /'Cf'/,n_core( 98) /86/,n_core2( 98) /86/
data Name( 99) /'Einsteinium  '/, Symbol( 99) /'Es'/,n_core( 99) /86/,n_core2( 99) /86/
data Name(100) /'Fermium      '/, Symbol(100) /'Fm'/,n_core(100) /86/,n_core2(100) /86/
data Name(101) /'Mendelevium  '/, Symbol(101) /'Md'/,n_core(101) /86/,n_core2(101) /86/
data Name(102) /'Nobelium     '/, Symbol(102) /'No'/,n_core(102) /86/,n_core2(102) /86/
data Name(103) /'Lawrencium   '/, Symbol(103) /'Lr'/,n_core(103) /86/,n_core2(103) /86/
data Name(104) /'Rutherfordium'/, Symbol(104) /'Rf'/,n_core(104) /86/,n_core2(104) /86/
data Name(105) /'Dubnium      '/, Symbol(105) /'Db'/,n_core(105) /86/,n_core2(105) /86/


! The list comes originally from Wikipedia https://en.wikipedia.org/wiki/List_of_nuclides
! However, all data were taken from the AME2012 database (masses) or from
! the NUBASE2012 one (half-lives, nuclear spins, abundances). Here are the references:
! Atomic mass data, from:
!
!`M. Wang G. Audi, A. H. Wapstra, F. G. Kondev, M. MacCormick, X. Xu and B. Pfeieffer,
!   `The AME2012 atomic mass evaluation (II). Tables, graphs and references',
!   Chinese Physics C, Vol. 36 No. 12 (Dec. 2012), pag. 1603-2014
!
!`G. Audi, F. G. Kondev, M. Wang,   B. Pfeieffer, X. Sun, J. Blachot and M. MacCormick,
!   `The NUBASE2012 evaluation of nuclear properties',
!   Chinese Physics C, Vol. 36 No. 12 (Dec. 2012), pag. 1157-1286

! Data copied from file 'mass.mas12' downloaded from  https://www-nds.iaea.org/amdc/
! The half-lives from Wikipedia agree very well (to 1% or so) with NUBASE2012.
! 
! To the Wikipedia list I added data for 33 metastable nuclear isomers with half-lives > 1h.
! However, masses for these are not reported in AME2012.
! For masses I used the same mass of the one for the corresponing (Z,A) nuclid in AME2012
! (normally the ground state); the uncertainty has been increased by 1e-3 Daltons,
! which is the mass-equivalent of 1 MeV (an estimation of the exitation energy of the metastable nucleus).
! Data are sorted according to:
! 1) Z quantum number
! 2) Natural isotopic nat_iso_ab
! 3) Half-life
! This is useful so that when optional argument are used the most abundant isotope is returned by default,
! and for radioactive ones the longest-lived one is returned by default.
! nuclide 0 has all values set to zero (unknown type)
data Z(  0) /  0/, A(  0) /  0/, mass(  0) /        0.0_dp/, mass_unc(  0) /      0.0_dp/
data Z(  1) /  1/, A(  1) /  1/, mass(  1) /  1.00782503223_dp/, mass_unc(  1) /0.00000000009_dp/  ! 1H
data Z(  2) /  1/, A(  2) /  2/, mass(  2) /  2.01410177812_dp/, mass_unc(  2) /0.00000000012_dp/  ! 2H
data Z(  3) /  1/, A(  3) /  3/, mass(  3) /  3.01604927791_dp/, mass_unc(  3) /0.00000000237_dp/  ! 3H
data Z(  4) /  2/, A(  4) /  4/, mass(  4) /  4.00260325413_dp/, mass_unc(  4) /0.00000000006_dp/  ! 4He
data Z(  5) /  2/, A(  5) /  3/, mass(  5) /  3.01602932008_dp/, mass_unc(  5) /0.00000000250_dp/  ! 3He
data Z(  6) /  3/, A(  6) /  7/, mass(  6) /  7.01600343659_dp/, mass_unc(  6) /0.00000000454_dp/  ! 7Li
data Z(  7) /  3/, A(  7) /  6/, mass(  7) /  6.01512288742_dp/, mass_unc(  7) /0.00000000155_dp/  ! 6Li
data Z(  8) /  4/, A(  8) /  9/, mass(  8) /  9.01218306500_dp/, mass_unc(  8) /0.00000008200_dp/  ! 9Be
data Z(  9) /  4/, A(  9) / 10/, mass(  9) / 10.01353469500_dp/, mass_unc(  9) /0.00000008600_dp/  ! 10Be
data Z( 10) /  4/, A( 10) /  7/, mass( 10) /  7.01692871700_dp/, mass_unc( 10) /0.00000007600_dp/  ! 7Be
data Z( 11) /  5/, A( 11) / 11/, mass( 11) / 11.00930535500_dp/, mass_unc( 11) /0.00000044600_dp/  ! 11B
data Z( 12) /  5/, A( 12) / 10/, mass( 12) / 10.01293694900_dp/, mass_unc( 12) /0.00000041200_dp/  ! 10B
data Z( 13) /  6/, A( 13) / 12/, mass( 13) / 12.00000000000_dp/, mass_unc( 13) /0.00000000000_dp/  ! 12C
data Z( 14) /  6/, A( 14) / 13/, mass( 14) / 13.00335483507_dp/, mass_unc( 14) /0.00000000023_dp/  ! 13C
data Z( 15) /  6/, A( 15) / 14/, mass( 15) / 14.00324198843_dp/, mass_unc( 15) /0.00000000403_dp/  ! 14C
data Z( 16) /  7/, A( 16) / 14/, mass( 16) / 14.00307400443_dp/, mass_unc( 16) /0.00000000020_dp/  ! 14N
data Z( 17) /  7/, A( 17) / 15/, mass( 17) / 15.00010889888_dp/, mass_unc( 17) /0.00000000064_dp/  ! 15N
data Z( 18) /  8/, A( 18) / 16/, mass( 18) / 15.99491461957_dp/, mass_unc( 18) /0.00000000017_dp/  ! 16O
data Z( 19) /  8/, A( 19) / 18/, mass( 19) / 17.99915961286_dp/, mass_unc( 19) /0.00000000076_dp/  ! 18O
data Z( 20) /  8/, A( 20) / 17/, mass( 20) / 16.99913175650_dp/, mass_unc( 20) /0.00000000069_dp/  ! 17O
data Z( 21) /  9/, A( 21) / 19/, mass( 21) / 18.99840316273_dp/, mass_unc( 21) /0.00000000092_dp/  ! 19F
data Z( 22) / 10/, A( 22) / 20/, mass( 22) / 19.99244017617_dp/, mass_unc( 22) /0.00000000168_dp/  ! 20Ne
data Z( 23) / 10/, A( 23) / 22/, mass( 23) / 21.99138511400_dp/, mass_unc( 23) /0.00000001800_dp/  ! 22Ne
data Z( 24) / 10/, A( 24) / 21/, mass( 24) / 20.99384668500_dp/, mass_unc( 24) /0.00000004100_dp/  ! 21Ne
data Z( 25) / 11/, A( 25) / 23/, mass( 25) / 22.98976928196_dp/, mass_unc( 25) /0.00000000194_dp/  ! 23Na
data Z( 26) / 11/, A( 26) / 22/, mass( 26) / 21.99443741100_dp/, mass_unc( 26) /0.00000018300_dp/  ! 22Na
data Z( 27) / 12/, A( 27) / 24/, mass( 27) / 23.98504169700_dp/, mass_unc( 27) /0.00000001400_dp/  ! 24Mg
data Z( 28) / 12/, A( 28) / 26/, mass( 28) / 25.98259296800_dp/, mass_unc( 28) /0.00000003100_dp/  ! 26Mg
data Z( 29) / 12/, A( 29) / 25/, mass( 29) / 24.98583697600_dp/, mass_unc( 29) /0.00000005000_dp/  ! 25Mg
data Z( 30) / 13/, A( 30) / 27/, mass( 30) / 26.98153853100_dp/, mass_unc( 30) /0.00000011100_dp/  ! 27Al
data Z( 31) / 13/, A( 31) / 26/, mass( 31) / 25.98689190400_dp/, mass_unc( 31) /0.00000006900_dp/  ! 26Al
data Z( 32) / 14/, A( 32) / 28/, mass( 32) / 27.97692653465_dp/, mass_unc( 32) /0.00000000044_dp/  ! 28Si
data Z( 33) / 14/, A( 33) / 29/, mass( 33) / 28.97649466490_dp/, mass_unc( 33) /0.00000000052_dp/  ! 29Si
data Z( 34) / 14/, A( 34) / 30/, mass( 34) / 29.97377013600_dp/, mass_unc( 34) /0.00000002300_dp/  ! 30Si
data Z( 35) / 14/, A( 35) / 32/, mass( 35) / 31.97415153900_dp/, mass_unc( 35) /0.00000032000_dp/  ! 32Si
data Z( 36) / 15/, A( 36) / 31/, mass( 36) / 30.97376199842_dp/, mass_unc( 36) /0.00000000070_dp/  ! 31P
data Z( 37) / 15/, A( 37) / 33/, mass( 37) / 32.97172569400_dp/, mass_unc( 37) /0.00000117000_dp/  ! 33P
data Z( 38) / 15/, A( 38) / 32/, mass( 38) / 31.97390764300_dp/, mass_unc( 38) /0.00000004200_dp/  ! 32P
data Z( 39) / 16/, A( 39) / 32/, mass( 39) / 31.97207117441_dp/, mass_unc( 39) /0.00000000141_dp/  ! 32S
data Z( 40) / 16/, A( 40) / 34/, mass( 40) / 33.96786700400_dp/, mass_unc( 40) /0.00000004700_dp/  ! 34S
data Z( 41) / 16/, A( 41) / 33/, mass( 41) / 32.97145890982_dp/, mass_unc( 41) /0.00000000145_dp/  ! 33S
data Z( 42) / 16/, A( 42) / 36/, mass( 42) / 35.96708070600_dp/, mass_unc( 42) /0.00000020100_dp/  ! 36S
data Z( 43) / 16/, A( 43) / 35/, mass( 43) / 34.96903231000_dp/, mass_unc( 43) /0.00000004300_dp/  ! 35S
data Z( 44) / 17/, A( 44) / 35/, mass( 44) / 34.96885268200_dp/, mass_unc( 44) /0.00000003700_dp/  ! 35Cl
data Z( 45) / 17/, A( 45) / 37/, mass( 45) / 36.96590260200_dp/, mass_unc( 45) /0.00000005500_dp/  ! 37Cl
data Z( 46) / 17/, A( 46) / 36/, mass( 46) / 35.96830680900_dp/, mass_unc( 46) /0.00000003800_dp/  ! 36Cl
data Z( 47) / 18/, A( 47) / 40/, mass( 47) / 39.96238312372_dp/, mass_unc( 47) /0.00000000240_dp/  ! 40Ar
data Z( 48) / 18/, A( 48) / 36/, mass( 48) / 35.96754510500_dp/, mass_unc( 48) /0.00000002800_dp/  ! 36Ar
data Z( 49) / 18/, A( 49) / 38/, mass( 49) / 37.96273210600_dp/, mass_unc( 49) /0.00000020900_dp/  ! 38Ar
data Z( 50) / 18/, A( 50) / 39/, mass( 50) / 38.96431303800_dp/, mass_unc( 50) /0.00000536700_dp/  ! 39Ar
data Z( 51) / 18/, A( 51) / 42/, mass( 51) / 41.96304573600_dp/, mass_unc( 51) /0.00000620000_dp/  ! 42Ar
data Z( 52) / 18/, A( 52) / 37/, mass( 52) / 36.96677633100_dp/, mass_unc( 52) /0.00000022100_dp/  ! 37Ar
data Z( 53) / 19/, A( 53) / 39/, mass( 53) / 38.96370648643_dp/, mass_unc( 53) /0.00000000492_dp/  ! 39K
data Z( 54) / 19/, A( 54) / 41/, mass( 54) / 40.96182525792_dp/, mass_unc( 54) /0.00000000408_dp/  ! 41K
data Z( 55) / 19/, A( 55) / 40/, mass( 55) / 39.96399816600_dp/, mass_unc( 55) /0.00000006000_dp/  ! 40K
data Z( 56) / 20/, A( 56) / 40/, mass( 56) / 39.96259086300_dp/, mass_unc( 56) /0.00000002200_dp/  ! 40Ca
data Z( 57) / 20/, A( 57) / 44/, mass( 57) / 43.95548156100_dp/, mass_unc( 57) /0.00000034800_dp/  ! 44Ca
data Z( 58) / 20/, A( 58) / 42/, mass( 58) / 41.95861783000_dp/, mass_unc( 58) /0.00000015900_dp/  ! 42Ca
data Z( 59) / 20/, A( 59) / 48/, mass( 59) / 47.95252276500_dp/, mass_unc( 59) /0.00000012900_dp/  ! 48Ca
data Z( 60) / 20/, A( 60) / 43/, mass( 60) / 42.95876643800_dp/, mass_unc( 60) /0.00000024400_dp/  ! 43Ca
data Z( 61) / 20/, A( 61) / 46/, mass( 61) / 45.95368902300_dp/, mass_unc( 61) /0.00000242600_dp/  ! 46Ca
data Z( 62) / 20/, A( 62) / 41/, mass( 62) / 40.96227792400_dp/, mass_unc( 62) /0.00000014700_dp/  ! 41Ca
data Z( 63) / 20/, A( 63) / 45/, mass( 63) / 44.95618635000_dp/, mass_unc( 63) /0.00000039200_dp/  ! 45Ca
data Z( 64) / 20/, A( 64) / 47/, mass( 64) / 46.95454243000_dp/, mass_unc( 64) /0.00000241300_dp/  ! 47Ca
data Z( 65) / 21/, A( 65) / 45/, mass( 65) / 44.95590827500_dp/, mass_unc( 65) /0.00000077300_dp/  ! 45Sc
data Z( 66) / 21/, A( 66) / 46/, mass( 66) / 45.95516825700_dp/, mass_unc( 66) /0.00000078100_dp/  ! 46Sc
data Z( 67) / 21/, A( 67) / 47/, mass( 67) / 46.95240374000_dp/, mass_unc( 67) /0.00000210600_dp/  ! 47Sc
data Z( 68) / 21/, A( 68) / 44/, mass( 68) / 43.95940287500_dp/, mass_unc( 68) /0.00000188400_dp/  ! 44m3Sc
data Z( 69) / 21/, A( 69) / 48/, mass( 69) / 47.95222361100_dp/, mass_unc( 69) /0.00000531700_dp/  ! 48Sc
data Z( 70) / 22/, A( 70) / 48/, mass( 70) / 47.94794197900_dp/, mass_unc( 70) /0.00000038400_dp/  ! 48Ti
data Z( 71) / 22/, A( 71) / 46/, mass( 71) / 45.95262771800_dp/, mass_unc( 71) /0.00000035100_dp/  ! 46Ti
data Z( 72) / 22/, A( 72) / 47/, mass( 72) / 46.95175878700_dp/, mass_unc( 72) /0.00000038200_dp/  ! 47Ti
data Z( 73) / 22/, A( 73) / 49/, mass( 73) / 48.94786567600_dp/, mass_unc( 73) /0.00000038600_dp/  ! 49Ti
data Z( 74) / 22/, A( 74) / 50/, mass( 74) / 49.94478688900_dp/, mass_unc( 74) /0.00000038900_dp/  ! 50Ti
data Z( 75) / 22/, A( 75) / 44/, mass( 75) / 43.95968994900_dp/, mass_unc( 75) /0.00000075100_dp/  ! 44Ti
data Z( 76) / 23/, A( 76) / 51/, mass( 76) / 50.94395703600_dp/, mass_unc( 76) /0.00000094200_dp/  ! 51V
data Z( 77) / 23/, A( 77) / 50/, mass( 77) / 49.94715601400_dp/, mass_unc( 77) /0.00000094500_dp/  ! 50V
data Z( 78) / 23/, A( 78) / 49/, mass( 78) / 48.94851179500_dp/, mass_unc( 78) /0.00000096100_dp/  ! 49V
data Z( 79) / 23/, A( 79) / 48/, mass( 79) / 47.95225222300_dp/, mass_unc( 79) /0.00000110200_dp/  ! 48V
data Z( 80) / 24/, A( 80) / 52/, mass( 80) / 51.94050623100_dp/, mass_unc( 80) /0.00000063100_dp/  ! 52Cr
data Z( 81) / 24/, A( 81) / 53/, mass( 81) / 52.94064814700_dp/, mass_unc( 81) /0.00000062000_dp/  ! 53Cr
data Z( 82) / 24/, A( 82) / 50/, mass( 82) / 49.94604183300_dp/, mass_unc( 82) /0.00000094300_dp/  ! 50Cr
data Z( 83) / 24/, A( 83) / 54/, mass( 83) / 53.93887915800_dp/, mass_unc( 83) /0.00000061100_dp/  ! 54Cr
data Z( 84) / 24/, A( 84) / 51/, mass( 84) / 50.94476501800_dp/, mass_unc( 84) /0.00000094300_dp/  ! 51Cr
data Z( 85) / 25/, A( 85) / 55/, mass( 85) / 54.93804391000_dp/, mass_unc( 85) /0.00000047600_dp/  ! 55Mn
data Z( 86) / 25/, A( 86) / 53/, mass( 86) / 52.94128889100_dp/, mass_unc( 86) /0.00000068400_dp/  ! 53Mn
data Z( 87) / 25/, A( 87) / 54/, mass( 87) / 53.94035761500_dp/, mass_unc( 87) /0.00000124000_dp/  ! 54Mn
data Z( 88) / 25/, A( 88) / 52/, mass( 88) / 51.94556394900_dp/, mass_unc( 88) /0.00000204000_dp/  ! 52Mn
data Z( 89) / 26/, A( 89) / 56/, mass( 89) / 55.93493632600_dp/, mass_unc( 89) /0.00000048900_dp/  ! 56Fe
data Z( 90) / 26/, A( 90) / 54/, mass( 90) / 53.93960898600_dp/, mass_unc( 90) /0.00000052900_dp/  ! 54Fe
data Z( 91) / 26/, A( 91) / 57/, mass( 91) / 56.93539284100_dp/, mass_unc( 91) /0.00000049000_dp/  ! 57Fe
data Z( 92) / 26/, A( 92) / 58/, mass( 92) / 57.93327443100_dp/, mass_unc( 92) /0.00000052700_dp/  ! 58Fe
data Z( 93) / 26/, A( 93) / 60/, mass( 93) / 59.93407110000_dp/, mass_unc( 93) /0.00000367800_dp/  ! 60Fe
data Z( 94) / 26/, A( 94) / 55/, mass( 94) / 54.93829199400_dp/, mass_unc( 94) /0.00000050600_dp/  ! 55Fe
data Z( 95) / 26/, A( 95) / 59/, mass( 95) / 58.93487433800_dp/, mass_unc( 95) /0.00000054000_dp/  ! 59Fe
data Z( 96) / 27/, A( 96) / 59/, mass( 96) / 58.93319428800_dp/, mass_unc( 96) /0.00000055600_dp/  ! 59Co
data Z( 97) / 27/, A( 97) / 60/, mass( 97) / 59.93381629900_dp/, mass_unc( 97) /0.00000056100_dp/  ! 60Co
data Z( 98) / 27/, A( 98) / 57/, mass( 98) / 56.93629057400_dp/, mass_unc( 98) /0.00000066400_dp/  ! 57Co
data Z( 99) / 27/, A( 99) / 56/, mass( 99) / 55.93983879800_dp/, mass_unc( 99) /0.00000062600_dp/  ! 56Co
data Z(100) / 27/, A(100) / 58/, mass(100) / 57.93575207300_dp/, mass_unc(100) /0.00000128900_dp/  ! 58Co
data Z(101) / 28/, A(101) / 58/, mass(101) / 57.93534241400_dp/, mass_unc(101) /0.00000051800_dp/  ! 58Ni
data Z(102) / 28/, A(102) / 60/, mass(102) / 59.93078588500_dp/, mass_unc(102) /0.00000051800_dp/  ! 60Ni
data Z(103) / 28/, A(103) / 62/, mass(103) / 61.92834536500_dp/, mass_unc(103) /0.00000055400_dp/  ! 62Ni
data Z(104) / 28/, A(104) / 61/, mass(104) / 60.93105557000_dp/, mass_unc(104) /0.00000051900_dp/  ! 61Ni
data Z(105) / 28/, A(105) / 64/, mass(105) / 63.92796681600_dp/, mass_unc(105) /0.00000058400_dp/  ! 64Ni
data Z(106) / 28/, A(106) / 59/, mass(106) / 58.93434620200_dp/, mass_unc(106) /0.00000051800_dp/  ! 59Ni
data Z(107) / 28/, A(107) / 63/, mass(107) / 62.92966962600_dp/, mass_unc(107) /0.00000055600_dp/  ! 63Ni
data Z(108) / 28/, A(108) / 56/, mass(108) / 55.94212854900_dp/, mass_unc(108) /0.00000057100_dp/  ! 56Ni
data Z(109) / 28/, A(109) / 66/, mass(109) / 65.92913933400_dp/, mass_unc(109) /0.00000150000_dp/  ! 66Ni
data Z(110) / 28/, A(110) / 57/, mass(110) / 56.93979218400_dp/, mass_unc(110) /0.00000071200_dp/  ! 57Ni
data Z(111) / 29/, A(111) / 63/, mass(111) / 62.92959772300_dp/, mass_unc(111) /0.00000055600_dp/  ! 63Cu
data Z(112) / 29/, A(112) / 65/, mass(112) / 64.92778970400_dp/, mass_unc(112) /0.00000070900_dp/  ! 65Cu
data Z(113) / 29/, A(113) / 67/, mass(113) / 66.92773031400_dp/, mass_unc(113) /0.00000130000_dp/  ! 67Cu
data Z(114) / 30/, A(114) / 64/, mass(114) / 63.92914201300_dp/, mass_unc(114) /0.00000070900_dp/  ! 64Zn
data Z(115) / 30/, A(115) / 66/, mass(115) / 65.92603380900_dp/, mass_unc(115) /0.00000093900_dp/  ! 66Zn
data Z(116) / 30/, A(116) / 68/, mass(116) / 67.92484455400_dp/, mass_unc(116) /0.00000098000_dp/  ! 68Zn
data Z(117) / 30/, A(117) / 67/, mass(117) / 66.92712774600_dp/, mass_unc(117) /0.00000095800_dp/  ! 67Zn
data Z(118) / 30/, A(118) / 70/, mass(118) / 69.92531920800_dp/, mass_unc(118) /0.00000205900_dp/  ! 70Zn
data Z(119) / 30/, A(119) / 65/, mass(119) / 64.92924077000_dp/, mass_unc(119) /0.00000071100_dp/  ! 65Zn
data Z(120) / 30/, A(120) / 72/, mass(120) / 71.92684280700_dp/, mass_unc(120) /0.00000230000_dp/  ! 72Zn
data Z(121) / 31/, A(121) / 69/, mass(121) / 68.92557354100_dp/, mass_unc(121) /0.00000128500_dp/  ! 69Ga
data Z(122) / 31/, A(122) / 71/, mass(122) / 70.92470257700_dp/, mass_unc(122) /0.00000087400_dp/  ! 71Ga
data Z(123) / 31/, A(123) / 67/, mass(123) / 66.92820254700_dp/, mass_unc(123) /0.00000130500_dp/  ! 67Ga
data Z(124) / 32/, A(124) / 74/, mass(124) / 73.92117776100_dp/, mass_unc(124) /0.00000001300_dp/  ! 74Ge
data Z(125) / 32/, A(125) / 72/, mass(125) / 71.92207582600_dp/, mass_unc(125) /0.00000008100_dp/  ! 72Ge
data Z(126) / 32/, A(126) / 70/, mass(126) / 69.92424875000_dp/, mass_unc(126) /0.00000090300_dp/  ! 70Ge
data Z(127) / 32/, A(127) / 73/, mass(127) / 72.92345895600_dp/, mass_unc(127) /0.00000006100_dp/  ! 73Ge
data Z(128) / 32/, A(128) / 76/, mass(128) / 75.92140272600_dp/, mass_unc(128) /0.00000001900_dp/  ! 76Ge
data Z(129) / 32/, A(129) / 68/, mass(129) / 67.92809530700_dp/, mass_unc(129) /0.00000201400_dp/  ! 68Ge
data Z(130) / 32/, A(130) / 71/, mass(130) / 70.92495232700_dp/, mass_unc(130) /0.00000089800_dp/  ! 71Ge
data Z(131) / 32/, A(131) / 69/, mass(131) / 68.92796448100_dp/, mass_unc(131) /0.00000141400_dp/  ! 69Ge
data Z(132) / 33/, A(132) / 75/, mass(132) / 74.92159456700_dp/, mass_unc(132) /0.00000094800_dp/  ! 75As
data Z(133) / 33/, A(133) / 73/, mass(133) / 72.92382908600_dp/, mass_unc(133) /0.00000413600_dp/  ! 73As
data Z(134) / 33/, A(134) / 74/, mass(134) / 73.92392859800_dp/, mass_unc(134) /0.00000181700_dp/  ! 74As
data Z(135) / 33/, A(135) / 71/, mass(135) / 70.92711380100_dp/, mass_unc(135) /0.00000447300_dp/  ! 71As
data Z(136) / 33/, A(136) / 77/, mass(136) / 76.92064756300_dp/, mass_unc(136) /0.00000184200_dp/  ! 77As
data Z(137) / 33/, A(137) / 72/, mass(137) / 71.92675229400_dp/, mass_unc(137) /0.00000438300_dp/  ! 72As
data Z(138) / 33/, A(138) / 76/, mass(138) / 75.92239201500_dp/, mass_unc(138) /0.00000095100_dp/  ! 76As
data Z(139) / 34/, A(139) / 80/, mass(139) / 79.91652176200_dp/, mass_unc(139) /0.00000133900_dp/  ! 80Se
data Z(140) / 34/, A(140) / 78/, mass(140) / 77.91730928000_dp/, mass_unc(140) /0.00000019500_dp/  ! 78Se
data Z(141) / 34/, A(141) / 76/, mass(141) / 75.91921370400_dp/, mass_unc(141) /0.00000001700_dp/  ! 76Se
data Z(142) / 34/, A(142) / 82/, mass(142) / 81.91669949700_dp/, mass_unc(142) /0.00000151600_dp/  ! 82Se
data Z(143) / 34/, A(143) / 77/, mass(143) / 76.91991415400_dp/, mass_unc(143) /0.00000006700_dp/  ! 77Se
data Z(144) / 34/, A(144) / 74/, mass(144) / 73.92247593400_dp/, mass_unc(144) /0.00000001500_dp/  ! 74Se
data Z(145) / 34/, A(145) / 79/, mass(145) / 78.91849928700_dp/, mass_unc(145) /0.00000024100_dp/  ! 79Se
data Z(146) / 34/, A(146) / 75/, mass(146) / 74.92252287000_dp/, mass_unc(146) /0.00000007800_dp/  ! 75Se
data Z(147) / 34/, A(147) / 72/, mass(147) / 71.92714050700_dp/, mass_unc(147) /0.00000210000_dp/  ! 72Se
data Z(148) / 35/, A(148) / 79/, mass(148) / 78.91833757900_dp/, mass_unc(148) /0.00000138700_dp/  ! 79Br
data Z(149) / 35/, A(149) / 81/, mass(149) / 80.91628969000_dp/, mass_unc(149) /0.00000137500_dp/  ! 81Br
data Z(150) / 35/, A(150) / 77/, mass(150) / 76.92137919800_dp/, mass_unc(150) /0.00000301700_dp/  ! 77Br
data Z(151) / 35/, A(151) / 82/, mass(151) / 81.91680324600_dp/, mass_unc(151) /0.00000137000_dp/  ! 82Br
data Z(152) / 36/, A(152) / 84/, mass(152) / 83.91149772816_dp/, mass_unc(152) /0.00000000442_dp/  ! 84Kr
data Z(153) / 36/, A(153) / 86/, mass(153) / 85.91061062693_dp/, mass_unc(153) /0.00000000407_dp/  ! 86Kr
data Z(154) / 36/, A(154) / 82/, mass(154) / 81.91348273000_dp/, mass_unc(154) /0.00000094000_dp/  ! 82Kr
data Z(155) / 36/, A(155) / 83/, mass(155) / 82.91412716400_dp/, mass_unc(155) /0.00000032100_dp/  ! 83Kr
data Z(156) / 36/, A(156) / 80/, mass(156) / 79.91637808400_dp/, mass_unc(156) /0.00000074700_dp/  ! 80Kr
data Z(157) / 36/, A(157) / 78/, mass(157) / 77.92036494400_dp/, mass_unc(157) /0.00000075600_dp/  ! 78Kr
data Z(158) / 36/, A(158) / 81/, mass(158) / 80.91659118100_dp/, mass_unc(158) /0.00000145100_dp/  ! 81Kr
data Z(159) / 36/, A(159) / 85/, mass(159) / 84.91252726200_dp/, mass_unc(159) /0.00000214700_dp/  ! 85Kr
data Z(160) / 36/, A(160) / 79/, mass(160) / 78.92008292300_dp/, mass_unc(160) /0.00000383800_dp/  ! 79Kr
data Z(161) / 37/, A(161) / 85/, mass(161) / 84.91178973790_dp/, mass_unc(161) /0.00000000535_dp/  ! 85Rb
data Z(162) / 37/, A(162) / 87/, mass(162) / 86.90918053100_dp/, mass_unc(162) /0.00000000600_dp/  ! 87Rb
data Z(163) / 37/, A(163) / 83/, mass(163) / 82.91511418300_dp/, mass_unc(163) /0.00000250000_dp/  ! 83Rb
data Z(164) / 37/, A(164) / 84/, mass(164) / 83.91437522900_dp/, mass_unc(164) /0.00000235500_dp/  ! 84Rb
data Z(165) / 37/, A(165) / 86/, mass(165) / 85.91116742500_dp/, mass_unc(165) /0.00000021300_dp/  ! 86Rb
data Z(166) / 38/, A(166) / 88/, mass(166) / 87.90561254200_dp/, mass_unc(166) /0.00000116400_dp/  ! 88Sr
data Z(167) / 38/, A(167) / 86/, mass(167) / 85.90926060800_dp/, mass_unc(167) /0.00000115400_dp/  ! 86Sr
data Z(168) / 38/, A(168) / 87/, mass(168) / 86.90887753100_dp/, mass_unc(168) /0.00000115400_dp/  ! 87Sr
data Z(169) / 38/, A(169) / 84/, mass(169) / 83.91341913600_dp/, mass_unc(169) /0.00000133400_dp/  ! 84Sr
data Z(170) / 38/, A(170) / 90/, mass(170) / 89.90773003700_dp/, mass_unc(170) /0.00000280900_dp/  ! 90Sr
data Z(171) / 38/, A(171) / 85/, mass(171) / 84.91293204600_dp/, mass_unc(171) /0.00000302000_dp/  ! 85Sr
data Z(172) / 38/, A(172) / 89/, mass(172) / 88.90745109500_dp/, mass_unc(172) /0.00000116800_dp/  ! 89Sr
data Z(173) / 38/, A(173) / 82/, mass(173) / 81.91839985500_dp/, mass_unc(173) /0.00000643200_dp/  ! 82Sr
data Z(174) / 38/, A(174) / 83/, mass(174) / 82.91755437400_dp/, mass_unc(174) /0.00000733600_dp/  ! 83Sr
data Z(175) / 39/, A(175) / 89/, mass(175) / 88.90584034800_dp/, mass_unc(175) /0.00000239700_dp/  ! 89Y
data Z(176) / 39/, A(176) / 88/, mass(176) / 87.90950156300_dp/, mass_unc(176) /0.00000198700_dp/  ! 88Y
data Z(177) / 39/, A(177) / 91/, mass(177) / 90.90729744200_dp/, mass_unc(177) /0.00000276100_dp/  ! 91Y
data Z(178) / 39/, A(178) / 87/, mass(178) / 86.91087613800_dp/, mass_unc(178) /0.00000167200_dp/  ! 87Y
data Z(179) / 39/, A(179) / 90/, mass(179) / 89.90714394200_dp/, mass_unc(179) /0.00000239700_dp/  ! 90Y
data Z(180) / 40/, A(180) / 90/, mass(180) / 89.90469765900_dp/, mass_unc(180) /0.00000198400_dp/  ! 90Zr
data Z(181) / 40/, A(181) / 94/, mass(181) / 93.90631082800_dp/, mass_unc(181) /0.00000201900_dp/  ! 94Zr
data Z(182) / 40/, A(182) / 92/, mass(182) / 91.90503467500_dp/, mass_unc(182) /0.00000196200_dp/  ! 92Zr
data Z(183) / 40/, A(183) / 91/, mass(183) / 90.90563958700_dp/, mass_unc(183) /0.00000196400_dp/  ! 91Zr
data Z(184) / 40/, A(184) / 96/, mass(184) / 95.90827143300_dp/, mass_unc(184) /0.00000214100_dp/  ! 96Zr
data Z(185) / 40/, A(185) / 93/, mass(185) / 92.90646994700_dp/, mass_unc(185) /0.00000195200_dp/  ! 93Zr
data Z(186) / 40/, A(186) / 88/, mass(186) / 87.91022129000_dp/, mass_unc(186) /0.00000584400_dp/  ! 88Zr
data Z(187) / 40/, A(187) / 95/, mass(187) / 94.90803853000_dp/, mass_unc(187) /0.00000194000_dp/  ! 95Zr
data Z(188) / 40/, A(188) / 89/, mass(188) / 88.90888144100_dp/, mass_unc(188) /0.00000372200_dp/  ! 89Zr
data Z(189) / 41/, A(189) / 93/, mass(189) / 92.90637300400_dp/, mass_unc(189) /0.00000195800_dp/  ! 93Nb
data Z(190) / 41/, A(190) / 92/, mass(190) / 91.90718808100_dp/, mass_unc(190) /0.00000257000_dp/  ! 92Nb
data Z(191) / 41/, A(191) / 94/, mass(191) / 93.90727883600_dp/, mass_unc(191) /0.00000195800_dp/  ! 94Nb
data Z(192) / 41/, A(192) / 91/, mass(192) / 90.90698965800_dp/, mass_unc(192) /0.00000367900_dp/  ! 91Nb
data Z(193) / 41/, A(193) / 93/, mass(193) / 92.90637300400_dp/, mass_unc(193) /0.00100195800_dp/  ! 93mNb
data Z(194) / 41/, A(194) / 95/, mass(194) / 94.90683240400_dp/, mass_unc(194) /0.00000070900_dp/  ! 95Nb
data Z(195) / 41/, A(195) / 95/, mass(195) / 94.90683240400_dp/, mass_unc(195) /0.00100070900_dp/  ! 95mNb
data Z(196) / 42/, A(196) / 98/, mass(196) / 97.90540482000_dp/, mass_unc(196) /0.00000049300_dp/  ! 98Mo
data Z(197) / 42/, A(197) / 96/, mass(197) / 95.90467611500_dp/, mass_unc(197) /0.00000047300_dp/  ! 96Mo
data Z(198) / 42/, A(198) / 95/, mass(198) / 94.90583876600_dp/, mass_unc(198) /0.00000047300_dp/  ! 95Mo
data Z(199) / 42/, A(199) / 92/, mass(199) / 91.90680795900_dp/, mass_unc(199) /0.00000083700_dp/  ! 92Mo
data Z(200) / 42/, A(200) /100/, mass(200) / 99.90747178200_dp/, mass_unc(200) /0.00000112400_dp/  ! 100Mo
data Z(201) / 42/, A(201) / 97/, mass(201) / 96.90601811800_dp/, mass_unc(201) /0.00000049000_dp/  ! 97Mo
data Z(202) / 42/, A(202) / 94/, mass(202) / 93.90508490300_dp/, mass_unc(202) /0.00000047900_dp/  ! 94Mo
data Z(203) / 42/, A(203) / 93/, mass(203) / 92.90680957700_dp/, mass_unc(203) /0.00000084300_dp/  ! 93Mo
data Z(204) / 42/, A(204) / 99/, mass(204) / 98.90770850900_dp/, mass_unc(204) /0.00000051700_dp/  ! 99Mo
data Z(205) / 43/, A(205) / 97/, mass(205) / 96.90636670600_dp/, mass_unc(205) /0.00000404500_dp/  ! 97Tc
data Z(206) / 43/, A(206) / 98/, mass(206) / 97.90721236500_dp/, mass_unc(206) /0.00000364900_dp/  ! 98Tc
data Z(207) / 43/, A(207) / 99/, mass(207) / 98.90625084400_dp/, mass_unc(207) /0.00000103700_dp/  ! 99Tc
data Z(208) / 43/, A(208) / 97/, mass(208) / 96.90636670600_dp/, mass_unc(208) /0.00100404500_dp/  ! 97mTc
data Z(209) / 43/, A(209) / 95/, mass(209) / 94.90765361200_dp/, mass_unc(209) /0.00000547200_dp/  ! 95mTc
data Z(210) / 43/, A(210) / 96/, mass(210) / 95.90786802200_dp/, mass_unc(210) /0.00000554300_dp/  ! 96Tc
data Z(211) / 44/, A(211) /102/, mass(211) /101.90434409600_dp/, mass_unc(211) /0.00000115500_dp/  ! 102Ru
data Z(212) / 44/, A(212) /104/, mass(212) /103.90542748100_dp/, mass_unc(212) /0.00000275500_dp/  ! 104Ru
data Z(213) / 44/, A(213) /101/, mass(213) /100.90557687200_dp/, mass_unc(213) /0.00000115400_dp/  ! 101Ru
data Z(214) / 44/, A(214) / 99/, mass(214) / 98.90593408200_dp/, mass_unc(214) /0.00000112900_dp/  ! 99Ru
data Z(215) / 44/, A(215) /100/, mass(215) / 99.90421425600_dp/, mass_unc(215) /0.00000112900_dp/  ! 100Ru
data Z(216) / 44/, A(216) / 96/, mass(216) / 95.90759025500_dp/, mass_unc(216) /0.00000049100_dp/  ! 96Ru
data Z(217) / 44/, A(217) / 98/, mass(217) / 97.90528681300_dp/, mass_unc(217) /0.00000693700_dp/  ! 98Ru
data Z(218) / 44/, A(218) /106/, mass(218) /105.90732910400_dp/, mass_unc(218) /0.00000579300_dp/  ! 106Ru
data Z(219) / 44/, A(219) /103/, mass(219) /102.90631862700_dp/, mass_unc(219) /0.00000116500_dp/  ! 103Ru
data Z(220) / 44/, A(220) / 97/, mass(220) / 96.90754712000_dp/, mass_unc(220) /0.00000300000_dp/  ! 97Ru
data Z(221) / 45/, A(221) /103/, mass(221) /102.90549799300_dp/, mass_unc(221) /0.00000261000_dp/  ! 103Rh
data Z(222) / 45/, A(222) /102/, mass(222) /101.90683737300_dp/, mass_unc(222) /0.00000502800_dp/  ! 102mRh
data Z(223) / 45/, A(223) /101/, mass(223) /100.90616061300_dp/, mass_unc(223) /0.00000628600_dp/  ! 101Rh
data Z(224) / 45/, A(224) /102/, mass(224) /101.90683737300_dp/, mass_unc(224) /0.00100502800_dp/  ! 102Rh
data Z(225) / 45/, A(225) / 99/, mass(225) / 98.90812823900_dp/, mass_unc(225) /0.00000725300_dp/  ! 99Rh
data Z(226) / 45/, A(226) /101/, mass(226) /100.90616061300_dp/, mass_unc(226) /0.00100628600_dp/  ! 101mRh
data Z(227) / 45/, A(227) /105/, mass(227) /104.90568854900_dp/, mass_unc(227) /0.00000269400_dp/  ! 105Rh
data Z(228) / 46/, A(228) /106/, mass(228) /105.90348042600_dp/, mass_unc(228) /0.00000118700_dp/  ! 106Pd
data Z(229) / 46/, A(229) /108/, mass(229) /107.90389164000_dp/, mass_unc(229) /0.00000119100_dp/  ! 108Pd
data Z(230) / 46/, A(230) /105/, mass(230) /104.90507962600_dp/, mass_unc(230) /0.00000122200_dp/  ! 105Pd
data Z(231) / 46/, A(231) /110/, mass(231) /109.90517219900_dp/, mass_unc(231) /0.00000074500_dp/  ! 110Pd
data Z(232) / 46/, A(232) /104/, mass(232) /103.90403054000_dp/, mass_unc(232) /0.00000143500_dp/  ! 104Pd
data Z(233) / 46/, A(233) /102/, mass(233) /101.90560218700_dp/, mass_unc(233) /0.00000283300_dp/  ! 102Pd
data Z(234) / 46/, A(234) /107/, mass(234) /106.90512819500_dp/, mass_unc(234) /0.00000129000_dp/  ! 107Pd
data Z(235) / 46/, A(235) /103/, mass(235) /102.90608094900_dp/, mass_unc(235) /0.00000272800_dp/  ! 103Pd
data Z(236) / 46/, A(236) /100/, mass(236) / 99.90850480500_dp/, mass_unc(236) /0.00001899600_dp/  ! 100Pd
data Z(237) / 47/, A(237) /107/, mass(237) /106.90509161100_dp/, mass_unc(237) /0.00000255700_dp/  ! 107Ag
data Z(238) / 47/, A(238) /109/, mass(238) /108.90475528200_dp/, mass_unc(238) /0.00000140500_dp/  ! 109Ag
data Z(239) / 47/, A(239) /108/, mass(239) /107.90595034600_dp/, mass_unc(239) /0.00000256300_dp/  ! 108mAg
data Z(240) / 47/, A(240) /110/, mass(240) /109.90611022600_dp/, mass_unc(240) /0.00000140400_dp/  ! 110m2Ag
data Z(241) / 47/, A(241) /105/, mass(241) /104.90652561500_dp/, mass_unc(241) /0.00000487700_dp/  ! 105Ag
data Z(242) / 47/, A(242) /106/, mass(242) /105.90666363700_dp/, mass_unc(242) /0.00000323700_dp/  ! 106mAg
data Z(243) / 47/, A(243) /111/, mass(243) /110.90529592300_dp/, mass_unc(243) /0.00000163400_dp/  ! 111Ag
data Z(244) / 48/, A(244) /114/, mass(244) /113.90336508600_dp/, mass_unc(244) /0.00000043000_dp/  ! 114Cd
data Z(245) / 48/, A(245) /112/, mass(245) /111.90276286800_dp/, mass_unc(245) /0.00000059500_dp/  ! 112Cd
data Z(246) / 48/, A(246) /111/, mass(246) /110.90418287200_dp/, mass_unc(246) /0.00000060500_dp/  ! 111Cd
data Z(247) / 48/, A(247) /110/, mass(247) /109.90300660600_dp/, mass_unc(247) /0.00000060600_dp/  ! 110Cd
data Z(248) / 48/, A(248) /113/, mass(248) /112.90440813300_dp/, mass_unc(248) /0.00000044500_dp/  ! 113Cd
data Z(249) / 48/, A(249) /116/, mass(249) /115.90476314800_dp/, mass_unc(249) /0.00000017100_dp/  ! 116Cd
data Z(250) / 48/, A(250) /106/, mass(250) /105.90645992800_dp/, mass_unc(250) /0.00000118500_dp/  ! 106Cd
data Z(251) / 48/, A(251) /108/, mass(251) /107.90418344000_dp/, mass_unc(251) /0.00000120600_dp/  ! 108Cd
data Z(252) / 48/, A(252) /113/, mass(252) /112.90440813300_dp/, mass_unc(252) /0.00100044500_dp/  ! 113mCd
data Z(253) / 48/, A(253) /109/, mass(253) /108.90498665300_dp/, mass_unc(253) /0.00000165000_dp/  ! 109Cd
data Z(254) / 48/, A(254) /115/, mass(254) /114.90543751300_dp/, mass_unc(254) /0.00000076500_dp/  ! 115mCd
data Z(255) / 48/, A(255) /115/, mass(255) /114.90543751300_dp/, mass_unc(255) /0.00100076500_dp/  ! 115Cd
data Z(256) / 49/, A(256) /115/, mass(256) /114.90387877600_dp/, mass_unc(256) /0.00000001200_dp/  ! 115In
data Z(257) / 49/, A(257) /113/, mass(257) /112.90406183900_dp/, mass_unc(257) /0.00000091400_dp/  ! 113In
data Z(258) / 49/, A(258) /114/, mass(258) /113.90491790900_dp/, mass_unc(258) /0.00000093800_dp/  ! 114mIn
data Z(259) / 49/, A(259) /111/, mass(259) /110.90510845800_dp/, mass_unc(259) /0.00000376800_dp/  ! 111In
data Z(260) / 50/, A(260) /120/, mass(260) /119.90220163400_dp/, mass_unc(260) /0.00000097000_dp/  ! 120Sn
data Z(261) / 50/, A(261) /118/, mass(261) /117.90160657400_dp/, mass_unc(261) /0.00000053600_dp/  ! 118Sn
data Z(262) / 50/, A(262) /116/, mass(262) /115.90174279700_dp/, mass_unc(262) /0.00000010200_dp/  ! 116Sn
data Z(263) / 50/, A(263) /119/, mass(263) /118.90331117200_dp/, mass_unc(263) /0.00000077900_dp/  ! 119Sn
data Z(264) / 50/, A(264) /117/, mass(264) /116.90295398300_dp/, mass_unc(264) /0.00000051900_dp/  ! 117Sn
data Z(265) / 50/, A(265) /124/, mass(265) /123.90527664500_dp/, mass_unc(265) /0.00000108900_dp/  ! 124Sn
data Z(266) / 50/, A(266) /122/, mass(266) /121.90344377400_dp/, mass_unc(266) /0.00000257600_dp/  ! 122Sn
data Z(267) / 50/, A(267) /112/, mass(267) /111.90482387400_dp/, mass_unc(267) /0.00000061400_dp/  ! 112Sn
data Z(268) / 50/, A(268) /114/, mass(268) /113.90278269500_dp/, mass_unc(268) /0.00000104200_dp/  ! 114Sn
data Z(269) / 50/, A(269) /115/, mass(269) /114.90334469900_dp/, mass_unc(269) /0.00000001600_dp/  ! 115Sn
data Z(270) / 50/, A(270) /126/, mass(270) /125.90765878600_dp/, mass_unc(270) /0.00001121500_dp/  ! 126Sn
data Z(271) / 50/, A(271) /121/, mass(271) /120.90424255400_dp/, mass_unc(271) /0.00000103300_dp/  ! 121mSn
data Z(272) / 50/, A(272) /119/, mass(272) /118.90331117200_dp/, mass_unc(272) /0.00100077900_dp/  ! 119mSn
data Z(273) / 50/, A(273) /123/, mass(273) /122.90572522100_dp/, mass_unc(273) /0.00000260100_dp/  ! 123Sn
data Z(274) / 50/, A(274) /113/, mass(274) /112.90517572800_dp/, mass_unc(274) /0.00000176200_dp/  ! 113Sn
data Z(275) / 50/, A(275) /117/, mass(275) /116.90295398300_dp/, mass_unc(275) /0.00100051900_dp/  ! 117mSn
data Z(276) / 50/, A(276) /125/, mass(276) /124.90778639500_dp/, mass_unc(276) /0.00000110900_dp/  ! 125Sn
data Z(277) / 50/, A(277) /121/, mass(277) /120.90424255400_dp/, mass_unc(277) /0.00100103300_dp/  ! 121Sn
data Z(278) / 51/, A(278) /121/, mass(278) /120.90381196700_dp/, mass_unc(278) /0.00000300900_dp/  ! 121Sb
data Z(279) / 51/, A(279) /123/, mass(279) /122.90421320400_dp/, mass_unc(279) /0.00000227800_dp/  ! 123Sb
data Z(280) / 51/, A(280) /125/, mass(280) /124.90525300700_dp/, mass_unc(280) /0.00000279100_dp/  ! 125Sb
data Z(281) / 51/, A(281) /124/, mass(281) /123.90593497800_dp/, mass_unc(281) /0.00000227700_dp/  ! 124Sb
data Z(282) / 51/, A(282) /126/, mass(282) /125.90725298700_dp/, mass_unc(282) /0.00003410300_dp/  ! 126Sb
data Z(283) / 51/, A(283) /120/, mass(283) /119.90507938500_dp/, mass_unc(283) /0.00000772600_dp/  ! 120mSb
data Z(284) / 51/, A(284) /127/, mass(284) /126.90692427800_dp/, mass_unc(284) /0.00000550300_dp/  ! 127Sb
data Z(285) / 51/, A(285) /122/, mass(285) /121.90516994800_dp/, mass_unc(285) /0.00000300600_dp/  ! 122Sb
data Z(286) / 51/, A(286) /119/, mass(286) /118.90394547100_dp/, mass_unc(286) /0.00000826700_dp/  ! 119Sb
data Z(287) / 52/, A(287) /130/, mass(287) /129.90622274800_dp/, mass_unc(287) /0.00000001200_dp/  ! 130Te
data Z(288) / 52/, A(288) /126/, mass(288) /125.90331088600_dp/, mass_unc(288) /0.00000161600_dp/  ! 126Te
data Z(289) / 52/, A(289) /125/, mass(289) /124.90442992000_dp/, mass_unc(289) /0.00000161400_dp/  ! 125Te
data Z(290) / 52/, A(290) /124/, mass(290) /123.90281708500_dp/, mass_unc(290) /0.00000161400_dp/  ! 124Te
data Z(291) / 52/, A(291) /122/, mass(291) /121.90304345500_dp/, mass_unc(291) /0.00000161900_dp/  ! 122Te
data Z(292) / 52/, A(292) /123/, mass(292) /122.90426976900_dp/, mass_unc(292) /0.00000161700_dp/  ! 123Te
data Z(293) / 52/, A(293) /120/, mass(293) /119.90405930000_dp/, mass_unc(293) /0.00000334100_dp/  ! 120Te
data Z(294) / 52/, A(294) /128/, mass(294) /127.90446127900_dp/, mass_unc(294) /0.00000092900_dp/  ! 128Te
data Z(295) / 52/, A(295) /121/, mass(295) /120.90494381200_dp/, mass_unc(295) /0.00002776400_dp/  ! 121mTe
data Z(296) / 52/, A(296) /123/, mass(296) /122.90426976900_dp/, mass_unc(296) /0.00100161700_dp/  ! 123mTe
data Z(297) / 52/, A(297) /127/, mass(297) /126.90522573500_dp/, mass_unc(297) /0.00000162700_dp/  ! 127mTe
data Z(298) / 52/, A(298) /125/, mass(298) /124.90442992000_dp/, mass_unc(298) /0.00100161400_dp/  ! 125mTe
data Z(299) / 52/, A(299) /129/, mass(299) /128.90659646000_dp/, mass_unc(299) /0.00000093300_dp/  ! 129mTe
data Z(300) / 52/, A(300) /121/, mass(300) /120.90494381200_dp/, mass_unc(300) /0.00102776400_dp/  ! 121Te
data Z(301) / 52/, A(301) /118/, mass(301) /117.90585362900_dp/, mass_unc(301) /0.00001984500_dp/  ! 118Te
data Z(302) / 52/, A(302) /119/, mass(302) /118.90640710800_dp/, mass_unc(302) /0.00000854100_dp/  ! 119mTe
data Z(303) / 52/, A(303) /132/, mass(303) /131.90854671600_dp/, mass_unc(303) /0.00000374200_dp/  ! 132Te
data Z(304) / 52/, A(304) /131/, mass(304) /130.90852221300_dp/, mass_unc(304) /0.00000006500_dp/  ! 131mTe
data Z(305) / 53/, A(305) /127/, mass(305) /126.90447185300_dp/, mass_unc(305) /0.00000391500_dp/  ! 127I
data Z(306) / 53/, A(306) /129/, mass(306) /128.90498366900_dp/, mass_unc(306) /0.00000340100_dp/  ! 129I
data Z(307) / 53/, A(307) /125/, mass(307) /124.90462935300_dp/, mass_unc(307) /0.00000161500_dp/  ! 125I
data Z(308) / 53/, A(308) /126/, mass(308) /125.90562332900_dp/, mass_unc(308) /0.00000409000_dp/  ! 126I
data Z(309) / 53/, A(309) /131/, mass(309) /130.90612630500_dp/, mass_unc(309) /0.00000069000_dp/  ! 131I
data Z(310) / 53/, A(310) /124/, mass(310) /123.90620904100_dp/, mass_unc(310) /0.00000256700_dp/  ! 124I
data Z(311) / 54/, A(311) /132/, mass(311) /131.90415508563_dp/, mass_unc(311) /0.00000000559_dp/  ! 132Xe
data Z(312) / 54/, A(312) /129/, mass(312) /128.90478086113_dp/, mass_unc(312) /0.00000000596_dp/  ! 129Xe
data Z(313) / 54/, A(313) /131/, mass(313) /130.90508405700_dp/, mass_unc(313) /0.00000023600_dp/  ! 131Xe
data Z(314) / 54/, A(314) /134/, mass(314) /133.90539466400_dp/, mass_unc(314) /0.00000090400_dp/  ! 134Xe
data Z(315) / 54/, A(315) /136/, mass(315) /135.90721448400_dp/, mass_unc(315) /0.00000001100_dp/  ! 136Xe
data Z(316) / 54/, A(316) /130/, mass(316) /129.90350934900_dp/, mass_unc(316) /0.00000001000_dp/  ! 130Xe
data Z(317) / 54/, A(317) /128/, mass(317) /127.90353101800_dp/, mass_unc(317) /0.00000113800_dp/  ! 128Xe
data Z(318) / 54/, A(318) /124/, mass(318) /123.90589198400_dp/, mass_unc(318) /0.00000193600_dp/  ! 124Xe
data Z(319) / 54/, A(319) /126/, mass(319) /125.90429829200_dp/, mass_unc(319) /0.00000384700_dp/  ! 126Xe
data Z(320) / 54/, A(320) /127/, mass(320) /126.90518291400_dp/, mass_unc(320) /0.00000441200_dp/  ! 127Xe
data Z(321) / 54/, A(321) /131/, mass(321) /130.90508405700_dp/, mass_unc(321) /0.00100023600_dp/  ! 131mXe
data Z(322) / 54/, A(322) /129/, mass(322) /128.90478086113_dp/, mass_unc(322) /0.00100000596_dp/  ! 129mXe
data Z(323) / 54/, A(323) /133/, mass(323) /132.90591075100_dp/, mass_unc(323) /0.00000257600_dp/  ! 133Xe
data Z(324) / 54/, A(324) /133/, mass(324) /132.90591075100_dp/, mass_unc(324) /0.00100257600_dp/  ! 133mXe
data Z(325) / 55/, A(325) /133/, mass(325) /132.90545196100_dp/, mass_unc(325) /0.00000000800_dp/  ! 133Cs
data Z(326) / 55/, A(326) /135/, mass(326) /134.90597704900_dp/, mass_unc(326) /0.00000106900_dp/  ! 135Cs
data Z(327) / 55/, A(327) /137/, mass(327) /136.90708923100_dp/, mass_unc(327) /0.00000035500_dp/  ! 137Cs
data Z(328) / 55/, A(328) /134/, mass(328) /133.90671850300_dp/, mass_unc(328) /0.00000001700_dp/  ! 134Cs
data Z(329) / 55/, A(329) /136/, mass(329) /135.90731135800_dp/, mass_unc(329) /0.00000201500_dp/  ! 136Cs
data Z(330) / 55/, A(330) /131/, mass(330) /130.90546489900_dp/, mass_unc(330) /0.00000534300_dp/  ! 131Cs
data Z(331) / 55/, A(331) /132/, mass(331) /131.90643391400_dp/, mass_unc(331) /0.00000214700_dp/  ! 132Cs
data Z(332) / 55/, A(332) /129/, mass(332) /128.90606568300_dp/, mass_unc(332) /0.00000488900_dp/  ! 129Cs
data Z(333) / 56/, A(333) /138/, mass(333) /137.90524699500_dp/, mass_unc(333) /0.00000030600_dp/  ! 138Ba
data Z(334) / 56/, A(334) /137/, mass(334) /136.90582714100_dp/, mass_unc(334) /0.00000030300_dp/  ! 137Ba
data Z(335) / 56/, A(335) /136/, mass(335) /135.90457572700_dp/, mass_unc(335) /0.00000029400_dp/  ! 136Ba
data Z(336) / 56/, A(336) /135/, mass(336) /134.90568837500_dp/, mass_unc(336) /0.00000029400_dp/  ! 135Ba
data Z(337) / 56/, A(337) /134/, mass(337) /133.90450818200_dp/, mass_unc(337) /0.00000029600_dp/  ! 134Ba
data Z(338) / 56/, A(338) /130/, mass(338) /129.90632066900_dp/, mass_unc(338) /0.00000276700_dp/  ! 130Ba
data Z(339) / 56/, A(339) /132/, mass(339) /131.90506112800_dp/, mass_unc(339) /0.00000113200_dp/  ! 132Ba
data Z(340) / 56/, A(340) /133/, mass(340) /132.90600735100_dp/, mass_unc(340) /0.00000106600_dp/  ! 133Ba
data Z(341) / 56/, A(341) /140/, mass(341) /139.91060573000_dp/, mass_unc(341) /0.00000852700_dp/  ! 140Ba
data Z(342) / 56/, A(342) /131/, mass(342) /130.90694097700_dp/, mass_unc(342) /0.00000278400_dp/  ! 131Ba
data Z(343) / 56/, A(343) /128/, mass(343) /127.90834196700_dp/, mass_unc(343) /0.00000563700_dp/  ! 128Ba
data Z(344) / 56/, A(344) /133/, mass(344) /132.90600735100_dp/, mass_unc(344) /0.00100106600_dp/  ! 133mBa
data Z(345) / 56/, A(345) /135/, mass(345) /134.90568837500_dp/, mass_unc(345) /0.00100029400_dp/  ! 135mBa
data Z(346) / 57/, A(346) /139/, mass(346) /138.90635625600_dp/, mass_unc(346) /0.00000243800_dp/  ! 139La
data Z(347) / 57/, A(347) /138/, mass(347) /137.90711491900_dp/, mass_unc(347) /0.00000368400_dp/  ! 138La
data Z(348) / 57/, A(348) /137/, mass(348) /136.90645038500_dp/, mass_unc(348) /0.00000177400_dp/  ! 137La
data Z(349) / 57/, A(349) /140/, mass(349) /139.90948063500_dp/, mass_unc(349) /0.00000243800_dp/  ! 140La
data Z(350) / 58/, A(350) /140/, mass(350) /139.90544310700_dp/, mass_unc(350) /0.00000234000_dp/  ! 140Ce
data Z(351) / 58/, A(351) /142/, mass(351) /141.90925037500_dp/, mass_unc(351) /0.00000291000_dp/  ! 142Ce
data Z(352) / 58/, A(352) /138/, mass(352) /137.90599108900_dp/, mass_unc(352) /0.00001061200_dp/  ! 138Ce
data Z(353) / 58/, A(353) /136/, mass(353) /135.90712920500_dp/, mass_unc(353) /0.00000041200_dp/  ! 136Ce
data Z(354) / 58/, A(354) /144/, mass(354) /143.91365293900_dp/, mass_unc(354) /0.00000335400_dp/  ! 144Ce
data Z(355) / 58/, A(355) /139/, mass(355) /138.90665511100_dp/, mass_unc(355) /0.00000784200_dp/  ! 139Ce
data Z(356) / 58/, A(356) /141/, mass(356) /140.90828067400_dp/, mass_unc(356) /0.00000234000_dp/  ! 141Ce
data Z(357) / 58/, A(357) /134/, mass(357) /133.90892814200_dp/, mass_unc(357) /0.00002188600_dp/  ! 134Ce
data Z(358) / 58/, A(358) /137/, mass(358) /136.90776236400_dp/, mass_unc(358) /0.00000044500_dp/  ! 137mCe
data Z(359) / 58/, A(359) /143/, mass(359) /142.91239212000_dp/, mass_unc(359) /0.00000290900_dp/  ! 143Ce
data Z(360) / 59/, A(360) /141/, mass(360) /140.90765756800_dp/, mass_unc(360) /0.00000229200_dp/  ! 141Pr
data Z(361) / 59/, A(361) /143/, mass(361) /142.91082279600_dp/, mass_unc(361) /0.00000239600_dp/  ! 143Pr
data Z(362) / 60/, A(362) /142/, mass(362) /141.90772899600_dp/, mass_unc(362) /0.00000195300_dp/  ! 142Nd
data Z(363) / 60/, A(363) /144/, mass(363) /143.91009297400_dp/, mass_unc(363) /0.00000195200_dp/  ! 144Nd
data Z(364) / 60/, A(364) /146/, mass(364) /145.91312262800_dp/, mass_unc(364) /0.00000195300_dp/  ! 146Nd
data Z(365) / 60/, A(365) /143/, mass(365) /142.90981998900_dp/, mass_unc(365) /0.00000195300_dp/  ! 143Nd
data Z(366) / 60/, A(366) /145/, mass(366) /144.91257932200_dp/, mass_unc(366) /0.00000195300_dp/  ! 145Nd
data Z(367) / 60/, A(367) /148/, mass(367) /147.91689929400_dp/, mass_unc(367) /0.00000259400_dp/  ! 148Nd
data Z(368) / 60/, A(368) /150/, mass(368) /149.92090224900_dp/, mass_unc(368) /0.00000182100_dp/  ! 150Nd
data Z(369) / 60/, A(369) /147/, mass(369) /146.91610613600_dp/, mass_unc(369) /0.00000195300_dp/  ! 147Nd
data Z(370) / 60/, A(370) /140/, mass(370) /139.90954984900_dp/, mass_unc(370) /0.00002792700_dp/  ! 140Nd
data Z(371) / 61/, A(371) /145/, mass(371) /144.91275593500_dp/, mass_unc(371) /0.00000330600_dp/  ! 145Pm
data Z(372) / 61/, A(372) /146/, mass(372) /145.91470239600_dp/, mass_unc(372) /0.00000477100_dp/  ! 146Pm
data Z(373) / 61/, A(373) /147/, mass(373) /146.91514498800_dp/, mass_unc(373) /0.00000192600_dp/  ! 147Pm
data Z(374) / 61/, A(374) /144/, mass(374) /143.91259639600_dp/, mass_unc(374) /0.00000340700_dp/  ! 144Pm
data Z(375) / 61/, A(375) /143/, mass(375) /142.91093826200_dp/, mass_unc(375) /0.00000343700_dp/  ! 143Pm
data Z(376) / 61/, A(376) /148/, mass(376) /147.91748194500_dp/, mass_unc(376) /0.00000625700_dp/  ! 148mPm
data Z(377) / 61/, A(377) /148/, mass(377) /147.91748194500_dp/, mass_unc(377) /0.00100625700_dp/  ! 148Pm
data Z(378) / 61/, A(378) /149/, mass(378) /148.91834227700_dp/, mass_unc(378) /0.00000269800_dp/  ! 149Pm
data Z(379) / 61/, A(379) /151/, mass(379) /150.92121753900_dp/, mass_unc(379) /0.00000513200_dp/  ! 151Pm
data Z(380) / 62/, A(380) /152/, mass(380) /151.91973972100_dp/, mass_unc(380) /0.00000175900_dp/  ! 152Sm
data Z(381) / 62/, A(381) /154/, mass(381) /153.92221686100_dp/, mass_unc(381) /0.00000196400_dp/  ! 154Sm
data Z(382) / 62/, A(382) /147/, mass(382) /146.91490443500_dp/, mass_unc(382) /0.00000191900_dp/  ! 147Sm
data Z(383) / 62/, A(383) /149/, mass(383) /148.91719206200_dp/, mass_unc(383) /0.00000183900_dp/  ! 149Sm
data Z(384) / 62/, A(384) /148/, mass(384) /147.91482922600_dp/, mass_unc(384) /0.00000192000_dp/  ! 148Sm
data Z(385) / 62/, A(385) /150/, mass(385) /149.91728291900_dp/, mass_unc(385) /0.00000180900_dp/  ! 150Sm
data Z(386) / 62/, A(386) /144/, mass(386) /143.91200646600_dp/, mass_unc(386) /0.00000209000_dp/  ! 144Sm
data Z(387) / 62/, A(387) /146/, mass(387) /145.91304699100_dp/, mass_unc(387) /0.00000352000_dp/  ! 146Sm
data Z(388) / 62/, A(388) /151/, mass(388) /150.91993979600_dp/, mass_unc(388) /0.00000180700_dp/  ! 151Sm
data Z(389) / 62/, A(389) /145/, mass(389) /144.91341733900_dp/, mass_unc(389) /0.00000211000_dp/  ! 145Sm
data Z(390) / 62/, A(390) /153/, mass(390) /152.92210465000_dp/, mass_unc(390) /0.00000176500_dp/  ! 153Sm
data Z(391) / 63/, A(391) /153/, mass(391) /152.92123800300_dp/, mass_unc(391) /0.00000184600_dp/  ! 153Eu
data Z(392) / 63/, A(392) /151/, mass(392) /150.91985780300_dp/, mass_unc(392) /0.00000184300_dp/  ! 151Eu
data Z(393) / 63/, A(393) /150/, mass(393) /149.91970767100_dp/, mass_unc(393) /0.00000678700_dp/  ! 150Eu
data Z(394) / 63/, A(394) /152/, mass(394) /151.92175218400_dp/, mass_unc(394) /0.00000184300_dp/  ! 152Eu
data Z(395) / 63/, A(395) /154/, mass(395) /153.92298696200_dp/, mass_unc(395) /0.00000185700_dp/  ! 154Eu
data Z(396) / 63/, A(396) /155/, mass(396) /154.92290110700_dp/, mass_unc(396) /0.00000190200_dp/  ! 155Eu
data Z(397) / 63/, A(397) /149/, mass(397) /148.91793776300_dp/, mass_unc(397) /0.00000439300_dp/  ! 149Eu
data Z(398) / 63/, A(398) /148/, mass(398) /147.91808924300_dp/, mass_unc(398) /0.00001076900_dp/  ! 148Eu
data Z(399) / 63/, A(399) /147/, mass(399) /146.91675265900_dp/, mass_unc(399) /0.00000305100_dp/  ! 147Eu
data Z(400) / 63/, A(400) /156/, mass(400) /155.92476049400_dp/, mass_unc(400) /0.00000592700_dp/  ! 156Eu
data Z(401) / 63/, A(401) /145/, mass(401) /144.91627262900_dp/, mass_unc(401) /0.00000356000_dp/  ! 145Eu
data Z(402) / 63/, A(402) /146/, mass(402) /145.91721103900_dp/, mass_unc(402) /0.00000654100_dp/  ! 146Eu
data Z(403) / 64/, A(403) /158/, mass(403) /157.92411234800_dp/, mass_unc(403) /0.00000174200_dp/  ! 158Gd
data Z(404) / 64/, A(404) /160/, mass(404) /159.92706241100_dp/, mass_unc(404) /0.00000184000_dp/  ! 160Gd
data Z(405) / 64/, A(405) /156/, mass(405) /155.92213124100_dp/, mass_unc(405) /0.00000174200_dp/  ! 156Gd
data Z(406) / 64/, A(406) /157/, mass(406) /156.92396856900_dp/, mass_unc(406) /0.00000174200_dp/  ! 157Gd
data Z(407) / 64/, A(407) /155/, mass(407) /154.92263047300_dp/, mass_unc(407) /0.00000174200_dp/  ! 155Gd
data Z(408) / 64/, A(408) /154/, mass(408) /153.92087406000_dp/, mass_unc(408) /0.00000174700_dp/  ! 154Gd
data Z(409) / 64/, A(409) /152/, mass(409) /151.91979949400_dp/, mass_unc(409) /0.00000175400_dp/  ! 152Gd
data Z(410) / 64/, A(410) /150/, mass(410) /149.91866442200_dp/, mass_unc(410) /0.00000659700_dp/  ! 150Gd
data Z(411) / 64/, A(411) /148/, mass(411) /147.91812151100_dp/, mass_unc(411) /0.00000209000_dp/  ! 148Gd
data Z(412) / 64/, A(412) /153/, mass(412) /152.92175802700_dp/, mass_unc(412) /0.00000175100_dp/  ! 153Gd
data Z(413) / 64/, A(413) /151/, mass(413) /150.92035595000_dp/, mass_unc(413) /0.00000348100_dp/  ! 151Gd
data Z(414) / 64/, A(414) /146/, mass(414) /145.91831881700_dp/, mass_unc(414) /0.00000456800_dp/  ! 146Gd
data Z(415) / 64/, A(415) /149/, mass(415) /148.91934811700_dp/, mass_unc(415) /0.00000379100_dp/  ! 149Gd
data Z(416) / 64/, A(416) /147/, mass(416) /146.91910138400_dp/, mass_unc(416) /0.00000245600_dp/  ! 147Gd
data Z(417) / 65/, A(417) /159/, mass(417) /158.92535471000_dp/, mass_unc(417) /0.00000187900_dp/  ! 159Tb
data Z(418) / 65/, A(418) /158/, mass(418) /157.92542094700_dp/, mass_unc(418) /0.00000199700_dp/  ! 158Tb
data Z(419) / 65/, A(419) /157/, mass(419) /156.92403302800_dp/, mass_unc(419) /0.00000177000_dp/  ! 157Tb
data Z(420) / 65/, A(420) /160/, mass(420) /159.92717555600_dp/, mass_unc(420) /0.00000188400_dp/  ! 160Tb
data Z(421) / 65/, A(421) /161/, mass(421) /160.92757782500_dp/, mass_unc(421) /0.00000196200_dp/  ! 161Tb
data Z(422) / 65/, A(422) /156/, mass(422) /155.92475518100_dp/, mass_unc(422) /0.00000429800_dp/  ! 156Tb
data Z(423) / 65/, A(423) /155/, mass(423) /154.92351054700_dp/, mass_unc(423) /0.00001062900_dp/  ! 155Tb
data Z(424) / 65/, A(424) /153/, mass(424) /152.92344240300_dp/, mass_unc(424) /0.00000444000_dp/  ! 153Tb
data Z(425) / 65/, A(425) /156/, mass(425) /155.92475518100_dp/, mass_unc(425) /0.00100429800_dp/  ! 156mTb
data Z(426) / 66/, A(426) /164/, mass(426) /163.92918187400_dp/, mass_unc(426) /0.00000202100_dp/  ! 164Dy
data Z(427) / 66/, A(427) /162/, mass(427) /161.92680557300_dp/, mass_unc(427) /0.00000202100_dp/  ! 162Dy
data Z(428) / 66/, A(428) /163/, mass(428) /162.92873828400_dp/, mass_unc(428) /0.00000202100_dp/  ! 163Dy
data Z(429) / 66/, A(429) /161/, mass(429) /160.92694049200_dp/, mass_unc(429) /0.00000202000_dp/  ! 161Dy
data Z(430) / 66/, A(430) /160/, mass(430) /159.92520464600_dp/, mass_unc(430) /0.00000202100_dp/  ! 160Dy
data Z(431) / 66/, A(431) /158/, mass(431) /157.92441587500_dp/, mass_unc(431) /0.00000307500_dp/  ! 158Dy
data Z(432) / 66/, A(432) /156/, mass(432) /155.92428471300_dp/, mass_unc(432) /0.00000174600_dp/  ! 156Dy
data Z(433) / 66/, A(433) /154/, mass(433) /153.92442927700_dp/, mass_unc(433) /0.00000804300_dp/  ! 154Dy
data Z(434) / 66/, A(434) /159/, mass(434) /158.92574695800_dp/, mass_unc(434) /0.00000218700_dp/  ! 159Dy
data Z(435) / 66/, A(435) /166/, mass(435) /165.93281386300_dp/, mass_unc(435) /0.00000206700_dp/  ! 166Dy
data Z(436) / 67/, A(436) /165/, mass(436) /164.93032883500_dp/, mass_unc(436) /0.00000212900_dp/  ! 165Ho
data Z(437) / 67/, A(437) /163/, mass(437) /162.92874102700_dp/, mass_unc(437) /0.00000202100_dp/  ! 163Ho
data Z(438) / 67/, A(438) /166/, mass(438) /165.93229092700_dp/, mass_unc(438) /0.00000212900_dp/  ! 166mHo
data Z(439) / 67/, A(439) /166/, mass(439) /165.93229092700_dp/, mass_unc(439) /0.00100212900_dp/  ! 166Ho
data Z(440) / 68/, A(440) /166/, mass(440) /165.93029953000_dp/, mass_unc(440) /0.00000217900_dp/  ! 166Er
data Z(441) / 68/, A(441) /168/, mass(441) /167.93237668800_dp/, mass_unc(441) /0.00000218200_dp/  ! 168Er
data Z(442) / 68/, A(442) /167/, mass(442) /166.93205461700_dp/, mass_unc(442) /0.00000218000_dp/  ! 167Er
data Z(443) / 68/, A(443) /170/, mass(443) /169.93547023000_dp/, mass_unc(443) /0.00000255600_dp/  ! 170Er
data Z(444) / 68/, A(444) /164/, mass(444) /163.92920879100_dp/, mass_unc(444) /0.00000202400_dp/  ! 164Er
data Z(445) / 68/, A(445) /162/, mass(445) /161.92878836400_dp/, mass_unc(445) /0.00000204600_dp/  ! 162Er
data Z(446) / 68/, A(446) /169/, mass(446) /168.93459685000_dp/, mass_unc(446) /0.00000218800_dp/  ! 169Er
data Z(447) / 68/, A(447) /172/, mass(447) /171.93936185800_dp/, mass_unc(447) /0.00000473100_dp/  ! 172Er
data Z(448) / 68/, A(448) /160/, mass(448) /159.92907713000_dp/, mass_unc(448) /0.00002603700_dp/  ! 160Er
data Z(449) / 69/, A(449) /169/, mass(449) /168.93421788900_dp/, mass_unc(449) /0.00000221700_dp/  ! 169Tm
data Z(450) / 69/, A(450) /171/, mass(450) /170.93643387100_dp/, mass_unc(450) /0.00000239600_dp/  ! 171Tm
data Z(451) / 69/, A(451) /170/, mass(451) /169.93580603200_dp/, mass_unc(451) /0.00000221800_dp/  ! 170Tm
data Z(452) / 69/, A(452) /168/, mass(452) /167.93417740900_dp/, mass_unc(452) /0.00000274100_dp/  ! 168Tm
data Z(453) / 69/, A(453) /167/, mass(453) /166.93285619200_dp/, mass_unc(453) /0.00000246400_dp/  ! 167Tm
data Z(454) / 69/, A(454) /172/, mass(454) /171.93840552100_dp/, mass_unc(454) /0.00000623900_dp/  ! 172Tm
data Z(455) / 69/, A(455) /165/, mass(455) /164.93244314000_dp/, mass_unc(455) /0.00000255200_dp/  ! 165Tm
data Z(456) / 70/, A(456) /174/, mass(456) /173.93886643700_dp/, mass_unc(456) /0.00000216000_dp/  ! 174Yb
data Z(457) / 70/, A(457) /172/, mass(457) /171.93638587200_dp/, mass_unc(457) /0.00000216600_dp/  ! 172Yb
data Z(458) / 70/, A(458) /173/, mass(458) /172.93821513600_dp/, mass_unc(458) /0.00000216000_dp/  ! 173Yb
data Z(459) / 70/, A(459) /171/, mass(459) /170.93633020800_dp/, mass_unc(459) /0.00000216700_dp/  ! 171Yb
data Z(460) / 70/, A(460) /176/, mass(460) /175.94257644700_dp/, mass_unc(460) /0.00000241400_dp/  ! 176Yb
data Z(461) / 70/, A(461) /170/, mass(461) /169.93476637600_dp/, mass_unc(461) /0.00000221000_dp/  ! 170Yb
data Z(462) / 70/, A(462) /168/, mass(462) /167.93388960200_dp/, mass_unc(462) /0.00000219800_dp/  ! 168Yb
data Z(463) / 70/, A(463) /169/, mass(463) /168.93518251200_dp/, mass_unc(463) /0.00000220400_dp/  ! 169Yb
data Z(464) / 70/, A(464) /175/, mass(464) /174.94128079700_dp/, mass_unc(464) /0.00000216100_dp/  ! 175Yb
data Z(465) / 70/, A(465) /166/, mass(465) /165.93387474900_dp/, mass_unc(465) /0.00000782900_dp/  ! 166Yb
data Z(466) / 71/, A(466) /175/, mass(466) /174.94077519100_dp/, mass_unc(466) /0.00000204100_dp/  ! 175Lu
data Z(467) / 71/, A(467) /176/, mass(467) /175.94268968000_dp/, mass_unc(467) /0.00000203900_dp/  ! 176Lu
data Z(468) / 71/, A(468) /174/, mass(468) /173.94034085400_dp/, mass_unc(468) /0.00000230700_dp/  ! 174Lu
data Z(469) / 71/, A(469) /173/, mass(469) /172.93893402900_dp/, mass_unc(469) /0.00000231900_dp/  ! 173Lu
data Z(470) / 71/, A(470) /177/, mass(470) /176.94376152500_dp/, mass_unc(470) /0.00000204000_dp/  ! 177m3Lu
data Z(471) / 71/, A(471) /174/, mass(471) /173.94034085400_dp/, mass_unc(471) /0.00100230700_dp/  ! 174mLu
data Z(472) / 71/, A(472) /171/, mass(472) /170.93791696000_dp/, mass_unc(472) /0.00000268800_dp/  ! 171Lu
data Z(473) / 71/, A(473) /172/, mass(473) /171.93908910300_dp/, mass_unc(473) /0.00000297400_dp/  ! 172Lu
data Z(474) / 71/, A(474) /177/, mass(474) /176.94376152500_dp/, mass_unc(474) /0.00100204000_dp/  ! 177Lu
data Z(475) / 71/, A(475) /170/, mass(475) /169.93847836500_dp/, mass_unc(475) /0.00001821600_dp/  ! 170Lu
data Z(476) / 71/, A(476) /169/, mass(476) /168.93764414900_dp/, mass_unc(476) /0.00000390300_dp/  ! 169Lu
data Z(477) / 72/, A(477) /180/, mass(477) /179.94655704200_dp/, mass_unc(477) /0.00000203200_dp/  ! 180Hf
data Z(478) / 72/, A(478) /178/, mass(478) /177.94370583300_dp/, mass_unc(478) /0.00000202900_dp/  ! 178Hf
data Z(479) / 72/, A(479) /177/, mass(479) /176.94322771700_dp/, mass_unc(479) /0.00000203100_dp/  ! 177Hf
data Z(480) / 72/, A(480) /179/, mass(480) /178.94582321200_dp/, mass_unc(480) /0.00000202800_dp/  ! 179Hf
data Z(481) / 72/, A(481) /176/, mass(481) /175.94140762800_dp/, mass_unc(481) /0.00000219000_dp/  ! 176Hf
data Z(482) / 72/, A(482) /174/, mass(482) /173.94004614100_dp/, mass_unc(482) /0.00000284500_dp/  ! 174Hf
data Z(483) / 72/, A(483) /182/, mass(483) /181.95056118500_dp/, mass_unc(483) /0.00000675400_dp/  ! 182Hf
data Z(484) / 72/, A(484) /178/, mass(484) /177.94370583300_dp/, mass_unc(484) /0.00100202900_dp/  ! 178m2Hf
data Z(485) / 72/, A(485) /172/, mass(485) /171.93944971600_dp/, mass_unc(485) /0.00002622400_dp/  ! 172Hf
data Z(486) / 72/, A(486) /175/, mass(486) /174.94150918700_dp/, mass_unc(486) /0.00000286600_dp/  ! 175Hf
data Z(487) / 72/, A(487) /181/, mass(487) /180.94910833800_dp/, mass_unc(487) /0.00000203300_dp/  ! 181Hf
data Z(488) / 72/, A(488) /179/, mass(488) /178.94582321200_dp/, mass_unc(488) /0.00100202800_dp/  ! 179m2Hf
data Z(489) / 73/, A(489) /181/, mass(489) /180.94799576900_dp/, mass_unc(489) /0.00000197900_dp/  ! 181Ta
data Z(490) / 73/, A(490) /180/, mass(490) /179.94746483200_dp/, mass_unc(490) /0.00000244500_dp/  ! 180mTa
data Z(491) / 73/, A(491) /179/, mass(491) /178.94593655500_dp/, mass_unc(491) /0.00000206600_dp/  ! 179Ta
data Z(492) / 73/, A(492) /182/, mass(492) /181.95015185300_dp/, mass_unc(492) /0.00000198000_dp/  ! 182Ta
data Z(493) / 73/, A(493) /183/, mass(493) /182.95137262000_dp/, mass_unc(493) /0.00000199100_dp/  ! 183Ta
data Z(494) / 73/, A(494) /177/, mass(494) /176.94447946900_dp/, mass_unc(494) /0.00000380700_dp/  ! 177Ta
data Z(495) / 74/, A(495) /184/, mass(495) /183.95093091600_dp/, mass_unc(495) /0.00000094200_dp/  ! 184W
data Z(496) / 74/, A(496) /186/, mass(496) /185.95436277100_dp/, mass_unc(496) /0.00000165500_dp/  ! 186W
data Z(497) / 74/, A(497) /182/, mass(497) /181.94820394500_dp/, mass_unc(497) /0.00000090500_dp/  ! 182W
data Z(498) / 74/, A(498) /183/, mass(498) /182.95022274800_dp/, mass_unc(498) /0.00000090300_dp/  ! 183W
data Z(499) / 74/, A(499) /180/, mass(499) /179.94671080500_dp/, mass_unc(499) /0.00000204300_dp/  ! 180W
data Z(500) / 74/, A(500) /181/, mass(500) /180.94819778300_dp/, mass_unc(500) /0.00000508100_dp/  ! 181W
data Z(501) / 74/, A(501) /185/, mass(501) /184.95341897400_dp/, mass_unc(501) /0.00000099300_dp/  ! 185W
data Z(502) / 74/, A(502) /188/, mass(502) /187.95848617700_dp/, mass_unc(502) /0.00000360000_dp/  ! 188W
data Z(503) / 74/, A(503) /178/, mass(503) /177.94588330300_dp/, mass_unc(503) /0.00001637100_dp/  ! 178W
data Z(504) / 75/, A(504) /187/, mass(504) /186.95575007100_dp/, mass_unc(504) /0.00000160900_dp/  ! 187Re
data Z(505) / 75/, A(505) /185/, mass(505) /184.95295448600_dp/, mass_unc(505) /0.00000132600_dp/  ! 185Re
data Z(506) / 75/, A(506) /186/, mass(506) /185.95498559500_dp/, mass_unc(506) /0.00000133600_dp/  ! 186mRe
data Z(507) / 75/, A(507) /184/, mass(507) /183.95252280900_dp/, mass_unc(507) /0.00000468800_dp/  ! 184mRe
data Z(508) / 75/, A(508) /183/, mass(508) /182.95081963800_dp/, mass_unc(508) /0.00000863500_dp/  ! 183Re
data Z(509) / 75/, A(509) /186/, mass(509) /185.95498559500_dp/, mass_unc(509) /0.00100133600_dp/  ! 186Re
data Z(510) / 75/, A(510) /182/, mass(510) /181.95120986900_dp/, mass_unc(510) /0.00010948400_dp/  ! 182Re
data Z(511) / 75/, A(511) /189/, mass(511) /188.95922602000_dp/, mass_unc(511) /0.00000892900_dp/  ! 189Re
data Z(512) / 76/, A(512) /192/, mass(512) /191.96147699800_dp/, mass_unc(512) /0.00000292600_dp/  ! 192Os
data Z(513) / 76/, A(513) /190/, mass(513) /189.95844370200_dp/, mass_unc(513) /0.00000170500_dp/  ! 190Os
data Z(514) / 76/, A(514) /189/, mass(514) /188.95814416200_dp/, mass_unc(514) /0.00000169400_dp/  ! 189Os
data Z(515) / 76/, A(515) /188/, mass(515) /187.95583517400_dp/, mass_unc(515) /0.00000161600_dp/  ! 188Os
data Z(516) / 76/, A(516) /187/, mass(516) /186.95574742200_dp/, mass_unc(516) /0.00000160900_dp/  ! 187Os
data Z(517) / 76/, A(517) /186/, mass(517) /185.95383504400_dp/, mass_unc(517) /0.00000159600_dp/  ! 186Os
data Z(518) / 76/, A(518) /184/, mass(518) /183.95248853600_dp/, mass_unc(518) /0.00000144100_dp/  ! 184Os
data Z(519) / 76/, A(519) /194/, mass(519) /193.96517724000_dp/, mass_unc(519) /0.00000298700_dp/  ! 194Os
data Z(520) / 76/, A(520) /185/, mass(520) /184.95404174100_dp/, mass_unc(520) /0.00000141000_dp/  ! 185Os
data Z(521) / 76/, A(521) /191/, mass(521) /190.96092636100_dp/, mass_unc(521) /0.00000170800_dp/  ! 191Os
data Z(522) / 76/, A(522) /193/, mass(522) /192.96414787000_dp/, mass_unc(522) /0.00000293100_dp/  ! 193Os
data Z(523) / 77/, A(523) /193/, mass(523) /192.96292158700_dp/, mass_unc(523) /0.00000207400_dp/  ! 193Ir
data Z(524) / 77/, A(524) /191/, mass(524) /190.96058929300_dp/, mass_unc(524) /0.00000206300_dp/  ! 191Ir
data Z(525) / 77/, A(525) /192/, mass(525) /191.96260024700_dp/, mass_unc(525) /0.00000206500_dp/  ! 192m2Ir
data Z(526) / 77/, A(526) /194/, mass(526) /193.96507353600_dp/, mass_unc(526) /0.00000207700_dp/  ! 194m2Ir
data Z(527) / 77/, A(527) /192/, mass(527) /191.96260024700_dp/, mass_unc(527) /0.00100206500_dp/  ! 192Ir
data Z(528) / 77/, A(528) /189/, mass(528) /188.95871502800_dp/, mass_unc(528) /0.00001371600_dp/  ! 189Ir
data Z(529) / 77/, A(529) /190/, mass(529) /189.96054121500_dp/, mass_unc(529) /0.00000210700_dp/  ! 190Ir
data Z(530) / 77/, A(530) /193/, mass(530) /192.96292158700_dp/, mass_unc(530) /0.00100207400_dp/  ! 193mIr
data Z(531) / 77/, A(531) /188/, mass(531) /187.95882809500_dp/, mass_unc(531) /0.00001026400_dp/  ! 188Ir
data Z(532) / 78/, A(532) /195/, mass(532) /194.96479171900_dp/, mass_unc(532) /0.00000099500_dp/  ! 195Pt
data Z(533) / 78/, A(533) /194/, mass(533) /193.96268085000_dp/, mass_unc(533) /0.00000100300_dp/  ! 194Pt
data Z(534) / 78/, A(534) /196/, mass(534) /195.96495209100_dp/, mass_unc(534) /0.00000098500_dp/  ! 196Pt
data Z(535) / 78/, A(535) /198/, mass(535) /197.96789492000_dp/, mass_unc(535) /0.00000232400_dp/  ! 198Pt
data Z(536) / 78/, A(536) /192/, mass(536) /191.96103874600_dp/, mass_unc(536) /0.00000322500_dp/  ! 192Pt
data Z(537) / 78/, A(537) /190/, mass(537) /189.95992970700_dp/, mass_unc(537) /0.00000626300_dp/  ! 190Pt
data Z(538) / 78/, A(538) /193/, mass(538) /192.96298238000_dp/, mass_unc(538) /0.00000209500_dp/  ! 193Pt
data Z(539) / 78/, A(539) /188/, mass(539) /187.95938888900_dp/, mass_unc(539) /0.00000605600_dp/  ! 188Pt
data Z(540) / 78/, A(540) /193/, mass(540) /192.96298238000_dp/, mass_unc(540) /0.00100209500_dp/  ! 193mPt
data Z(541) / 78/, A(541) /195/, mass(541) /194.96479171900_dp/, mass_unc(541) /0.00100099500_dp/  ! 195mPt
data Z(542) / 78/, A(542) /191/, mass(542) /190.96167291200_dp/, mass_unc(542) /0.00000532500_dp/  ! 191Pt
data Z(543) / 78/, A(543) /202/, mass(543) /201.97563900000_dp/, mass_unc(543) /0.00002700000_dp/  ! 202Pt
data Z(544) / 79/, A(544) /197/, mass(544) /196.96656878600_dp/, mass_unc(544) /0.00000070500_dp/  ! 197Au
data Z(545) / 79/, A(545) /195/, mass(545) /194.96503522500_dp/, mass_unc(545) /0.00000146300_dp/  ! 195Au
data Z(546) / 79/, A(546) /196/, mass(546) /195.96656990800_dp/, mass_unc(546) /0.00000320400_dp/  ! 196Au
data Z(547) / 79/, A(547) /199/, mass(547) /198.96876528200_dp/, mass_unc(547) /0.00000069900_dp/  ! 199Au
data Z(548) / 79/, A(548) /198/, mass(548) /197.96824242000_dp/, mass_unc(548) /0.00000069800_dp/  ! 198Au
data Z(549) / 79/, A(549) /198/, mass(549) /197.96824242000_dp/, mass_unc(549) /0.00100069800_dp/  ! 198m2Au
data Z(550) / 79/, A(550) /194/, mass(550) /193.96541775400_dp/, mass_unc(550) /0.00000230700_dp/  ! 194Au
data Z(551) / 80/, A(551) /202/, mass(551) /201.97064340000_dp/, mass_unc(551) /0.00000069200_dp/  ! 202Hg
data Z(552) / 80/, A(552) /200/, mass(552) /199.96832659000_dp/, mass_unc(552) /0.00000046700_dp/  ! 200Hg
data Z(553) / 80/, A(553) /199/, mass(553) /198.96828064300_dp/, mass_unc(553) /0.00000045700_dp/  ! 199Hg
data Z(554) / 80/, A(554) /201/, mass(554) /200.97030283900_dp/, mass_unc(554) /0.00000069300_dp/  ! 201Hg
data Z(555) / 80/, A(555) /198/, mass(555) /197.96676860200_dp/, mass_unc(555) /0.00000051500_dp/  ! 198Hg
data Z(556) / 80/, A(556) /204/, mass(556) /203.97349398100_dp/, mass_unc(556) /0.00000052700_dp/  ! 204Hg
data Z(557) / 80/, A(557) /196/, mass(557) /195.96583256000_dp/, mass_unc(557) /0.00000317200_dp/  ! 196Hg
data Z(558) / 80/, A(558) /194/, mass(558) /193.96544911200_dp/, mass_unc(558) /0.00000310000_dp/  ! 194Hg
data Z(559) / 80/, A(559) /203/, mass(559) /202.97287283700_dp/, mass_unc(559) /0.00000181400_dp/  ! 203Hg
data Z(560) / 80/, A(560) /197/, mass(560) /196.96721284700_dp/, mass_unc(560) /0.00000345000_dp/  ! 197Hg
data Z(561) / 80/, A(561) /195/, mass(561) /194.96672054000_dp/, mass_unc(561) /0.00002484700_dp/  ! 195mHg
data Z(562) / 81/, A(562) /205/, mass(562) /204.97442780100_dp/, mass_unc(562) /0.00000142600_dp/  ! 205Tl
data Z(563) / 81/, A(563) /203/, mass(563) /202.97234459300_dp/, mass_unc(563) /0.00000136300_dp/  ! 203Tl
data Z(564) / 81/, A(564) /204/, mass(564) /203.97386390000_dp/, mass_unc(564) /0.00000134100_dp/  ! 204Tl
data Z(565) / 81/, A(565) /202/, mass(565) /201.97210244100_dp/, mass_unc(565) /0.00001533500_dp/  ! 202Tl
data Z(566) / 81/, A(566) /201/, mass(566) /200.97082221200_dp/, mass_unc(566) /0.00001539700_dp/  ! 201Tl
data Z(567) / 81/, A(567) /200/, mass(567) /199.97096325800_dp/, mass_unc(567) /0.00000617400_dp/  ! 200Tl
data Z(568) / 82/, A(568) /208/, mass(568) /207.97665248100_dp/, mass_unc(568) /0.00000133300_dp/  ! 208Pb
data Z(569) / 82/, A(569) /206/, mass(569) /205.97446568300_dp/, mass_unc(569) /0.00000132900_dp/  ! 206Pb
data Z(570) / 82/, A(570) /207/, mass(570) /206.97589729700_dp/, mass_unc(570) /0.00000133200_dp/  ! 207Pb
data Z(571) / 82/, A(571) /204/, mass(571) /203.97304398100_dp/, mass_unc(571) /0.00000133300_dp/  ! 204Pb
data Z(572) / 82/, A(572) /205/, mass(572) /204.97448215700_dp/, mass_unc(572) /0.00000133000_dp/  ! 205Pb
data Z(573) / 82/, A(573) /202/, mass(573) /201.97215202600_dp/, mass_unc(573) /0.00000404000_dp/  ! 202Pb
data Z(574) / 82/, A(574) /210/, mass(574) /209.98418886100_dp/, mass_unc(574) /0.00000163800_dp/  ! 210Pb
data Z(575) / 82/, A(575) /203/, mass(575) /202.97339112600_dp/, mass_unc(575) /0.00000705500_dp/  ! 203Pb
data Z(576) / 83/, A(576) /209/, mass(576) /208.98039906800_dp/, mass_unc(576) /0.00000155300_dp/  ! 209Bi
data Z(577) / 83/, A(577) /210/, mass(577) /209.98412070500_dp/, mass_unc(577) /0.00000155200_dp/  ! 210mBi
data Z(578) / 83/, A(578) /208/, mass(578) /207.97974253000_dp/, mass_unc(578) /0.00000252700_dp/  ! 208Bi
data Z(579) / 83/, A(579) /207/, mass(579) /206.97847102200_dp/, mass_unc(579) /0.00000262400_dp/  ! 207Bi
data Z(580) / 83/, A(580) /205/, mass(580) /204.97738669400_dp/, mass_unc(580) /0.00000549300_dp/  ! 205Bi
data Z(581) / 83/, A(581) /206/, mass(581) /205.97849931700_dp/, mass_unc(581) /0.00000820900_dp/  ! 206Bi
data Z(582) / 83/, A(582) /210/, mass(582) /209.98412070500_dp/, mass_unc(582) /0.00100155200_dp/  ! 210Bi
data Z(583) / 84/, A(583) /209/, mass(583) /208.98243083600_dp/, mass_unc(583) /0.00000197500_dp/  ! 209Po
data Z(584) / 84/, A(584) /208/, mass(584) /207.98124609200_dp/, mass_unc(584) /0.00000193300_dp/  ! 208Po
data Z(585) / 84/, A(585) /210/, mass(585) /209.98287407600_dp/, mass_unc(585) /0.00000133100_dp/  ! 210Po
data Z(586) / 84/, A(586) /206/, mass(586) /205.98047399100_dp/, mass_unc(586) /0.00000427400_dp/  ! 206Po
data Z(587) / 86/, A(587) /222/, mass(587) /222.01757824600_dp/, mass_unc(587) /0.00000249200_dp/  ! 222Rn
data Z(588) / 88/, A(588) /226/, mass(588) /226.02541033000_dp/, mass_unc(588) /0.00000247700_dp/  ! 226Ra
data Z(589) / 88/, A(589) /228/, mass(589) /228.03107072800_dp/, mass_unc(589) /0.00000257600_dp/  ! 228Ra
data Z(590) / 88/, A(590) /225/, mass(590) /225.02361185700_dp/, mass_unc(590) /0.00000316600_dp/  ! 225Ra
data Z(591) / 88/, A(591) /223/, mass(591) /223.01850232700_dp/, mass_unc(591) /0.00000268400_dp/  ! 223Ra
data Z(592) / 88/, A(592) /224/, mass(592) /224.02021196800_dp/, mass_unc(592) /0.00000231600_dp/  ! 224Ra
data Z(593) / 89/, A(593) /227/, mass(593) /227.02775228300_dp/, mass_unc(593) /0.00000254400_dp/  ! 227Ac
data Z(594) / 89/, A(594) /225/, mass(594) /225.02322998700_dp/, mass_unc(594) /0.00000525600_dp/  ! 225Ac
data Z(595) / 89/, A(595) /226/, mass(595) /226.02609838300_dp/, mass_unc(595) /0.00000356100_dp/  ! 226Ac
data Z(596) / 90/, A(596) /232/, mass(596) /232.03805576000_dp/, mass_unc(596) /0.00000209200_dp/  ! 232Th
data Z(597) / 90/, A(597) /230/, mass(597) /230.03313413000_dp/, mass_unc(597) /0.00000188600_dp/  ! 230Th
data Z(598) / 90/, A(598) /229/, mass(598) /229.03176271300_dp/, mass_unc(598) /0.00000298900_dp/  ! 229Th
data Z(599) / 90/, A(599) /228/, mass(599) /228.02874127200_dp/, mass_unc(599) /0.00000231500_dp/  ! 228Th
data Z(600) / 90/, A(600) /234/, mass(600) /234.04360140700_dp/, mass_unc(600) /0.00000374300_dp/  ! 234Th
data Z(601) / 90/, A(601) /227/, mass(601) /227.02770422700_dp/, mass_unc(601) /0.00000268200_dp/  ! 227Th
data Z(602) / 90/, A(602) /231/, mass(602) /231.03630462800_dp/, mass_unc(602) /0.00000189500_dp/  ! 231Th
data Z(603) / 91/, A(603) /231/, mass(603) /231.03588424300_dp/, mass_unc(603) /0.00000239500_dp/  ! 231Pa
data Z(604) / 91/, A(604) /233/, mass(604) /233.04024722200_dp/, mass_unc(604) /0.00000224300_dp/  ! 233Pa
data Z(605) / 91/, A(605) /230/, mass(605) /230.03454104700_dp/, mass_unc(605) /0.00000350100_dp/  ! 230Pa
data Z(606) / 91/, A(606) /229/, mass(606) /229.03209723600_dp/, mass_unc(606) /0.00000378900_dp/  ! 229Pa
data Z(607) / 91/, A(607) /232/, mass(607) /232.03859173700_dp/, mass_unc(607) /0.00000830300_dp/  ! 232Pa
data Z(608) / 92/, A(608) /238/, mass(608) /238.05078842300_dp/, mass_unc(608) /0.00000201900_dp/  ! 238U
data Z(609) / 92/, A(609) /235/, mass(609) /235.04393013100_dp/, mass_unc(609) /0.00000192000_dp/  ! 235U
data Z(610) / 92/, A(610) /234/, mass(610) /234.04095230600_dp/, mass_unc(610) /0.00000192300_dp/  ! 234U
data Z(611) / 92/, A(611) /236/, mass(611) /236.04556821000_dp/, mass_unc(611) /0.00000192300_dp/  ! 236U
data Z(612) / 92/, A(612) /233/, mass(612) /233.03963552500_dp/, mass_unc(612) /0.00000286300_dp/  ! 233U
data Z(613) / 92/, A(613) /232/, mass(613) /232.03715629700_dp/, mass_unc(613) /0.00000231700_dp/  ! 232U
data Z(614) / 92/, A(614) /230/, mass(614) /230.03394009600_dp/, mass_unc(614) /0.00000510600_dp/  ! 230U
data Z(615) / 92/, A(615) /237/, mass(615) /237.04873037800_dp/, mass_unc(615) /0.00000198500_dp/  ! 237U
data Z(616) / 92/, A(616) /231/, mass(616) /231.03629386100_dp/, mass_unc(616) /0.00000322300_dp/  ! 231U
data Z(617) / 93/, A(617) /237/, mass(617) /237.04817364900_dp/, mass_unc(617) /0.00000193000_dp/  ! 237Np
data Z(618) / 93/, A(618) /236/, mass(618) /236.04656974400_dp/, mass_unc(618) /0.00005414400_dp/  ! 236Np
data Z(619) / 93/, A(619) /235/, mass(619) /235.04406348700_dp/, mass_unc(619) /0.00000210600_dp/  ! 235Np
data Z(620) / 93/, A(620) /234/, mass(620) /234.04289525600_dp/, mass_unc(620) /0.00000913700_dp/  ! 234Np
data Z(621) / 93/, A(621) /239/, mass(621) /239.05293924100_dp/, mass_unc(621) /0.00000219700_dp/  ! 239Np
data Z(622) / 93/, A(622) /238/, mass(622) /238.05094661100_dp/, mass_unc(622) /0.00000194200_dp/  ! 238Np
data Z(623) / 94/, A(623) /244/, mass(623) /244.06420526000_dp/, mass_unc(623) /0.00000557100_dp/  ! 244Pu
data Z(624) / 94/, A(624) /242/, mass(624) /242.05874280900_dp/, mass_unc(624) /0.00000196300_dp/  ! 242Pu
data Z(625) / 94/, A(625) /239/, mass(625) /239.05216359100_dp/, mass_unc(625) /0.00000192400_dp/  ! 239Pu
data Z(626) / 94/, A(626) /240/, mass(626) /240.05381375000_dp/, mass_unc(626) /0.00000192100_dp/  ! 240Pu
data Z(627) / 94/, A(627) /238/, mass(627) /238.04956011100_dp/, mass_unc(627) /0.00000193100_dp/  ! 238Pu
data Z(628) / 94/, A(628) /241/, mass(628) /241.05685166100_dp/, mass_unc(628) /0.00000192100_dp/  ! 241Pu
data Z(629) / 94/, A(629) /236/, mass(629) /236.04605810900_dp/, mass_unc(629) /0.00000231800_dp/  ! 236Pu
data Z(630) / 94/, A(630) /237/, mass(630) /237.04840983000_dp/, mass_unc(630) /0.00000235800_dp/  ! 237Pu
data Z(631) / 94/, A(631) /246/, mass(631) /246.07020545800_dp/, mass_unc(631) /0.00001638800_dp/  ! 246Pu
data Z(632) / 94/, A(632) /247/, mass(632) /247.07419000000_dp/, mass_unc(632) /0.00021000000_dp/  ! 247Pu
data Z(633) / 95/, A(633) /243/, mass(633) /243.06138130200_dp/, mass_unc(633) /0.00000243200_dp/  ! 243Am
data Z(634) / 95/, A(634) /241/, mass(634) /241.05682934900_dp/, mass_unc(634) /0.00000192600_dp/  ! 241Am
data Z(635) / 95/, A(635) /242/, mass(635) /242.05954936400_dp/, mass_unc(635) /0.00000192900_dp/  ! 242mAm
data Z(636) / 95/, A(636) /240/, mass(636) /240.05530038400_dp/, mass_unc(636) /0.00001492500_dp/  ! 240Am
data Z(637) / 96/, A(637) /247/, mass(637) /247.07035413100_dp/, mass_unc(637) /0.00000469000_dp/  ! 247Cm
data Z(638) / 96/, A(638) /248/, mass(638) /248.07234986200_dp/, mass_unc(638) /0.00000557500_dp/  ! 248Cm
data Z(639) / 96/, A(639) /245/, mass(639) /245.06549145400_dp/, mass_unc(639) /0.00000220100_dp/  ! 245Cm
data Z(640) / 96/, A(640) /250/, mass(640) /250.07835831300_dp/, mass_unc(640) /0.00001209600_dp/  ! 250Cm
data Z(641) / 96/, A(641) /246/, mass(641) /246.06722384100_dp/, mass_unc(641) /0.00000218300_dp/  ! 246Cm
data Z(642) / 96/, A(642) /243/, mass(642) /243.06138932500_dp/, mass_unc(642) /0.00000220300_dp/  ! 243Cm
data Z(643) / 96/, A(643) /244/, mass(643) /244.06275278300_dp/, mass_unc(643) /0.00000192200_dp/  ! 244Cm
data Z(644) / 96/, A(644) /242/, mass(644) /242.05883603900_dp/, mass_unc(644) /0.00000193300_dp/  ! 242Cm
data Z(645) / 96/, A(645) /241/, mass(645) /241.05765317000_dp/, mass_unc(645) /0.00000228500_dp/  ! 241Cm
data Z(646) / 96/, A(646) /240/, mass(646) /240.05552968100_dp/, mass_unc(646) /0.00000240400_dp/  ! 240Cm
data Z(647) / 96/, A(647) /252/, mass(647) /252.08487000000_dp/, mass_unc(647) /0.00032000000_dp/  ! 252Cm
data Z(648) / 97/, A(648) /247/, mass(648) /247.07030730200_dp/, mass_unc(648) /0.00000589300_dp/  ! 247Bk
data Z(649) / 97/, A(649) /248/, mass(649) /248.07308800000_dp/, mass_unc(649) /0.00007600000_dp/  ! 248Bk
data Z(650) / 97/, A(650) /249/, mass(650) /249.07498767600_dp/, mass_unc(650) /0.00000273800_dp/  ! 249Bk
data Z(651) / 97/, A(651) /245/, mass(651) /245.06636182100_dp/, mass_unc(651) /0.00000244600_dp/  ! 245Bk
data Z(652) / 97/, A(652) /246/, mass(652) /246.06867312600_dp/, mass_unc(652) /0.00006444900_dp/  ! 246Bk
data Z(653) / 98/, A(653) /251/, mass(653) /251.07958862500_dp/, mass_unc(653) /0.00000478800_dp/  ! 251Cf
data Z(654) / 98/, A(654) /249/, mass(654) /249.07485390300_dp/, mass_unc(654) /0.00000232400_dp/  ! 249Cf
data Z(655) / 98/, A(655) /250/, mass(655) /250.07640624400_dp/, mass_unc(655) /0.00000219300_dp/  ! 250Cf
data Z(656) / 98/, A(656) /252/, mass(656) /252.08162719900_dp/, mass_unc(656) /0.00000557500_dp/  ! 252Cf
data Z(657) / 98/, A(657) /248/, mass(657) /248.07218506600_dp/, mass_unc(657) /0.00000570100_dp/  ! 248Cf
data Z(658) / 98/, A(658) /254/, mass(658) /254.08732426300_dp/, mass_unc(658) /0.00001326900_dp/  ! 254Cf
data Z(659) / 98/, A(659) /253/, mass(659) /253.08513449900_dp/, mass_unc(659) /0.00000674900_dp/  ! 253Cf
data Z(660) / 98/, A(660) /246/, mass(660) /246.06880553100_dp/, mass_unc(660) /0.00000220600_dp/  ! 246Cf
data Z(661) / 99/, A(661) /252/, mass(661) /252.08297986500_dp/, mass_unc(661) /0.00005396500_dp/  ! 252Es
data Z(662) / 99/, A(662) /254/, mass(662) /254.08802219900_dp/, mass_unc(662) /0.00000454000_dp/  ! 254Es
data Z(663) / 99/, A(663) /255/, mass(663) /255.09027495800_dp/, mass_unc(663) /0.00001184200_dp/  ! 255Es
data Z(664) / 99/, A(664) /253/, mass(664) /253.08482571500_dp/, mass_unc(664) /0.00000273800_dp/  ! 253Es
data Z(665) / 99/, A(665) /257/, mass(665) /257.09597900000_dp/, mass_unc(665) /0.00044100000_dp/  ! 257Es
data Z(666) / 99/, A(666) /254/, mass(666) /254.08802219900_dp/, mass_unc(666) /0.00100454000_dp/  ! 254mEs
data Z(667) / 99/, A(667) /251/, mass(667) /251.07999358600_dp/, mass_unc(667) /0.00000671500_dp/  ! 251Es
data Z(668) /100/, A(668) /257/, mass(668) /257.09510607800_dp/, mass_unc(668) /0.00000691800_dp/  ! 257Fm
data Z(669) /100/, A(669) /253/, mass(669) /253.08518457100_dp/, mass_unc(669) /0.00000370100_dp/  ! 253Fm
data Z(670) /100/, A(670) /252/, mass(670) /252.08246706000_dp/, mass_unc(670) /0.00000609200_dp/  ! 252Fm
data Z(671) /101/, A(671) /258/, mass(671) /258.09843149600_dp/, mass_unc(671) /0.00000495800_dp/  ! 258Md
data Z(672) /101/, A(672) /260/, mass(672) /260.10365300000_dp/, mass_unc(672) /0.00034000000_dp/  ! 260Md
data Z(673) /105/, A(673) /268/, mass(673) /268.12567100000_dp/, mass_unc(673) /0.00056800000_dp/  ! 268Db

! nuclear stability data from the NUBASE2012 database
data flag(  0) / 0/, h_f(  0) /long_time /
data flag(  1) / 1/, h_f(  1) /long_time /         ! 1H
data flag(  2) / 1/, h_f(  2) /long_time /         ! 2H
data flag(  3) / 4/, h_f(  3) /3.8878E+08_dp/      ! 3H
data flag(  4) / 1/, h_f(  4) /long_time /         ! 4He
data flag(  5) / 1/, h_f(  5) /long_time /         ! 3He
data flag(  6) / 1/, h_f(  6) /long_time /         ! 7Li
data flag(  7) / 1/, h_f(  7) /long_time /         ! 6Li
data flag(  8) / 1/, h_f(  8) /long_time /         ! 9Be
data flag(  9) / 4/, h_f(  9) /4.7651E+13_dp/      ! 10Be
data flag( 10) / 4/, h_f( 10) /4.5982E+06_dp/      ! 7Be
data flag( 11) / 1/, h_f( 11) /long_time /         ! 11B
data flag( 12) / 1/, h_f( 12) /long_time /         ! 10B
data flag( 13) / 1/, h_f( 13) /long_time /         ! 12C
data flag( 14) / 1/, h_f( 14) /long_time /         ! 13C
data flag( 15) / 4/, h_f( 15) /1.7987E+11_dp/      ! 14C
data flag( 16) / 1/, h_f( 16) /long_time /         ! 14N
data flag( 17) / 1/, h_f( 17) /long_time /         ! 15N
data flag( 18) / 1/, h_f( 18) /long_time /         ! 16O
data flag( 19) / 1/, h_f( 19) /long_time /         ! 18O
data flag( 20) / 1/, h_f( 20) /long_time /         ! 17O
data flag( 21) / 1/, h_f( 21) /long_time /         ! 19F
data flag( 22) / 1/, h_f( 22) /long_time /         ! 20Ne
data flag( 23) / 1/, h_f( 23) /long_time /         ! 22Ne
data flag( 24) / 1/, h_f( 24) /long_time /         ! 21Ne
data flag( 25) / 1/, h_f( 25) /long_time /         ! 23Na
data flag( 26) / 4/, h_f( 26) /8.2133E+07_dp/      ! 22Na
data flag( 27) / 1/, h_f( 27) /long_time /         ! 24Mg
data flag( 28) / 1/, h_f( 28) /long_time /         ! 26Mg
data flag( 29) / 1/, h_f( 29) /long_time /         ! 25Mg
data flag( 30) / 1/, h_f( 30) /long_time /         ! 27Al
data flag( 31) / 4/, h_f( 31) /2.2626E+13_dp/      ! 26Al
data flag( 32) / 1/, h_f( 32) /long_time /         ! 28Si
data flag( 33) / 1/, h_f( 33) /long_time /         ! 29Si
data flag( 34) / 1/, h_f( 34) /long_time /         ! 30Si
data flag( 35) / 4/, h_f( 35) /4.8282E+09_dp/      ! 32Si
data flag( 36) / 1/, h_f( 36) /long_time /         ! 31P
data flag( 37) / 4/, h_f( 37) /2.1902E+06_dp/      ! 33P
data flag( 38) / 4/, h_f( 38) /1.2323E+06_dp/      ! 32P
data flag( 39) / 1/, h_f( 39) /long_time /         ! 32S
data flag( 40) / 1/, h_f( 40) /long_time /         ! 34S
data flag( 41) / 1/, h_f( 41) /long_time /         ! 33S
data flag( 42) / 1/, h_f( 42) /long_time /         ! 36S
data flag( 43) / 4/, h_f( 43) /7.5488E+06_dp/      ! 35S
data flag( 44) / 1/, h_f( 44) /long_time /         ! 35Cl
data flag( 45) / 1/, h_f( 45) /long_time /         ! 37Cl
data flag( 46) / 4/, h_f( 46) /9.5081E+12_dp/      ! 36Cl
data flag( 47) / 1/, h_f( 47) /long_time /         ! 40Ar
data flag( 48) / 1/, h_f( 48) /long_time /         ! 36Ar
data flag( 49) / 1/, h_f( 49) /long_time /         ! 38Ar
data flag( 50) / 4/, h_f( 50) /8.4888E+09_dp/      ! 39Ar
data flag( 51) / 4/, h_f( 51) /1.0382E+09_dp/      ! 42Ar
data flag( 52) / 4/, h_f( 52) /3.0250E+06_dp/      ! 37Ar
data flag( 53) / 1/, h_f( 53) /long_time /         ! 39K
data flag( 54) / 1/, h_f( 54) /long_time /         ! 41K
data flag( 55) / 3/, h_f( 55) /3.9383E+16_dp/      ! 40K
data flag( 56) / 2/, h_f( 56) /long_time /         ! 40Ca
data flag( 57) / 1/, h_f( 57) /long_time /         ! 44Ca
data flag( 58) / 1/, h_f( 58) /long_time /         ! 42Ca
data flag( 59) / 3/, h_f( 59) /1.6725E+27_dp/      ! 48Ca
data flag( 60) / 1/, h_f( 60) /long_time /         ! 43Ca
data flag( 61) / 2/, h_f( 61) /long_time /         ! 46Ca
data flag( 62) / 4/, h_f( 62) /3.1368E+12_dp/      ! 41Ca
data flag( 63) / 4/, h_f( 63) /1.4050E+07_dp/      ! 45Ca
data flag( 64) / 4/, h_f( 64) /3.9191E+05_dp/      ! 47Ca
data flag( 65) / 1/, h_f( 65) /long_time /         ! 45Sc
data flag( 66) / 4/, h_f( 66) /7.2395E+06_dp/      ! 46Sc
data flag( 67) / 4/, h_f( 67) /2.8937E+05_dp/      ! 47Sc
data flag( 68) /34/, h_f( 68) /2.1100E+05_dp/      ! 44m3Sc
data flag( 69) / 4/, h_f( 69) /1.5721E+05_dp/      ! 48Sc
data flag( 70) / 1/, h_f( 70) /long_time /         ! 48Ti
data flag( 71) / 1/, h_f( 71) /long_time /         ! 46Ti
data flag( 72) / 1/, h_f( 72) /long_time /         ! 47Ti
data flag( 73) / 1/, h_f( 73) /long_time /         ! 49Ti
data flag( 74) / 1/, h_f( 74) /long_time /         ! 50Ti
data flag( 75) / 4/, h_f( 75) /1.8650E+09_dp/      ! 44Ti
data flag( 76) / 1/, h_f( 76) /long_time /         ! 51V
data flag( 77) / 3/, h_f( 77) /4.7335E+24_dp/      ! 50V
data flag( 78) / 4/, h_f( 78) /2.8512E+07_dp/      ! 49V
data flag( 79) / 4/, h_f( 79) /1.3801E+06_dp/      ! 48V
data flag( 80) / 1/, h_f( 80) /long_time /         ! 52Cr
data flag( 81) / 1/, h_f( 81) /long_time /         ! 53Cr
data flag( 82) / 2/, h_f( 82) /long_time /         ! 50Cr
data flag( 83) / 1/, h_f( 83) /long_time /         ! 54Cr
data flag( 84) / 4/, h_f( 84) /2.3934E+06_dp/      ! 51Cr
data flag( 85) / 1/, h_f( 85) /long_time /         ! 55Mn
data flag( 86) / 4/, h_f( 86) /1.1676E+14_dp/      ! 53Mn
data flag( 87) / 4/, h_f( 87) /2.6961E+07_dp/      ! 54Mn
data flag( 88) / 4/, h_f( 88) /4.8306E+05_dp/      ! 52Mn
data flag( 89) / 1/, h_f( 89) /long_time /         ! 56Fe
data flag( 90) / 1/, h_f( 90) /long_time /         ! 54Fe
data flag( 91) / 1/, h_f( 91) /long_time /         ! 57Fe
data flag( 92) / 1/, h_f( 92) /long_time /         ! 58Fe
data flag( 93) / 4/, h_f( 93) /8.2679E+13_dp/      ! 60Fe
data flag( 94) / 4/, h_f( 94) /8.6592E+07_dp/      ! 55Fe
data flag( 95) / 4/, h_f( 95) /3.8444E+06_dp/      ! 59Fe
data flag( 96) / 1/, h_f( 96) /long_time /         ! 59Co
data flag( 97) / 4/, h_f( 97) /1.6634E+08_dp/      ! 60Co
data flag( 98) / 4/, h_f( 98) /2.3478E+07_dp/      ! 57Co
data flag( 99) / 4/, h_f( 99) /6.6732E+06_dp/      ! 56Co
data flag(100) / 4/, h_f(100) /6.1223E+06_dp/      ! 58Co
data flag(101) / 1/, h_f(101) /long_time /         ! 58Ni
data flag(102) / 1/, h_f(102) /long_time /         ! 60Ni
data flag(103) / 1/, h_f(103) /long_time /         ! 62Ni
data flag(104) / 1/, h_f(104) /long_time /         ! 61Ni
data flag(105) / 1/, h_f(105) /long_time /         ! 64Ni
data flag(106) / 4/, h_f(106) /3.1872E+12_dp/      ! 59Ni
data flag(107) / 4/, h_f(107) /3.1936E+09_dp/      ! 63Ni
data flag(108) / 4/, h_f(108) /5.2488E+05_dp/      ! 56Ni
data flag(109) / 4/, h_f(109) /1.9656E+05_dp/      ! 66Ni
data flag(110) / 4/, h_f(110) /1.2816E+05_dp/      ! 57Ni
data flag(111) / 1/, h_f(111) /long_time /         ! 63Cu
data flag(112) / 1/, h_f(112) /long_time /         ! 65Cu
data flag(113) / 4/, h_f(113) /2.2259E+05_dp/      ! 67Cu
data flag(114) / 2/, h_f(114) /long_time /         ! 64Zn
data flag(115) / 1/, h_f(115) /long_time /         ! 66Zn
data flag(116) / 1/, h_f(116) /long_time /         ! 68Zn
data flag(117) / 1/, h_f(117) /long_time /         ! 67Zn
data flag(118) / 2/, h_f(118) /long_time /         ! 70Zn
data flag(119) / 4/, h_f(119) /2.1076E+07_dp/      ! 65Zn
data flag(120) / 4/, h_f(120) /1.6740E+05_dp/      ! 72Zn
data flag(121) / 1/, h_f(121) /long_time /         ! 69Ga
data flag(122) / 1/, h_f(122) /long_time /         ! 71Ga
data flag(123) / 4/, h_f(123) /2.8181E+05_dp/      ! 67Ga
data flag(124) / 1/, h_f(124) /long_time /         ! 74Ge
data flag(125) / 1/, h_f(125) /long_time /         ! 72Ge
data flag(126) / 1/, h_f(126) /long_time /         ! 70Ge
data flag(127) / 1/, h_f(127) /long_time /         ! 73Ge
data flag(128) / 3/, h_f(128) /4.9860E+28_dp/      ! 76Ge
data flag(129) / 4/, h_f(129) /2.3408E+07_dp/      ! 68Ge
data flag(130) / 4/, h_f(130) /9.8755E+05_dp/      ! 71Ge
data flag(131) / 4/, h_f(131) /1.4058E+05_dp/      ! 69Ge
data flag(132) / 1/, h_f(132) /long_time /         ! 75As
data flag(133) / 4/, h_f(133) /6.9379E+06_dp/      ! 73As
data flag(134) / 4/, h_f(134) /1.5353E+06_dp/      ! 74As
data flag(135) / 4/, h_f(135) /2.3508E+05_dp/      ! 71As
data flag(136) / 4/, h_f(136) /1.3964E+05_dp/      ! 77As
data flag(137) / 4/, h_f(137) /9.3600E+04_dp/      ! 72As
data flag(138) / 4/, h_f(138) /9.3122E+04_dp/      ! 76As
data flag(139) / 1/, h_f(139) /long_time /         ! 80Se
data flag(140) / 1/, h_f(140) /long_time /         ! 78Se
data flag(141) / 1/, h_f(141) /long_time /         ! 76Se
data flag(142) / 3/, h_f(142) /3.0610E+27_dp/      ! 82Se
data flag(143) / 1/, h_f(143) /long_time /         ! 77Se
data flag(144) / 1/, h_f(144) /long_time /         ! 74Se
data flag(145) / 4/, h_f(145) /1.0572E+13_dp/      ! 79Se
data flag(146) / 4/, h_f(146) /1.0349E+07_dp/      ! 75Se
data flag(147) / 4/, h_f(147) /7.2576E+05_dp/      ! 72Se
data flag(148) / 1/, h_f(148) /long_time /         ! 79Br
data flag(149) / 1/, h_f(149) /long_time /         ! 81Br
data flag(150) / 4/, h_f(150) /2.0534E+05_dp/      ! 77Br
data flag(151) / 4/, h_f(151) /1.2702E+05_dp/      ! 82Br
data flag(152) / 1/, h_f(152) /long_time /         ! 84Kr
data flag(153) / 1/, h_f(153) /long_time /         ! 86Kr
data flag(154) / 1/, h_f(154) /long_time /         ! 82Kr
data flag(155) / 1/, h_f(155) /long_time /         ! 83Kr
data flag(156) / 1/, h_f(156) /long_time /         ! 80Kr
data flag(157) / 2/, h_f(157) /long_time /         ! 78Kr
data flag(158) / 4/, h_f(158) /7.2265E+12_dp/      ! 81Kr
data flag(159) / 4/, h_f(159) /3.4006E+08_dp/      ! 85Kr
data flag(160) / 4/, h_f(160) /1.2614E+05_dp/      ! 79Kr
data flag(161) / 1/, h_f(161) /long_time /         ! 85Rb
data flag(162) / 3/, h_f(162) /1.5535E+18_dp/      ! 87Rb
data flag(163) / 4/, h_f(163) /7.4477E+06_dp/      ! 83Rb
data flag(164) / 4/, h_f(164) /2.8356E+06_dp/      ! 84Rb
data flag(165) / 4/, h_f(165) /1.6107E+06_dp/      ! 86Rb
data flag(166) / 1/, h_f(166) /long_time /         ! 88Sr
data flag(167) / 1/, h_f(167) /long_time /         ! 86Sr
data flag(168) / 1/, h_f(168) /long_time /         ! 87Sr
data flag(169) / 1/, h_f(169) /long_time /         ! 84Sr
data flag(170) / 4/, h_f(170) /9.0852E+08_dp/      ! 90Sr
data flag(171) / 4/, h_f(171) /5.6033E+06_dp/      ! 85Sr
data flag(172) / 4/, h_f(172) /4.3658E+06_dp/      ! 89Sr
data flag(173) / 4/, h_f(173) /2.1911E+06_dp/      ! 82Sr
data flag(174) / 4/, h_f(174) /1.1668E+05_dp/      ! 83Sr
data flag(175) / 1/, h_f(175) /long_time /         ! 89Y
data flag(176) / 4/, h_f(176) /9.2125E+06_dp/      ! 88Y
data flag(177) / 4/, h_f(177) /5.0553E+06_dp/      ! 91Y
data flag(178) / 4/, h_f(178) /2.8728E+05_dp/      ! 87Y
data flag(179) / 4/, h_f(179) /2.3040E+05_dp/      ! 90Y
data flag(180) / 1/, h_f(180) /long_time /         ! 90Zr
data flag(181) / 1/, h_f(181) /long_time /         ! 94Zr
data flag(182) / 1/, h_f(182) /long_time /         ! 92Zr
data flag(183) / 1/, h_f(183) /long_time /         ! 91Zr
data flag(184) / 3/, h_f(184) /6.3114E+26_dp/      ! 96Zr
data flag(185) / 4/, h_f(185) /5.0807E+13_dp/      ! 93Zr
data flag(186) / 4/, h_f(186) /7.2058E+06_dp/      ! 88Zr
data flag(187) / 4/, h_f(187) /5.5324E+06_dp/      ! 95Zr
data flag(188) / 4/, h_f(188) /2.8228E+05_dp/      ! 89Zr
data flag(189) / 1/, h_f(189) /long_time /         ! 93Nb
data flag(190) / 4/, h_f(190) /1.0950E+15_dp/      ! 92Nb
data flag(191) / 4/, h_f(191) /6.4061E+11_dp/      ! 94Nb
data flag(192) / 4/, h_f(192) /2.1459E+10_dp/      ! 91Nb
data flag(193) /14/, h_f(193) /5.0870E+08_dp/      ! 93mNb
data flag(194) / 4/, h_f(194) /3.0232E+06_dp/      ! 95Nb
data flag(195) /14/, h_f(195) /3.1190E+05_dp/      ! 95mNb
data flag(196) / 1/, h_f(196) /long_time /         ! 98Mo
data flag(197) / 1/, h_f(197) /long_time /         ! 96Mo
data flag(198) / 1/, h_f(198) /long_time /         ! 95Mo
data flag(199) / 1/, h_f(199) /long_time /         ! 92Mo
data flag(200) / 3/, h_f(200) /2.3037E+26_dp/      ! 100Mo
data flag(201) / 1/, h_f(201) /long_time /         ! 97Mo
data flag(202) / 1/, h_f(202) /long_time /         ! 94Mo
data flag(203) / 4/, h_f(203) /1.2623E+11_dp/      ! 93Mo
data flag(204) / 4/, h_f(204) /2.3751E+05_dp/      ! 99Mo
data flag(205) / 4/, h_f(205) /1.3285E+14_dp/      ! 97Tc
data flag(206) / 4/, h_f(206) /1.3254E+14_dp/      ! 98Tc
data flag(207) / 4/, h_f(207) /6.6617E+12_dp/      ! 99Tc
data flag(208) /14/, h_f(208) /7.8624E+06_dp/      ! 97mTc
data flag(209) /14/, h_f(209) /5.2704E+06_dp/      ! 95mTc
data flag(210) / 4/, h_f(210) /3.6979E+05_dp/      ! 96Tc
data flag(211) / 1/, h_f(211) /long_time /         ! 102Ru
data flag(212) / 1/, h_f(212) /long_time /         ! 104Ru
data flag(213) / 1/, h_f(213) /long_time /         ! 101Ru
data flag(214) / 1/, h_f(214) /long_time /         ! 99Ru
data flag(215) / 1/, h_f(215) /long_time /         ! 100Ru
data flag(216) / 1/, h_f(216) /long_time /         ! 96Ru
data flag(217) / 1/, h_f(217) /long_time /         ! 98Ru
data flag(218) / 4/, h_f(218) /3.2124E+07_dp/      ! 106Ru
data flag(219) / 4/, h_f(219) /3.3909E+06_dp/      ! 103Ru
data flag(220) / 4/, h_f(220) /2.4512E+05_dp/      ! 97Ru
data flag(221) / 1/, h_f(221) /long_time /         ! 103Rh
data flag(222) /14/, h_f(222) /1.1809E+08_dp/      ! 102mRh
data flag(223) / 4/, h_f(223) /1.0414E+08_dp/      ! 101Rh
data flag(224) / 4/, h_f(224) /1.7885E+07_dp/      ! 102Rh
data flag(225) / 4/, h_f(225) /1.3910E+06_dp/      ! 99Rh
data flag(226) /14/, h_f(226) /3.7498E+05_dp/      ! 101mRh
data flag(227) / 4/, h_f(227) /1.2729E+05_dp/      ! 105Rh
data flag(228) / 1/, h_f(228) /long_time /         ! 106Pd
data flag(229) / 1/, h_f(229) /long_time /         ! 108Pd
data flag(230) / 1/, h_f(230) /long_time /         ! 105Pd
data flag(231) / 1/, h_f(231) /long_time /         ! 110Pd
data flag(232) / 1/, h_f(232) /long_time /         ! 104Pd
data flag(233) / 1/, h_f(233) /long_time /         ! 102Pd
data flag(234) / 4/, h_f(234) /2.0512E+14_dp/      ! 107Pd
data flag(235) / 4/, h_f(235) /1.4680E+06_dp/      ! 103Pd
data flag(236) / 4/, h_f(236) /3.1363E+05_dp/      ! 100Pd
data flag(237) / 1/, h_f(237) /long_time /         ! 107Ag
data flag(238) / 1/, h_f(238) /long_time /         ! 109Ag
data flag(239) /14/, h_f(239) /1.3822E+10_dp/      ! 108mAg
data flag(240) /24/, h_f(240) /2.1585E+07_dp/      ! 110m2Ag
data flag(241) / 4/, h_f(241) /3.5675E+06_dp/      ! 105Ag
data flag(242) /14/, h_f(242) /7.1539E+05_dp/      ! 106mAg
data flag(243) / 4/, h_f(243) /6.4368E+05_dp/      ! 111Ag
data flag(244) / 2/, h_f(244) /long_time /         ! 114Cd
data flag(245) / 1/, h_f(245) /long_time /         ! 112Cd
data flag(246) / 1/, h_f(246) /long_time /         ! 111Cd
data flag(247) / 1/, h_f(247) /long_time /         ! 110Cd
data flag(248) / 3/, h_f(248) /2.5372E+23_dp/      ! 113Cd
data flag(249) / 3/, h_f(249) /9.4671E+26_dp/      ! 116Cd
data flag(250) / 2/, h_f(250) /long_time /         ! 106Cd
data flag(251) / 2/, h_f(251) /long_time /         ! 108Cd
data flag(252) /14/, h_f(252) /4.3833E+08_dp/      ! 113mCd
data flag(253) / 4/, h_f(253) /3.9865E+07_dp/      ! 109Cd
data flag(254) /14/, h_f(254) /3.8500E+06_dp/      ! 115mCd
data flag(255) / 4/, h_f(255) /1.9246E+05_dp/      ! 115Cd
data flag(256) / 3/, h_f(256) /1.3917E+22_dp/      ! 115In
data flag(257) / 1/, h_f(257) /long_time /         ! 113In
data flag(258) /14/, h_f(258) /4.2777E+06_dp/      ! 114mIn
data flag(259) / 4/, h_f(259) /2.4233E+05_dp/      ! 111In
data flag(260) / 1/, h_f(260) /long_time /         ! 120Sn
data flag(261) / 1/, h_f(261) /long_time /         ! 118Sn
data flag(262) / 1/, h_f(262) /long_time /         ! 116Sn
data flag(263) / 1/, h_f(263) /long_time /         ! 119Sn
data flag(264) / 1/, h_f(264) /long_time /         ! 117Sn
data flag(265) / 1/, h_f(265) /long_time /         ! 124Sn
data flag(266) / 1/, h_f(266) /long_time /         ! 122Sn
data flag(267) / 1/, h_f(267) /long_time /         ! 112Sn
data flag(268) / 1/, h_f(268) /long_time /         ! 114Sn
data flag(269) / 1/, h_f(269) /long_time /         ! 115Sn
data flag(270) / 4/, h_f(270) /7.2581E+12_dp/      ! 126Sn
data flag(271) /14/, h_f(271) /1.3853E+09_dp/      ! 121mSn
data flag(272) /14/, h_f(272) /2.5324E+07_dp/      ! 119mSn
data flag(273) / 4/, h_f(273) /1.1163E+07_dp/      ! 123Sn
data flag(274) / 4/, h_f(274) /9.9438E+06_dp/      ! 113Sn
data flag(275) /14/, h_f(275) /1.2096E+06_dp/      ! 117mSn
data flag(276) / 4/, h_f(276) /8.3290E+05_dp/      ! 125Sn
data flag(277) / 4/, h_f(277) /9.7308E+04_dp/      ! 121Sn
data flag(278) / 1/, h_f(278) /long_time /         ! 121Sb
data flag(279) / 1/, h_f(279) /long_time /         ! 123Sb
data flag(280) / 4/, h_f(280) /8.7053E+07_dp/      ! 125Sb
data flag(281) / 4/, h_f(281) /5.2013E+06_dp/      ! 124Sb
data flag(282) / 4/, h_f(282) /1.0670E+06_dp/      ! 126Sb
data flag(283) /14/, h_f(283) /4.9766E+05_dp/      ! 120mSb
data flag(284) / 4/, h_f(284) /3.3264E+05_dp/      ! 127Sb
data flag(285) / 4/, h_f(285) /2.3534E+05_dp/      ! 122Sb
data flag(286) / 4/, h_f(286) /1.3748E+05_dp/      ! 119Sb
data flag(287) / 3/, h_f(287) /2.4930E+28_dp/      ! 130Te
data flag(288) / 1/, h_f(288) /long_time /         ! 126Te
data flag(289) / 1/, h_f(289) /long_time /         ! 125Te
data flag(290) / 1/, h_f(290) /long_time /         ! 124Te
data flag(291) / 1/, h_f(291) /long_time /         ! 122Te
data flag(292) / 2/, h_f(292) /long_time /         ! 123Te
data flag(293) / 2/, h_f(293) /long_time /         ! 120Te
data flag(294) / 3/, h_f(294) /6.9425E+31_dp/      ! 128Te
data flag(295) /14/, h_f(295) /1.4187E+07_dp/      ! 121mTe
data flag(296) /14/, h_f(296) /1.0299E+07_dp/      ! 123mTe
data flag(297) /14/, h_f(297) /9.1670E+06_dp/      ! 127mTe
data flag(298) /14/, h_f(298) /4.9594E+06_dp/      ! 125mTe
data flag(299) /14/, h_f(299) /2.9030E+06_dp/      ! 129mTe
data flag(300) / 4/, h_f(300) /1.6563E+06_dp/      ! 121Te
data flag(301) / 4/, h_f(301) /5.1840E+05_dp/      ! 118Te
data flag(302) /14/, h_f(302) /4.0608E+05_dp/      ! 119mTe
data flag(303) / 4/, h_f(303) /2.7683E+05_dp/      ! 132Te
data flag(304) /14/, h_f(304) /1.1693E+05_dp/      ! 131mTe
data flag(305) / 1/, h_f(305) /long_time /         ! 127I
data flag(306) / 4/, h_f(306) /4.9544E+14_dp/      ! 129I
data flag(307) / 4/, h_f(307) /5.1328E+06_dp/      ! 125I
data flag(308) / 4/, h_f(308) /1.1172E+06_dp/      ! 126I
data flag(309) / 4/, h_f(309) /6.9338E+05_dp/      ! 131I
data flag(310) / 4/, h_f(310) /3.6081E+05_dp/      ! 124I
data flag(311) / 1/, h_f(311) /long_time /         ! 132Xe
data flag(312) / 1/, h_f(312) /long_time /         ! 129Xe
data flag(313) / 1/, h_f(313) /long_time /         ! 131Xe
data flag(314) / 2/, h_f(314) /long_time /         ! 134Xe
data flag(315) / 3/, h_f(315) /7.5000E+28_dp/      ! 136Xe
data flag(316) / 1/, h_f(316) /long_time /         ! 130Xe
data flag(317) / 1/, h_f(317) /long_time /         ! 128Xe
data flag(318) / 2/, h_f(318) /long_time /         ! 124Xe
data flag(319) / 1/, h_f(319) /long_time /         ! 126Xe
data flag(320) / 4/, h_f(320) /3.1403E+06_dp/      ! 127Xe
data flag(321) /14/, h_f(321) /1.0230E+06_dp/      ! 131mXe
data flag(322) /14/, h_f(322) /7.6723E+05_dp/      ! 129mXe
data flag(323) / 4/, h_f(323) /4.5338E+05_dp/      ! 133Xe
data flag(324) /14/, h_f(324) /1.8991E+05_dp/      ! 133mXe
data flag(325) / 1/, h_f(325) /long_time /         ! 133Cs
data flag(326) / 4/, h_f(326) /7.2581E+13_dp/      ! 135Cs
data flag(327) / 4/, h_f(327) /9.4923E+08_dp/      ! 137Cs
data flag(328) / 4/, h_f(328) /6.5171E+07_dp/      ! 134Cs
data flag(329) / 4/, h_f(329) /1.1370E+06_dp/      ! 136Cs
data flag(330) / 4/, h_f(330) /8.3713E+05_dp/      ! 131Cs
data flag(331) / 4/, h_f(331) /5.5987E+05_dp/      ! 132Cs
data flag(332) / 4/, h_f(332) /1.1542E+05_dp/      ! 129Cs
data flag(333) / 1/, h_f(333) /long_time /         ! 138Ba
data flag(334) / 1/, h_f(334) /long_time /         ! 137Ba
data flag(335) / 1/, h_f(335) /long_time /         ! 136Ba
data flag(336) / 1/, h_f(336) /long_time /         ! 135Ba
data flag(337) / 1/, h_f(337) /long_time /         ! 134Ba
data flag(338) / 3/, h_f(338) /3.7900E+28_dp/      ! 130Ba
data flag(339) / 2/, h_f(339) /long_time /         ! 132Ba
data flag(340) / 4/, h_f(340) /3.3296E+08_dp/      ! 133Ba
data flag(341) / 4/, h_f(341) /1.1018E+06_dp/      ! 140Ba
data flag(342) / 4/, h_f(342) /9.9533E+05_dp/      ! 131Ba
data flag(343) / 4/, h_f(343) /2.0995E+05_dp/      ! 128Ba
data flag(344) /14/, h_f(344) /1.4004E+05_dp/      ! 133mBa
data flag(345) /14/, h_f(345) /1.0120E+05_dp/      ! 135mBa
data flag(346) / 1/, h_f(346) /long_time /         ! 139La
data flag(347) / 3/, h_f(347) /3.2188E+18_dp/      ! 138La
data flag(348) / 4/, h_f(348) /1.8934E+12_dp/      ! 137La
data flag(349) / 4/, h_f(349) /1.4503E+05_dp/      ! 140La
data flag(350) / 1/, h_f(350) /long_time /         ! 140Ce
data flag(351) / 2/, h_f(351) /long_time /         ! 142Ce
data flag(352) / 2/, h_f(352) /long_time /         ! 138Ce
data flag(353) / 2/, h_f(353) /long_time /         ! 136Ce
data flag(354) / 4/, h_f(354) /2.4616E+07_dp/      ! 144Ce
data flag(355) / 4/, h_f(355) /1.1892E+07_dp/      ! 139Ce
data flag(356) / 4/, h_f(356) /2.8087E+06_dp/      ! 141Ce
data flag(357) / 4/, h_f(357) /2.7302E+05_dp/      ! 134Ce
data flag(358) /14/, h_f(358) /1.2384E+05_dp/      ! 137mCe
data flag(359) / 4/, h_f(359) /1.1894E+05_dp/      ! 143Ce
data flag(360) / 1/, h_f(360) /long_time /         ! 141Pr
data flag(361) / 4/, h_f(361) /1.1724E+06_dp/      ! 143Pr
data flag(362) / 1/, h_f(362) /long_time /         ! 142Nd
data flag(363) / 3/, h_f(363) /7.2265E+22_dp/      ! 144Nd
data flag(364) / 1/, h_f(364) /long_time /         ! 146Nd
data flag(365) / 1/, h_f(365) /long_time /         ! 143Nd
data flag(366) / 1/, h_f(366) /long_time /         ! 145Nd
data flag(367) / 1/, h_f(367) /long_time /         ! 148Nd
data flag(368) / 3/, h_f(368) /2.1143E+26_dp/      ! 150Nd
data flag(369) / 4/, h_f(369) /9.4867E+05_dp/      ! 147Nd
data flag(370) / 4/, h_f(370) /2.9117E+05_dp/      ! 140Nd
data flag(371) / 4/, h_f(371) /5.5856E+08_dp/      ! 145Pm
data flag(372) / 4/, h_f(372) /1.7451E+08_dp/      ! 146Pm
data flag(373) / 4/, h_f(373) /8.2786E+07_dp/      ! 147Pm
data flag(374) / 4/, h_f(374) /3.1363E+07_dp/      ! 144Pm
data flag(375) / 4/, h_f(375) /2.2896E+07_dp/      ! 143Pm
data flag(376) /14/, h_f(376) /3.5675E+06_dp/      ! 148mPm
data flag(377) / 4/, h_f(377) /4.6380E+05_dp/      ! 148Pm
data flag(378) / 4/, h_f(378) /1.9109E+05_dp/      ! 149Pm
data flag(379) / 4/, h_f(379) /1.0224E+05_dp/      ! 151Pm
data flag(380) / 1/, h_f(380) /long_time /         ! 152Sm
data flag(381) / 1/, h_f(381) /long_time /         ! 154Sm
data flag(382) / 3/, h_f(382) /3.3640E+18_dp/      ! 147Sm
data flag(383) / 2/, h_f(383) /long_time /         ! 149Sm
data flag(384) / 3/, h_f(384) /2.2090E+23_dp/      ! 148Sm
data flag(385) / 1/, h_f(385) /long_time /         ! 150Sm
data flag(386) / 1/, h_f(386) /long_time /         ! 144Sm
data flag(387) / 3/, h_f(387) /2.1459E+15_dp/      ! 146Sm
data flag(388) / 4/, h_f(388) /2.8401E+09_dp/      ! 151Sm
data flag(389) / 4/, h_f(389) /2.9376E+07_dp/      ! 145Sm
data flag(390) / 4/, h_f(390) /1.6662E+05_dp/      ! 153Sm
data flag(391) / 1/, h_f(391) /long_time /         ! 153Eu
data flag(392) / 3/, h_f(392) /1.5780E+26_dp/      ! 151Eu
data flag(393) / 4/, h_f(393) /1.1645E+09_dp/      ! 150Eu
data flag(394) / 4/, h_f(394) /4.2719E+08_dp/      ! 152Eu
data flag(395) / 4/, h_f(395) /2.7142E+08_dp/      ! 154Eu
data flag(396) / 4/, h_f(396) /1.4999E+08_dp/      ! 155Eu
data flag(397) / 4/, h_f(397) /8.0438E+06_dp/      ! 149Eu
data flag(398) / 4/, h_f(398) /4.7088E+06_dp/      ! 148Eu
data flag(399) / 4/, h_f(399) /2.0822E+06_dp/      ! 147Eu
data flag(400) / 4/, h_f(400) /1.3124E+06_dp/      ! 156Eu
data flag(401) / 4/, h_f(401) /5.1235E+05_dp/      ! 145Eu
data flag(402) / 4/, h_f(402) /3.9830E+05_dp/      ! 146Eu
data flag(403) / 1/, h_f(403) /long_time /         ! 158Gd
data flag(404) / 2/, h_f(404) /long_time /         ! 160Gd
data flag(405) / 1/, h_f(405) /long_time /         ! 156Gd
data flag(406) / 1/, h_f(406) /long_time /         ! 157Gd
data flag(407) / 1/, h_f(407) /long_time /         ! 155Gd
data flag(408) / 1/, h_f(408) /long_time /         ! 154Gd
data flag(409) / 3/, h_f(409) /3.4081E+21_dp/      ! 152Gd
data flag(410) / 4/, h_f(410) /5.6487E+13_dp/      ! 150Gd
data flag(411) / 4/, h_f(411) /2.2374E+09_dp/      ! 148Gd
data flag(412) / 4/, h_f(412) /2.0771E+07_dp/      ! 153Gd
data flag(413) / 4/, h_f(413) /1.0705E+07_dp/      ! 151Gd
data flag(414) / 4/, h_f(414) /4.1705E+06_dp/      ! 146Gd
data flag(415) / 4/, h_f(415) /8.0179E+05_dp/      ! 149Gd
data flag(416) / 4/, h_f(416) /1.3702E+05_dp/      ! 147Gd
data flag(417) / 1/, h_f(417) /long_time /         ! 159Tb
data flag(418) / 4/, h_f(418) /5.6802E+09_dp/      ! 158Tb
data flag(419) / 4/, h_f(419) /2.2405E+09_dp/      ! 157Tb
data flag(420) / 4/, h_f(420) /6.2467E+06_dp/      ! 160Tb
data flag(421) / 4/, h_f(421) /5.9530E+05_dp/      ! 161Tb
data flag(422) / 4/, h_f(422) /4.6224E+05_dp/      ! 156Tb
data flag(423) / 4/, h_f(423) /4.5965E+05_dp/      ! 155Tb
data flag(424) / 4/, h_f(424) /2.0218E+05_dp/      ! 153Tb
data flag(425) /14/, h_f(425) /8.7840E+04_dp/      ! 156mTb
data flag(426) / 1/, h_f(426) /long_time /         ! 164Dy
data flag(427) / 1/, h_f(427) /long_time /         ! 162Dy
data flag(428) / 1/, h_f(428) /long_time /         ! 163Dy
data flag(429) / 1/, h_f(429) /long_time /         ! 161Dy
data flag(430) / 1/, h_f(430) /long_time /         ! 160Dy
data flag(431) / 1/, h_f(431) /long_time /         ! 158Dy
data flag(432) / 1/, h_f(432) /long_time /         ! 156Dy
data flag(433) / 4/, h_f(433) /9.4671E+13_dp/      ! 154Dy
data flag(434) / 4/, h_f(434) /1.2476E+07_dp/      ! 159Dy
data flag(435) / 4/, h_f(435) /2.9376E+05_dp/      ! 166Dy
data flag(436) / 1/, h_f(436) /long_time /         ! 165Ho
data flag(437) / 4/, h_f(437) /1.4422E+11_dp/      ! 163Ho
data flag(438) /14/, h_f(438) /3.5754E+10_dp/      ! 166mHo
data flag(439) / 4/, h_f(439) /9.6566E+04_dp/      ! 166Ho
data flag(440) / 1/, h_f(440) /long_time /         ! 166Er
data flag(441) / 1/, h_f(441) /long_time /         ! 168Er
data flag(442) / 1/, h_f(442) /long_time /         ! 167Er
data flag(443) / 1/, h_f(443) /long_time /         ! 170Er
data flag(444) / 1/, h_f(444) /long_time /         ! 164Er
data flag(445) / 1/, h_f(445) /long_time /         ! 162Er
data flag(446) / 4/, h_f(446) /8.1147E+05_dp/      ! 169Er
data flag(447) / 4/, h_f(447) /1.7748E+05_dp/      ! 172Er
data flag(448) / 4/, h_f(448) /1.0289E+05_dp/      ! 160Er
data flag(449) / 1/, h_f(449) /long_time /         ! 169Tm
data flag(450) / 4/, h_f(450) /6.0589E+07_dp/      ! 171Tm
data flag(451) / 4/, h_f(451) /1.1111E+07_dp/      ! 170Tm
data flag(452) / 4/, h_f(452) /8.0438E+06_dp/      ! 168Tm
data flag(453) / 4/, h_f(453) /7.9920E+05_dp/      ! 167Tm
data flag(454) / 4/, h_f(454) /2.2896E+05_dp/      ! 172Tm
data flag(455) / 4/, h_f(455) /1.0822E+05_dp/      ! 165Tm
data flag(456) / 1/, h_f(456) /long_time /         ! 174Yb
data flag(457) / 1/, h_f(457) /long_time /         ! 172Yb
data flag(458) / 1/, h_f(458) /long_time /         ! 173Yb
data flag(459) / 1/, h_f(459) /long_time /         ! 171Yb
data flag(460) / 2/, h_f(460) /long_time /         ! 176Yb
data flag(461) / 1/, h_f(461) /long_time /         ! 170Yb
data flag(462) / 1/, h_f(462) /long_time /         ! 168Yb
data flag(463) / 4/, h_f(463) /2.7664E+06_dp/      ! 169Yb
data flag(464) / 4/, h_f(464) /3.6158E+05_dp/      ! 175Yb
data flag(465) / 4/, h_f(465) /2.0412E+05_dp/      ! 166Yb
data flag(466) / 1/, h_f(466) /long_time /         ! 175Lu
data flag(467) / 3/, h_f(467) /1.1865E+18_dp/      ! 176Lu
data flag(468) / 4/, h_f(468) /1.0445E+08_dp/      ! 174Lu
data flag(469) / 4/, h_f(469) /4.3233E+07_dp/      ! 173Lu
data flag(470) /34/, h_f(470) /1.3862E+07_dp/      ! 177m3Lu
data flag(471) /14/, h_f(471) /1.2269E+07_dp/      ! 174mLu
data flag(472) / 4/, h_f(472) /7.1194E+05_dp/      ! 171Lu
data flag(473) / 4/, h_f(473) /5.7888E+05_dp/      ! 172Lu
data flag(474) / 4/, h_f(474) /5.7430E+05_dp/      ! 177Lu
data flag(475) / 4/, h_f(475) /1.7384E+05_dp/      ! 170Lu
data flag(476) / 4/, h_f(476) /1.2262E+05_dp/      ! 169Lu
data flag(477) / 1/, h_f(477) /long_time /         ! 180Hf
data flag(478) / 1/, h_f(478) /long_time /         ! 178Hf
data flag(479) / 1/, h_f(479) /long_time /         ! 177Hf
data flag(480) / 1/, h_f(480) /long_time /         ! 179Hf
data flag(481) / 1/, h_f(481) /long_time /         ! 176Hf
data flag(482) / 3/, h_f(482) /6.3114E+22_dp/      ! 174Hf
data flag(483) / 4/, h_f(483) /2.8086E+14_dp/      ! 182Hf
data flag(484) /24/, h_f(484) /9.7826E+08_dp/      ! 178m2Hf
data flag(485) / 4/, h_f(485) /5.9011E+07_dp/      ! 172Hf
data flag(486) / 4/, h_f(486) /6.0480E+06_dp/      ! 175Hf
data flag(487) / 4/, h_f(487) /3.6625E+06_dp/      ! 181Hf
data flag(488) /24/, h_f(488) /2.1643E+06_dp/      ! 179m2Hf
data flag(489) / 1/, h_f(489) /long_time /         ! 181Ta
data flag(490) /12/, h_f(490) /long_time /         ! 180mTa
data flag(491) / 4/, h_f(491) /5.7434E+07_dp/      ! 179Ta
data flag(492) / 4/, h_f(492) /9.9135E+06_dp/      ! 182Ta
data flag(493) / 4/, h_f(493) /4.4064E+05_dp/      ! 183Ta
data flag(494) / 4/, h_f(494) /2.0362E+05_dp/      ! 177Ta
data flag(495) / 2/, h_f(495) /long_time /         ! 184W
data flag(496) / 2/, h_f(496) /long_time /         ! 186W
data flag(497) / 2/, h_f(497) /long_time /         ! 182W
data flag(498) / 2/, h_f(498) /long_time /         ! 183W
data flag(499) / 3/, h_f(499) /5.6802E+25_dp/      ! 180W
data flag(500) / 4/, h_f(500) /1.0472E+07_dp/      ! 181W
data flag(501) / 4/, h_f(501) /6.4886E+06_dp/      ! 185W
data flag(502) / 4/, h_f(502) /6.0290E+06_dp/      ! 188W
data flag(503) / 4/, h_f(503) /1.8662E+06_dp/      ! 178W
data flag(504) / 3/, h_f(504) /1.3664E+18_dp/      ! 187Re
data flag(505) / 1/, h_f(505) /long_time /         ! 185Re
data flag(506) /14/, h_f(506) /6.3114E+12_dp/      ! 186mRe
data flag(507) /14/, h_f(507) /1.4602E+07_dp/      ! 184mRe
data flag(508) / 4/, h_f(508) /6.0480E+06_dp/      ! 183Re
data flag(509) / 4/, h_f(509) /3.2126E+05_dp/      ! 186Re
data flag(510) / 4/, h_f(510) /2.3040E+05_dp/      ! 182Re
data flag(511) / 4/, h_f(511) /8.7480E+04_dp/      ! 189Re
data flag(512) / 1/, h_f(512) /long_time /         ! 192Os
data flag(513) / 1/, h_f(513) /long_time /         ! 190Os
data flag(514) / 1/, h_f(514) /long_time /         ! 189Os
data flag(515) / 1/, h_f(515) /long_time /         ! 188Os
data flag(516) / 1/, h_f(516) /long_time /         ! 187Os
data flag(517) / 3/, h_f(517) /6.3114E+22_dp/      ! 186Os
data flag(518) / 2/, h_f(518) /long_time /         ! 184Os
data flag(519) / 4/, h_f(519) /1.8934E+08_dp/      ! 194Os
data flag(520) / 4/, h_f(520) /8.0309E+06_dp/      ! 185Os
data flag(521) / 4/, h_f(521) /1.2951E+06_dp/      ! 191Os
data flag(522) / 4/, h_f(522) /1.0739E+05_dp/      ! 193Os
data flag(523) / 1/, h_f(523) /long_time /         ! 193Ir
data flag(524) / 1/, h_f(524) /long_time /         ! 191Ir
data flag(525) /24/, h_f(525) /7.6052E+09_dp/      ! 192m2Ir
data flag(526) /24/, h_f(526) /1.4774E+07_dp/      ! 194m2Ir
data flag(527) / 4/, h_f(527) /6.3787E+06_dp/      ! 192Ir
data flag(528) / 4/, h_f(528) /1.1405E+06_dp/      ! 189Ir
data flag(529) / 4/, h_f(529) /1.0178E+06_dp/      ! 190Ir
data flag(530) /14/, h_f(530) /9.0979E+05_dp/      ! 193mIr
data flag(531) / 4/, h_f(531) /1.4940E+05_dp/      ! 188Ir
data flag(532) / 1/, h_f(532) /long_time /         ! 195Pt
data flag(533) / 1/, h_f(533) /long_time /         ! 194Pt
data flag(534) / 1/, h_f(534) /long_time /         ! 196Pt
data flag(535) / 1/, h_f(535) /long_time /         ! 198Pt
data flag(536) / 1/, h_f(536) /long_time /         ! 192Pt
data flag(537) / 3/, h_f(537) /2.0512E+19_dp/      ! 190Pt
data flag(538) / 4/, h_f(538) /1.5778E+09_dp/      ! 193Pt
data flag(539) / 4/, h_f(539) /8.8128E+05_dp/      ! 188Pt
data flag(540) /14/, h_f(540) /3.7411E+05_dp/      ! 193mPt
data flag(541) /14/, h_f(541) /3.4646E+05_dp/      ! 195mPt
data flag(542) / 4/, h_f(542) /2.4451E+05_dp/      ! 191Pt
data flag(543) / 4/, h_f(543) /1.5840E+05_dp/      ! 202Pt
data flag(544) / 1/, h_f(544) /long_time /         ! 197Au
data flag(545) / 4/, h_f(545) /1.6079E+07_dp/      ! 195Au
data flag(546) / 4/, h_f(546) /5.3282E+05_dp/      ! 196Au
data flag(547) / 4/, h_f(547) /2.7121E+05_dp/      ! 199Au
data flag(548) / 4/, h_f(548) /2.3283E+05_dp/      ! 198Au
data flag(549) /24/, h_f(549) /1.9630E+05_dp/      ! 198m2Au
data flag(550) / 4/, h_f(550) /1.3687E+05_dp/      ! 194Au
data flag(551) / 1/, h_f(551) /long_time /         ! 202Hg
data flag(552) / 1/, h_f(552) /long_time /         ! 200Hg
data flag(553) / 1/, h_f(553) /long_time /         ! 199Hg
data flag(554) / 1/, h_f(554) /long_time /         ! 201Hg
data flag(555) / 1/, h_f(555) /long_time /         ! 198Hg
data flag(556) / 1/, h_f(556) /long_time /         ! 204Hg
data flag(557) / 1/, h_f(557) /long_time /         ! 196Hg
data flag(558) / 4/, h_f(558) /1.3885E+10_dp/      ! 194Hg
data flag(559) / 4/, h_f(559) /4.0257E+06_dp/      ! 203Hg
data flag(560) / 4/, h_f(560) /2.3378E+05_dp/      ! 197Hg
data flag(561) /14/, h_f(561) /1.4976E+05_dp/      ! 195mHg
data flag(562) / 1/, h_f(562) /long_time /         ! 205Tl
data flag(563) / 1/, h_f(563) /long_time /         ! 203Tl
data flag(564) / 4/, h_f(564) /1.1938E+08_dp/      ! 204Tl
data flag(565) / 4/, h_f(565) /1.0636E+06_dp/      ! 202Tl
data flag(566) / 4/, h_f(566) /2.6284E+05_dp/      ! 201Tl
data flag(567) / 4/, h_f(567) /9.3960E+04_dp/      ! 200Tl
data flag(568) / 1/, h_f(568) /long_time /         ! 208Pb
data flag(569) / 1/, h_f(569) /long_time /         ! 206Pb
data flag(570) / 1/, h_f(570) /long_time /         ! 207Pb
data flag(571) / 2/, h_f(571) /long_time /         ! 204Pb
data flag(572) / 4/, h_f(572) /5.4593E+14_dp/      ! 205Pb
data flag(573) / 4/, h_f(573) /1.6567E+12_dp/      ! 202Pb
data flag(574) / 4/, h_f(574) /7.0056E+08_dp/      ! 210Pb
data flag(575) / 4/, h_f(575) /1.8691E+05_dp/      ! 203Pb
data flag(576) / 3/, h_f(576) /6.2798E+26_dp/      ! 209Bi
data flag(577) /14/, h_f(577) /9.5933E+13_dp/      ! 210mBi
data flag(578) / 4/, h_f(578) /1.1613E+13_dp/      ! 208Bi
data flag(579) / 4/, h_f(579) /9.9562E+08_dp/      ! 207Bi
data flag(580) / 4/, h_f(580) /1.3228E+06_dp/      ! 205Bi
data flag(581) / 4/, h_f(581) /5.3940E+05_dp/      ! 206Bi
data flag(582) / 4/, h_f(582) /4.3304E+05_dp/      ! 210Bi
data flag(583) / 4/, h_f(583) /3.2188E+09_dp/      ! 209Po
data flag(584) / 4/, h_f(584) /9.1452E+07_dp/      ! 208Po
data flag(585) / 4/, h_f(585) /1.1956E+07_dp/      ! 210Po
data flag(586) / 4/, h_f(586) /7.6032E+05_dp/      ! 206Po
data flag(587) / 4/, h_f(587) /3.3035E+05_dp/      ! 222Rn
data flag(588) / 4/, h_f(588) /5.0491E+10_dp/      ! 226Ra
data flag(589) / 4/, h_f(589) /1.8145E+08_dp/      ! 228Ra
data flag(590) / 4/, h_f(590) /1.2874E+06_dp/      ! 225Ra
data flag(591) / 4/, h_f(591) /9.8755E+05_dp/      ! 223Ra
data flag(592) / 4/, h_f(592) /3.1622E+05_dp/      ! 224Ra
data flag(593) / 4/, h_f(593) /6.8706E+08_dp/      ! 227Ac
data flag(594) / 4/, h_f(594) /8.5709E+05_dp/      ! 225Ac
data flag(595) / 4/, h_f(595) /1.0573E+05_dp/      ! 226Ac
data flag(596) / 3/, h_f(596) /4.4180E+17_dp/      ! 232Th
data flag(597) / 4/, h_f(597) /2.3794E+12_dp/      ! 230Th
data flag(598) / 4/, h_f(598) /2.5031E+11_dp/      ! 229Th
data flag(599) / 4/, h_f(599) /6.0324E+07_dp/      ! 228Th
data flag(600) / 4/, h_f(600) /2.0822E+06_dp/      ! 234Th
data flag(601) / 4/, h_f(601) /1.6140E+06_dp/      ! 227Th
data flag(602) / 4/, h_f(602) /9.1872E+04_dp/      ! 231Th
data flag(603) / 4/, h_f(603) /1.0338E+12_dp/      ! 231Pa
data flag(604) / 4/, h_f(604) /2.3306E+06_dp/      ! 233Pa
data flag(605) / 4/, h_f(605) /1.5034E+06_dp/      ! 230Pa
data flag(606) / 4/, h_f(606) /1.2960E+05_dp/      ! 229Pa
data flag(607) / 4/, h_f(607) /1.1405E+05_dp/      ! 232Pa
data flag(608) / 3/, h_f(608) /1.4100E+17_dp/      ! 238U
data flag(609) / 3/, h_f(609) /2.2216E+16_dp/      ! 235U
data flag(610) / 4/, h_f(610) /7.7472E+12_dp/      ! 234U
data flag(611) / 4/, h_f(611) /7.3906E+14_dp/      ! 236U
data flag(612) / 4/, h_f(612) /5.0239E+12_dp/      ! 233U
data flag(613) / 4/, h_f(613) /2.1743E+09_dp/      ! 232U
data flag(614) / 4/, h_f(614) /1.7479E+06_dp/      ! 230U
data flag(615) / 4/, h_f(615) /5.8337E+05_dp/      ! 237U
data flag(616) / 4/, h_f(616) /3.6288E+05_dp/      ! 231U
data flag(617) / 4/, h_f(617) /6.7658E+13_dp/      ! 237Np
data flag(618) / 4/, h_f(618) /4.8282E+12_dp/      ! 236Np
data flag(619) / 4/, h_f(619) /3.4223E+07_dp/      ! 235Np
data flag(620) / 4/, h_f(620) /3.8016E+05_dp/      ! 234Np
data flag(621) / 4/, h_f(621) /2.0356E+05_dp/      ! 239Np
data flag(622) / 4/, h_f(622) /1.8291E+05_dp/      ! 238Np
data flag(623) / 3/, h_f(623) /2.5246E+15_dp/      ! 244Pu
data flag(624) / 4/, h_f(624) /1.1834E+13_dp/      ! 242Pu
data flag(625) / 4/, h_f(625) /7.6084E+11_dp/      ! 239Pu
data flag(626) / 4/, h_f(626) /2.0704E+11_dp/      ! 240Pu
data flag(627) / 4/, h_f(627) /2.7675E+09_dp/      ! 238Pu
data flag(628) / 4/, h_f(628) /4.5095E+08_dp/      ! 241Pu
data flag(629) / 4/, h_f(629) /9.0190E+07_dp/      ! 236Pu
data flag(630) / 4/, h_f(630) /3.9433E+06_dp/      ! 237Pu
data flag(631) / 4/, h_f(631) /9.3658E+05_dp/      ! 246Pu
data flag(632) / 4/, h_f(632) /1.9613E+05_dp/      ! 247Pu
data flag(633) / 4/, h_f(633) /2.3257E+11_dp/      ! 243Am
data flag(634) / 4/, h_f(634) /1.3652E+10_dp/      ! 241Am
data flag(635) /14/, h_f(635) /4.4495E+09_dp/      ! 242mAm
data flag(636) / 4/, h_f(636) /1.8288E+05_dp/      ! 240Am
data flag(637) / 4/, h_f(637) /4.9229E+14_dp/      ! 247Cm
data flag(638) / 4/, h_f(638) /1.0982E+13_dp/      ! 248Cm
data flag(639) / 4/, h_f(639) /2.6580E+11_dp/      ! 245Cm
data flag(640) / 4/, h_f(640) /2.6192E+11_dp/      ! 250Cm
data flag(641) / 4/, h_f(641) /1.4851E+11_dp/      ! 246Cm
data flag(642) / 4/, h_f(642) /9.1831E+08_dp/      ! 243Cm
data flag(643) / 4/, h_f(643) /5.7118E+08_dp/      ! 244Cm
data flag(644) / 4/, h_f(644) /1.4066E+07_dp/      ! 242Cm
data flag(645) / 4/, h_f(645) /2.8339E+06_dp/      ! 241Cm
data flag(646) / 4/, h_f(646) /2.3328E+06_dp/      ! 240Cm
data flag(647) / 4/, h_f(647) /1.7280E+05_dp/      ! 252Cm
data flag(648) / 4/, h_f(648) /4.3549E+10_dp/      ! 247Bk
data flag(649) / 4/, h_f(649) /2.8401E+08_dp/      ! 248Bk
data flag(650) / 4/, h_f(650) /2.8512E+07_dp/      ! 249Bk
data flag(651) / 4/, h_f(651) /4.2768E+05_dp/      ! 245Bk
data flag(652) / 4/, h_f(652) /1.5552E+05_dp/      ! 246Bk
data flag(653) / 4/, h_f(653) /2.8401E+10_dp/      ! 251Cf
data flag(654) / 4/, h_f(654) /1.1076E+10_dp/      ! 249Cf
data flag(655) / 4/, h_f(655) /4.1276E+08_dp/      ! 250Cf
data flag(656) / 4/, h_f(656) /8.3468E+07_dp/      ! 252Cf
data flag(657) / 4/, h_f(657) /2.8858E+07_dp/      ! 248Cf
data flag(658) / 4/, h_f(658) /5.2272E+06_dp/      ! 254Cf
data flag(659) / 4/, h_f(659) /1.5388E+06_dp/      ! 253Cf
data flag(660) / 4/, h_f(660) /1.2852E+05_dp/      ! 246Cf
data flag(661) / 4/, h_f(661) /4.0755E+07_dp/      ! 252Es
data flag(662) / 4/, h_f(662) /2.3820E+07_dp/      ! 254Es
data flag(663) / 4/, h_f(663) /3.4387E+06_dp/      ! 255Es
data flag(664) / 4/, h_f(664) /1.7686E+06_dp/      ! 253Es
data flag(665) / 4/, h_f(665) /6.6528E+05_dp/      ! 257Es
data flag(666) /14/, h_f(666) /1.4148E+05_dp/      ! 254mEs
data flag(667) / 4/, h_f(667) /1.1880E+05_dp/      ! 251Es
data flag(668) / 4/, h_f(668) /8.6832E+06_dp/      ! 257Fm
data flag(669) / 4/, h_f(669) /2.5920E+05_dp/      ! 253Fm
data flag(670) / 4/, h_f(670) /9.1404E+04_dp/      ! 252Fm
data flag(671) / 4/, h_f(671) /4.4496E+06_dp/      ! 258Md
data flag(672) / 4/, h_f(672) /2.4019E+06_dp/      ! 260Md
data flag(673) / 4/, h_f(673) /1.1088E+05_dp/      ! 268Db

data spin(  0) /spin_undef/, nat_iso_ab(  0) /       0.0_dp/
data spin(  1) /0.5_dp  /  , nat_iso_ab(  1) / 99.988500_dp/ ! 1H
data spin(  2) /1.0_dp  /  , nat_iso_ab(  2) /  0.011500_dp/ ! 2H
data spin(  3) /0.5_dp  /  , nat_iso_ab(  3) /  0.000000_dp/ ! 3H
data spin(  4) /0.0_dp  /  , nat_iso_ab(  4) / 99.999866_dp/ ! 4He
data spin(  5) /0.5_dp  /  , nat_iso_ab(  5) /  0.000134_dp/ ! 3He
data spin(  6) /1.5_dp  /  , nat_iso_ab(  6) / 92.410000_dp/ ! 7Li
data spin(  7) /1.0_dp  /  , nat_iso_ab(  7) /  7.590000_dp/ ! 6Li
data spin(  8) /1.5_dp  /  , nat_iso_ab(  8) /100.000000_dp/ ! 9Be
data spin(  9) /0.0_dp  /  , nat_iso_ab(  9) /  0.000000_dp/ ! 10Be
data spin( 10) /1.5_dp  /  , nat_iso_ab( 10) /  0.000000_dp/ ! 7Be
data spin( 11) /1.5_dp  /  , nat_iso_ab( 11) / 80.100000_dp/ ! 11B
data spin( 12) /3.0_dp  /  , nat_iso_ab( 12) / 19.900000_dp/ ! 10B
data spin( 13) /0.0_dp  /  , nat_iso_ab( 13) / 98.930000_dp/ ! 12C
data spin( 14) /0.5_dp  /  , nat_iso_ab( 14) /  1.070000_dp/ ! 13C
data spin( 15) /0.0_dp  /  , nat_iso_ab( 15) /  0.000000_dp/ ! 14C
data spin( 16) /1.0_dp  /  , nat_iso_ab( 16) / 99.636000_dp/ ! 14N
data spin( 17) /0.5_dp  /  , nat_iso_ab( 17) /  0.364000_dp/ ! 15N
data spin( 18) /0.0_dp  /  , nat_iso_ab( 18) / 99.757000_dp/ ! 16O
data spin( 19) /0.0_dp  /  , nat_iso_ab( 19) /  0.205000_dp/ ! 18O
data spin( 20) /2.5_dp  /  , nat_iso_ab( 20) /  0.038000_dp/ ! 17O
data spin( 21) /0.5_dp  /  , nat_iso_ab( 21) /100.000000_dp/ ! 19F
data spin( 22) /0.0_dp  /  , nat_iso_ab( 22) / 90.480000_dp/ ! 20Ne
data spin( 23) /0.0_dp  /  , nat_iso_ab( 23) /  9.250000_dp/ ! 22Ne
data spin( 24) /1.5_dp  /  , nat_iso_ab( 24) /  0.270000_dp/ ! 21Ne
data spin( 25) /1.5_dp  /  , nat_iso_ab( 25) /100.000000_dp/ ! 23Na
data spin( 26) /3.0_dp  /  , nat_iso_ab( 26) /  0.000000_dp/ ! 22Na
data spin( 27) /0.0_dp  /  , nat_iso_ab( 27) / 78.990000_dp/ ! 24Mg
data spin( 28) /0.0_dp  /  , nat_iso_ab( 28) / 11.010000_dp/ ! 26Mg
data spin( 29) /2.5_dp  /  , nat_iso_ab( 29) / 10.000000_dp/ ! 25Mg
data spin( 30) /2.5_dp  /  , nat_iso_ab( 30) /100.000000_dp/ ! 27Al
data spin( 31) /5.0_dp  /  , nat_iso_ab( 31) /  0.000000_dp/ ! 26Al
data spin( 32) /0.0_dp  /  , nat_iso_ab( 32) / 92.223000_dp/ ! 28Si
data spin( 33) /0.5_dp  /  , nat_iso_ab( 33) /  4.685000_dp/ ! 29Si
data spin( 34) /0.0_dp  /  , nat_iso_ab( 34) /  3.092000_dp/ ! 30Si
data spin( 35) /0.0_dp  /  , nat_iso_ab( 35) /  0.000000_dp/ ! 32Si
data spin( 36) /0.5_dp  /  , nat_iso_ab( 36) /100.000000_dp/ ! 31P
data spin( 37) /0.5_dp  /  , nat_iso_ab( 37) /  0.000000_dp/ ! 33P
data spin( 38) /1.0_dp  /  , nat_iso_ab( 38) /  0.000000_dp/ ! 32P
data spin( 39) /0.0_dp  /  , nat_iso_ab( 39) / 94.990000_dp/ ! 32S
data spin( 40) /0.0_dp  /  , nat_iso_ab( 40) /  4.250000_dp/ ! 34S
data spin( 41) /1.5_dp  /  , nat_iso_ab( 41) /  0.750000_dp/ ! 33S
data spin( 42) /0.0_dp  /  , nat_iso_ab( 42) /  0.010000_dp/ ! 36S
data spin( 43) /1.5_dp  /  , nat_iso_ab( 43) /  0.000000_dp/ ! 35S
data spin( 44) /1.5_dp  /  , nat_iso_ab( 44) / 75.760000_dp/ ! 35Cl
data spin( 45) /1.5_dp  /  , nat_iso_ab( 45) / 24.240000_dp/ ! 37Cl
data spin( 46) /2.0_dp  /  , nat_iso_ab( 46) /  0.000000_dp/ ! 36Cl
data spin( 47) /0.0_dp  /  , nat_iso_ab( 47) / 99.603500_dp/ ! 40Ar
data spin( 48) /0.0_dp  /  , nat_iso_ab( 48) /  0.333600_dp/ ! 36Ar
data spin( 49) /0.0_dp  /  , nat_iso_ab( 49) /  0.062900_dp/ ! 38Ar
data spin( 50) /3.5_dp  /  , nat_iso_ab( 50) /  0.000000_dp/ ! 39Ar
data spin( 51) /0.0_dp  /  , nat_iso_ab( 51) /  0.000000_dp/ ! 42Ar
data spin( 52) /1.5_dp  /  , nat_iso_ab( 52) /  0.000000_dp/ ! 37Ar
data spin( 53) /1.5_dp  /  , nat_iso_ab( 53) / 93.258100_dp/ ! 39K
data spin( 54) /1.5_dp  /  , nat_iso_ab( 54) /  6.730200_dp/ ! 41K
data spin( 55) /4.0_dp  /  , nat_iso_ab( 55) /  0.011700_dp/ ! 40K
data spin( 56) /0.0_dp  /  , nat_iso_ab( 56) / 96.940000_dp/ ! 40Ca
data spin( 57) /0.0_dp  /  , nat_iso_ab( 57) /  2.090000_dp/ ! 44Ca
data spin( 58) /0.0_dp  /  , nat_iso_ab( 58) /  0.647000_dp/ ! 42Ca
data spin( 59) /0.0_dp  /  , nat_iso_ab( 59) /  0.187000_dp/ ! 48Ca
data spin( 60) /3.5_dp  /  , nat_iso_ab( 60) /  0.135000_dp/ ! 43Ca
data spin( 61) /0.0_dp  /  , nat_iso_ab( 61) /  0.004000_dp/ ! 46Ca
data spin( 62) /3.5_dp  /  , nat_iso_ab( 62) /  0.000000_dp/ ! 41Ca
data spin( 63) /3.5_dp  /  , nat_iso_ab( 63) /  0.000000_dp/ ! 45Ca
data spin( 64) /3.5_dp  /  , nat_iso_ab( 64) /  0.000000_dp/ ! 47Ca
data spin( 65) /3.5_dp  /  , nat_iso_ab( 65) /100.000000_dp/ ! 45Sc
data spin( 66) /4.0_dp  /  , nat_iso_ab( 66) /  0.000000_dp/ ! 46Sc
data spin( 67) /3.5_dp  /  , nat_iso_ab( 67) /  0.000000_dp/ ! 47Sc
data spin( 68) /6.0_dp  /  , nat_iso_ab( 68) /  0.000000_dp/ ! 44m3Sc
data spin( 69) /6.0_dp  /  , nat_iso_ab( 69) /  0.000000_dp/ ! 48Sc
data spin( 70) /0.0_dp  /  , nat_iso_ab( 70) / 73.720000_dp/ ! 48Ti
data spin( 71) /0.0_dp  /  , nat_iso_ab( 71) /  8.250000_dp/ ! 46Ti
data spin( 72) /2.5_dp  /  , nat_iso_ab( 72) /  7.440000_dp/ ! 47Ti
data spin( 73) /3.5_dp  /  , nat_iso_ab( 73) /  5.410000_dp/ ! 49Ti
data spin( 74) /0.0_dp  /  , nat_iso_ab( 74) /  5.180000_dp/ ! 50Ti
data spin( 75) /0.0_dp  /  , nat_iso_ab( 75) /  0.000000_dp/ ! 44Ti
data spin( 76) /3.5_dp  /  , nat_iso_ab( 76) / 99.750000_dp/ ! 51V
data spin( 77) /6.0_dp  /  , nat_iso_ab( 77) /  0.250000_dp/ ! 50V
data spin( 78) /3.5_dp  /  , nat_iso_ab( 78) /  0.000000_dp/ ! 49V
data spin( 79) /4.0_dp  /  , nat_iso_ab( 79) /  0.000000_dp/ ! 48V
data spin( 80) /0.0_dp  /  , nat_iso_ab( 80) / 83.789000_dp/ ! 52Cr
data spin( 81) /1.5_dp  /  , nat_iso_ab( 81) /  9.501000_dp/ ! 53Cr
data spin( 82) /0.0_dp  /  , nat_iso_ab( 82) /  4.345000_dp/ ! 50Cr
data spin( 83) /0.0_dp  /  , nat_iso_ab( 83) /  2.365000_dp/ ! 54Cr
data spin( 84) /3.5_dp  /  , nat_iso_ab( 84) /  0.000000_dp/ ! 51Cr
data spin( 85) /2.5_dp  /  , nat_iso_ab( 85) /100.000000_dp/ ! 55Mn
data spin( 86) /3.5_dp  /  , nat_iso_ab( 86) /  0.000000_dp/ ! 53Mn
data spin( 87) /3.0_dp  /  , nat_iso_ab( 87) /  0.000000_dp/ ! 54Mn
data spin( 88) /6.0_dp  /  , nat_iso_ab( 88) /  0.000000_dp/ ! 52Mn
data spin( 89) /0.0_dp  /  , nat_iso_ab( 89) / 91.754000_dp/ ! 56Fe
data spin( 90) /0.0_dp  /  , nat_iso_ab( 90) /  5.845000_dp/ ! 54Fe
data spin( 91) /0.5_dp  /  , nat_iso_ab( 91) /  2.119000_dp/ ! 57Fe
data spin( 92) /0.0_dp  /  , nat_iso_ab( 92) /  0.282000_dp/ ! 58Fe
data spin( 93) /0.0_dp  /  , nat_iso_ab( 93) /  0.000000_dp/ ! 60Fe
data spin( 94) /1.5_dp  /  , nat_iso_ab( 94) /  0.000000_dp/ ! 55Fe
data spin( 95) /1.5_dp  /  , nat_iso_ab( 95) /  0.000000_dp/ ! 59Fe
data spin( 96) /3.5_dp  /  , nat_iso_ab( 96) /100.000000_dp/ ! 59Co
data spin( 97) /5.0_dp  /  , nat_iso_ab( 97) /  0.000000_dp/ ! 60Co
data spin( 98) /3.5_dp  /  , nat_iso_ab( 98) /  0.000000_dp/ ! 57Co
data spin( 99) /4.0_dp  /  , nat_iso_ab( 99) /  0.000000_dp/ ! 56Co
data spin(100) /2.0_dp  /  , nat_iso_ab(100) /  0.000000_dp/ ! 58Co
data spin(101) /0.0_dp  /  , nat_iso_ab(101) / 68.077000_dp/ ! 58Ni
data spin(102) /0.0_dp  /  , nat_iso_ab(102) / 26.223000_dp/ ! 60Ni
data spin(103) /0.0_dp  /  , nat_iso_ab(103) /  3.634600_dp/ ! 62Ni
data spin(104) /1.5_dp  /  , nat_iso_ab(104) /  1.139900_dp/ ! 61Ni
data spin(105) /0.0_dp  /  , nat_iso_ab(105) /  0.925500_dp/ ! 64Ni
data spin(106) /1.5_dp  /  , nat_iso_ab(106) /  0.000000_dp/ ! 59Ni
data spin(107) /0.5_dp  /  , nat_iso_ab(107) /  0.000000_dp/ ! 63Ni
data spin(108) /0.0_dp  /  , nat_iso_ab(108) /  0.000000_dp/ ! 56Ni
data spin(109) /0.0_dp  /  , nat_iso_ab(109) /  0.000000_dp/ ! 66Ni
data spin(110) /1.5_dp  /  , nat_iso_ab(110) /  0.000000_dp/ ! 57Ni
data spin(111) /1.5_dp  /  , nat_iso_ab(111) / 69.150000_dp/ ! 63Cu
data spin(112) /1.5_dp  /  , nat_iso_ab(112) / 30.850000_dp/ ! 65Cu
data spin(113) /1.5_dp  /  , nat_iso_ab(113) /  0.000000_dp/ ! 67Cu
data spin(114) /0.0_dp  /  , nat_iso_ab(114) / 49.170000_dp/ ! 64Zn
data spin(115) /0.0_dp  /  , nat_iso_ab(115) / 27.730000_dp/ ! 66Zn
data spin(116) /0.0_dp  /  , nat_iso_ab(116) / 18.450000_dp/ ! 68Zn
data spin(117) /2.5_dp  /  , nat_iso_ab(117) /  4.040000_dp/ ! 67Zn
data spin(118) /0.0_dp  /  , nat_iso_ab(118) /  0.610000_dp/ ! 70Zn
data spin(119) /2.5_dp  /  , nat_iso_ab(119) /  0.000000_dp/ ! 65Zn
data spin(120) /0.0_dp  /  , nat_iso_ab(120) /  0.000000_dp/ ! 72Zn
data spin(121) /1.5_dp  /  , nat_iso_ab(121) / 60.108000_dp/ ! 69Ga
data spin(122) /1.5_dp  /  , nat_iso_ab(122) / 39.892000_dp/ ! 71Ga
data spin(123) /1.5_dp  /  , nat_iso_ab(123) /  0.000000_dp/ ! 67Ga
data spin(124) /0.0_dp  /  , nat_iso_ab(124) / 36.500000_dp/ ! 74Ge
data spin(125) /0.0_dp  /  , nat_iso_ab(125) / 27.450000_dp/ ! 72Ge
data spin(126) /0.0_dp  /  , nat_iso_ab(126) / 20.570000_dp/ ! 70Ge
data spin(127) /4.5_dp  /  , nat_iso_ab(127) /  7.750000_dp/ ! 73Ge
data spin(128) /0.0_dp  /  , nat_iso_ab(128) /  7.730000_dp/ ! 76Ge
data spin(129) /0.0_dp  /  , nat_iso_ab(129) /  0.000000_dp/ ! 68Ge
data spin(130) /0.5_dp  /  , nat_iso_ab(130) /  0.000000_dp/ ! 71Ge
data spin(131) /2.5_dp  /  , nat_iso_ab(131) /  0.000000_dp/ ! 69Ge
data spin(132) /1.5_dp  /  , nat_iso_ab(132) /100.000000_dp/ ! 75As
data spin(133) /1.5_dp  /  , nat_iso_ab(133) /  0.000000_dp/ ! 73As
data spin(134) /2.0_dp  /  , nat_iso_ab(134) /  0.000000_dp/ ! 74As
data spin(135) /2.5_dp  /  , nat_iso_ab(135) /  0.000000_dp/ ! 71As
data spin(136) /1.5_dp  /  , nat_iso_ab(136) /  0.000000_dp/ ! 77As
data spin(137) /2.0_dp  /  , nat_iso_ab(137) /  0.000000_dp/ ! 72As
data spin(138) /2.0_dp  /  , nat_iso_ab(138) /  0.000000_dp/ ! 76As
data spin(139) /0.0_dp  /  , nat_iso_ab(139) / 49.610000_dp/ ! 80Se
data spin(140) /0.0_dp  /  , nat_iso_ab(140) / 23.770000_dp/ ! 78Se
data spin(141) /0.0_dp  /  , nat_iso_ab(141) /  9.370000_dp/ ! 76Se
data spin(142) /0.0_dp  /  , nat_iso_ab(142) /  8.730000_dp/ ! 82Se
data spin(143) /0.5_dp  /  , nat_iso_ab(143) /  7.630000_dp/ ! 77Se
data spin(144) /0.0_dp  /  , nat_iso_ab(144) /  0.890000_dp/ ! 74Se
data spin(145) /3.5_dp  /  , nat_iso_ab(145) /  0.000000_dp/ ! 79Se
data spin(146) /2.5_dp  /  , nat_iso_ab(146) /  0.000000_dp/ ! 75Se
data spin(147) /0.0_dp  /  , nat_iso_ab(147) /  0.000000_dp/ ! 72Se
data spin(148) /1.5_dp  /  , nat_iso_ab(148) / 50.690000_dp/ ! 79Br
data spin(149) /1.5_dp  /  , nat_iso_ab(149) / 49.310000_dp/ ! 81Br
data spin(150) /1.5_dp  /  , nat_iso_ab(150) /  0.000000_dp/ ! 77Br
data spin(151) /5.0_dp  /  , nat_iso_ab(151) /  0.000000_dp/ ! 82Br
data spin(152) /0.0_dp  /  , nat_iso_ab(152) / 56.987000_dp/ ! 84Kr
data spin(153) /0.0_dp  /  , nat_iso_ab(153) / 17.279000_dp/ ! 86Kr
data spin(154) /0.0_dp  /  , nat_iso_ab(154) / 11.593000_dp/ ! 82Kr
data spin(155) /4.5_dp  /  , nat_iso_ab(155) / 11.500000_dp/ ! 83Kr
data spin(156) /0.0_dp  /  , nat_iso_ab(156) /  2.286000_dp/ ! 80Kr
data spin(157) /0.0_dp  /  , nat_iso_ab(157) /  0.355000_dp/ ! 78Kr
data spin(158) /3.5_dp  /  , nat_iso_ab(158) /  0.000000_dp/ ! 81Kr
data spin(159) /4.5_dp  /  , nat_iso_ab(159) /  0.000000_dp/ ! 85Kr
data spin(160) /0.5_dp  /  , nat_iso_ab(160) /  0.000000_dp/ ! 79Kr
data spin(161) /2.5_dp  /  , nat_iso_ab(161) / 72.170000_dp/ ! 85Rb
data spin(162) /1.5_dp  /  , nat_iso_ab(162) / 27.830000_dp/ ! 87Rb
data spin(163) /2.5_dp  /  , nat_iso_ab(163) /  0.000000_dp/ ! 83Rb
data spin(164) /2.0_dp  /  , nat_iso_ab(164) /  0.000000_dp/ ! 84Rb
data spin(165) /2.0_dp  /  , nat_iso_ab(165) /  0.000000_dp/ ! 86Rb
data spin(166) /0.0_dp  /  , nat_iso_ab(166) / 82.580000_dp/ ! 88Sr
data spin(167) /0.0_dp  /  , nat_iso_ab(167) /  9.860000_dp/ ! 86Sr
data spin(168) /4.5_dp  /  , nat_iso_ab(168) /  7.000000_dp/ ! 87Sr
data spin(169) /0.0_dp  /  , nat_iso_ab(169) /  0.560000_dp/ ! 84Sr
data spin(170) /0.0_dp  /  , nat_iso_ab(170) /  0.000000_dp/ ! 90Sr
data spin(171) /4.5_dp  /  , nat_iso_ab(171) /  0.000000_dp/ ! 85Sr
data spin(172) /2.5_dp  /  , nat_iso_ab(172) /  0.000000_dp/ ! 89Sr
data spin(173) /0.0_dp  /  , nat_iso_ab(173) /  0.000000_dp/ ! 82Sr
data spin(174) /3.5_dp  /  , nat_iso_ab(174) /  0.000000_dp/ ! 83Sr
data spin(175) /0.5_dp  /  , nat_iso_ab(175) /100.000000_dp/ ! 89Y
data spin(176) /4.0_dp  /  , nat_iso_ab(176) /  0.000000_dp/ ! 88Y
data spin(177) /0.5_dp  /  , nat_iso_ab(177) /  0.000000_dp/ ! 91Y
data spin(178) /0.5_dp  /  , nat_iso_ab(178) /  0.000000_dp/ ! 87Y
data spin(179) /2.0_dp  /  , nat_iso_ab(179) /  0.000000_dp/ ! 90Y
data spin(180) /0.0_dp  /  , nat_iso_ab(180) / 51.450000_dp/ ! 90Zr
data spin(181) /0.0_dp  /  , nat_iso_ab(181) / 17.380000_dp/ ! 94Zr
data spin(182) /0.0_dp  /  , nat_iso_ab(182) / 17.150000_dp/ ! 92Zr
data spin(183) /2.5_dp  /  , nat_iso_ab(183) / 11.220000_dp/ ! 91Zr
data spin(184) /0.0_dp  /  , nat_iso_ab(184) /  2.800000_dp/ ! 96Zr
data spin(185) /2.5_dp  /  , nat_iso_ab(185) /  0.000000_dp/ ! 93Zr
data spin(186) /0.0_dp  /  , nat_iso_ab(186) /  0.000000_dp/ ! 88Zr
data spin(187) /2.5_dp  /  , nat_iso_ab(187) /  0.000000_dp/ ! 95Zr
data spin(188) /4.5_dp  /  , nat_iso_ab(188) /  0.000000_dp/ ! 89Zr
data spin(189) /4.5_dp  /  , nat_iso_ab(189) /100.000000_dp/ ! 93Nb
data spin(190) /7.0_dp  /  , nat_iso_ab(190) /  0.000000_dp/ ! 92Nb
data spin(191) /6.0_dp  /  , nat_iso_ab(191) /  0.000000_dp/ ! 94Nb
data spin(192) /4.5_dp  /  , nat_iso_ab(192) /  0.000000_dp/ ! 91Nb
data spin(193) /0.5_dp  /  , nat_iso_ab(193) /  0.000000_dp/ ! 93mNb
data spin(194) /4.5_dp  /  , nat_iso_ab(194) /  0.000000_dp/ ! 95Nb
data spin(195) /0.5_dp  /  , nat_iso_ab(195) /  0.000000_dp/ ! 95mNb
data spin(196) /0.0_dp  /  , nat_iso_ab(196) / 24.390000_dp/ ! 98Mo
data spin(197) /0.0_dp  /  , nat_iso_ab(197) / 16.670000_dp/ ! 96Mo
data spin(198) /2.5_dp  /  , nat_iso_ab(198) / 15.840000_dp/ ! 95Mo
data spin(199) /0.0_dp  /  , nat_iso_ab(199) / 14.530000_dp/ ! 92Mo
data spin(200) /0.0_dp  /  , nat_iso_ab(200) /  9.820000_dp/ ! 100Mo
data spin(201) /2.5_dp  /  , nat_iso_ab(201) /  9.600000_dp/ ! 97Mo
data spin(202) /0.0_dp  /  , nat_iso_ab(202) /  9.150000_dp/ ! 94Mo
data spin(203) /2.5_dp  /  , nat_iso_ab(203) /  0.000000_dp/ ! 93Mo
data spin(204) /0.5_dp  /  , nat_iso_ab(204) /  0.000000_dp/ ! 99Mo
data spin(205) /4.5_dp  /  , nat_iso_ab(205) /  0.000000_dp/ ! 97Tc
data spin(206) /6.0_dp  /  , nat_iso_ab(206) /  0.000000_dp/ ! 98Tc
data spin(207) /4.5_dp  /  , nat_iso_ab(207) /  0.000000_dp/ ! 99Tc
data spin(208) /0.5_dp  /  , nat_iso_ab(208) /  0.000000_dp/ ! 97mTc
data spin(209) /0.5_dp  /  , nat_iso_ab(209) /  0.000000_dp/ ! 95mTc
data spin(210) /7.0_dp  /  , nat_iso_ab(210) /  0.000000_dp/ ! 96Tc
data spin(211) /0.0_dp  /  , nat_iso_ab(211) / 31.550000_dp/ ! 102Ru
data spin(212) /0.0_dp  /  , nat_iso_ab(212) / 18.620000_dp/ ! 104Ru
data spin(213) /2.5_dp  /  , nat_iso_ab(213) / 17.060000_dp/ ! 101Ru
data spin(214) /2.5_dp  /  , nat_iso_ab(214) / 12.760000_dp/ ! 99Ru
data spin(215) /0.0_dp  /  , nat_iso_ab(215) / 12.600000_dp/ ! 100Ru
data spin(216) /0.0_dp  /  , nat_iso_ab(216) /  5.540000_dp/ ! 96Ru
data spin(217) /0.0_dp  /  , nat_iso_ab(217) /  1.870000_dp/ ! 98Ru
data spin(218) /0.0_dp  /  , nat_iso_ab(218) /  0.000000_dp/ ! 106Ru
data spin(219) /1.5_dp  /  , nat_iso_ab(219) /  0.000000_dp/ ! 103Ru
data spin(220) /2.5_dp  /  , nat_iso_ab(220) /  0.000000_dp/ ! 97Ru
data spin(221) /0.5_dp  /  , nat_iso_ab(221) /100.000000_dp/ ! 103Rh
data spin(222) /6.0_dp  /  , nat_iso_ab(222) /  0.000000_dp/ ! 102mRh
data spin(223) /0.5_dp  /  , nat_iso_ab(223) /  0.000000_dp/ ! 101Rh
data spin(224) /1.0_dp  /  , nat_iso_ab(224) /  0.000000_dp/ ! 102Rh
data spin(225) /0.5_dp  /  , nat_iso_ab(225) /  0.000000_dp/ ! 99Rh
data spin(226) /4.5_dp  /  , nat_iso_ab(226) /  0.000000_dp/ ! 101mRh
data spin(227) /3.5_dp  /  , nat_iso_ab(227) /  0.000000_dp/ ! 105Rh
data spin(228) /0.0_dp  /  , nat_iso_ab(228) / 27.330000_dp/ ! 106Pd
data spin(229) /0.0_dp  /  , nat_iso_ab(229) / 26.460000_dp/ ! 108Pd
data spin(230) /2.5_dp  /  , nat_iso_ab(230) / 22.330000_dp/ ! 105Pd
data spin(231) /0.0_dp  /  , nat_iso_ab(231) / 11.720000_dp/ ! 110Pd
data spin(232) /0.0_dp  /  , nat_iso_ab(232) / 11.140000_dp/ ! 104Pd
data spin(233) /0.0_dp  /  , nat_iso_ab(233) /  1.020000_dp/ ! 102Pd
data spin(234) /2.5_dp  /  , nat_iso_ab(234) /  0.000000_dp/ ! 107Pd
data spin(235) /2.5_dp  /  , nat_iso_ab(235) /  0.000000_dp/ ! 103Pd
data spin(236) /0.0_dp  /  , nat_iso_ab(236) /  0.000000_dp/ ! 100Pd
data spin(237) /0.5_dp  /  , nat_iso_ab(237) / 51.839000_dp/ ! 107Ag
data spin(238) /0.5_dp  /  , nat_iso_ab(238) / 48.161000_dp/ ! 109Ag
data spin(239) /6.0_dp  /  , nat_iso_ab(239) /  0.000000_dp/ ! 108mAg
data spin(240) /6.0_dp  /  , nat_iso_ab(240) /  0.000000_dp/ ! 110m2Ag
data spin(241) /0.5_dp  /  , nat_iso_ab(241) /  0.000000_dp/ ! 105Ag
data spin(242) /6.0_dp  /  , nat_iso_ab(242) /  0.000000_dp/ ! 106mAg
data spin(243) /0.5_dp  /  , nat_iso_ab(243) /  0.000000_dp/ ! 111Ag
data spin(244) /0.0_dp  /  , nat_iso_ab(244) / 28.730000_dp/ ! 114Cd
data spin(245) /0.0_dp  /  , nat_iso_ab(245) / 24.130000_dp/ ! 112Cd
data spin(246) /0.5_dp  /  , nat_iso_ab(246) / 12.800000_dp/ ! 111Cd
data spin(247) /0.0_dp  /  , nat_iso_ab(247) / 12.490000_dp/ ! 110Cd
data spin(248) /0.5_dp  /  , nat_iso_ab(248) / 12.220000_dp/ ! 113Cd
data spin(249) /0.0_dp  /  , nat_iso_ab(249) /  7.490000_dp/ ! 116Cd
data spin(250) /0.0_dp  /  , nat_iso_ab(250) /  1.250000_dp/ ! 106Cd
data spin(251) /0.0_dp  /  , nat_iso_ab(251) /  0.890000_dp/ ! 108Cd
data spin(252) /5.5_dp  /  , nat_iso_ab(252) /  0.000000_dp/ ! 113mCd
data spin(253) /2.5_dp  /  , nat_iso_ab(253) /  0.000000_dp/ ! 109Cd
data spin(254) /5.5_dp  /  , nat_iso_ab(254) /  0.000000_dp/ ! 115mCd
data spin(255) /0.5_dp  /  , nat_iso_ab(255) /  0.000000_dp/ ! 115Cd
data spin(256) /4.5_dp  /  , nat_iso_ab(256) / 95.710000_dp/ ! 115In
data spin(257) /4.5_dp  /  , nat_iso_ab(257) /  4.290000_dp/ ! 113In
data spin(258) /5.0_dp  /  , nat_iso_ab(258) /  0.000000_dp/ ! 114mIn
data spin(259) /4.5_dp  /  , nat_iso_ab(259) /  0.000000_dp/ ! 111In
data spin(260) /0.0_dp  /  , nat_iso_ab(260) / 32.580000_dp/ ! 120Sn
data spin(261) /0.0_dp  /  , nat_iso_ab(261) / 24.220000_dp/ ! 118Sn
data spin(262) /0.0_dp  /  , nat_iso_ab(262) / 14.540000_dp/ ! 116Sn
data spin(263) /0.5_dp  /  , nat_iso_ab(263) /  8.590000_dp/ ! 119Sn
data spin(264) /0.5_dp  /  , nat_iso_ab(264) /  7.680000_dp/ ! 117Sn
data spin(265) /0.0_dp  /  , nat_iso_ab(265) /  5.790000_dp/ ! 124Sn
data spin(266) /0.0_dp  /  , nat_iso_ab(266) /  4.630000_dp/ ! 122Sn
data spin(267) /0.0_dp  /  , nat_iso_ab(267) /  0.970000_dp/ ! 112Sn
data spin(268) /0.0_dp  /  , nat_iso_ab(268) /  0.660000_dp/ ! 114Sn
data spin(269) /0.5_dp  /  , nat_iso_ab(269) /  0.340000_dp/ ! 115Sn
data spin(270) /0.0_dp  /  , nat_iso_ab(270) /  0.000000_dp/ ! 126Sn
data spin(271) /5.5_dp  /  , nat_iso_ab(271) /  0.000000_dp/ ! 121mSn
data spin(272) /5.5_dp  /  , nat_iso_ab(272) /  0.000000_dp/ ! 119mSn
data spin(273) /5.5_dp  /  , nat_iso_ab(273) /  0.000000_dp/ ! 123Sn
data spin(274) /0.5_dp  /  , nat_iso_ab(274) /  0.000000_dp/ ! 113Sn
data spin(275) /5.5_dp  /  , nat_iso_ab(275) /  0.000000_dp/ ! 117mSn
data spin(276) /5.5_dp  /  , nat_iso_ab(276) /  0.000000_dp/ ! 125Sn
data spin(277) /1.5_dp  /  , nat_iso_ab(277) /  0.000000_dp/ ! 121Sn
data spin(278) /2.5_dp  /  , nat_iso_ab(278) / 57.210000_dp/ ! 121Sb
data spin(279) /3.5_dp  /  , nat_iso_ab(279) / 42.790000_dp/ ! 123Sb
data spin(280) /3.5_dp  /  , nat_iso_ab(280) /  0.000000_dp/ ! 125Sb
data spin(281) /3.0_dp  /  , nat_iso_ab(281) /  0.000000_dp/ ! 124Sb
data spin(282) /8.0_dp  /  , nat_iso_ab(282) /  0.000000_dp/ ! 126Sb
data spin(283) /8.0_dp  /  , nat_iso_ab(283) /  0.000000_dp/ ! 120mSb
data spin(284) /3.5_dp  /  , nat_iso_ab(284) /  0.000000_dp/ ! 127Sb
data spin(285) /2.0_dp  /  , nat_iso_ab(285) /  0.000000_dp/ ! 122Sb
data spin(286) /2.5_dp  /  , nat_iso_ab(286) /  0.000000_dp/ ! 119Sb
data spin(287) /0.0_dp  /  , nat_iso_ab(287) / 34.080000_dp/ ! 130Te
data spin(288) /0.0_dp  /  , nat_iso_ab(288) / 18.840000_dp/ ! 126Te
data spin(289) /0.5_dp  /  , nat_iso_ab(289) /  7.070000_dp/ ! 125Te
data spin(290) /0.0_dp  /  , nat_iso_ab(290) /  4.740000_dp/ ! 124Te
data spin(291) /0.0_dp  /  , nat_iso_ab(291) /  2.550000_dp/ ! 122Te
data spin(292) /0.5_dp  /  , nat_iso_ab(292) /  0.890000_dp/ ! 123Te
data spin(293) /0.0_dp  /  , nat_iso_ab(293) /  0.090000_dp/ ! 120Te
data spin(294) /0.0_dp  /  , nat_iso_ab(294) / 31.740000_dp/ ! 128Te
data spin(295) /5.5_dp  /  , nat_iso_ab(295) /  0.000000_dp/ ! 121mTe
data spin(296) /5.5_dp  /  , nat_iso_ab(296) /  0.000000_dp/ ! 123mTe
data spin(297) /5.5_dp  /  , nat_iso_ab(297) /  0.000000_dp/ ! 127mTe
data spin(298) /5.5_dp  /  , nat_iso_ab(298) /  0.000000_dp/ ! 125mTe
data spin(299) /5.5_dp  /  , nat_iso_ab(299) /  0.000000_dp/ ! 129mTe
data spin(300) /0.5_dp  /  , nat_iso_ab(300) /  0.000000_dp/ ! 121Te
data spin(301) /0.0_dp  /  , nat_iso_ab(301) /  0.000000_dp/ ! 118Te
data spin(302) /5.5_dp  /  , nat_iso_ab(302) /  0.000000_dp/ ! 119mTe
data spin(303) /0.0_dp  /  , nat_iso_ab(303) /  0.000000_dp/ ! 132Te
data spin(304) /5.5_dp  /  , nat_iso_ab(304) /  0.000000_dp/ ! 131mTe
data spin(305) /2.5_dp  /  , nat_iso_ab(305) /100.000000_dp/ ! 127I
data spin(306) /3.5_dp  /  , nat_iso_ab(306) /  0.000000_dp/ ! 129I
data spin(307) /2.5_dp  /  , nat_iso_ab(307) /  0.000000_dp/ ! 125I
data spin(308) /2.0_dp  /  , nat_iso_ab(308) /  0.000000_dp/ ! 126I
data spin(309) /3.5_dp  /  , nat_iso_ab(309) /  0.000000_dp/ ! 131I
data spin(310) /2.0_dp  /  , nat_iso_ab(310) /  0.000000_dp/ ! 124I
data spin(311) /0.0_dp  /  , nat_iso_ab(311) / 26.908600_dp/ ! 132Xe
data spin(312) /0.5_dp  /  , nat_iso_ab(312) / 26.400600_dp/ ! 129Xe
data spin(313) /1.5_dp  /  , nat_iso_ab(313) / 21.232400_dp/ ! 131Xe
data spin(314) /0.0_dp  /  , nat_iso_ab(314) / 10.435700_dp/ ! 134Xe
data spin(315) /0.0_dp  /  , nat_iso_ab(315) /  8.857300_dp/ ! 136Xe
data spin(316) /0.0_dp  /  , nat_iso_ab(316) /  4.071000_dp/ ! 130Xe
data spin(317) /0.0_dp  /  , nat_iso_ab(317) /  1.910200_dp/ ! 128Xe
data spin(318) /0.0_dp  /  , nat_iso_ab(318) /  0.095200_dp/ ! 124Xe
data spin(319) /0.0_dp  /  , nat_iso_ab(319) /  0.089000_dp/ ! 126Xe
data spin(320) /0.5_dp  /  , nat_iso_ab(320) /  0.000000_dp/ ! 127Xe
data spin(321) /5.5_dp  /  , nat_iso_ab(321) /  0.000000_dp/ ! 131mXe
data spin(322) /5.5_dp  /  , nat_iso_ab(322) /  0.000000_dp/ ! 129mXe
data spin(323) /1.5_dp  /  , nat_iso_ab(323) /  0.000000_dp/ ! 133Xe
data spin(324) /5.5_dp  /  , nat_iso_ab(324) /  0.000000_dp/ ! 133mXe
data spin(325) /3.5_dp  /  , nat_iso_ab(325) /100.000000_dp/ ! 133Cs
data spin(326) /3.5_dp  /  , nat_iso_ab(326) /  0.000000_dp/ ! 135Cs
data spin(327) /3.5_dp  /  , nat_iso_ab(327) /  0.000000_dp/ ! 137Cs
data spin(328) /4.0_dp  /  , nat_iso_ab(328) /  0.000000_dp/ ! 134Cs
data spin(329) /5.0_dp  /  , nat_iso_ab(329) /  0.000000_dp/ ! 136Cs
data spin(330) /2.5_dp  /  , nat_iso_ab(330) /  0.000000_dp/ ! 131Cs
data spin(331) /2.0_dp  /  , nat_iso_ab(331) /  0.000000_dp/ ! 132Cs
data spin(332) /0.5_dp  /  , nat_iso_ab(332) /  0.000000_dp/ ! 129Cs
data spin(333) /0.0_dp  /  , nat_iso_ab(333) / 71.698000_dp/ ! 138Ba
data spin(334) /1.5_dp  /  , nat_iso_ab(334) / 11.232000_dp/ ! 137Ba
data spin(335) /0.0_dp  /  , nat_iso_ab(335) /  7.854000_dp/ ! 136Ba
data spin(336) /1.5_dp  /  , nat_iso_ab(336) /  6.592000_dp/ ! 135Ba
data spin(337) /0.0_dp  /  , nat_iso_ab(337) /  2.417000_dp/ ! 134Ba
data spin(338) /0.0_dp  /  , nat_iso_ab(338) /  0.106000_dp/ ! 130Ba
data spin(339) /0.0_dp  /  , nat_iso_ab(339) /  0.101000_dp/ ! 132Ba
data spin(340) /0.5_dp  /  , nat_iso_ab(340) /  0.000000_dp/ ! 133Ba
data spin(341) /0.0_dp  /  , nat_iso_ab(341) /  0.000000_dp/ ! 140Ba
data spin(342) /0.5_dp  /  , nat_iso_ab(342) /  0.000000_dp/ ! 131Ba
data spin(343) /0.0_dp  /  , nat_iso_ab(343) /  0.000000_dp/ ! 128Ba
data spin(344) /5.5_dp  /  , nat_iso_ab(344) /  0.000000_dp/ ! 133mBa
data spin(345) /5.5_dp  /  , nat_iso_ab(345) /  0.000000_dp/ ! 135mBa
data spin(346) /3.5_dp  /  , nat_iso_ab(346) / 99.911190_dp/ ! 139La
data spin(347) /5.0_dp  /  , nat_iso_ab(347) /  0.088810_dp/ ! 138La
data spin(348) /3.5_dp  /  , nat_iso_ab(348) /  0.000000_dp/ ! 137La
data spin(349) /3.0_dp  /  , nat_iso_ab(349) /  0.000000_dp/ ! 140La
data spin(350) /0.0_dp  /  , nat_iso_ab(350) / 88.450000_dp/ ! 140Ce
data spin(351) /0.0_dp  /  , nat_iso_ab(351) / 11.114000_dp/ ! 142Ce
data spin(352) /0.0_dp  /  , nat_iso_ab(352) /  0.251000_dp/ ! 138Ce
data spin(353) /0.0_dp  /  , nat_iso_ab(353) /  0.185000_dp/ ! 136Ce
data spin(354) /0.0_dp  /  , nat_iso_ab(354) /  0.000000_dp/ ! 144Ce
data spin(355) /1.5_dp  /  , nat_iso_ab(355) /  0.000000_dp/ ! 139Ce
data spin(356) /3.5_dp  /  , nat_iso_ab(356) /  0.000000_dp/ ! 141Ce
data spin(357) /0.0_dp  /  , nat_iso_ab(357) /  0.000000_dp/ ! 134Ce
data spin(358) /5.5_dp  /  , nat_iso_ab(358) /  0.000000_dp/ ! 137mCe
data spin(359) /1.5_dp  /  , nat_iso_ab(359) /  0.000000_dp/ ! 143Ce
data spin(360) /2.5_dp  /  , nat_iso_ab(360) /100.000000_dp/ ! 141Pr
data spin(361) /3.5_dp  /  , nat_iso_ab(361) /  0.000000_dp/ ! 143Pr
data spin(362) /0.0_dp  /  , nat_iso_ab(362) / 27.152000_dp/ ! 142Nd
data spin(363) /0.0_dp  /  , nat_iso_ab(363) / 23.798000_dp/ ! 144Nd
data spin(364) /0.0_dp  /  , nat_iso_ab(364) / 17.189000_dp/ ! 146Nd
data spin(365) /3.5_dp  /  , nat_iso_ab(365) / 12.174000_dp/ ! 143Nd
data spin(366) /3.5_dp  /  , nat_iso_ab(366) /  8.293000_dp/ ! 145Nd
data spin(367) /0.0_dp  /  , nat_iso_ab(367) /  5.756000_dp/ ! 148Nd
data spin(368) /0.0_dp  /  , nat_iso_ab(368) /  5.638000_dp/ ! 150Nd
data spin(369) /2.5_dp  /  , nat_iso_ab(369) /  0.000000_dp/ ! 147Nd
data spin(370) /0.0_dp  /  , nat_iso_ab(370) /  0.000000_dp/ ! 140Nd
data spin(371) /2.5_dp  /  , nat_iso_ab(371) /  0.000000_dp/ ! 145Pm
data spin(372) /3.0_dp  /  , nat_iso_ab(372) /  0.000000_dp/ ! 146Pm
data spin(373) /3.5_dp  /  , nat_iso_ab(373) /  0.000000_dp/ ! 147Pm
data spin(374) /5.0_dp  /  , nat_iso_ab(374) /  0.000000_dp/ ! 144Pm
data spin(375) /2.5_dp  /  , nat_iso_ab(375) /  0.000000_dp/ ! 143Pm
data spin(376) /5.0_dp  /  , nat_iso_ab(376) /  0.000000_dp/ ! 148mPm
data spin(377) /1.0_dp  /  , nat_iso_ab(377) /  0.000000_dp/ ! 148Pm
data spin(378) /3.5_dp  /  , nat_iso_ab(378) /  0.000000_dp/ ! 149Pm
data spin(379) /2.5_dp  /  , nat_iso_ab(379) /  0.000000_dp/ ! 151Pm
data spin(380) /0.0_dp  /  , nat_iso_ab(380) / 26.750000_dp/ ! 152Sm
data spin(381) /0.0_dp  /  , nat_iso_ab(381) / 22.750000_dp/ ! 154Sm
data spin(382) /3.5_dp  /  , nat_iso_ab(382) / 14.990000_dp/ ! 147Sm
data spin(383) /3.5_dp  /  , nat_iso_ab(383) / 13.820000_dp/ ! 149Sm
data spin(384) /0.0_dp  /  , nat_iso_ab(384) / 11.240000_dp/ ! 148Sm
data spin(385) /0.0_dp  /  , nat_iso_ab(385) /  7.380000_dp/ ! 150Sm
data spin(386) /0.0_dp  /  , nat_iso_ab(386) /  3.070000_dp/ ! 144Sm
data spin(387) /0.0_dp  /  , nat_iso_ab(387) /  0.000000_dp/ ! 146Sm
data spin(388) /2.5_dp  /  , nat_iso_ab(388) /  0.000000_dp/ ! 151Sm
data spin(389) /3.5_dp  /  , nat_iso_ab(389) /  0.000000_dp/ ! 145Sm
data spin(390) /1.5_dp  /  , nat_iso_ab(390) /  0.000000_dp/ ! 153Sm
data spin(391) /2.5_dp  /  , nat_iso_ab(391) / 52.190000_dp/ ! 153Eu
data spin(392) /2.5_dp  /  , nat_iso_ab(392) / 47.810000_dp/ ! 151Eu
data spin(393) /5.0_dp  /  , nat_iso_ab(393) /  0.000000_dp/ ! 150Eu
data spin(394) /3.0_dp  /  , nat_iso_ab(394) /  0.000000_dp/ ! 152Eu
data spin(395) /3.0_dp  /  , nat_iso_ab(395) /  0.000000_dp/ ! 154Eu
data spin(396) /2.5_dp  /  , nat_iso_ab(396) /  0.000000_dp/ ! 155Eu
data spin(397) /2.5_dp  /  , nat_iso_ab(397) /  0.000000_dp/ ! 149Eu
data spin(398) /5.0_dp  /  , nat_iso_ab(398) /  0.000000_dp/ ! 148Eu
data spin(399) /2.5_dp  /  , nat_iso_ab(399) /  0.000000_dp/ ! 147Eu
data spin(400) /0.0_dp  /  , nat_iso_ab(400) /  0.000000_dp/ ! 156Eu
data spin(401) /2.5_dp  /  , nat_iso_ab(401) /  0.000000_dp/ ! 145Eu
data spin(402) /4.0_dp  /  , nat_iso_ab(402) /  0.000000_dp/ ! 146Eu
data spin(403) /0.0_dp  /  , nat_iso_ab(403) / 24.840000_dp/ ! 158Gd
data spin(404) /0.0_dp  /  , nat_iso_ab(404) / 21.860000_dp/ ! 160Gd
data spin(405) /0.0_dp  /  , nat_iso_ab(405) / 20.470000_dp/ ! 156Gd
data spin(406) /1.5_dp  /  , nat_iso_ab(406) / 15.650000_dp/ ! 157Gd
data spin(407) /1.5_dp  /  , nat_iso_ab(407) / 14.800000_dp/ ! 155Gd
data spin(408) /0.0_dp  /  , nat_iso_ab(408) /  2.180000_dp/ ! 154Gd
data spin(409) /0.0_dp  /  , nat_iso_ab(409) /  0.200000_dp/ ! 152Gd
data spin(410) /0.0_dp  /  , nat_iso_ab(410) /  0.000000_dp/ ! 150Gd
data spin(411) /0.0_dp  /  , nat_iso_ab(411) /  0.000000_dp/ ! 148Gd
data spin(412) /1.5_dp  /  , nat_iso_ab(412) /  0.000000_dp/ ! 153Gd
data spin(413) /3.5_dp  /  , nat_iso_ab(413) /  0.000000_dp/ ! 151Gd
data spin(414) /0.0_dp  /  , nat_iso_ab(414) /  0.000000_dp/ ! 146Gd
data spin(415) /3.5_dp  /  , nat_iso_ab(415) /  0.000000_dp/ ! 149Gd
data spin(416) /3.5_dp  /  , nat_iso_ab(416) /  0.000000_dp/ ! 147Gd
data spin(417) /1.5_dp  /  , nat_iso_ab(417) /100.000000_dp/ ! 159Tb
data spin(418) /3.0_dp  /  , nat_iso_ab(418) /  0.000000_dp/ ! 158Tb
data spin(419) /1.5_dp  /  , nat_iso_ab(419) /  0.000000_dp/ ! 157Tb
data spin(420) /3.0_dp  /  , nat_iso_ab(420) /  0.000000_dp/ ! 160Tb
data spin(421) /1.5_dp  /  , nat_iso_ab(421) /  0.000000_dp/ ! 161Tb
data spin(422) /3.0_dp  /  , nat_iso_ab(422) /  0.000000_dp/ ! 156Tb
data spin(423) /1.5_dp  /  , nat_iso_ab(423) /  0.000000_dp/ ! 155Tb
data spin(424) /2.5_dp  /  , nat_iso_ab(424) /  0.000000_dp/ ! 153Tb
data spin(425) /7.0_dp  /  , nat_iso_ab(425) /  0.000000_dp/ ! 156mTb
data spin(426) /0.0_dp  /  , nat_iso_ab(426) / 28.260000_dp/ ! 164Dy
data spin(427) /0.0_dp  /  , nat_iso_ab(427) / 25.475000_dp/ ! 162Dy
data spin(428) /2.5_dp  /  , nat_iso_ab(428) / 24.896000_dp/ ! 163Dy
data spin(429) /2.5_dp  /  , nat_iso_ab(429) / 18.889000_dp/ ! 161Dy
data spin(430) /0.0_dp  /  , nat_iso_ab(430) /  2.329000_dp/ ! 160Dy
data spin(431) /0.0_dp  /  , nat_iso_ab(431) /  0.095000_dp/ ! 158Dy
data spin(432) /0.0_dp  /  , nat_iso_ab(432) /  0.056000_dp/ ! 156Dy
data spin(433) /0.0_dp  /  , nat_iso_ab(433) /  0.000000_dp/ ! 154Dy
data spin(434) /1.5_dp  /  , nat_iso_ab(434) /  0.000000_dp/ ! 159Dy
data spin(435) /0.0_dp  /  , nat_iso_ab(435) /  0.000000_dp/ ! 166Dy
data spin(436) /3.5_dp  /  , nat_iso_ab(436) /100.000000_dp/ ! 165Ho
data spin(437) /3.5_dp  /  , nat_iso_ab(437) /  0.000000_dp/ ! 163Ho
data spin(438) /7.0_dp  /  , nat_iso_ab(438) /  0.000000_dp/ ! 166mHo
data spin(439) /0.0_dp  /  , nat_iso_ab(439) /  0.000000_dp/ ! 166Ho
data spin(440) /0.0_dp  /  , nat_iso_ab(440) / 33.503000_dp/ ! 166Er
data spin(441) /0.0_dp  /  , nat_iso_ab(441) / 26.978000_dp/ ! 168Er
data spin(442) /3.5_dp  /  , nat_iso_ab(442) / 22.869000_dp/ ! 167Er
data spin(443) /0.0_dp  /  , nat_iso_ab(443) / 14.910000_dp/ ! 170Er
data spin(444) /0.0_dp  /  , nat_iso_ab(444) /  1.601000_dp/ ! 164Er
data spin(445) /0.0_dp  /  , nat_iso_ab(445) /  0.139000_dp/ ! 162Er
data spin(446) /0.5_dp  /  , nat_iso_ab(446) /  0.000000_dp/ ! 169Er
data spin(447) /0.0_dp  /  , nat_iso_ab(447) /  0.000000_dp/ ! 172Er
data spin(448) /0.0_dp  /  , nat_iso_ab(448) /  0.000000_dp/ ! 160Er
data spin(449) /0.5_dp  /  , nat_iso_ab(449) /100.000000_dp/ ! 169Tm
data spin(450) /0.5_dp  /  , nat_iso_ab(450) /  0.000000_dp/ ! 171Tm
data spin(451) /1.0_dp  /  , nat_iso_ab(451) /  0.000000_dp/ ! 170Tm
data spin(452) /3.0_dp  /  , nat_iso_ab(452) /  0.000000_dp/ ! 168Tm
data spin(453) /0.5_dp  /  , nat_iso_ab(453) /  0.000000_dp/ ! 167Tm
data spin(454) /2.0_dp  /  , nat_iso_ab(454) /  0.000000_dp/ ! 172Tm
data spin(455) /0.5_dp  /  , nat_iso_ab(455) /  0.000000_dp/ ! 165Tm
data spin(456) /0.0_dp  /  , nat_iso_ab(456) / 32.026000_dp/ ! 174Yb
data spin(457) /0.0_dp  /  , nat_iso_ab(457) / 21.680000_dp/ ! 172Yb
data spin(458) /2.5_dp  /  , nat_iso_ab(458) / 16.103000_dp/ ! 173Yb
data spin(459) /0.5_dp  /  , nat_iso_ab(459) / 14.090000_dp/ ! 171Yb
data spin(460) /0.0_dp  /  , nat_iso_ab(460) / 12.996000_dp/ ! 176Yb
data spin(461) /0.0_dp  /  , nat_iso_ab(461) /  2.982000_dp/ ! 170Yb
data spin(462) /0.0_dp  /  , nat_iso_ab(462) /  0.123000_dp/ ! 168Yb
data spin(463) /3.5_dp  /  , nat_iso_ab(463) /  0.000000_dp/ ! 169Yb
data spin(464) /3.5_dp  /  , nat_iso_ab(464) /  0.000000_dp/ ! 175Yb
data spin(465) /0.0_dp  /  , nat_iso_ab(465) /  0.000000_dp/ ! 166Yb
data spin(466) /3.5_dp  /  , nat_iso_ab(466) / 97.401000_dp/ ! 175Lu
data spin(467) /7.0_dp  /  , nat_iso_ab(467) /  2.599000_dp/ ! 176Lu
data spin(468) /1.0_dp  /  , nat_iso_ab(468) /  0.000000_dp/ ! 174Lu
data spin(469) /3.5_dp  /  , nat_iso_ab(469) /  0.000000_dp/ ! 173Lu
data spin(470) /11.5_dp /  , nat_iso_ab(470) /  0.000000_dp/   ! 177m3Lu
data spin(471) /6.0_dp  /  , nat_iso_ab(471) /  0.000000_dp/ ! 174mLu
data spin(472) /3.5_dp  /  , nat_iso_ab(472) /  0.000000_dp/ ! 171Lu
data spin(473) /4.0_dp  /  , nat_iso_ab(473) /  0.000000_dp/ ! 172Lu
data spin(474) /3.5_dp  /  , nat_iso_ab(474) /  0.000000_dp/ ! 177Lu
data spin(475) /0.0_dp  /  , nat_iso_ab(475) /  0.000000_dp/ ! 170Lu
data spin(476) /3.5_dp  /  , nat_iso_ab(476) /  0.000000_dp/ ! 169Lu
data spin(477) /0.0_dp  /  , nat_iso_ab(477) / 35.080000_dp/ ! 180Hf
data spin(478) /0.0_dp  /  , nat_iso_ab(478) / 27.280000_dp/ ! 178Hf
data spin(479) /3.5_dp  /  , nat_iso_ab(479) / 18.600000_dp/ ! 177Hf
data spin(480) /4.5_dp  /  , nat_iso_ab(480) / 13.620000_dp/ ! 179Hf
data spin(481) /0.0_dp  /  , nat_iso_ab(481) /  5.260000_dp/ ! 176Hf
data spin(482) /0.0_dp  /  , nat_iso_ab(482) /  0.160000_dp/ ! 174Hf
data spin(483) /0.0_dp  /  , nat_iso_ab(483) /  0.000000_dp/ ! 182Hf
data spin(484) /16.0_dp /  , nat_iso_ab(484) /  0.000000_dp/   ! 178m2Hf
data spin(485) /0.0_dp  /  , nat_iso_ab(485) /  0.000000_dp/ ! 172Hf
data spin(486) /2.5_dp  /  , nat_iso_ab(486) /  0.000000_dp/ ! 175Hf
data spin(487) /0.5_dp  /  , nat_iso_ab(487) /  0.000000_dp/ ! 181Hf
data spin(488) /12.5_dp /  , nat_iso_ab(488) /  0.000000_dp/   ! 179m2Hf
data spin(489) /3.5_dp  /  , nat_iso_ab(489) / 99.987990_dp/ ! 181Ta
data spin(490) /9.0_dp  /  , nat_iso_ab(490) /  0.012010_dp/ ! 180mTa
data spin(491) /3.5_dp  /  , nat_iso_ab(491) /  0.000000_dp/ ! 179Ta
data spin(492) /3.0_dp  /  , nat_iso_ab(492) /  0.000000_dp/ ! 182Ta
data spin(493) /3.5_dp  /  , nat_iso_ab(493) /  0.000000_dp/ ! 183Ta
data spin(494) /3.5_dp  /  , nat_iso_ab(494) /  0.000000_dp/ ! 177Ta
data spin(495) /0.0_dp  /  , nat_iso_ab(495) / 30.640000_dp/ ! 184W
data spin(496) /0.0_dp  /  , nat_iso_ab(496) / 28.430000_dp/ ! 186W
data spin(497) /0.0_dp  /  , nat_iso_ab(497) / 26.500000_dp/ ! 182W
data spin(498) /0.5_dp  /  , nat_iso_ab(498) / 14.310000_dp/ ! 183W
data spin(499) /0.0_dp  /  , nat_iso_ab(499) /  0.120000_dp/ ! 180W
data spin(500) /4.5_dp  /  , nat_iso_ab(500) /  0.000000_dp/ ! 181W
data spin(501) /1.5_dp  /  , nat_iso_ab(501) /  0.000000_dp/ ! 185W
data spin(502) /0.0_dp  /  , nat_iso_ab(502) /  0.000000_dp/ ! 188W
data spin(503) /0.0_dp  /  , nat_iso_ab(503) /  0.000000_dp/ ! 178W
data spin(504) /2.5_dp  /  , nat_iso_ab(504) / 62.600000_dp/ ! 187Re
data spin(505) /2.5_dp  /  , nat_iso_ab(505) / 37.400000_dp/ ! 185Re
data spin(506) /8.0_dp  /  , nat_iso_ab(506) /  0.000000_dp/ ! 186mRe
data spin(507) /8.0_dp  /  , nat_iso_ab(507) /  0.000000_dp/ ! 184mRe
data spin(508) /2.5_dp  /  , nat_iso_ab(508) /  0.000000_dp/ ! 183Re
data spin(509) /1.0_dp  /  , nat_iso_ab(509) /  0.000000_dp/ ! 186Re
data spin(510) /7.0_dp  /  , nat_iso_ab(510) /  0.000000_dp/ ! 182Re
data spin(511) /2.5_dp  /  , nat_iso_ab(511) /  0.000000_dp/ ! 189Re
data spin(512) /0.0_dp  /  , nat_iso_ab(512) / 40.780000_dp/ ! 192Os
data spin(513) /0.0_dp  /  , nat_iso_ab(513) / 26.260000_dp/ ! 190Os
data spin(514) /1.5_dp  /  , nat_iso_ab(514) / 16.150000_dp/ ! 189Os
data spin(515) /0.0_dp  /  , nat_iso_ab(515) / 13.240000_dp/ ! 188Os
data spin(516) /0.5_dp  /  , nat_iso_ab(516) /  1.960000_dp/ ! 187Os
data spin(517) /0.0_dp  /  , nat_iso_ab(517) /  1.590000_dp/ ! 186Os
data spin(518) /0.0_dp  /  , nat_iso_ab(518) /  0.020000_dp/ ! 184Os
data spin(519) /0.0_dp  /  , nat_iso_ab(519) /  0.000000_dp/ ! 194Os
data spin(520) /0.5_dp  /  , nat_iso_ab(520) /  0.000000_dp/ ! 185Os
data spin(521) /4.5_dp  /  , nat_iso_ab(521) /  0.000000_dp/ ! 191Os
data spin(522) /1.5_dp  /  , nat_iso_ab(522) /  0.000000_dp/ ! 193Os
data spin(523) /1.5_dp  /  , nat_iso_ab(523) / 62.700000_dp/ ! 193Ir
data spin(524) /1.5_dp  /  , nat_iso_ab(524) / 37.300000_dp/ ! 191Ir
data spin(525) /11.0_dp /  , nat_iso_ab(525) /  0.000000_dp/   ! 192m2Ir
data spin(526) /10.0_dp /  , nat_iso_ab(526) /  0.000000_dp/   ! 194m2Ir
data spin(527) /4.0_dp  /  , nat_iso_ab(527) /  0.000000_dp/ ! 192Ir
data spin(528) /1.5_dp  /  , nat_iso_ab(528) /  0.000000_dp/ ! 189Ir
data spin(529) /4.0_dp  /  , nat_iso_ab(529) /  0.000000_dp/ ! 190Ir
data spin(530) /5.5_dp  /  , nat_iso_ab(530) /  0.000000_dp/ ! 193mIr
data spin(531) /1.0_dp  /  , nat_iso_ab(531) /  0.000000_dp/ ! 188Ir
data spin(532) /0.5_dp  /  , nat_iso_ab(532) / 33.780000_dp/ ! 195Pt
data spin(533) /0.0_dp  /  , nat_iso_ab(533) / 32.860000_dp/ ! 194Pt
data spin(534) /0.0_dp  /  , nat_iso_ab(534) / 25.210000_dp/ ! 196Pt
data spin(535) /0.0_dp  /  , nat_iso_ab(535) /  7.360000_dp/ ! 198Pt
data spin(536) /0.0_dp  /  , nat_iso_ab(536) /  0.782000_dp/ ! 192Pt
data spin(537) /0.0_dp  /  , nat_iso_ab(537) /  0.012000_dp/ ! 190Pt
data spin(538) /0.5_dp  /  , nat_iso_ab(538) /  0.000000_dp/ ! 193Pt
data spin(539) /0.0_dp  /  , nat_iso_ab(539) /  0.000000_dp/ ! 188Pt
data spin(540) /6.5_dp  /  , nat_iso_ab(540) /  0.000000_dp/ ! 193mPt
data spin(541) /6.5_dp  /  , nat_iso_ab(541) /  0.000000_dp/ ! 195mPt
data spin(542) /1.5_dp  /  , nat_iso_ab(542) /  0.000000_dp/ ! 191Pt
data spin(543) /0.0_dp  /  , nat_iso_ab(543) /  0.000000_dp/ ! 202Pt
data spin(544) /1.5_dp  /  , nat_iso_ab(544) /100.000000_dp/ ! 197Au
data spin(545) /1.5_dp  /  , nat_iso_ab(545) /  0.000000_dp/ ! 195Au
data spin(546) /2.0_dp  /  , nat_iso_ab(546) /  0.000000_dp/ ! 196Au
data spin(547) /1.5_dp  /  , nat_iso_ab(547) /  0.000000_dp/ ! 199Au
data spin(548) /2.0_dp  /  , nat_iso_ab(548) /  0.000000_dp/ ! 198Au
data spin(549) /12.0_dp /  , nat_iso_ab(549) /  0.000000_dp/   ! 198m2Au
data spin(550) /1.0_dp  /  , nat_iso_ab(550) /  0.000000_dp/ ! 194Au
data spin(551) /0.0_dp  /  , nat_iso_ab(551) / 29.860000_dp/ ! 202Hg
data spin(552) /0.0_dp  /  , nat_iso_ab(552) / 23.100000_dp/ ! 200Hg
data spin(553) /0.5_dp  /  , nat_iso_ab(553) / 16.870000_dp/ ! 199Hg
data spin(554) /1.5_dp  /  , nat_iso_ab(554) / 13.180000_dp/ ! 201Hg
data spin(555) /0.0_dp  /  , nat_iso_ab(555) /  9.970000_dp/ ! 198Hg
data spin(556) /0.0_dp  /  , nat_iso_ab(556) /  6.870000_dp/ ! 204Hg
data spin(557) /0.0_dp  /  , nat_iso_ab(557) /  0.150000_dp/ ! 196Hg
data spin(558) /0.0_dp  /  , nat_iso_ab(558) /  0.000000_dp/ ! 194Hg
data spin(559) /2.5_dp  /  , nat_iso_ab(559) /  0.000000_dp/ ! 203Hg
data spin(560) /0.5_dp  /  , nat_iso_ab(560) /  0.000000_dp/ ! 197Hg
data spin(561) /6.5_dp  /  , nat_iso_ab(561) /  0.000000_dp/ ! 195mHg
data spin(562) /0.5_dp  /  , nat_iso_ab(562) / 70.480000_dp/ ! 205Tl
data spin(563) /0.5_dp  /  , nat_iso_ab(563) / 29.520000_dp/ ! 203Tl
data spin(564) /2.0_dp  /  , nat_iso_ab(564) /  0.000000_dp/ ! 204Tl
data spin(565) /2.0_dp  /  , nat_iso_ab(565) /  0.000000_dp/ ! 202Tl
data spin(566) /0.5_dp  /  , nat_iso_ab(566) /  0.000000_dp/ ! 201Tl
data spin(567) /2.0_dp  /  , nat_iso_ab(567) /  0.000000_dp/ ! 200Tl
data spin(568) /0.0_dp  /  , nat_iso_ab(568) / 52.400000_dp/ ! 208Pb
data spin(569) /0.0_dp  /  , nat_iso_ab(569) / 24.100000_dp/ ! 206Pb
data spin(570) /0.5_dp  /  , nat_iso_ab(570) / 22.100000_dp/ ! 207Pb
data spin(571) /0.0_dp  /  , nat_iso_ab(571) /  1.400000_dp/ ! 204Pb
data spin(572) /2.5_dp  /  , nat_iso_ab(572) /  0.000000_dp/ ! 205Pb
data spin(573) /0.0_dp  /  , nat_iso_ab(573) /  0.000000_dp/ ! 202Pb
data spin(574) /0.0_dp  /  , nat_iso_ab(574) /  0.000000_dp/ ! 210Pb
data spin(575) /2.5_dp  /  , nat_iso_ab(575) /  0.000000_dp/ ! 203Pb
data spin(576) /4.5_dp  /  , nat_iso_ab(576) /100.000000_dp/ ! 209Bi
data spin(577) /9.0_dp  /  , nat_iso_ab(577) /  0.000000_dp/ ! 210mBi
data spin(578) /5.0_dp  /  , nat_iso_ab(578) /  0.000000_dp/ ! 208Bi
data spin(579) /4.5_dp  /  , nat_iso_ab(579) /  0.000000_dp/ ! 207Bi
data spin(580) /4.5_dp  /  , nat_iso_ab(580) /  0.000000_dp/ ! 205Bi
data spin(581) /6.0_dp  /  , nat_iso_ab(581) /  0.000000_dp/ ! 206Bi
data spin(582) /1.0_dp  /  , nat_iso_ab(582) /  0.000000_dp/ ! 210Bi
data spin(583) /0.5_dp  /  , nat_iso_ab(583) /  0.000000_dp/ ! 209Po
data spin(584) /0.0_dp  /  , nat_iso_ab(584) /  0.000000_dp/ ! 208Po
data spin(585) /0.0_dp  /  , nat_iso_ab(585) /  0.000000_dp/ ! 210Po
data spin(586) /0.0_dp  /  , nat_iso_ab(586) /  0.000000_dp/ ! 206Po
data spin(587) /0.0_dp  /  , nat_iso_ab(587) /  0.000000_dp/ ! 222Rn
data spin(588) /0.0_dp  /  , nat_iso_ab(588) /  0.000000_dp/ ! 226Ra
data spin(589) /0.0_dp  /  , nat_iso_ab(589) /  0.000000_dp/ ! 228Ra
data spin(590) /0.5_dp  /  , nat_iso_ab(590) /  0.000000_dp/ ! 225Ra
data spin(591) /1.5_dp  /  , nat_iso_ab(591) /  0.000000_dp/ ! 223Ra
data spin(592) /0.0_dp  /  , nat_iso_ab(592) /  0.000000_dp/ ! 224Ra
data spin(593) /1.5_dp  /  , nat_iso_ab(593) /  0.000000_dp/ ! 227Ac
data spin(594) /1.5_dp  /  , nat_iso_ab(594) /  0.000000_dp/ ! 225Ac
data spin(595) /1.0_dp  /  , nat_iso_ab(595) /  0.000000_dp/ ! 226Ac
data spin(596) /0.0_dp  /  , nat_iso_ab(596) /100.000000_dp/ ! 232Th
data spin(597) /0.0_dp  /  , nat_iso_ab(597) /  0.000000_dp/ ! 230Th
data spin(598) /2.5_dp  /  , nat_iso_ab(598) /  0.000000_dp/ ! 229Th
data spin(599) /0.0_dp  /  , nat_iso_ab(599) /  0.000000_dp/ ! 228Th
data spin(600) /0.0_dp  /  , nat_iso_ab(600) /  0.000000_dp/ ! 234Th
data spin(601) /0.5_dp  /  , nat_iso_ab(601) /  0.000000_dp/ ! 227Th
data spin(602) /2.5_dp  /  , nat_iso_ab(602) /  0.000000_dp/ ! 231Th
data spin(603) /1.5_dp  /  , nat_iso_ab(603) /  0.000000_dp/ ! 231Pa
data spin(604) /1.5_dp  /  , nat_iso_ab(604) /  0.000000_dp/ ! 233Pa
data spin(605) /2.0_dp  /  , nat_iso_ab(605) /  0.000000_dp/ ! 230Pa
data spin(606) /2.5_dp  /  , nat_iso_ab(606) /  0.000000_dp/ ! 229Pa
data spin(607) /2.0_dp  /  , nat_iso_ab(607) /  0.000000_dp/ ! 232Pa
data spin(608) /0.0_dp  /  , nat_iso_ab(608) / 99.274200_dp/ ! 238U
data spin(609) /3.5_dp  /  , nat_iso_ab(609) /  0.720400_dp/ ! 235U
data spin(610) /0.0_dp  /  , nat_iso_ab(610) /  0.005400_dp/ ! 234U
data spin(611) /0.0_dp  /  , nat_iso_ab(611) /  0.000000_dp/ ! 236U
data spin(612) /2.5_dp  /  , nat_iso_ab(612) /  0.000000_dp/ ! 233U
data spin(613) /0.0_dp  /  , nat_iso_ab(613) /  0.000000_dp/ ! 232U
data spin(614) /0.0_dp  /  , nat_iso_ab(614) /  0.000000_dp/ ! 230U
data spin(615) /0.5_dp  /  , nat_iso_ab(615) /  0.000000_dp/ ! 237U
data spin(616) /2.5_dp  /  , nat_iso_ab(616) /  0.000000_dp/ ! 231U
data spin(617) /2.5_dp  /  , nat_iso_ab(617) /  0.000000_dp/ ! 237Np
data spin(618) /6.0_dp  /  , nat_iso_ab(618) /  0.000000_dp/ ! 236Np
data spin(619) /2.5_dp  /  , nat_iso_ab(619) /  0.000000_dp/ ! 235Np
data spin(620) /0.0_dp  /  , nat_iso_ab(620) /  0.000000_dp/ ! 234Np
data spin(621) /2.5_dp  /  , nat_iso_ab(621) /  0.000000_dp/ ! 239Np
data spin(622) /2.0_dp  /  , nat_iso_ab(622) /  0.000000_dp/ ! 238Np
data spin(623) /0.0_dp  /  , nat_iso_ab(623) /  0.000000_dp/ ! 244Pu
data spin(624) /0.0_dp  /  , nat_iso_ab(624) /  0.000000_dp/ ! 242Pu
data spin(625) /0.5_dp  /  , nat_iso_ab(625) /  0.000000_dp/ ! 239Pu
data spin(626) /0.0_dp  /  , nat_iso_ab(626) /  0.000000_dp/ ! 240Pu
data spin(627) /0.0_dp  /  , nat_iso_ab(627) /  0.000000_dp/ ! 238Pu
data spin(628) /2.5_dp  /  , nat_iso_ab(628) /  0.000000_dp/ ! 241Pu
data spin(629) /0.0_dp  /  , nat_iso_ab(629) /  0.000000_dp/ ! 236Pu
data spin(630) /3.5_dp  /  , nat_iso_ab(630) /  0.000000_dp/ ! 237Pu
data spin(631) /0.0_dp  /  , nat_iso_ab(631) /  0.000000_dp/ ! 246Pu
data spin(632) /0.5_dp  /  , nat_iso_ab(632) /  0.000000_dp/ ! 247Pu
data spin(633) /2.5_dp  /  , nat_iso_ab(633) /  0.000000_dp/ ! 243Am
data spin(634) /2.5_dp  /  , nat_iso_ab(634) /  0.000000_dp/ ! 241Am
data spin(635) /5.0_dp  /  , nat_iso_ab(635) /  0.000000_dp/ ! 242mAm
data spin(636) /3.0_dp  /  , nat_iso_ab(636) /  0.000000_dp/ ! 240Am
data spin(637) /4.5_dp  /  , nat_iso_ab(637) /  0.000000_dp/ ! 247Cm
data spin(638) /0.0_dp  /  , nat_iso_ab(638) /  0.000000_dp/ ! 248Cm
data spin(639) /3.5_dp  /  , nat_iso_ab(639) /  0.000000_dp/ ! 245Cm
data spin(640) /0.0_dp  /  , nat_iso_ab(640) /  0.000000_dp/ ! 250Cm
data spin(641) /0.0_dp  /  , nat_iso_ab(641) /  0.000000_dp/ ! 246Cm
data spin(642) /2.5_dp  /  , nat_iso_ab(642) /  0.000000_dp/ ! 243Cm
data spin(643) /0.0_dp  /  , nat_iso_ab(643) /  0.000000_dp/ ! 244Cm
data spin(644) /0.0_dp  /  , nat_iso_ab(644) /  0.000000_dp/ ! 242Cm
data spin(645) /0.5_dp  /  , nat_iso_ab(645) /  0.000000_dp/ ! 241Cm
data spin(646) /0.0_dp  /  , nat_iso_ab(646) /  0.000000_dp/ ! 240Cm
data spin(647) /0.0_dp  /  , nat_iso_ab(647) /  0.000000_dp/ ! 252Cm
data spin(648) /1.5_dp  /  , nat_iso_ab(648) /  0.000000_dp/ ! 247Bk
data spin(649) /6.0_dp  /  , nat_iso_ab(649) /  0.000000_dp/ ! 248Bk
data spin(650) /3.5_dp  /  , nat_iso_ab(650) /  0.000000_dp/ ! 249Bk
data spin(651) /1.5_dp  /  , nat_iso_ab(651) /  0.000000_dp/ ! 245Bk
data spin(652) /2.0_dp  /  , nat_iso_ab(652) /  0.000000_dp/ ! 246Bk
data spin(653) /0.5_dp  /  , nat_iso_ab(653) /  0.000000_dp/ ! 251Cf
data spin(654) /4.5_dp  /  , nat_iso_ab(654) /  0.000000_dp/ ! 249Cf
data spin(655) /0.0_dp  /  , nat_iso_ab(655) /  0.000000_dp/ ! 250Cf
data spin(656) /0.0_dp  /  , nat_iso_ab(656) /  0.000000_dp/ ! 252Cf
data spin(657) /0.0_dp  /  , nat_iso_ab(657) /  0.000000_dp/ ! 248Cf
data spin(658) /0.0_dp  /  , nat_iso_ab(658) /  0.000000_dp/ ! 254Cf
data spin(659) /3.5_dp  /  , nat_iso_ab(659) /  0.000000_dp/ ! 253Cf
data spin(660) /0.0_dp  /  , nat_iso_ab(660) /  0.000000_dp/ ! 246Cf
data spin(661) /4.0_dp  /  , nat_iso_ab(661) /  0.000000_dp/ ! 252Es
data spin(662) /7.0_dp  /  , nat_iso_ab(662) /  0.000000_dp/ ! 254Es
data spin(663) /3.5_dp  /  , nat_iso_ab(663) /  0.000000_dp/ ! 255Es
data spin(664) /3.5_dp  /  , nat_iso_ab(664) /  0.000000_dp/ ! 253Es
data spin(665) /3.5_dp  /  , nat_iso_ab(665) /  0.000000_dp/ ! 257Es
data spin(666) /2.0_dp  /  , nat_iso_ab(666) /  0.000000_dp/ ! 254mEs
data spin(667) /1.5_dp  /  , nat_iso_ab(667) /  0.000000_dp/ ! 251Es
data spin(668) /4.5_dp  /  , nat_iso_ab(668) /  0.000000_dp/ ! 257Fm
data spin(669) /0.5_dp  /  , nat_iso_ab(669) /  0.000000_dp/ ! 253Fm
data spin(670) /0.0_dp  /  , nat_iso_ab(670) /  0.000000_dp/ ! 252Fm
data spin(671) /8.0_dp  /  , nat_iso_ab(671) /  0.000000_dp/ ! 258Md
data spin(672) /spin_undef/, nat_iso_ab(672) /  0.000000_dp/ ! 260Md
data spin(673) /spin_undef/, nat_iso_ab(673) /  0.000000_dp/ ! 268Db

contains

! takes as input an atomic mass and finds the atom in the database
! which is closest
subroutine get_name_from_mass(xmass, atom_name)
  real(kind=dp) :: xmass
  character(len=15), intent(out) :: atom_name
  real(kind=dp) :: min_diff, diff 
  real(kind=dp) :: max_allowed_diff=0.1_dp ! maximum difference in mass for allowing a match
  integer       :: i, imin
  integer       :: m
  character(len=2) :: m_tail ! contains m1, m2, m3 (nuclear isomer flag) if necessary

  imin = 0 ; min_diff = 1d30
  do i=1, n_nuclides
    diff =abs(mass(i) - xmass)
    if( diff < min_diff) then
       min_diff = diff
       imin = i
    endif
  enddo

  if( min_diff >= max_allowed_diff) imin=0

  m=flag(imin)/10

  select case(m)
    case default
          m_tail=''
    case(1)
          m_tail='m'
    case(2)
          m_tail='m2'
    case(3)
          m_tail='m3'
  end select

!   atom_name = nuclide_name(Z(imin),A(imin), flag(imin)/10 ) ! flag is 10-19 for m, 20-29 for m1 etc (nuclear isomers)
   if (iverbose>=4) write (atom_name,'(A,a1,i0,a)') trim(element_symbol(Z(imin))), '-',A(imin), m_tail

end subroutine get_name_from_mass

! gives a character variable with the name of the element
pure function element_name(zin, ain)
  integer, intent(in)           :: zin
  integer, intent(in), optional :: ain
  character(len=15) :: element_name

  select case(zin)
   case default
     element_name = 'Unknown'
   case(1) ! special case for hydrogen
     element_name = Name(zin)
     if( present(ain) .and. ain == 2) element_name = 'Deuterium'
     if( present(ain) .and. ain == 3) element_name = 'Tritium'
   case(2:zmax)
    element_name = Name(zin)
  end select

end function element_name


! gives a character variable with the full name of the nuclide (eg., Carbon-14 )
pure function nuclide_name(zin, ain, m)
  integer, intent(in)  :: zin
  integer, optional, intent(in)  :: ain
  integer, optional, intent(in)  :: m ! m=0,1,2,3 counts over excited nuclear states
  character(len=25) :: nuclide_name
  integer :: my_flag
  character(len=2) :: label

  select case(zin)
   case default
     nuclide_name = 'Unknown'

   case(1:zmax)
    if( present(ain)         .and. present(m)         ) my_flag =flag(ii(zin,ain,m))
    if( present(ain)         .and. (.not. present(m)) ) my_flag =flag(ii(zin,ain))
    if( (.not. present(ain)) .and. present(m)         ) my_flag =flag(ii(zin,m=m))
    if( (.not. present(ain)) .and. (.not. present(m)) ) my_flag =flag(ii(zin ))

    if( present(ain)                                   ) write(nuclide_name, '(A,I0)')  trim(Name(zin)) // "-", ain
    if( (.not. present(ain)) .and. present(m)          ) write(nuclide_name, '(A,I0)')  trim(Name(zin)) // "-", A(ii(zin,m=m))
    if( (.not. present(ain)) .and. (.not. present(m))  ) write(nuclide_name, '(A,I0)')  trim(Name(zin)) // "-", A(ii(zin ))


    select case(my_flag)
     case default  ; label = " "
     case( 1:9)    ; label = " "
     case(11:19)   ; label = "m"
     case(21:29)   ; label = "m2"
     case(31:39)   ; label = "m3"
     end select
    nuclide_name = trim(nuclide_name) // trim(label)
   end select


end function nuclide_name


! gives a character variable with the symbol of the element
pure function element_symbol(zin, ain)
  integer, intent(in)           :: zin
  integer, intent(in), optional :: ain
  character(len=15) :: element_symbol

  select case(zin)
   case default
     element_symbol = 'Unknown'
   case(1) ! special case for hydrogen
     element_symbol = Symbol(zin)
     if( present(ain) ) then 
      if (ain == 2) element_symbol = 'D'
     endif
     if( present(ain) ) then 
      if (ain == 3) element_symbol = 'T'
     endif
   case(2:zmax)
    element_symbol = Symbol(zin)
  end select

end function element_symbol


! takes in Z and A and returns the counter ii identifying each nuclide
! used only internally in this module
! gives ii = 0 if not found
pure function ii(zin, ain, m)
  integer, intent(in)  :: zin
  integer, optional, intent(in)  :: ain
  integer, optional, intent(in)  :: m
  integer :: ii, j, ioption

  if( present(ain)         .and. present(m)         ) ioption=1 ! match zin, ain, m
  if( present(ain)         .and. (.not. present(m)) ) ioption=2 ! match zin, ain
  if( (.not. present(ain)) .and. present(m)         ) ioption=3 ! match zin, m
  if( (.not. present(ain)) .and. (.not. present(m)) ) ioption=4 ! match zin


  ii = 0 ! means: not found
  search_loop: do j=1, n_nuclides
    if( Z(j)  > zin ) exit  search_loop
    if( Z(j)  < zin ) cycle search_loop

    select case(ioption)
      case default
        return
      case(1)
        if( ain == A(j) .and. m == flag(j)/10) then
         ii = j
         return
         endif
      case(2)
        if( ain == A(j)) then
         ii = j
         return
        endif
        case(3)
        if( m == flag(j)/10) then
        ii = j
        return
        endif

      case(4)
        ii = j
        return
      end select

    enddo search_loop

end function ii


! gives atomic mass of the element
pure function atomic_mass(zin,ain,m)
   integer, intent(in)  :: zin
   integer, optional, intent(in)  :: ain
   integer, optional, intent(in)  :: m
   real(kind=dp) :: atomic_mass
 
    if( present(ain)         .and. present(m)         ) atomic_mass = mass( ii(zin, ain,m) )
    if( present(ain)         .and. (.not. present(m)) ) atomic_mass = mass( ii(zin, ain) )
    if( (.not. present(ain)) .and. present(m)         ) atomic_mass = mass( ii(zin, m=m) )
    if( (.not. present(ain)) .and. (.not. present(m)) ) atomic_mass = mass( ii(zin) )
   
 
end function atomic_mass

! gives atomic mass uncertainty of the element
pure function atomic_mass_unc(zin,ain,m)
   integer, intent(in)  :: zin
   integer, intent(in), optional  :: ain
   integer, intent(in), optional  :: m
   real(kind=dp) :: atomic_mass_unc

    if( present(ain)         .and. present(m)         ) atomic_mass_unc = mass_unc( ii(zin, ain,m) )
    if( present(ain)         .and. (.not. present(m)) ) atomic_mass_unc = mass_unc( ii(zin, ain) )
    if( (.not. present(ain)) .and. present(m)         ) atomic_mass_unc = mass_unc( ii(zin, m=m) )
    if( (.not. present(ain)) .and. (.not. present(m)) ) atomic_mass_unc = mass_unc( ii(zin) )
      
end function atomic_mass_unc


pure function stability_info(zin, ain, m)
   integer, intent(in)              :: zin
   integer, intent(in),optional  :: ain
   integer, intent(in), optional  :: m   
   character(len=60)    :: stability_info
   integer :: my_flag
   
    if( present(ain)         .and. present(m)         ) my_flag = flag( ii(zin, ain,m) )
    if( present(ain)         .and. (.not. present(m)) ) my_flag = flag( ii(zin, ain) )
    if( (.not. present(ain)) .and. present(m)         ) my_flag = flag( ii(zin, m=m) )
    if( (.not. present(ain)) .and. (.not. present(m)) ) my_flag = flag( ii(zin) )
   
   
   select case( my_flag )
     case default
       stability_info = 'Unknown'

     case(1)
       stability_info = 'Stable nuclide'
     case(2)
       stability_info = 'Observationally stable nuclide'
     case(3)
       stability_info = 'Radioactive nuclide (primordial)'
     case(4)
       stability_info = 'Radioactive nuclide'
       
     case(11)
       stability_info = 'Stable nuclide (metastable isomer, m)' !This case should NEVER happen!
     case(12)
       stability_info = 'Observationally stable nuclide (metastable isomer, m)'
     case(13)
       stability_info = 'Radioactive nuclide (primordial; metastable isomer, m)'
     case(14)
       stability_info = 'Radioactive nuclide (metastable isomer, m)'

     case(21)
       stability_info = 'Stable nuclide (metastable isomer, m2)' !This case should NEVER happen!
     case(22)
       stability_info = 'Observationally stable nuclide (metastable isomer, m2)'
     case(23)
       stability_info = 'Radioactive nuclide (primordial; metastable isomer, m2)'
     case(24)
       stability_info = 'Radioactive nuclide (metastable isomer, m2)'

     case(31)
       stability_info = 'Stable nuclide (metastable isomer, m3)' !This case should NEVER happen!
     case(32)
       stability_info = 'Observationally stable nuclide (metastable isomer, m3)'
     case(33)
       stability_info = 'Radioactive nuclide (primordial; metastable isomer, m3)'
     case(34)
       stability_info = 'Radioactive nuclide (metastable isomer, m3)'       
       
   end select


end function stability_info



pure function stability_info_flag(zin, ain,m)
   integer, intent(in)  :: zin
   integer, intent(in), optional  :: ain
   integer, intent(in), optional  :: m
   integer    :: stability_info_flag

   stability_info_flag = flag(ii(zin, ain,m))
  
end function stability_info_flag



pure function half_life(zin, ain,m)
   integer, intent(in)   :: zin
   integer, intent(in), optional  :: ain
   integer, intent(in), optional  :: m

   real(kind=dp)        :: half_life

   half_life = h_f( ii(zin, ain, m) )

end function half_life

pure function half_life_txt(zin, ain,m)
   character(len=30) :: half_life_txt
   integer, intent(in)   :: zin
   integer, intent(in), optional  :: ain
   integer, intent(in), optional  :: m
   real(kind=dp)        :: half_life
   real(kind=dp),parameter    :: hour=3600._dp 
   real(kind=dp),parameter    :: day=86400._dp
   real(kind=dp),parameter    :: year=31.556926e6_dp ! definition used in Nubase2012
   real(kind=dp),parameter    :: kyear=1e3_dp*year
   real(kind=dp),parameter    :: Myear=1e6_dp*year
   real(kind=dp),parameter    :: Gyear=1e9_dp*year
   real(kind=dp),parameter    :: Tyear=1000._dp*Gyear
   
    if( present(ain)         .and. present(m)         ) half_life = h_f( ii(zin, ain, m) )
    if( present(ain)         .and. (.not. present(m)) ) half_life = h_f( ii(zin, ain) )
    if( (.not. present(ain)) .and. present(m)         ) half_life = h_f( ii(zin, m=m) )
    if( (.not. present(ain)) .and. (.not. present(m)) ) half_life = h_f( ii(zin) )

   if(half_life <= day)                     write(half_life_txt,'(f10.3, A)')  half_life/hour  , ' hours'
   if(half_life  > day   .and. half_life <= year ) write(half_life_txt,'(f10.3, A)')  half_life/day   , ' days'
   if(half_life  > year  .and. half_life <= kyear) write(half_life_txt,'(f10.3, A)')  half_life/year  , ' years'
   if(half_life  > kyear .and. half_life <= Myear) write(half_life_txt,'(f10.3, A)')  half_life/kyear , ' thousand years'
   if(half_life  > Myear .and. half_life <= Gyear) write(half_life_txt,'(f10.3, A)')  half_life/Myear , ' million years'
   if(half_life  > Gyear .and. half_life <= Tyear) write(half_life_txt,'(f10.3, A)')  half_life/Gyear , ' billion years'
   if(half_life  > Tyear)                          write(half_life_txt,'(es10.3, A)') half_life/year  , ' years'   
   
   
end function half_life_txt



pure function nuclear_spin(zin, ain, m)
   real(kind=dp) :: nuclear_spin
   integer, intent(in)   :: zin
   integer, intent(in), optional  :: ain
   integer, intent(in), optional  :: m


    if( present(ain)         .and. present(m)         ) nuclear_spin = spin( ii(zin, ain,m) )
    if( present(ain)         .and. (.not. present(m)) ) nuclear_spin = spin( ii(zin, ain) )
    if( (.not. present(ain)) .and. present(m)         ) nuclear_spin = spin( ii(zin, m=m) )
    if( (.not. present(ain)) .and. (.not. present(m)) ) nuclear_spin = spin( ii(zin) )

end function nuclear_spin


pure function nat_iso_abundance(zin, ain, m)
   real(kind=dp) :: nat_iso_abundance
   integer, intent(in)   :: zin
   integer, intent(in), optional  :: ain
   integer, intent(in), optional  :: m

    if( present(ain)         .and. present(m)         ) nat_iso_abundance = nat_iso_ab( ii(zin, ain,m) )
    if( present(ain)         .and. (.not. present(m)) ) nat_iso_abundance = nat_iso_ab( ii(zin, ain) )
    if( (.not. present(ain)) .and. present(m)         ) nat_iso_abundance = nat_iso_ab( ii(zin, m=m) )
    if( (.not. present(ain)) .and. (.not. present(m)) ) nat_iso_abundance = nat_iso_ab( ii(zin) )


end function nat_iso_abundance

! Takes in a chemical symbol and gives the corresponding Z atomic number (=nuclear charge)
! Accepts the following formats:
! 1) chemical symbol only:  C, O, Sc,...
! 2) chemical symbol with post-fixed atomic number: C-14, O-16, Sc-44m3
! Gives back zero if not found.
pure function get_z(string)
integer, parameter :: max_length=8
character(len=max_length), intent(in) :: string
character(len=2) :: my_symbol
integer :: get_z
integer :: j

get_z = 0

my_symbol = string(1:2)
if( my_symbol(2:2) == '-') my_symbol(2:2)=' '

do j=1, zmax
  if( my_symbol(1:2) == Symbol(j)) then
    get_z = j
    return
  elseif ( my_symbol(1:2) == 'D '.or. my_symbol(1:2) == 'T ') then    ! Deal with deuterium and tritium separately
    get_z = 1
    return
  endif
enddo




end function get_z


! Takes in a chemical symbol and gives the corresponding A atomic number (=number of nucleons)
! Accepts the following formats:
! 1) chemical symbol with post-fixed atomic number: C-14, O-16, Sc-44m3
! Gives back zero if not found.
pure function get_a(string)
integer, parameter :: max_length=8
character(len=max_length), intent(in) :: string
integer :: get_a
integer :: j, jstart, jend, ierr

get_a = 0

! get starting position for atomic number
jstart=1
do j=1, max_length-1
  if( string(j:j)== '-') jstart=j+1
enddo

! get end position for atomic number
jend=max_length
do j=jstart+1, max_length
  if( string(j:j) == 'm') then
     jend=j-1
     exit
  endif
enddo

read( string(jstart:jend), *, iostat=ierr ) get_a
if(ierr /=0)   get_a = 0
if( get_a < 0) get_a = 0

! deal here with deuterium and tritium
if( string(1:2) == "D ") get_a=2
if( string(1:2) == "T ") get_a=3



end function get_a


! Takes in a chemical symbol and gives the corresponding m nuclear number, if given (=excitation of nucleus)
! Accepts the following formats:
! 1) chemical symbol with post-fixed atomic number: C-14, O-16, Sc-44m3
! Gives back -1 if not found.
pure function get_m(string)
integer, parameter :: max_length=8
character(len=max_length), intent(in) :: string
integer :: get_m
integer :: j, ierr

get_m = -1

do j=1, max_length-1
  if( string(j:j)/= 'm') cycle
  read(string(j+1:j+1), '(I1)', iostat=ierr) get_m
  if( ierr /=0) get_m = -1
  if( string(j+1:j+1) == " ") get_m = 1
  return
enddo

end function get_m


subroutine print_atomic_and_nuclear_info(my_string,verbose,unit)
 integer, parameter :: max_length=8
 character(len=max_length), intent(in) :: my_string
 integer, intent(in)  :: verbose
 integer, optional, intent(in)  :: unit
 integer :: u1, ioption, zi, ai, m
 character(len=20) my_element_name
 character(len=15) my_element_symbol
 character(len=30) my_nuclide_name
 real(kind=dp)  :: my_atomic_mass, my_atomic_mass_unc, my_nuclear_mass, my_nuclear_plus_core_mass, &
                   my_nuclear_plus_core_mass2, my_nuclear_spin, my_half_life, my_iso_abundance
 character(len=60) :: my_stability_info, my_half_life_txt
 character(len=60) :: my_fmt
 
 u1 = 6 ! use standard (`screen') output if unit number is missing
 if( present(unit) ) u1=unit
 !
 iverbose = verbose
 !
 zi = get_z(my_string)
 ai = get_a(my_string)
 m  = get_m(my_string)

  if( zi <= 0) then
    write(u1,'(A)') 'Error: could not properly read Z from element symbol: ' // trim(my_string)
    return
  endif

  !
  ! Both ai and m are optional, so I deal separately with each of the four possibilities
  !  (both given, only one given, none given)
  !
  if( ai >0   .and. m >= 0 ) ioption=1 ! ain and m both present      (rare case)
  if( ai >0   .and. m <  0 ) ioption=2 ! ain present,  m not present (common case)
  if( ai <= 0 .and. m >= 0 ) ioption=3 ! m   present, ai not present (very rare!!!)
  if( ai <= 0 .and. m  < 0 ) ioption=4 ! neither ai nor m present    (common case)

  select case(ioption)
  case default
     write(u1,'(A,I0)') 'Error while reading chemical symbol. ioption= ', ioption
  case(1)
    my_element_name   = element_name(     zi,ai)
    my_element_symbol = element_symbol(   zi,ai)
    my_nuclide_name   = nuclide_name(     zi,ai,m)
    my_atomic_mass    = atomic_mass(      zi,ai,m)
    my_atomic_mass_unc= atomic_mass_unc(  zi,ai,m)
    my_nuclear_spin   = nuclear_spin(     zi,ai,m)
    my_stability_info = stability_info(   zi,ai,m) 
    my_half_life      = half_life(        zi,ai,m)
    my_half_life_txt  = half_life_txt(    zi,ai,m)
    my_iso_abundance  = nat_iso_abundance(zi,ai,m)
  case(2)
    my_element_name   = element_name(     zi,ai)
    my_element_symbol = element_symbol(   zi,ai)
    my_nuclide_name   = nuclide_name(     zi,ai)
    my_atomic_mass    = atomic_mass(      zi,ai)
    my_atomic_mass_unc= atomic_mass_unc(  zi,ai)
    my_nuclear_spin   = nuclear_spin(     zi,ai)
    my_stability_info = stability_info(   zi,ai)
    my_half_life      = half_life(        zi,ai)
    my_half_life_txt  = half_life_txt(    zi,ai)
    my_iso_abundance  = nat_iso_abundance(zi,ai)
  case(3)
    my_element_name   = element_name(     zi)
    my_element_symbol = element_symbol(   zi)
    my_nuclide_name   = nuclide_name(     zi,m=m)
    my_atomic_mass    = atomic_mass(      zi,m=m)
    my_atomic_mass_unc= atomic_mass_unc(  zi,m=m)
    my_nuclear_spin   = nuclear_spin(     zi,m=m)
    my_stability_info = stability_info(   zi,m=m)
    my_half_life      = half_life(        zi,m=m)
    my_half_life_txt  = half_life_txt(    zi,m=m)
    my_iso_abundance  = nat_iso_abundance(zi,m=m)
  case(4)
    my_element_name   = element_name(     zi)
    my_element_symbol = element_symbol(   zi)
    my_nuclide_name   = nuclide_name(     zi)
    my_atomic_mass    = atomic_mass(      zi)
    my_atomic_mass_unc= atomic_mass_unc(  zi)
    my_nuclear_spin   = nuclear_spin(     zi)
    my_stability_info = stability_info(   zi)
    my_half_life      = half_life(        zi)
    my_half_life_txt  = half_life_txt(    zi)
    my_iso_abundance  = nat_iso_abundance(zi)
  end select

  ! NB in principle for consistency I should use the value of the electron mass from the physical_constants
  ! module (but it really doesn't matter).
  my_nuclear_mass            =  my_atomic_mass  - real(zi,dp)          / 1822.8884861185961_dp
  my_nuclear_plus_core_mass  =  my_nuclear_mass + real(n_core(zi),dp)  / 1822.8884861185961_dp
  my_nuclear_plus_core_mass2 =  my_nuclear_mass + real(n_core2(zi),dp) / 1822.8884861185961_dp

  if (iverbose>=4) write(u1,'(A,a)'       )  'Nuclide full name          = ', trim(my_nuclide_name)
  if (iverbose>=4) write(u1,'(A,A)'       )  'Element name               = ', trim(my_element_name)
  if (iverbose>=4) write(u1,'(A,A)'       )  'Element symbol             = ', trim(my_element_symbol)
  if (iverbose>=4) write(u1,'(A,I5)'      )  'Nuclear charge Z           = ', zi
  if( my_atomic_mass <= 0._dp) then
     write(u1,'(A)') 'This element was not found in the internal database. '
     write(u1,'(A)') 'Probably radioactive with a half-life < 1 hour.'
     return
  endif
  my_fmt = '(A,F20.12,A, es9.2,a)'
  if (iverbose>=4) write(u1,my_fmt)          'Atomic mass                = ', my_atomic_mass, ' +/- ', my_atomic_mass_unc , &
                                             ' Daltons'
  if (iverbose>=4) write(u1,'(A,F20.12,A)')  'Nuclear mass               = ', my_nuclear_mass,'   (=atomic mass - Z*me)'

  if (iverbose>=4) write(u1,'(A,F20.12,A,i3,a)')     'Nuclear mass + core elec.  = ', my_nuclear_plus_core_mass, &
                                                                   '   (=nuclear mass + ',n_core(zi),'*me)' !&
!    // " [core elec. = those in filled noble-gas shells]"

!   write(u1,'(A,F20.12,A,i3,a)')     'Nuclear mass + core elec.  = ', my_nuclear_plus_core_mass2, &
!                                                                    '   (nuclear mass + ',n_core2(zi),'*me)' &
!    // " [core elec. = those in orbitals with energy < -1.5 Eh]"


  if (iverbose>=4) write(u1,'(A,F5.1,A)'  )  'Nuclear spin               = ', my_nuclear_spin, ' hbar'
  if (iverbose>=4) write(u1,'(A, A)'      )  'Nuclear stability          = ', trim(my_stability_info)

  if(iverbose>=4 .and. trim(my_stability_info) /= "Stable nuclide") &
  write(u1,'(A, es14.3,A)') 'Half-life                  = ',my_half_life,' seconds' // ' ( ' // trim(my_half_life_txt) // ')'
  if (iverbose>=4) write(u1,'(A, F12.6,a)' ) 'Natural isotopic abundance = ', my_iso_abundance, ' %'

end subroutine print_atomic_and_nuclear_info

end module atomic_and_nuclear_data

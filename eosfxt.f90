      program teos
      include 'implno.dek'
      include 'vector_eos.dek'

! tests the eos routine
! 
! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

      integer          ionmax
      parameter        (ionmax=3)
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),temp,den,abar,zbar


! set the mass fractions, z's and a's of the composition
! hydrogen, heliu, and carbon
      xmass(1) = 0.75d0 ; aion(1)  = 1.0d0  ; zion(1)  = 1.0d0
      xmass(2) = 0.23d0 ; aion(2)  = 4.0d0  ; zion(2)  = 2.0d0
      xmass(3) = 0.02d0 ; aion(3)  = 12.0d0 ; zion(3)  = 6.0d0

! average atomic weight and charge
      abar   = 1.0d0/sum(xmass(1:ionmax)/aion(1:ionmax))
      zbar   = abar * sum(xmass(1:ionmax) * zion(1:ionmax)/aion(1:ionmax))

! set the input vector. pipeline is only 1 element long
      temp_row(1) = 1.0d8 ; den_row(1)  = 1.0d6 ; abar_row(1) = abar ; zbar_row(1) = zbar
      jlo_eos = 1 ; jhi_eos = 1

! call the eos
      call eosfxt

! write out the results
      call pretty_eos_out('eosfxt:  ')

      end   




!
! this file contains
! timmes eos: eosfxt
! routine xneroot does electron-positron thermodynamics
! routine coulomb implments coulomb corrections
! routine etages makes a good guess for the electron chemical potential






      subroutine eosfxt
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'

! given a temperature temp [K], density den [g/cm**3], and a composition
! characterized by abar (average weight) and zbar (average charge),
! this routine returns all the other thermodynamic quantities.

! of interest is the pressure [erg/cm**3], specific thermal energy [erg/gr],
! the entropy [erg/g/K], with their derivatives with respect to temperature,
! density, abar, and zbar.

! other quantites such the normalized chemical potential eta (plus its
! derivatives), number density of electrons and positron pair (along
! with their derivatives), adiabatic indices, specific heats, and
! relativistically correct sound speed are also returned.

! this routine assumes planckian photons, an ideal gas of ions,
! and an electron-positron gas with an arbitrary degree of relativity
! and degeneracy. the full fermi-dirac integrals and their derivatives
! with respect to eta and beta are computed to machine precision, and
! all other derivatives are analytic.

! references: cox & giuli (c&g) chapter 24,
!             timmes & arnett, apj supp. 125, 277, 1999
!             timmes & swesty, apj supp. 126, 501, 2000



! bring in the dictionary

! a dictionary of terms used:
! this routine has now been pipelined.
! all the input and output variables are in the file vector_eos.dek.
! the vector name is the scaler name appended with an "_row",
! for example, temp_row(i), den_row(i), and so on.



! input:
! temp     = temperature
! den      = density
! abar     = average number of nucleons per nuclei
! zbar     = average number of protons per nuclei


! output:

! pres     = total pressure
! dpresdd  = derivative of total pressure with respect to density
! dpresdt  = derivative of total pressure with respect to temperature
! dpresda  = derivative of total pressure with respect to abar
! dpresdz  = derivative of total pressure with respect to zbar

! ener     = total internal energy
! denerdd  = derivative of total energy with respect to density
! denerdt  = derivative of total energy with respect to temperature
! denerda  = derivative of total energy with respect to abar
! denerdz  = derivative of total energy with respect to zbar

! entr     = total entropy
! dentrdd  = derivative of total entropy with respect to density
! dentrdt  = derivative of total entropy with respect to temperature
! dentrda  = derivative of total entropy with respect to abar
! dentrdz  = derivative of total entropy with respect to zbar



! prad     = radiation pressure
! dpraddd  = derivative of the radiation pressure with density
! dpraddt  = derivative of the radiation pressure with temperature
! dpradda  = derivative of the radiation pressure with abar
! dpraddz  = derivative of the radiation pressure with zbar

! erad     = radiation energy
! deraddd  = derivative of the radiation energy with density
! deraddt  = derivative of the radiation energy with temperature
! deradda  = derivative of the radiation energy with abar
! deraddz  = derivative of the radiation energy with zbar

! srad     = radiation entropy
! dsraddd  = derivative of the radiation entropy with density
! dsraddt  = derivative of the radiation entropy with temperature
! dsradda  = derivative of the radiation entropy with abar
! dsraddz  = derivative of the radiation entropy with zbar

! radmult  = radiation multiplier (useful for turning radiation off/on)




! xni      = number density of ions
! dxnidd   = derivative of the ion number density with density
! dxnidt   = derivative of the ion number density with temperature
! dxnida   = derivative of the ion number density with abar
! dxnidz   = derivative of the ion number density with zbar

! pion     = ion pressure
! dpiondd  = derivative of the ion pressure with density
! dpiondt  = derivative of the ion pressure with temperature
! dpionda  = derivative of the ion pressure with abar
! dpiondz  = derivative of the ion pressure with zbar

! eion     = ion energy
! deiondd  = derivative of the ion energy with density
! deiondt  = derivative of the ion energy with temperature
! deionda  = derivative of the ion energy with abar
! deiondz  = derivative of the ion energy with zbar

! sion     = ion entropy
! dsiondd  = derivative of the ion entropy with density
! dsiondt  = derivative of the ion entropy with temperature
! dsionda  = derivative of the ion entropy with abar
! dsiondz  = derivative of the ion entropy with zbar

! ionmult  = ion multiplier (useful for turning ions off/on)


! etaele   = electron chemical potential
! detadd   = derivative of the electron chem potential with density
! detadt   = derivative of the electron chem potential with temperature
! detada   = derivative of the electron chem potential with abar
! detadz   = derivative of the electron chem potential with zbar

! etapos   = positron degeneracy parameter

! xne       = number density of electrons
! dxnedd    = derivative of the electron number density with density
! dxnedt    = derivative of the electron number density with temperature
! dxneda    = derivative of the electron number density with abar
! dxnedz    = derivative of the electron number density with zbar

! xnefer    = fermi integral electron number density
! dxneferdd = derivative of the fermi electron number density with density
! dxneferdt = derivative of the fermi electron number density with temperature
! dxneferda = derivative of the fermi electron number density with abar
! dxneferdz = derivative of the fermi electron number density with zbar

! xnpfer    = fermi integral positron number density
! dxnpferdd = derivative of the fermi positron number density with density
! dxnpferdt = derivative of the fermi positron number density with temperature
! dxnpferda = derivative of the fermi positron number density with abar
! dxnpferdz = derivative of the fermi positron number density with zbar

! pele      = electron pressure
! dpeledd   = derivative of the electron pressure with density
! dpeledt   = derivative of the electron pressure with temperature
! dpeleda   = derivative of the electron pressure with abar
! dpeledz   = derivative of the electron pressure with zbar

! eele     = electron energy
! deeledd   = derivative of the electron energy with density
! deeledt   = derivative of the electron energy with temperature
! deeleda   = derivative of the electron energy with abar
! deeledz   = derivative of the electron energy with zbar

! sele     = electron entropy
! dseledd   = derivative of the electron entropy with density
! dseledt   = derivative of the electron entropy with temperature
! dseleda   = derivative of the electron entropy with abar
! dseledz   = derivative of the electron entropy with zbar


! ppos     = positron pressure
! dpposdd   = derivative of the positron pressure with density
! dpposdt   = derivative of the positron pressure with temperature
! dpposda   = derivative of the positron pressure with abar
! dpposdz   = derivative of the positron pressure with zbar

! epos     = electron energy
! deposdd   = derivative of the positron energy with density
! deposdt   = derivative of the positron energy with temperature
! deposda   = derivative of the positron energy with abar
! deposdz   = derivative of the positron energy with zbar

! spos     = electron entropy
! dsposdd   = derivative of the positron entropy with density
! dsposdt   = derivative of the positron entropy with temperature
! dsposda   = derivative of the positron entropy with abar
! dsposdz   = derivative of the positron entropy with zbar

! pep      = electron + positron pressure
! dpepdd   = derivative of the electron+positron pressure with density
! dpepdt   = derivative of the electron+positron pressure with temperature
! dpepda   = derivative of the electron+positron pressure with abar
! dpepdz   = derivative of the electron+positron pressure with zbar

! eep      = electron + positron energy
! deepdd   = derivative of the electron+positron energy with density
! deepdt   = derivative of the electron+positron energy with temperature
! deepda   = derivative of the electron+positron energy with abar
! deepdz   = derivative of the electron+positron energy with zbar

! sep      = electron + positron entropy
! dsepdd   = derivative of the electron+positron entropy with density
! dsepdt   = derivative of the electron+positron entropy with temperature
! dsepda   = derivative of the electron+positron entropy with abar
! dsepdz   = derivative of the electron+positron entropy with zbar

! elemult  = electron multiplier (useful for turning e-e+ off/on)


! eip      = ionization potential ennergy
! deipdd   = derivative of ionization energy with density
! deipdt   = derivative of ionization energy with temperature
! deipda   = derivative of ionization energy with abar
! deipdz   = derivative of ionization energy with zbar


! sip      = ionization potential ennergy
! dsipdd   = derivative of ionization energy with density
! dsipdt   = derivative of ionization energy with temperature
! dsipda   = derivative of ionization energy with abar
! dsipdz   = derivative of ionization energy with zbar

! potmult  = ionization energy multiplier (useful for turning off ionization additions)



! pcoul    = coulomb pressure correction
! coulmult = coulomb component multiplier
! dpcouldd = derivative of the coulomb pressure with density
! dpcouldt = derivative of the coulomb pressure with temperature
! dpcoulda = derivative of the coulomb pressure with abar
! dpcouldz = derivative of the coulomb pressure with zbar

! ecoul    = coulomb energy correction
! decouldd = derivative of the coulomb energy with density
! decouldt = derivative of the coulomb energy with temperature
! decoulda = derivative of the coulomb energy with abar
! decouldz = derivative of the coulomb energy with zbar

! scoul    = coulomb entropy correction
! dscouldd = derivative of the coulomb entropy with density
! dscouldt = derivative of the coulomb entropy with temperature
! dscoulda = derivative of the coulomb entropy with abar
! dscouldz = derivative of the coulomb entropy with zbar


! kt       = kerg * temperature
! beta     = dimensionless ratio of kerg*temp/me*c^2

! chit     = temperature exponent in the pressure equation of state
! chid     = density exponent in the pressure equation of state
! cv       = specific heat at constant volume
! cp       = specific heat at constant pressure
! gam1     = first adiabatic exponent
! gam2     = second adiabatic exponent
! gam3     = third adiabatic exponent
! nabad    = adiabatic gradient
! sound    = relativistically correct adiabatic sound speed
! plasg    = ratio of electrostatic to thermal energy


! dse      = thermodynamic consistency check de/dt = t*ds/dt
! dpe      = thermodynamic consistency check p = d**2 de/dd + t*dpdt
! dsp      = thermodynamic consistency check dp/dt = - d**2 ds/dd







! declare the input
      double precision temp,den,zbar,abar


! declare local variables
! totals
      double precision pres, &
                       dpresdd,dpresdt,dpresda,dpresdz, &
                       dpresddd,dpresddt,dpresdda,dpresddz, &
                       dpresdtt,dpresdta,dpresdtz, &
                       dpresdaa,dpresdaz,dpresdzz, &
                       ener, &
                       denerdd,denerdt,denerda,denerdz, &
                       denerddd,denerddt,denerdda,denerddz, &
                       denerdtt,denerdta,denerdtz, &
                       denerdaa,denerdaz,denerdzz, &
                       entr, &
                       dentrdd,dentrdt,dentrda,dentrdz, &
                       dentrddd,dentrddt,dentrdda,dentrddz, &
                       dentrdtt,dentrdta,dentrdtz, &
                       dentrdaa,dentrdaz,dentrdzz


! for the gas
      double precision pgas, &
                       dpgasdd,dpgasdt,dpgasda,dpgasdz, &
                       dpgasddd,dpgasddt,dpgasdda,dpgasddz, &
                       dpgasdtt,dpgasdta,dpgasdtz, &
                       dpgasdaa,dpgasdaz,dpgasdzz, &
                       egas, &
                       degasdd,degasdt,degasda,degasdz, &
                       degasddd,degasddt,degasdda,degasddz, &
                       degasdtt,degasdta,degasdtz, &
                       degasdaa,degasdaz,degasdzz, &
                       sgas, &
                       dsgasdd,dsgasdt,dsgasda,dsgasdz, &
                       dsgasddd,dsgasddt,dsgasdda,dsgasddz, &
                       dsgasdtt,dsgasdta,dsgasdtz, &
                       dsgasdaa,dsgasdaz,dsgasdzz


! radiation
      integer          radmult
      double precision prad, &
                       dpraddd,dpraddt,dpradda,dpraddz, &
                       dpradddd,dpradddt,dpraddda,dpradddz, &
                       dpraddtt,dpraddta,dpraddtz, &
                       dpraddaa,dpraddaz,dpraddzz, &
                       erad, &
                       deraddd,deraddt,deradda,deraddz, &
                       deradddd,deradddt,deraddda,deradddz, &
                       deraddtt,deraddta,deraddtz, &
                       deraddaa,deraddaz,deraddzz, &
                       srad, &
                       dsraddd,dsraddt,dsradda,dsraddz, &
                       dsradddd,dsradddt,dsraddda,dsradddz, &
                       dsraddtt,dsraddta,dsraddtz, &
                       dsraddaa,dsraddaz,dsraddzz


! ions; now in xniroot_common.dek
      integer          ionmult
      double precision eta_ion_try


! electron-positrons; now in the file xneroot_common.dek
      integer          elemult

! ionization contributions; now in the file xneroot_common.dek
      integer          ionized,potmult


! coulomb corrections; now in the file coulomb_common.dek
      integer          coulmult

      double precision s,sinv,dsdd,dsdt,dsda,dsdz, &
                       dsddd,dsddt,dsdda,dsddz,dsdtt,dsdta,dsdtz, &
                       dsdaa,dsdaz,dsdzz

      double precision aele,aeleinv,daeledd,daeledt,daeleda,daeledz, &
                       daeleddd,daeleddt,daeledda,daeleddz,daeledtt, &
                       daeledta,daeledtz,daeledaa,daeledaz,daeledzz

      double precision eplasg, &
                       deplasgdd,deplasgdt,deplasgda,deplasgdz, &
                       deplasgddd,deplasgddt,deplasgdda,deplasgddz, &
                       deplasgdtt,deplasgdta,deplasgdtz,deplasgdaa, &
                       deplasgdaz,deplasgdzz

      double precision &
                       dplasgdd,dplasgdt,dplasgda,dplasgdz, &
                       dplasgddd,dplasgddt,dplasgdda,dplasgddz, &
                       dplasgdtt,dplasgdta,dplasgdtz,dplasgdaa, &
                       dplasgdaz,dplasgdzz

      double precision u0,du0,ddu0,p1,p2,p3,p4,p5,p6,ion_radius
      double precision a1,b1,c1,d1,e1,a2,b2,c2
      parameter        (a1 = -0.898004d0, &
                        b1 =  0.96786d0, &
                        c1 =  0.220703d0, &
                        d1 = -0.86097d0, &
                        e1 =  2.5269d0, &
                        a2 =  0.29561d0, &
                        b2 =  1.9885d0, &
                        c2 =  0.288675d0)



! various physical quantities based on derivatives
      double precision &
                       chit, &
                       dchitdd,dchitdt,dchitda,dchitdz, &
                       chid, &
                       dchiddd,dchiddt,dchidda,dchiddz, &
                       cv, &
                       dcvdd,dcvdt,dcvda,dcvdz, &
                       cp, &
                       dcpdd,dcpdt,dcpda,dcpdz, &
                       gam1, &
                       dgam1dd,dgam1dt,dgam1da,dgam1dz, &
                       gam2, &
                       dgam2dd,dgam2dt,dgam2da,dgam2dz, &
                       gam3, &
                       dgam3dd,dgam3dt,dgam3da,dgam3dz, &
                       nabad, &
                       dnabdd,dnabdt,dnabda,dnabdz, &
                       sound, &
                       dcsdd,dcsdt,dcsda,dcsdz

      double precision &
                       chit_gas, &
                       dchit_gasdd,dchit_gasdt,dchit_gasda,dchit_gasdz, &
                       chid_gas, &
                       dchid_gasdd,dchid_gasdt,dchid_gasda,dchid_gasdz, &
                       cv_gas, &
                       dcv_gasdd,dcv_gasdt,dcv_gasda,dcv_gasdz, &
                       cp_gas, &
                       dcp_gasdd,dcp_gasdt,dcp_gasda,dcp_gasdz, &
                       gam1_gas, &
                       dgam1_gasdd,dgam1_gasdt,dgam1_gasda,dgam1_gasdz, &
                       gam2_gas, &
                       dgam2_gasdd,dgam2_gasdt,dgam2_gasda,dgam2_gasdz, &
                       gam3_gas, &
                       dgam3_gasdd,dgam3_gasdt,dgam3_gasda,dgam3_gasdz, &
                       nabad_gas, &
                       dnab_gasdd,dnab_gasdt,dnab_gasda,dnab_gasdz, &
                       sound_gas, &
                       dcs_gasdd,dcs_gasdt,dcs_gasda,dcs_gasdz


! for the maxwell relations
      double precision dse,dpe,dsp


! miscelaneous local variables
      integer          i,j,k,kend,niter,mode
      double precision kt,ktinv,x,y,z,ww,xx,yy,zz,zzi,ytot1,ye, &
                       ages,agesav,agesnew,ratio,fk,dfk, &
                       deninv,tempinv,presinv,plasginv,zbarxx


! various derived constants
! f90 allows expressions 
      double precision third,sioncon,sifac,kergavo,asoli3, &
                       clight2,eostol,fpmin
      parameter        (third  = 1.0d0/3.0d0, &
                        sioncon = (2.0d0 * pi * amu * kerg)/(h*h), &
                        sifac  =  h**3/(2.0d0 * pi* amu)**1.5d0, &
                        kergavo = kerg * avo, &
                        asoli3  = asol/3.0d0, &
                        clight2 = clight*clight, &
                        eostol = 1.0d-13, &
                        fpmin  = 1.0d-14)

      double precision forth,fiveth,teninth,esqu,forthpi
      parameter        (forth   = 4.0d0/3.0d0, &
                        fiveth  = 5.0d0/3.0d0, &
                        teninth = 10.0d0/9.0d0, &
                        esqu    = qe*qe, &
                        forthpi = forth * pi)




! bring in common block variables

! common block communication with routine xneroot

! electron-positrons

      double precision etaele,detadd,detadt,detada,detadz, &
                       detaddd,detaddt,detadda,detaddz,detadtt, &
                       detadta,detadtz,detadaa,detadaz,detadzz

      common /xneta/   etaele,detadd,detadt,detada,detadz, &
                       detaddd,detaddt,detadda,detaddz,detadtt, &
                       detadta,detadtz,detadaa,detadaz,detadzz


      double precision &
                       pep,dpepdd,dpepdt,dpepda,dpepdz, &
                       dpepddd,dpepddt,dpepdda,dpepddz, &
                       dpepdtt,dpepdta,dpepdtz,dpepdaa, &
                       dpepdaz,dpepdzz, &
                       eep,deepdd,deepdt,deepda,deepdz, &
                       deepddd,deepddt,deepdda,deepddz, &
                       deepdtt,deepdta,deepdtz,deepdaa, &
                       deepdaz,deepdzz, &
                       sep,dsepdd,dsepdt,dsepda,dsepdz, &
                       dsepddd,dsepddt,dsepdda,dsepddz, &
                       dsepdtt,dsepdta,dsepdtz,dsepdaa, &
                       dsepdaz,dsepdzz

      common /epc1/ &
                       pep,dpepdd,dpepdt,dpepda,dpepdz, &
                       dpepddd,dpepddt,dpepdda,dpepddz, &
                       dpepdtt,dpepdta,dpepdtz,dpepdaa, &
                       dpepdaz,dpepdzz, &
                       eep,deepdd,deepdt,deepda,deepdz, &
                       deepddd,deepddt,deepdda,deepddz, &
                       deepdtt,deepdta,deepdtz,deepdaa, &
                       deepdaz,deepdzz, &
                       sep,dsepdd,dsepdt,dsepda,dsepdz, &
                       dsepddd,dsepddt,dsepdda,dsepddz, &
                       dsepdtt,dsepdta,dsepdtz,dsepdaa, &
                       dsepdaz,dsepdzz


      double precision &
                       etapos,zeff
      common /xnec1/ &
                       etapos,zeff


      double precision &
                       pele,dpeledd,dpeledt,dpeleda,dpeledz, &
                       dpeleddd,dpeleddt,dpeledda,dpeleddz, &
                       dpeledtt,dpeledta,dpeledtz,dpeledaa, &
                       dpeledaz,dpeledzz, &
                       eele,deeledd,deeledt,deeleda,deeledz, &
                       deeleddd,deeleddt,deeledda,deeleddz, &
                       deeledtt,deeledta,deeledtz,deeledaa, &
                       deeledaz,deeledzz, &
                       sele,dseledd,dseledt,dseleda,dseledz, &
                       dseleddd,dseleddt,dseledda,dseleddz, &
                       dseledtt,dseledta,dseledtz,dseledaa, &
                       dseledaz,dseledzz
      common /eleth1/ &
                       pele,dpeledd,dpeledt,dpeleda,dpeledz, &
                       dpeleddd,dpeleddt,dpeledda,dpeleddz, &
                       dpeledtt,dpeledta,dpeledtz,dpeledaa, &
                       dpeledaz,dpeledzz, &
                       eele,deeledd,deeledt,deeleda,deeledz, &
                       deeleddd,deeleddt,deeledda,deeleddz, &
                       deeledtt,deeledta,deeledtz,deeledaa, &
                       deeledaz,deeledzz, &
                       sele,dseledd,dseledt,dseleda,dseledz, &
                       dseleddd,dseleddt,dseledda,dseleddz, &
                       dseledtt,dseledta,dseledtz,dseledaa, &
                       dseledaz,dseledzz

      double precision &
                       ppos,dpposdd,dpposdt,dpposda,dpposdz, &
                       dpposddd,dpposddt,dpposdda,dpposddz, &
                       dpposdtt,dpposdta,dpposdtz,dpposdaa, &
                       dpposdaz,dpposdzz, &
                       epos,deposdd,deposdt,deposda,deposdz, &
                       deposddd,deposddt,deposdda,deposddz, &
                       deposdtt,deposdta,deposdtz,deposdaa, &
                       deposdaz,deposdzz, &
                       spos,dsposdd,dsposdt,dsposda,dsposdz, &
                       dsposddd,dsposddt,dsposdda,dsposddz, &
                       dsposdtt,dsposdta,dsposdtz,dsposdaa, &
                       dsposdaz,dsposdzz
      common /posth1/ &
                       ppos,dpposdd,dpposdt,dpposda,dpposdz, &
                       dpposddd,dpposddt,dpposdda,dpposddz, &
                       dpposdtt,dpposdta,dpposdtz,dpposdaa, &
                       dpposdaz,dpposdzz, &
                       epos,deposdd,deposdt,deposda,deposdz, &
                       deposddd,deposddt,deposdda,deposddz, &
                       deposdtt,deposdta,deposdtz,deposdaa, &
                       deposdaz,deposdzz, &
                       spos,dsposdd,dsposdt,dsposda,dsposdz, &
                       dsposddd,dsposddt,dsposdda,dsposddz, &
                       dsposdtt,dsposdta,dsposdtz,dsposdaa, &
                       dsposdaz,dsposdzz


      double precision xne, &
                       dxnedd,dxnedt,dxneda,dxnedz, &
                       dxneddd,dxneddt,dxnedda,dxneddz, &
                       dxnedtt,dxnedta,dxnedtz,dxnedaa, &
                       dxnedaz,dxnedzz

      common /xnec2/   xne, &
                       dxnedd,dxnedt,dxneda,dxnedz, &
                       dxneddd,dxneddt,dxnedda,dxneddz, &
                       dxnedtt,dxnedta,dxnedtz,dxnedaa, &
                       dxnedaz,dxnedzz

      double precision &
                       xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz, &
                       dxneferddd,dxneferddt,dxneferdda,dxneferddz, &
                       dxneferdtt,dxneferdta,dxneferdtz,dxneferdaa, &
                       dxneferdaz,dxneferdzz, &
                       xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz, &
                       dxnpferddd,dxnpferddt,dxnpferdda,dxnpferddz, &
                       dxnpferdtt,dxnpferdta,dxnpferdtz,dxnpferdaa, &
                       dxnpferdaz,dxnpferdzz
      common /xnec3/ &
                       xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz, &
                       dxneferddd,dxneferddt,dxneferdda,dxneferddz, &
                       dxneferdtt,dxneferdta,dxneferdtz,dxneferdaa, &
                       dxneferdaz,dxneferdzz, &
                       xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz, &
                       dxnpferddd,dxnpferddt,dxnpferdda,dxnpferddz, &
                       dxnpferdtt,dxnpferdta,dxnpferdtz,dxnpferdaa, &
                       dxnpferdaz,dxnpferdzz





! ionization contributions
      double precision eip,deipdd,deipdt,deipda,deipdz, &
                       sip,dsipdd,dsipdt,dsipda,dsipdz, &
                       pip,dpipdd,dpipdt,dpipda,dpipdz

      common /xnec4/   eip,deipdd,deipdt,deipda,deipdz, &
                       sip,dsipdd,dsipdt,dsipda,dsipdz, &
                       pip,dpipdd,dpipdt,dpipda,dpipdz




! common block communication with routine xniroot

! ions

! pressure, energy, entropy
      double precision pion,eion,sion, &
                       dpiondd,dpiondt,dpionda,dpiondz, &
                       deiondd,deiondt,deionda,deiondz, &
                       dsiondd,dsiondt,dsionda,dsiondz, &
                       dpionddd,dpionddt,dpiondda,dpionddz, &
                       dpiondtt,dpiondta,dpiondtz, &
                       dpiondaa,dpiondaz,dpiondzz, &
                       deionddd,deionddt,deiondda,deionddz, &
                       deiondtt,deiondta,deiondtz, &
                       deiondaa,deiondaz,deiondzz, &
                       dsionddd,dsionddt,dsiondda,dsionddz, &
                       dsiondtt,dsiondta,dsiondtz, &
                       dsiondaa,dsiondaz,dsiondzz

      common /ionc1/   pion,eion,sion, &
                       dpiondd,dpiondt,dpionda,dpiondz, &
                       deiondd,deiondt,deionda,deiondz, &
                       dsiondd,dsiondt,dsionda,dsiondz, &
                       dpionddd,dpionddt,dpiondda,dpionddz, &
                       dpiondtt,dpiondta,dpiondtz, &
                       dpiondaa,dpiondaz,dpiondzz, &
                       deionddd,deionddt,deiondda,deionddz, &
                       deiondtt,deiondta,deiondtz, &
                       deiondaa,deiondaz,deiondzz, &
                       dsionddd,dsionddt,dsiondda,dsionddz, &
                       dsiondtt,dsiondta,dsiondtz, &
                       dsiondaa,dsiondaz,dsiondzz

! number densities
      double precision xni, &
                       dxnidd,dxnidt,dxnida,dxnidz, &
                       dxniddd,dxniddt,dxnidda,dxniddz, &
                       dxnidtt,dxnidta,dxnidtz, &
                       dxnidaa,dxnidaz,dxnidzz

      common /ionc2/   xni, &
                       dxnidd,dxnidt,dxnida,dxnidz, &
                       dxniddd,dxniddt,dxnidda,dxniddz, &
                       dxnidtt,dxnidta,dxnidtz, &
                       dxnidaa,dxnidaz,dxnidzz


      double precision &
                       xnifer,dxniferdd,dxniferdt,dxniferda,dxniferdz, &
                       dxniferddd,dxniferddt,dxniferdda,dxniferddz, &
                       dxniferdtt,dxniferdta,dxniferdtz,dxniferdaa, &
                       dxniferdaz,dxniferdzz
      common /xnic3/ &
                       xnifer,dxniferdd,dxniferdt,dxniferda,dxniferdz, &
                       dxniferddd,dxniferddt,dxniferdda,dxniferddz, &
                       dxniferdtt,dxniferdta,dxniferdtz,dxniferdaa, &
                       dxniferdaz,dxniferdzz


! chemical potential
      double precision etaion, &
                       detaidd,detaidt,detaida,detaidz, &
                       detaiddd,detaiddt,detaidda,detaiddz, &
                       detaidtt,detaidta,detaidtz, &
                       detaidaa,detaidaz,detaidzz

      common /ionc3/   etaion, &
                       detaidd,detaidt,detaida,detaidz, &
                       detaiddd,detaiddt,detaidda,detaiddz, &
                       detaidtt,detaidta,detaidtz, &
                       detaidaa,detaidaz,detaidzz



! common block communication with routine coulomb

! coulomb corrections

! pressure, energy, entropy

      double precision plasg
      common /coulc1/  plasg

      double precision pcoul, &
                       dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                       dpcoulddd,dpcoulddt,dpcouldda,dpcoulddz, &
                       dpcouldtt,dpcouldta,dpcouldtz,dpcouldaa, &
                       dpcouldaz,dpcouldzz

      common /coulc2/  pcoul, &
                       dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                       dpcoulddd,dpcoulddt,dpcouldda,dpcoulddz, &
                       dpcouldtt,dpcouldta,dpcouldtz,dpcouldaa, &
                       dpcouldaz,dpcouldzz


      double precision ecoul, &
                       decouldd,decouldt,decoulda,decouldz, &
                       decoulddd,decoulddt,decouldda,decoulddz, &
                       decouldtt,decouldta,decouldtz,decouldaa, &
                       decouldaz,decouldzz

      common /coulc4/  ecoul, &
                       decouldd,decouldt,decoulda,decouldz, &
                       decoulddd,decoulddt,decouldda,decoulddz, &
                       decouldtt,decouldta,decouldtz,decouldaa, &
                       decouldaz,decouldzz


      double precision scoul, &
                       dscouldd,dscouldt,dscoulda,dscouldz, &
                       dscoulddd,dscoulddt,dscouldda,dscoulddz, &
                       dscouldtt,dscouldta,dscouldtz,dscouldaa, &
                       dscouldaz,dscouldzz

      common /coulc4/  scoul, &
                       dscouldd,dscouldt,dscoulda,dscouldz, &
                       dscoulddd,dscoulddt,dscouldda,dscoulddz, &
                       dscouldtt,dscouldta,dscouldtz,dscouldaa, &
                       dscouldaz,dscouldzz


! popular format statements for debugging
01    format(1x,5(a,1pe24.16))
02    format(1x,5(a,1pe16.8))
03    format(1x,1p5e16.8)



! set the on/off switches
      radmult  = 1
      ionmult  = 1
      ionized  = 1
      elemult  = 1
      coulmult = 1
      potmult  = 0


! start pipeline loop
      do j=jlo_eos,jhi_eos

       if ((temp_row(j) .le. 0.0) .or. (den_row(j) .le. 0.0)) then
        call zero_eos_vector(j,j)

       else
        temp  = temp_row(j)
        den   = den_row(j)
        abar  = abar_row(j)
        zbar  = zbar_row(j)
        ytot1 = 1.0d0/abar
 

! initialize local variables to zero

       prad     = 0.0d0
       dpraddd  = 0.0d0
       dpraddt  = 0.0d0
       dpradda  = 0.0d0
       dpraddz  = 0.0d0
       dpradddd = 0.0d0
       dpradddt = 0.0d0
       dpraddda = 0.0d0
       dpradddz = 0.0d0
       dpraddtt = 0.0d0
       dpraddta = 0.0d0
       dpraddtz = 0.0d0
       dpraddaa = 0.0d0
       dpraddaz = 0.0d0
       dpraddzz = 0.0d0

       erad     = 0.0d0
       deraddd  = 0.0d0
       deraddt  = 0.0d0
       deradda  = 0.0d0
       deraddz  = 0.0d0
       deradddd = 0.0d0
       deradddt = 0.0d0
       deraddda = 0.0d0
       deradddz = 0.0d0
       deraddtt = 0.0d0
       deraddta = 0.0d0
       deraddtz = 0.0d0
       deraddaa = 0.0d0
       deraddaz = 0.0d0
       deraddzz = 0.0d0

       srad     = 0.0d0
       dsraddd  = 0.0d0
       dsraddt  = 0.0d0
       dsradda  = 0.0d0
       dsraddz  = 0.0d0
       dsradddd = 0.0d0
       dsradddt = 0.0d0
       dsraddda = 0.0d0
       dsradddz = 0.0d0
       dsraddtt = 0.0d0
       dsraddta = 0.0d0
       dsraddtz = 0.0d0
       dsraddaa = 0.0d0
       dsraddaz = 0.0d0
       dsraddzz = 0.0d0

       xni     = 0.0d0
       dxnidd  = 0.0d0
       dxnidt  = 0.0d0
       dxnida  = 0.0d0
       dxnidz  = 0.0d0
       dxniddd = 0.0d0
       dxniddt = 0.0d0
       dxnidda = 0.0d0
       dxniddz = 0.0d0
       dxnidtt = 0.0d0
       dxnidta = 0.0d0
       dxnidtz = 0.0d0
       dxnidaa = 0.0d0
       dxnidaz = 0.0d0
       dxnidzz = 0.0d0

       pion     = 0.0d0
       dpiondd  = 0.0d0
       dpiondt  = 0.0d0
       dpionda  = 0.0d0
       dpiondz  = 0.0d0
       dpionddd = 0.0d0
       dpionddt = 0.0d0
       dpiondda = 0.0d0
       dpionddz = 0.0d0
       dpiondtt = 0.0d0
       dpiondta = 0.0d0
       dpiondtz = 0.0d0
       dpiondaa = 0.0d0
       dpiondaz = 0.0d0
       dpiondzz = 0.0d0

       eion     = 0.0d0
       deiondd  = 0.0d0
       deiondt  = 0.0d0
       deionda  = 0.0d0
       deiondz  = 0.0d0
       deionddd = 0.0d0
       deionddt = 0.0d0
       deiondda = 0.0d0
       deionddz = 0.0d0
       deiondtt = 0.0d0
       deiondta = 0.0d0
       deiondtz = 0.0d0
       deiondaa = 0.0d0
       deiondaz = 0.0d0
       deiondzz = 0.0d0

       sion     = 0.0d0
       dsiondd  = 0.0d0
       dsiondt  = 0.0d0
       dsionda  = 0.0d0
       dsiondz  = 0.0d0
       dsionddd = 0.0d0
       dsionddt = 0.0d0
       dsiondda = 0.0d0
       dsionddz = 0.0d0
       dsiondtt = 0.0d0
       dsiondta = 0.0d0
       dsiondtz = 0.0d0
       dsiondaa = 0.0d0
       dsiondaz = 0.0d0
       dsiondzz = 0.0d0

       etaion   = 0.0d0
       detaidd  = 0.0d0
       detaidt  = 0.0d0
       detaida  = 0.0d0
       detaidz  = 0.0d0
       detaiddd = 0.0d0
       detaiddt = 0.0d0
       detaidda = 0.0d0
       detaiddz = 0.0d0
       detaidtt = 0.0d0
       detaidta = 0.0d0
       detaidtz = 0.0d0
       detaidaa = 0.0d0
       detaidaz = 0.0d0
       detaidzz = 0.0d0

       xne     = 0.0d0
       dxnedd  = 0.0d0
       dxnedt  = 0.0d0
       dxneda  = 0.0d0
       dxnedz  = 0.0d0
       dxneddd = 0.0d0
       dxneddt = 0.0d0
       dxnedda = 0.0d0
       dxneddz = 0.0d0
       dxnedtt = 0.0d0
       dxnedta = 0.0d0
       dxnedtz = 0.0d0
       dxnedaa = 0.0d0
       dxnedaz = 0.0d0
       dxnedzz = 0.0d0


       etaele  = 0.0d0
       detadd  = 0.0d0
       detadt  = 0.0d0
       detada  = 0.0d0
       detadz  = 0.0d0
       detaddd = 0.0d0
       detaddt = 0.0d0
       detadda = 0.0d0
       detaddz = 0.0d0
       detadtt = 0.0d0
       detadta = 0.0d0
       detadtz = 0.0d0
       detadaa = 0.0d0
       detadaz = 0.0d0
       detadzz = 0.0d0

       etapos   = 0.0d0

       xnefer     = 0.0d0
       dxneferdd  = 0.0d0
       dxneferdt  = 0.0d0
       dxneferda  = 0.0d0
       dxneferdz  = 0.0d0
       dxneferddd = 0.0d0
       dxneferddt = 0.0d0
       dxneferdda = 0.0d0
       dxneferddz = 0.0d0
       dxneferdtt = 0.0d0
       dxneferdta = 0.0d0
       dxneferdtz = 0.0d0
       dxneferdaa = 0.0d0
       dxneferdaz = 0.0d0
       dxneferdzz = 0.0d0

       xnpfer     = 0.0d0
       dxnpferdd  = 0.0d0
       dxnpferdt  = 0.0d0
       dxnpferda  = 0.0d0
       dxnpferdz  = 0.0d0
       dxnpferddd = 0.0d0
       dxnpferddt = 0.0d0
       dxnpferdda = 0.0d0
       dxnpferddz = 0.0d0
       dxnpferdtt = 0.0d0
       dxnpferdta = 0.0d0
       dxnpferdtz = 0.0d0
       dxnpferdaa = 0.0d0
       dxnpferdaz = 0.0d0
       dxnpferdzz = 0.0d0

       pele     = 0.0d0
       dpeledd  = 0.0d0
       dpeledt  = 0.0d0
       dpeleda  = 0.0d0
       dpeledz  = 0.0d0
       dpeleddd = 0.0d0
       dpeleddt = 0.0d0
       dpeledda = 0.0d0
       dpeleddz = 0.0d0
       dpeledtt = 0.0d0
       dpeledta = 0.0d0
       dpeledtz = 0.0d0
       dpeledaa = 0.0d0
       dpeledaz = 0.0d0
       dpeledzz = 0.0d0

       eele     = 0.0d0
       deeledd  = 0.0d0
       deeledt  = 0.0d0
       deeleda  = 0.0d0
       deeledz  = 0.0d0
       deeleddd = 0.0d0
       deeleddt = 0.0d0
       deeledda = 0.0d0
       deeleddz = 0.0d0
       deeledtt = 0.0d0
       deeledta = 0.0d0
       deeledtz = 0.0d0
       deeledaa = 0.0d0
       deeledaz = 0.0d0
       deeledzz = 0.0d0

       sele     = 0.0d0
       dseledd  = 0.0d0
       dseledt  = 0.0d0
       dseleda  = 0.0d0
       dseledz  = 0.0d0
       dseleddd = 0.0d0
       dseleddt = 0.0d0
       dseledda = 0.0d0
       dseleddz = 0.0d0
       dseledtt = 0.0d0
       dseledta = 0.0d0
       dseledtz = 0.0d0
       dseledaa = 0.0d0
       dseledaz = 0.0d0
       dseledzz = 0.0d0

       ppos     = 0.0d0
       dpposdd  = 0.0d0
       dpposdt  = 0.0d0
       dpposda  = 0.0d0
       dpeledz  = 0.0d0
       dpposddd = 0.0d0
       dpposddt = 0.0d0
       dpposdda = 0.0d0
       dpposddz = 0.0d0
       dpposdtt = 0.0d0
       dpposdta = 0.0d0
       dpposdtz = 0.0d0
       dpposdaa = 0.0d0
       dpposdaz = 0.0d0
       dpposdzz = 0.0d0

       epos     = 0.0d0
       deposdd  = 0.0d0
       deposdt  = 0.0d0
       deposda  = 0.0d0
       deeledz  = 0.0d0
       deposddd = 0.0d0
       deposddt = 0.0d0
       deposdda = 0.0d0
       deposddz = 0.0d0
       deposdtt = 0.0d0
       deposdta = 0.0d0
       deposdtz = 0.0d0
       deposdaa = 0.0d0
       deposdaz = 0.0d0
       deposdzz = 0.0d0

       spos     = 0.0d0
       dsposdd  = 0.0d0
       dsposdt  = 0.0d0
       dsposda  = 0.0d0
       dseledz  = 0.0d0
       dseleddd = 0.0d0
       dseleddt = 0.0d0
       dseledda = 0.0d0
       dseleddz = 0.0d0
       dseledtt = 0.0d0
       dseledta = 0.0d0
       dseledtz = 0.0d0
       dseledaa = 0.0d0
       dseledaz = 0.0d0
       dseledzz = 0.0d0

       pep     = 0.0d0
       dpepdd  = 0.0d0
       dpepdt  = 0.0d0
       dpepda  = 0.0d0
       dpepdz  = 0.0d0
       dpepddd = 0.0d0
       dpepddt = 0.0d0
       dpepdda = 0.0d0
       dpepddz = 0.0d0
       dpepdtt = 0.0d0
       dpepdta = 0.0d0
       dpepdtz = 0.0d0
       dpepdaa = 0.0d0
       dpepdaz = 0.0d0
       dpepdzz = 0.0d0

       eep     = 0.0d0
       deepdd  = 0.0d0
       deepdt  = 0.0d0
       deepda  = 0.0d0
       deepdz  = 0.0d0
       deepddd = 0.0d0
       deepddt = 0.0d0
       deepdda = 0.0d0
       deepddz = 0.0d0
       deepdtt = 0.0d0
       deepdta = 0.0d0
       deepdtz = 0.0d0
       deepdaa = 0.0d0
       deepdaz = 0.0d0
       deepdzz = 0.0d0

       sep     = 0.0d0
       dsepdd  = 0.0d0
       dsepdt  = 0.0d0
       dsepda  = 0.0d0
       dsepdz  = 0.0d0
       dsepddd = 0.0d0
       dsepddt = 0.0d0
       dsepdda = 0.0d0
       dsepddz = 0.0d0
       dsepdtt = 0.0d0
       dsepdta = 0.0d0
       dsepdtz = 0.0d0
       dsepdaa = 0.0d0
       dsepdaz = 0.0d0
       dsepdzz = 0.0d0

       eip      = 0.0d0
       deipdd   = 0.0d0
       deipdt   = 0.0d0
       deipda   = 0.0d0
       deipdz   = 0.0d0

       sip      = 0.0d0
       dsipdd   = 0.0d0
       dsipdt   = 0.0d0
       dsipda   = 0.0d0
       dsipdz   = 0.0d0

       pcoul     = 0.0d0
       dpcouldd  = 0.0d0
       dpcouldt  = 0.0d0
       dpcoulda  = 0.0d0
       dpcouldz  = 0.0d0
       dpcoulddd = 0.0d0
       dpcoulddt = 0.0d0
       dpcouldda = 0.0d0
       dpcoulddz = 0.0d0
       dpcouldtt = 0.0d0
       dpcouldta = 0.0d0
       dpcouldtz = 0.0d0
       dpcouldaa = 0.0d0
       dpcouldaz = 0.0d0
       dpcouldzz = 0.0d0

       ecoul     = 0.0d0
       decouldd  = 0.0d0
       decouldt  = 0.0d0
       decoulda  = 0.0d0
       decouldz  = 0.0d0
       decoulddd = 0.0d0
       decoulddt = 0.0d0
       decouldda = 0.0d0
       decoulddz = 0.0d0
       decouldtt = 0.0d0
       decouldta = 0.0d0
       decouldtz = 0.0d0
       decouldaa = 0.0d0
       decouldaz = 0.0d0
       decouldzz = 0.0d0

       scoul     = 0.0d0
       dscouldd  = 0.0d0
       dscouldt  = 0.0d0
       dscoulda  = 0.0d0
       dscouldz  = 0.0d0
       dscoulddd = 0.0d0
       dscoulddt = 0.0d0
       dscouldda = 0.0d0
       dscoulddz = 0.0d0
       dscouldtt = 0.0d0
       dscouldta = 0.0d0
       dscouldtz = 0.0d0
       dscouldaa = 0.0d0
       dscouldaz = 0.0d0
       dscouldzz = 0.0d0

       chit     = 0.0d0
       dchitdd  = 0.0d0
       dchitdt  = 0.0d0
       dchitda  = 0.0d0
       dchitdz  = 0.0d0

       chid     = 0.0d0
       dchiddd  = 0.0d0
       dchiddt  = 0.0d0
       dchidda  = 0.0d0
       dchiddz  = 0.0d0

       cv     = 0.0d0
       dcvdd  = 0.0d0
       dcvdt  = 0.0d0
       dcvda  = 0.0d0
       dcvdz  = 0.0d0

       cp     = 0.0d0
       dcpdd  = 0.0d0
       dcpdt  = 0.0d0
       dcpda  = 0.0d0
       dcpdz  = 0.0d0

       gam1     = 0.0d0
       dgam1dd  = 0.0d0
       dgam1dt  = 0.0d0
       dgam1da  = 0.0d0
       dgam1dz  = 0.0d0

       gam2     = 0.0d0
       dgam2dd  = 0.0d0
       dgam2dt  = 0.0d0
       dgam2da  = 0.0d0
       dgam2dz  = 0.0d0

       gam3     = 0.0d0
       dgam3dd  = 0.0d0
       dgam3dt  = 0.0d0
       dgam3da  = 0.0d0
       dgam3dz  = 0.0d0

       nabad     = 0.0d0
       dnabdd  = 0.0d0
       dnabdt  = 0.0d0
       dnabda  = 0.0d0
       dnabdz  = 0.0d0

       sound  = 0.0d0
       dcsdd  = 0.0d0
       dcsdt  = 0.0d0
       dcsda  = 0.0d0
       dcsdz  = 0.0d0

       chit_gas     = 0.0d0
       dchit_gasdd  = 0.0d0
       dchit_gasdt  = 0.0d0
       dchit_gasda  = 0.0d0
       dchit_gasdz  = 0.0d0

       chid_gas     = 0.0d0
       dchid_gasdd  = 0.0d0
       dchid_gasdt  = 0.0d0
       dchid_gasda  = 0.0d0
       dchid_gasdz  = 0.0d0

       cv_gas     = 0.0d0
       dcv_gasdd  = 0.0d0
       dcv_gasdt  = 0.0d0
       dcv_gasda  = 0.0d0
       dcv_gasdz  = 0.0d0

       cp_gas     = 0.0d0
       dcp_gasdd  = 0.0d0
       dcp_gasdt  = 0.0d0
       dcp_gasda  = 0.0d0
       dcp_gasdz  = 0.0d0

       gam1_gas     = 0.0d0
       dgam1_gasdd  = 0.0d0
       dgam1_gasdt  = 0.0d0
       dgam1_gasda  = 0.0d0
       dgam1_gasdz  = 0.0d0

       gam2_gas     = 0.0d0
       dgam2_gasdd  = 0.0d0
       dgam2_gasdt  = 0.0d0
       dgam2_gasda  = 0.0d0
       dgam2_gasdz  = 0.0d0

       gam3_gas     = 0.0d0
       dgam3_gasdd  = 0.0d0
       dgam3_gasdt  = 0.0d0
       dgam3_gasda  = 0.0d0
       dgam3_gasdz  = 0.0d0

       nabad_gas     = 0.0d0
       dnab_gasdd  = 0.0d0
       dnab_gasdt  = 0.0d0
       dnab_gasda  = 0.0d0
       dnab_gasdz  = 0.0d0

       sound_gas  = 0.0d0
       dcs_gasdd  = 0.0d0
       dcs_gasdt  = 0.0d0
       dcs_gasda  = 0.0d0
       dcs_gasdz  = 0.0d0

       dse = 0.0d0
       dpe = 0.0d0
       dsp = 0.0d0


! frequent combinations
       kt      = kerg * temp
       ktinv   = 1.0d0/kt
       deninv  = 1.0d0/den
       tempinv = 1.0d0/temp





! radiation section:
       if (radmult .ne. 0) then

! pressure in erg/cm**3
        prad    = asol * third * temp * temp * temp * temp
        dpraddd = 0.0d0
        dpraddt = 4.0d0 * prad/temp
        dpradda = 0.0d0
        dpraddz = 0.0d0

! energy in erg/gr
        erad    = 3.0d0 * prad * deninv
        deraddd = -erad * deninv
        deraddt = 3.0d0 * dpraddt * deninv
        deradda = 0.0d0
        deraddz = 0.0d0

! entropy in erg/g/kelvin
        srad    = (prad*deninv + erad) * tempinv
        dsraddd = (dpraddd*deninv - prad*deninv**2 + deraddd) * tempinv
        dsraddt = (dpraddt*deninv + deraddt - srad)  * tempinv
        dsradda = 0.0d0
        dsraddz = 0.0d0
       end if




! ion section:

! number density in 1/cm**3,
        xni     = avo * ytot1 * den
        dxnidd  = avo * ytot1
        dxnidt  = 0.0d0
        dxnida  = -xni * ytot1
        dxnidz  = 0.0d0

       if (ionmult .ne. 0) then

! pressure in erg/cm**3
        pion    = xni * kt
        dpiondd = dxnidd * kt
        dpiondt = xni * kerg
        dpionda = -pion * ytot1
        dpiondz = 0.0d0

!  energy in erg/gr
        eion    = 1.5d0 * pion*deninv
        deiondd = (1.5d0 * dpiondd - eion)*deninv
        deiondt = 1.5d0 * dpiondt*deninv
        deionda = 1.5d0 * dpionda*deninv
        deiondz = 0.0d0


! ion degeneracy parameter (c&g 9.60)
        y       = 1.0d0/(abar*kt)
        yy      = y * sqrt(y)
        z       = xni * sifac * yy
        etaion  = log(z)

        xx      = 1.0d0/xni
        detaidd = dxnidd*xx
        detaidt = dxnidt*xx - 1.5d0*tempinv
        detaida = dxnida*xx - 1.5d0*ytot1
        detaidz = dxnidz*xx


! entropy in erg/gr/kelvin
! the last term is the usual  etaion * kerg * xni/den
! sometimes called the sacker-tetrode equation

        sion    = (eion + pion*deninv)*tempinv - etaion * kerg*avo*ytot1

        dsiondd = (deiondd + dpiondd*deninv - pion*deninv**2)*tempinv &
                  - detaidd * kerg * avo*ytot1

        dsiondt = (deiondt + dpiondt*deninv)*tempinv &
                  - (eion + pion*deninv)*tempinv**2 &
                  - detaidt * kerg * avo*ytot1

        dsionda = (deionda + dpionda*deninv)*tempinv &
                  - detaida * kerg * avo*ytot1 &
                  + etaion * kerg * avo * ytot1**2

        dsiondz = 0.0d0
       end if






! electron-positron section:
       if (elemult .ne. 0) then

! make a good guess at the the electron degeneracy parameter eta
        call etages(xni,zbar,temp,ages)
        agesav = ages


! newton-raphson to get the electron/positron quantities
        eosfail = .false.
        do i=1,100

         call xneroot(0,den,temp,abar,zbar,ionized,potmult,ages,fk,dfk)

         if (dfk .eq. 0.0) goto 11
         ratio   = fk/dfk
         agesnew = ages - ratio
         z       = abs((agesnew - ages)/ages)
         ages    = agesnew
         niter   = i
         if (z .lt. eostol .or. abs(ratio) .le. fpmin) goto 20
        enddo

11      write(6,*)
        write(6,*) 'newton-raphson failed in routine eosfxt'
        write(6,01) 'temp  =',temp,' den =',den
        write(6,01) 'z     =',z,' ages=',ages, ' agesav=',agesav
        write(6,01) 'eostol=',eostol
        write(6,01) 'f/df  =',fk/dfk,' f   =',fk,    ' df    =',dfk
        write(6,01) 'fpmin =',fpmin
        write(6,*)
        eosfail = .true.
        return
20      continue



! with the converged values, get the energy, pressure, and entropy
         call xneroot(1,den,temp,abar,zbar,ionized,potmult,ages,fk,dfk)
       end if





! coulomb corrections section:
       if (coulmult .ne. 0) then



! this code fragment implments coulomb corrections
! see yakovlev & shalybkov 1989, uniform background corrections

! input:
! den  = density g/cc
! temp = temperature k
! abar = avarage atomic weight
! zbar = avarge charge
! pion and all its derivatives = ion pressure through common block
! xne and all its derivatives = electron number density through common block

! output:
! plasg = ion coupling parameter
! pcoul and all its derivatives = coulomb pressure through common block
! ecoul and all its derivatives = coulomb energy through common block
! scoul and all its derivatives = coulomb entropy through common block



! yakovlev & shalybkov eqs 5, 9 and 10
! use the ion number density instead of the free electron number desnity
! to avoid issues with positrons being included or not (eosfxt vs helmeos)

      y        = forthpi * zbar
      s        = y * xni
      sinv     = 1.0d0/s

! first derivatives
      dsdd     = y * dxnidd
      dsdt     = y * dxnidt
      dsda     = y * dxnida
      dsdz     = y * dxnidz + forthpi*xni

! second derivatives
      dsddd   = y * dxniddd
      dsddt   = y * dxniddt
      dsdda   = y * dxnidda
      dsddz   = y * dxniddz + forthpi*dxnidd
      dsdtt   = y * dxnidtt
      dsdta   = y * dxnidta
      dsdtz   = y * dxnidtz + forthpi*dxnidt
      dsdaa   = y * dxnidaa
      dsdaz   = y * dxnidaz + forthpi*dxnida
      dsdzz   = y * dxnidzz + 2.0d0*forthpi*dxnidz


! electron-sphere radius aele
      aele     = sinv**third
      aeleinv  = 1.0d0/aele
      z        = -third * aele * sinv
      y        = -forth * z * sinv

! first derivatives
      daeledd  = z * dsdd
      daeledt  = z * dsdt
      daeleda  = z * dsda
      daeledz  = z * dsdz

! second derivatives
      daeleddd = y*dsdd*dsdd + z*dsddd
      daeleddt = y*dsdt*dsdd + z*dsddt
      daeledda = y*dsda*dsdd + z*dsdda
      daeleddz = y*dsdz*dsdd + z*dsddz
      daeledtt = y*dsdt*dsdt + z*dsdtt
      daeledta = y*dsda*dsdt + z*dsdta
      daeledtz = y*dsdz*dsdt + z*dsdtz
      daeledaa = y*dsda*dsda + z*dsdaa
      daeledaz = y*dsdz*dsda + z*dsdaz
      daeledzz = y*dsdz*dsdz + z*dsdzz


! electron coupling parameter eplasg
      eplasg   = esqu * ktinv * aeleinv
      z        = -eplasg * aeleinv
      y        = -2.0d0 * z * aeleinv

! first derivatives
      deplasgdd = z * daeledd
      deplasgdt = z * daeledt - eplasg*tempinv
      deplasgda = z * daeleda
      deplasgdz = z * daeledz

! second derivatives
      deplasgddd = y*daeledd*daeledd + z*daeleddd
      deplasgddt = y*daeledt*daeledd - deplasgdd*tempinv + z*daeleddt
      deplasgdda = y*daeleda*daeledd + z*daeledda
      deplasgddz = y*daeledz*daeledd + z*daeleddz
      deplasgdtt = y*daeledt*daeledt + z*daeledtt &
                   + (2.0d0*z*daeledt + 2.0d0*eplasg*tempinv)*tempinv
      deplasgdta = y*daeleda*daeledt + z*daeledta &
                   - deplasgda*tempinv
      deplasgdtz = y*daeledz*daeledt + z*daeledtz &
                   - deplasgdz*tempinv
      deplasgdaa = y*daeleda*daeleda + z*daeledaa
      deplasgdaz = y*daeledz*daeleda + z*daeledaz
      deplasgdzz = y*daeledz*daeledz + z*daeledzz


! ion-sphere radius aion
      x          = zbar**third
      z          = x*x*x*x*x
      ww         = fiveth * x * x
      ion_radius = x * aele


! ion coupling parameter plasg
      plasg    = z * eplasg
      plasginv = 1.0d0/plasg


! first derivatives
      dplasgdd  = z * deplasgdd
      dplasgdt  = z * deplasgdt
      dplasgda  = z * deplasgda
      dplasgdz  = z * deplasgdz + ww*eplasg

! second derivatives
      dplasgddd  = z * deplasgddd
      dplasgddt  = z * deplasgddt
      dplasgdda  = z * deplasgdda
      dplasgddz  = z * deplasgddz + ww*deplasgdd
      dplasgdtt  = z * deplasgdtt
      dplasgdta  = z * deplasgdta
      dplasgdtz  = z * deplasgdtz + ww*deplasgdt
      dplasgdaa  = z * deplasgdaa
      dplasgdaz  = z * deplasgdaz + ww*deplasgda
      dplasgdzz  = z * deplasgdzz + 2.0d0*ww*deplasgdz +teninth/x*eplasg


! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
      if (plasg .ge. 1.0) then
       x        = sqrt(sqrt(plasg))
       p1       = x
       p2       = 1.0d0/x
       p3       = p1*plasginv
       p4       = p2*plasginv
       p5       = p3*plasginv
       p6       = p4*plasginv

       u0       = a1*plasg + b1*p1 + c1*p2 + d1
       du0      = a1 + 0.25d0*b1*p3 - 0.25d0*c1*p4
       ddu0     = -0.1875d0*b1*p5 + 0.3125d0*c1*p6


! energy in erg/gr
       z        = pion * deninv
       ecoul    = z * u0

       x        = deninv*u0
       y        = deninv*du0
       ww       = z*du0
       dfk       = z*ddu0

! first derivatives
       decouldd = dpiondd*x - z*x + ww*dplasgdd
       decouldt = dpiondt*x + ww*dplasgdt
       decoulda = dpionda*x + ww*dplasgda
       decouldz = dpiondz*x + ww*dplasgdz

! second derivatives
       decoulddd = dpionddd*x + 2.0d0*dpiondd*(y*dplasgdd - deninv*x) &
                   + z*(2.0d0*deninv*x &
                   + (ddu0*dplasgdd - 2.0d0*y)*dplasgdd + du0*dplasgddd)
       decoulddt = dpionddt*x + dpiondd*y*dplasgdt &
                   - dpiondt*deninv*x  - z*y*dplasgdt &
                   + (dpiondt*y + dfk*dplasgdt)*dplasgdd + ww*dplasgddt
       decouldda = dpiondda*x + dpiondd*y*dplasgda &
                   - dpionda*deninv*x - z*y*dplasgda &
                   + (dpionda*y + dfk*dplasgda)*dplasgdd + ww*dplasgdda
       decoulddz = dpionddz*x + dpiondd*y*dplasgdz &
                   - dpiondz*deninv*x - z*y*dplasgdz &
                   + (dpiondz*y + dfk*dplasgdz)*dplasgdd + ww*dplasgddz
       decouldtt = dpiondtt*x + (2.0d0*dpiondt*y &
                  + dfk*dplasgdt)*dplasgdt + ww*dplasgdtt
       decouldta = dpiondta*x + dpiondt*y*dplasgda &
                  + (dpionda*y + dfk*dplasgda)*dplasgdt + ww*dplasgdta
       decouldtz = dpiondtz*x + dpiondt*y*dplasgdz &
                  + (dpiondz*y + dfk*dplasgdz)*dplasgdt + ww*dplasgdtz
       decouldaa = dpiondaa*x + (2.0d0*dpionda*y &
                  + dfk*dplasgda)*dplasgda + ww*dplasgdaa
       decouldaz = dpiondaz*x + dpionda*y*dplasgdz &
                  + (dpiondz*y + dfk*dplasgdz)*dplasgda + ww*dplasgdaz
       decouldzz = dpiondzz*x + (2.0d0*dpiondz*y &
                  + dfk*dplasgdz)*dplasgdz + ww*dplasgdzz


! pressure in erg/cc
       y        = third * den
       pcoul    = y * ecoul

! first derivatives
       dpcouldd = third*ecoul + y*decouldd
       dpcouldt = y * decouldt
       dpcoulda = y * decoulda
       dpcouldz = y * decouldz

! second derivatives
       dpcoulddd = 2.0d0*third*decouldd + y*decoulddd
       dpcoulddt = third*decouldt + y*decoulddt
       dpcouldda = third*decoulda + y*decouldda
       dpcoulddz = third*decouldz + y*decoulddz
       dpcouldtt = y * decouldtt
       dpcouldta = y * decouldta
       dpcouldtz = y * decouldtz
       dpcouldaa = y * decouldaa
       dpcouldaz = y * decouldaz
       dpcouldzz = y * decouldzz


! entropy in erg/g/kelvin
       u0   = 3.0d0*b1*p1 - 5.0d0*c1*p2 + d1*(log(plasg) - 1.0d0) - e1
       du0  = 0.75d0*b1*p3 + 1.25d0*c1*p4 + d1*plasginv
       ddu0 = -0.5625d0*b1*p5 - 1.5625d0*c1*p6 - d1*plasginv*plasginv

       z    = -avo*ytot1*kerg

       scoul = z*u0
       ww    = z*du0
       x     = z*ddu0

! first derivatives
       dscouldd = ww*dplasgdd
       dscouldt = ww*dplasgdt
       dscoulda = ww*dplasgda - scoul*ytot1
       dscouldz = ww*dplasgdz

! second derivatives
       dscoulddd = x*dplasgdd*dplasgdd + ww*dplasgddd
       dscoulddt = x*dplasgdt*dplasgdd + ww*dplasgddt
       dscouldda = x*dplasgda*dplasgdd + ww*dplasgdda - x*ytot1*dplasgdd
       dscoulddz = x*dplasgdz*dplasgdd + ww*dplasgddz
       dscouldtt = x*dplasgdt*dplasgdt + ww*dplasgdtt
       dscouldta = x*dplasgda*dplasgdt + ww*dplasgdta - x*ytot1*dplasgdt
       dscouldtz = x*dplasgdz*dplasgdt + ww*dplasgdtz
       dscouldaa = x*dplasgda*dplasgda + ww*dplasgdaa - x*ytot1*dplasgda &
                   - ww*dplasgda*ytot1 + 2.0d0*scoul*ytot1*ytot1
       dscouldaz = x*dplasgdz*dplasgda + ww*dplasgdaz &
                   - ww*dplasgdz*ytot1
       dscouldzz = x*dplasgdz*dplasgdz + ww*dplasgdzz





! yakovlev & shalybkov 1989 equations 102, 103, 104
      else if (plasg .lt. 1.0) then

       x        = sqrt(plasg)
       p1       = plasg*x
       p2       = plasg**b2
       p3       = x
       p4       = p2*plasginv
       p5       = p3*plasginv
       p6       = p4*plasginv

       u0   = c2*p1 - third*a2*p2
       du0  = 1.5d0*c2*p3 - third*a2*b2*p4
       ddu0 = 0.75d0*c2*p5 - third*a2*b2*(b2-1.0d0)*p6


! pressure
       pcoul    = -pion * u0

       x        = pion*du0
       y        = pion*ddu0

! first derivatives
       dpcouldd = -dpiondd*u0 - x*dplasgdd
       dpcouldt = -dpiondt*u0 - x*dplasgdt
       dpcoulda = -dpionda*u0 - x*dplasgda
       dpcouldz = -dpiondz*u0 - x*dplasgdz


! second derivatives
       dpcoulddd = -dpionddd*u0 - (2.0d0*dpiondd*du0 &
                  + y*dplasgdd)*dplasgdd - x*dplasgddd
       dpcoulddt = -dpionddt*u0 - dpiondd*du0*dplasgdt &
                  - (dpiondt*du0 + y*dplasgdt)*dplasgdd - x*dplasgddt
       dpcouldda = -dpiondda*u0 - dpiondd*du0*dplasgda &
                  - (dpionda*du0 + y*dplasgda)*dplasgdd - x*dplasgdda
       dpcoulddz = -dpionddz*u0 - dpiondd*du0*dplasgdz &
                  - (dpiondz*du0 + y*dplasgdz)*dplasgdd - x*dplasgddz
       dpcouldtt = -dpiondtt*u0 - (2.0d0*dpiondt*du0 &
                  + y*dplasgdt)*dplasgdt - x*dplasgdtt
       dpcouldta = -dpiondta*u0 - dpiondt*du0*dplasgda &
                  - (dpionda*du0 + y*dplasgda)*dplasgdt - x*dplasgdta
       dpcouldtz = -dpiondtz*u0 - dpiondt*du0*dplasgdz &
                  - (dpiondz*du0 + y*dplasgdz)*dplasgdt - x*dplasgdtz
       dpcouldaa = -dpiondaa*u0 - (2.0d0*dpionda*du0 &
                  + y*dplasgda)*dplasgda - x*dplasgdaa
       dpcouldaz = -dpiondaz*u0 - dpionda*du0*dplasgdz &
                  - (dpiondz*du0 +  y*dplasgdz)*dplasgda - x*dplasgdaz
       dpcouldzz = -dpiondzz*u0 - (2.0d0*dpiondz*du0 &
                  + y*dplasgdz)*dplasgdz - x*dplasgdzz


! energy in erg/gr
       z        = 3.0d0*deninv
       y        = -z*deninv

       ecoul    = z*pcoul

       x = deninv*deninv

! first derivatives
       decouldd = z*dpcouldd - ecoul*deninv
       decouldt = z*dpcouldt
       decoulda = z*dpcoulda
       decouldz = z*dpcouldz

! second derivatives
       decoulddd = z*dpcoulddd + y*dpcouldd &
                  + (2.0d0*ecoul - 3.0d0*dpcouldd)*x
       decoulddt = z*dpcoulddt - 3.0d0*dpcouldt*x
       decouldda = z*dpcouldda - 3.0d0*dpcoulda*x
       decoulddz = z*dpcoulddz - 3.0d0*dpcouldz*x
       decouldtt = z*dpcouldtt
       decouldta = z*dpcouldta
       decouldtz = z*dpcouldtz
       decouldaa = z*dpcouldaa
       decouldaz = z*dpcouldaz
       decouldzz = z*dpcouldzz



! entropy in erg/g/kelvin
       u0    = c2*p1 - a2/b2*(b2-1.0d0)*p2
       du0   = 1.5d0*c2*p3 - a2*(b2-1.0d0)*p4
       ddu0  = 0.75d0*c2*p5 - a2*(b2-1.0d0)*(b2-1.0d0)*p6
       z     = -avo*ytot1*kerg
       y     = -z*ytot1

       scoul = z*u0
       x     = z*du0
       y     = z*ddu0

! first derivatives
       dscouldd = x*dplasgdd
       dscouldt = x*dplasgdt
       dscoulda = x*dplasgda - scoul*ytot1
       dscouldz = x*dplasgdz

! second derivatives
       dscoulddd = y*dplasgdd*dplasgdd + x*dplasgddd
       dscoulddt = y*dplasgdt*dplasgdd + x*dplasgddt
       dscouldda = y*dplasgda*dplasgdd + x*dplasgdda - x*ytot1*dplasgdd
       dscoulddz = y*dplasgdz*dplasgdd + x*dplasgddz
       dscouldtt = y*dplasgdt*dplasgdt + x*dplasgdtt
       dscouldta = y*dplasgda*dplasgdt + x*dplasgdta - x*ytot1*dplasgdt
       dscouldtz = y*dplasgdz*dplasgdt + x*dplasgdtz
       dscouldaa = y*dplasgda*dplasgda + x*dplasgdaa - x*ytot1*dplasgda &
                   - x*dplasgda*ytot1 + 2.0d0*scoul*ytot1*ytot1
       dscouldaz = y*dplasgdz*dplasgda + x*dplasgdaz - x*dplasgdz*ytot1
       dscouldzz = y*dplasgdz*dplasgdz + x*dplasgdzz

! end of plasg if block
      end if



! bomb proof
! wish we didn't have to do this if statement
        x   = prad + pion + pele + pcoul
        y   = erad + eion + eele + ecoul
        z   = srad + sion + sele + scoul
        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then

!         write(6,*)
!         write(6,*) 'coulomb corrections are causing a negative pressure'
!         write(6,*) 'setting all coulomb corrections to zero'
!         write(6,*)

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
        end if



! bomb proof the coulomb corrections
        x   = prad + pion + pele + ppos + pcoul
        y   = erad + eion + eele + epos + ecoul
        z   = srad + sion + sele + spos + scoul
        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then

!         write(6,*)
!         write(6,*) 'coulomb corrections are causing a negative pressure'
!         write(6,*) 'setting all coulomb corrections to zero'
!         write(6,*)

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
        end if
       end if





! sum all the gas components
       pgas    = pion + pele + ppos + pcoul
       egas    = eion + eele + epos + ecoul
       sgas    = sion + sele + spos + scoul

       dpgasdd = dpiondd + dpepdd + dpcouldd
       dpgasdt = dpiondt + dpepdt + dpcouldt
       dpgasda = dpionda + dpepda + dpcoulda
       dpgasdz = dpiondz + dpepdz + dpcouldz

       degasdd = deiondd + deepdd + decouldd
       degasdt = deiondt + deepdt + decouldt
       degasda = deionda + deepda + decoulda
       degasdz = deiondz + deepdz + decouldz

       dsgasdd = dsiondd + dsepdd + dscouldd
       dsgasdt = dsiondt + dsepdt + dscouldt
       dsgasda = dsionda + dsepda + dscoulda
       dsgasdz = dsiondz + dsepdz + dscouldz




! add in radiation to get the total
       pres    = prad + pgas
       ener    = erad + egas
       entr    = srad + sgas

       dpresdd = dpraddd + dpgasdd
       dpresdt = dpraddt + dpgasdt
       dpresda = dpradda + dpgasda
       dpresdz = dpraddz + dpgasdz

       denerdd = deraddd + degasdd
       denerdt = deraddt + degasdt
       denerda = deradda + degasda
       denerdz = deraddz + degasdz

       dentrdd = dsraddd + dsgasdd
       dentrdt = dsraddt + dsgasdt
       dentrda = dsradda + dsgasda
       dentrdz = dsraddz + dsgasdz



! for the gas
! the temperature and density exponents (c&g 9.81 9.82)
! the specific heat at constant volume (c&g 9.92)
! the third adiabatic exponent (c&g 9.93)
! the first adiabatic exponent (c&g 9.97)
! the second adiabatic exponent (c&g 9.105)
! the specific heat at constant pressure (c&g 9.98)
! and relativistic formula for the sound speed (c&g 14.29)

       zz        = pgas/den
       chit_gas  = temp/pgas * dpgasdt
       chid_gas  = dpgasdd/zz
       cv_gas    = degasdt
       x         = zz * chit_gas/(temp * cv_gas)
       gam3_gas  = x + 1.0d0
       gam1_gas  = chit_gas*x + chid_gas
       nabad_gas = x/gam1_gas
       gam2_gas  = 1.0d0/(1.0d0 - nabad_gas)
       cp_gas    = cv_gas * gam1_gas/chid_gas
       z         = 1.0d0 + (egas + clight*clight)/zz
       sound_gas = clight * sqrt(gam1_gas/z)


! for the totals
       zz    = pres/den
       chit  = temp/pres * dpresdt
       chid  = dpresdd/zz
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + clight*clight)/zz
       sound = clight * sqrt(gam1/z)




! maxwell relations; each is zero if the consistency is perfect
! delicate subtraction in very degenerate regions causes roundoff error

       dse = temp*dentrdt/denerdt - 1.0d0

       dpe = (denerdd*den**2 + temp*dpresdt)/pres - 1.0d0

       dsp = -dentrdd*(den**2/dpresdt) - 1.0d0


! store this row

! store this row
        ptot_row(j) = pres
        dpt_row(j)  = dpresdt
        dpd_row(j)  = dpresdd
        dpa_row(j)  = dpresda
        dpz_row(j)  = dpresdz
        dpdd_row(j) = dpresddd
        dpdt_row(j) = dpresddt
        dpda_row(j) = dpresdda
        dpdz_row(j) = dpresddz
        dptt_row(j) = dpresdtt
        dpta_row(j) = dpresdta
        dptz_row(j) = dpresdtz
        dpaa_row(j) = dpresdaa
        dpaz_row(j) = dpresdaz
        dpzz_row(j) = dpresdzz

        etot_row(j) = ener
        det_row(j)  = denerdt
        ded_row(j)  = denerdd
        dea_row(j)  = denerda
        dez_row(j)  = denerdz
        dedd_row(j) = denerddd
        dedt_row(j) = denerddt
        deda_row(j) = denerdda
        dedz_row(j) = denerddz
        dett_row(j) = denerdtt
        deta_row(j) = denerdta
        detz_row(j) = denerdtz
        deaa_row(j) = denerdaa
        deaz_row(j) = denerdaz
        dezz_row(j) = denerdzz

        stot_row(j) = entr
        dst_row(j)  = dentrdt
        dsd_row(j)  = dentrdd
        dsa_row(j)  = dentrda
        dsz_row(j)  = dentrdz
        dsdd_row(j) = dentrddd
        dsdt_row(j) = dentrddt
        dsda_row(j) = dentrdda
        dsdz_row(j) = dentrddz
        dstt_row(j) = dentrdtt
        dsta_row(j) = dentrdta
        dstz_row(j) = dentrdtz
        dsaa_row(j) = dentrdaa
        dsaz_row(j) = dentrdaz
        dszz_row(j) = dentrdzz


        pgas_row(j)    = pgas
        dpgast_row(j)  = dpgasdt
        dpgasd_row(j)  = dpgasdd
        dpgasa_row(j)  = dpgasda
        dpgasz_row(j)  = dpgasdz
        dpgasdd_row(j) = dpgasddd
        dpgasdt_row(j) = dpgasddt
        dpgasda_row(j) = dpgasdda
        dpgasdz_row(j) = dpgasddz
        dpgastt_row(j) = dpgasdtt
        dpgasta_row(j) = dpgasdta
        dpgastz_row(j) = dpgasdtz
        dpgasaa_row(j) = dpgasdaa
        dpgasaz_row(j) = dpgasdaz
        dpgaszz_row(j) = dpgasdzz

        egas_row(j)    = egas
        degast_row(j)  = degasdt
        degasd_row(j)  = degasdd
        degasa_row(j)  = degasda
        degasz_row(j)  = degasdz
        degasdd_row(j) = degasddd
        degasdt_row(j) = degasddt
        degasda_row(j) = degasdda
        degasdz_row(j) = degasddz
        degastt_row(j) = degasdtt
        degasta_row(j) = degasdta
        degastz_row(j) = degasdtz
        degasaa_row(j) = degasdaa
        degasaz_row(j) = degasdaz
        degaszz_row(j) = degasdzz

        sgas_row(j)    = sgas
        dsgast_row(j)  = dsgasdt
        dsgasd_row(j)  = dsgasdd
        dsgasa_row(j)  = dsgasda
        dsgasz_row(j)  = dsgasdz
        dsgasdd_row(j) = dsgasddd
        dsgasdt_row(j) = dsgasddt
        dsgasda_row(j) = dsgasdda
        dsgasdz_row(j) = dsgasddz
        dsgastt_row(j) = dsgasdtt
        dsgasta_row(j) = dsgasdta
        dsgastz_row(j) = dsgasdtz
        dsgasaa_row(j) = dsgasdaa
        dsgasaz_row(j) = dsgasdaz
        dsgaszz_row(j) = dsgasdzz

        prad_row(j)    = prad
        dpradt_row(j)  = dpraddt
        dpradd_row(j)  = dpraddd
        dprada_row(j)  = dpradda
        dpradz_row(j)  = dpraddz
        dpraddd_row(j) = dpradddd
        dpraddt_row(j) = dpradddt
        dpradda_row(j) = dpraddda
        dpraddz_row(j) = dpradddz
        dpradtt_row(j) = dpraddtt
        dpradta_row(j) = dpraddta
        dpradtz_row(j) = dpraddtz
        dpradaa_row(j) = dpraddaa
        dpradaz_row(j) = dpraddaz
        dpradzz_row(j) = dpraddzz

        erad_row(j)    = erad
        deradt_row(j)  = deraddt
        deradd_row(j)  = deraddd
        derada_row(j)  = deradda
        deradz_row(j)  = deraddz
        deraddd_row(j) = deradddd
        deraddt_row(j) = deradddt
        deradda_row(j) = deraddda
        deraddz_row(j) = deradddz
        deradtt_row(j) = deraddtt
        deradta_row(j) = deraddta
        deradtz_row(j) = deraddtz
        deradaa_row(j) = deraddaa
        deradaz_row(j) = deraddaz
        deradzz_row(j) = deraddzz

        srad_row(j)    = srad
        dsradt_row(j)  = dsraddt
        dsradd_row(j)  = dsraddd
        dsrada_row(j)  = dsradda
        dsradz_row(j)  = dsraddz
        dsraddd_row(j) = dsradddd
        dsraddt_row(j) = dsradddt
        dsradda_row(j) = dsraddda
        dsraddz_row(j) = dsradddz
        dsradtt_row(j) = dsraddtt
        dsradta_row(j) = dsraddta
        dsradtz_row(j) = dsraddtz
        dsradaa_row(j) = dsraddaa
        dsradaz_row(j) = dsraddaz
        dsradzz_row(j) = dsraddzz

        pion_row(j)   = pion
        dpiont_row(j) = dpiondt
        dpiond_row(j) = dpiondd
        dpiona_row(j) = dpionda
        dpionz_row(j) = dpiondz
        dpiondd_row(j) = dpionddd
        dpiondt_row(j) = dpionddt
        dpionda_row(j) = dpiondda
        dpiondz_row(j) = dpionddz
        dpiontt_row(j) = dpiondtt
        dpionta_row(j) = dpiondta
        dpiontz_row(j) = dpiondtz
        dpionaa_row(j) = dpiondaa
        dpionaz_row(j) = dpiondaz
        dpionzz_row(j) = dpiondzz

        eion_row(j)   = eion
        deiont_row(j) = deiondt
        deiond_row(j) = deiondd
        deiona_row(j) = deionda
        deionz_row(j) = deiondz
        deiondd_row(j) = deionddd
        deiondt_row(j) = deionddt
        deionda_row(j) = deiondda
        deiondz_row(j) = deionddz
        deiontt_row(j) = deiondtt
        deionta_row(j) = deiondta
        deiontz_row(j) = deiondtz
        deionaa_row(j) = deiondaa
        deionaz_row(j) = deiondaz
        deionzz_row(j) = deiondz

        sion_row(j)   = sion
        dsiont_row(j) = dsiondt
        dsiond_row(j) = dsiondd
        dsiona_row(j) = dsionda
        dsionz_row(j) = dsiondz
        dsiondd_row(j) = dsionddd
        dsiondt_row(j) = dsionddt
        dsionda_row(j) = dsiondda
        dsiondz_row(j) = dsionddz
        dsiontt_row(j) = dsiondtt
        dsionta_row(j) = dsiondta
        dsiontz_row(j) = dsiondtz
        dsionaa_row(j) = dsiondaa
        dsionaz_row(j) = dsiondaz
        dsionzz_row(j) = dsiondzz

        xnim_row(j)   = xni
        xni_row(j)    = xnifer
        dxnit_row(j)  = dxniferdt
        dxnid_row(j)  = dxniferdd
        dxnia_row(j)  = dxniferda
        dxniz_row(j)  = dxniferdz
        dxnidd_row(j) = dxniferddd
        dxnidt_row(j) = dxniferddt
        dxnida_row(j) = dxniferdda
        dxnidz_row(j) = dxniferddz
        dxnitt_row(j) = dxniferdtt
        dxnita_row(j) = dxniferdta
        dxnitz_row(j) = dxniferdtz
        dxniaa_row(j) = dxniferdaa
        dxniaz_row(j) = dxniferdaz
        dxnizz_row(j) = dxniferdzz

        etaion_row(j)  = etaion
        detait_row(j)  = detaidt
        detaid_row(j)  = detaidd
        detaia_row(j)  = detaida
        detaiz_row(j)  = detaidz
        detaidd_row(j) = detaiddd
        detaidt_row(j) = detaiddt
        detaida_row(j) = detaidda
        detaidz_row(j) = detaiddz
        detaitt_row(j) = detaidtt
        detaita_row(j) = detaidta
        detaitz_row(j) = detaidtz
        detaiaa_row(j) = detaidaa
        detaiaz_row(j) = detaidaz
        detaizz_row(j) = detaidzz

        pele_row(j)   = pele
        ppos_row(j)   = ppos
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = dpepda
        dpepz_row(j)  = dpepdz
        dpepdd_row(j) = dpepddd
        dpepdt_row(j) = dpepddt
        dpepda_row(j) = dpepdda
        dpepdz_row(j) = dpepddz
        dpeptt_row(j) = dpepdtt
        dpepta_row(j) = dpepdta
        dpeptz_row(j) = dpepdtz
        dpepaa_row(j) = dpepdaa
        dpepaz_row(j) = dpepdaz
        dpepzz_row(j) = dpepdzz

        eele_row(j)   = eele
        epos_row(j)   = epos
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = deepda
        deepz_row(j)  = deepdz
        deepdd_row(j) = deepddd
        deepdt_row(j) = deepddt
        deepda_row(j) = deepdda
        deepdz_row(j) = deepddz
        deeptt_row(j) = deepdtt
        deepta_row(j) = deepdta
        deeptz_row(j) = deepdtz
        deepaa_row(j) = deepdaa
        deepaz_row(j) = deepdaz
        deepzz_row(j) = deepdzz

        sele_row(j)   = sele
        spos_row(j)   = spos
        dsept_row(j)  = dsepdt
        dsepd_row(j)  = dsepdd
        dsepa_row(j)  = dsepda
        dsepz_row(j)  = dsepdz
        dsepdd_row(j) = dsepddd
        dsepdt_row(j) = dsepddt
        dsepda_row(j) = dsepdda
        dsepdz_row(j) = dsepddz
        dseptt_row(j) = dsepdtt
        dsepta_row(j) = dsepdta
        dseptz_row(j) = dsepdtz
        dsepaa_row(j) = dsepdaa
        dsepaz_row(j) = dsepdaz
        dsepzz_row(j) = dsepdzz

        xnem_row(j)   = xne
        xne_row(j)    = xnefer
        xnp_row(j)    = xnpfer
        dxnet_row(j)  = dxneferdt + dxnpferdt
        dxned_row(j)  = dxneferdd + dxnpferdd
        dxnea_row(j)  = dxneferda + dxnpferda
        dxnez_row(j)  = dxneferdz + dxnpferdz
        dxnedd_row(j) = dxneferddd + dxnpferddd
        dxnedt_row(j) = dxneferddt + dxnpferddt
        dxneda_row(j) = dxneferdda + dxnpferdda
        dxnedz_row(j) = dxneferddz + dxnpferddz
        dxnett_row(j) = dxneferdtt + dxnpferdtt
        dxneta_row(j) = dxneferdta + dxnpferdta
        dxnetz_row(j) = dxneferdtz + dxnpferdtz
        dxneaa_row(j) = dxneferdaa + dxnpferdaa
        dxneaz_row(j) = dxneferdaz + dxnpferdaz
        dxnezz_row(j) = dxneferdzz + dxnpferdzz


        zeff_row(j)   = zeff

        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd
        detaa_row(j)  = detada
        detaz_row(j)  = detadz
        detadd_row(j) = detaddd
        detadt_row(j) = detaddt
        detada_row(j) = detadda
        detadz_row(j) = detaddz
        detatt_row(j) = detadtt
        detata_row(j) = detadta
        detatz_row(j) = detadtz
        detaaa_row(j) = detadaa
        detaaz_row(j) = detadaz
        detazz_row(j) = detadzz


        etapos_row(j) = etapos

        pip_row(j)    = pip
        eip_row(j)    = eip
        sip_row(j)    = sip


        pcou_row(j)   = pcoul
        dpcout_row(j) = dpcouldt
        dpcoud_row(j) = dpcouldd
        dpcoua_row(j) = dpcoulda
        dpcouz_row(j) = dpcouldz

        ecou_row(j)   = ecoul
        decout_row(j) = decouldt
        decoud_row(j) = decouldd
        decoua_row(j) = decoulda
        decouz_row(j) = decouldz

        scou_row(j)   = scoul
        dscout_row(j) = dscouldt
        dscoud_row(j) = dscouldd
        dscoua_row(j) = dscoulda
        dscouz_row(j) = dscouldz

        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_gas_row(j)    = cv_gas
        dcv_gasdd_row(j) = dcv_gasdd
        dcv_gasdt_row(j) = dcv_gasdt
        dcv_gasda_row(j) = dcv_gasda
        dcv_gasdz_row(j) = dcv_gasdz

        cp_gas_row(j)    = cp_gas
        dcp_gasdd_row(j) = dcp_gasdd
        dcp_gasdt_row(j) = dcp_gasdt
        dcp_gasda_row(j) = dcp_gasda
        dcp_gasdz_row(j) = dcp_gasdz

        gam1_gas_row(j)    = gam1_gas
        dgam1_gasdd_row(j) = dgam1_gasdd
        dgam1_gasdt_row(j) = dgam1_gasdt
        dgam1_gasda_row(j) = dgam1_gasda
        dgam1_gasdz_row(j) = dgam1_gasdz

        gam2_gas_row(j)    = gam2_gas
        dgam2_gasdd_row(j) = dgam2_gasdd
        dgam2_gasdt_row(j) = dgam2_gasdt
        dgam2_gasda_row(j) = dgam2_gasda
        dgam2_gasdz_row(j) = dgam2_gasdz

        gam3_gas_row(j)    = gam3_gas
        dgam3_gasdd_row(j) = dgam3_gasdd
        dgam3_gasdt_row(j) = dgam3_gasdt
        dgam3_gasda_row(j) = dgam3_gasda
        dgam3_gasdz_row(j) = dgam3_gasdz

        nabad_gas_row(j)  = nabad_gas
        dnab_gasdd_row(j) = dnab_gasdd
        dnab_gasdt_row(j) = dnab_gasdt
        dnab_gasda_row(j) = dnab_gasda
        dnab_gasdz_row(j) = dnab_gasdz

        cs_gas_row(j)    = sound_gas
        dcs_gasdd_row(j) = dcs_gasdd
        dcs_gasdt_row(j) = dcs_gasdt
        dcs_gasda_row(j) = dcs_gasda
        dcs_gasdz_row(j) = dcs_gasdz

        cv_row(j)     = cv
        dcvdd_row(j)  = dcvdd
        dcvdt_row(j)  = dcvdt
        dcvda_row(j)  = dcvda
        dcvdz_row(j)  = dcvdz

        cp_row(j)     = cp
        dcpdd_row(j)  = dcpdd
        dcpdt_row(j)  = dcpdt
        dcpda_row(j)  = dcpda
        dcpdz_row(j)  = dcpdz

        gam1_row(j)    = gam1
        dgam1dd_row(j) = dgam1dd
        dgam1dt_row(j) = dgam1dt
        dgam1da_row(j) = dgam1da
        dgam1dz_row(j) = dgam1dz

        gam2_row(j)    = gam2
        dgam2dd_row(j) = dgam2dd
        dgam2dt_row(j) = dgam2dt
        dgam2da_row(j) = dgam2da
        dgam2dz_row(j) = dgam2dz

        gam3_row(j)    = gam3
        dgam3dd_row(j) = dgam3dd
        dgam3dt_row(j) = dgam3dt
        dgam3da_row(j) = dgam3da
        dgam3dz_row(j) = dgam3dz

        nabad_row(j)  = nabad
        dnabdd_row(j) = dnabdd
        dnabdt_row(j) = dnabdt
        dnabda_row(j) = dnabda
        dnabdz_row(j) = dnabdz

        cs_row(j)    = sound
        dcsdd_row(j) = dcsdd
        dcsdt_row(j) = dcsdt
        dcsda_row(j) = dcsda
        dcsdz_row(j) = dcsdz


! end of pipeline loop
       end if
      enddo
      return
      end











      subroutine xneroot(mode,den,temp,abar,zbar,ionized,potmult,aa, &
                         f,df)

      include 'implno.dek'
      include 'const.dek'

! this routine is called by a root finder to find the degeneracy parameter aa
! where the number density from a saha equation equals the number density as
! computed by the fermi-dirac integrals.

! input:
! mode    = 0 = root find on electron number density, = 1 = full calculation
! temp    = temperature
! den     = density
! abar    = average weight
! zbar    = average charge
! aa      = degeneracy parameter (chemical potential/kerg*temp)
! ionized = flag to turn off/on ionization contributions
! potmult = flag to turn off/on ionization ptential contributions

! output
! f = value of the function to be zeroed ; xne - xnefer + xnepos = 0
! df = derivative with respect to aa of a


! declare the pass
      integer          mode,ionized,potmult
      double precision den,temp,abar,zbar,aa,f,df


! bring in common block variables
! common block communication with routine xneroot

! electron-positrons

      double precision etaele,detadd,detadt,detada,detadz, &
                       detaddd,detaddt,detadda,detaddz,detadtt, &
                       detadta,detadtz,detadaa,detadaz,detadzz

      common /xneta/   etaele,detadd,detadt,detada,detadz, &
                       detaddd,detaddt,detadda,detaddz,detadtt, &
                       detadta,detadtz,detadaa,detadaz,detadzz


      double precision &
                       pep,dpepdd,dpepdt,dpepda,dpepdz, &
                       dpepddd,dpepddt,dpepdda,dpepddz, &
                       dpepdtt,dpepdta,dpepdtz,dpepdaa, &
                       dpepdaz,dpepdzz, &
                       eep,deepdd,deepdt,deepda,deepdz, &
                       deepddd,deepddt,deepdda,deepddz, &
                       deepdtt,deepdta,deepdtz,deepdaa, &
                       deepdaz,deepdzz, &
                       sep,dsepdd,dsepdt,dsepda,dsepdz, &
                       dsepddd,dsepddt,dsepdda,dsepddz, &
                       dsepdtt,dsepdta,dsepdtz,dsepdaa, &
                       dsepdaz,dsepdzz

      common /epc1/ &
                       pep,dpepdd,dpepdt,dpepda,dpepdz, &
                       dpepddd,dpepddt,dpepdda,dpepddz, &
                       dpepdtt,dpepdta,dpepdtz,dpepdaa, &
                       dpepdaz,dpepdzz, &
                       eep,deepdd,deepdt,deepda,deepdz, &
                       deepddd,deepddt,deepdda,deepddz, &
                       deepdtt,deepdta,deepdtz,deepdaa, &
                       deepdaz,deepdzz, &
                       sep,dsepdd,dsepdt,dsepda,dsepdz, &
                       dsepddd,dsepddt,dsepdda,dsepddz, &
                       dsepdtt,dsepdta,dsepdtz,dsepdaa, &
                       dsepdaz,dsepdzz


      double precision &
                       etapos,zeff
      common /xnec1/ &
                       etapos,zeff


      double precision &
                       pele,dpeledd,dpeledt,dpeleda,dpeledz, &
                       dpeleddd,dpeleddt,dpeledda,dpeleddz, &
                       dpeledtt,dpeledta,dpeledtz,dpeledaa, &
                       dpeledaz,dpeledzz, &
                       eele,deeledd,deeledt,deeleda,deeledz, &
                       deeleddd,deeleddt,deeledda,deeleddz, &
                       deeledtt,deeledta,deeledtz,deeledaa, &
                       deeledaz,deeledzz, &
                       sele,dseledd,dseledt,dseleda,dseledz, &
                       dseleddd,dseleddt,dseledda,dseleddz, &
                       dseledtt,dseledta,dseledtz,dseledaa, &
                       dseledaz,dseledzz
      common /eleth1/ &
                       pele,dpeledd,dpeledt,dpeleda,dpeledz, &
                       dpeleddd,dpeleddt,dpeledda,dpeleddz, &
                       dpeledtt,dpeledta,dpeledtz,dpeledaa, &
                       dpeledaz,dpeledzz, &
                       eele,deeledd,deeledt,deeleda,deeledz, &
                       deeleddd,deeleddt,deeledda,deeleddz, &
                       deeledtt,deeledta,deeledtz,deeledaa, &
                       deeledaz,deeledzz, &
                       sele,dseledd,dseledt,dseleda,dseledz, &
                       dseleddd,dseleddt,dseledda,dseleddz, &
                       dseledtt,dseledta,dseledtz,dseledaa, &
                       dseledaz,dseledzz

      double precision &
                       ppos,dpposdd,dpposdt,dpposda,dpposdz, &
                       dpposddd,dpposddt,dpposdda,dpposddz, &
                       dpposdtt,dpposdta,dpposdtz,dpposdaa, &
                       dpposdaz,dpposdzz, &
                       epos,deposdd,deposdt,deposda,deposdz, &
                       deposddd,deposddt,deposdda,deposddz, &
                       deposdtt,deposdta,deposdtz,deposdaa, &
                       deposdaz,deposdzz, &
                       spos,dsposdd,dsposdt,dsposda,dsposdz, &
                       dsposddd,dsposddt,dsposdda,dsposddz, &
                       dsposdtt,dsposdta,dsposdtz,dsposdaa, &
                       dsposdaz,dsposdzz
      common /posth1/ &
                       ppos,dpposdd,dpposdt,dpposda,dpposdz, &
                       dpposddd,dpposddt,dpposdda,dpposddz, &
                       dpposdtt,dpposdta,dpposdtz,dpposdaa, &
                       dpposdaz,dpposdzz, &
                       epos,deposdd,deposdt,deposda,deposdz, &
                       deposddd,deposddt,deposdda,deposddz, &
                       deposdtt,deposdta,deposdtz,deposdaa, &
                       deposdaz,deposdzz, &
                       spos,dsposdd,dsposdt,dsposda,dsposdz, &
                       dsposddd,dsposddt,dsposdda,dsposddz, &
                       dsposdtt,dsposdta,dsposdtz,dsposdaa, &
                       dsposdaz,dsposdzz


      double precision xne, &
                       dxnedd,dxnedt,dxneda,dxnedz, &
                       dxneddd,dxneddt,dxnedda,dxneddz, &
                       dxnedtt,dxnedta,dxnedtz,dxnedaa, &
                       dxnedaz,dxnedzz

      common /xnec2/   xne, &
                       dxnedd,dxnedt,dxneda,dxnedz, &
                       dxneddd,dxneddt,dxnedda,dxneddz, &
                       dxnedtt,dxnedta,dxnedtz,dxnedaa, &
                       dxnedaz,dxnedzz

      double precision &
                       xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz, &
                       dxneferddd,dxneferddt,dxneferdda,dxneferddz, &
                       dxneferdtt,dxneferdta,dxneferdtz,dxneferdaa, &
                       dxneferdaz,dxneferdzz, &
                       xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz, &
                       dxnpferddd,dxnpferddt,dxnpferdda,dxnpferddz, &
                       dxnpferdtt,dxnpferdta,dxnpferdtz,dxnpferdaa, &
                       dxnpferdaz,dxnpferdzz
      common /xnec3/ &
                       xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz, &
                       dxneferddd,dxneferddt,dxneferdda,dxneferddz, &
                       dxneferdtt,dxneferdta,dxneferdtz,dxneferdaa, &
                       dxneferdaz,dxneferdzz, &
                       xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz, &
                       dxnpferddd,dxnpferddt,dxnpferdda,dxnpferddz, &
                       dxnpferdtt,dxnpferdta,dxnpferdtz,dxnpferdaa, &
                       dxnpferdaz,dxnpferdzz





! ionization contributions
      double precision eip,deipdd,deipdt,deipda,deipdz, &
                       sip,dsipdd,dsipdt,dsipda,dsipdz, &
                       pip,dpipdd,dpipdt,dpipda,dpipdz

      common /xnec4/   eip,deipdd,deipdt,deipda,deipdz, &
                       sip,dsipdd,dsipdt,dsipda,dsipdz, &
                       pip,dpipdd,dpipdt,dpipda,dpipdz




! local variables
      double precision deni,kt,kti,beta,beta12,beta32,beta52, &
                       f12,f12eta,f12beta,f12eta2,f12beta2,f12etabeta, &
                       f32,f32eta,f32beta,f32eta2,f32beta2,f32etabeta, &
                       f52,f52eta,f52beta,f52eta2,f52beta2,f52etabeta, &
                       ytot1,zz,y,yy,ww,dum1,dum2,dum3,denion

      double precision xni,dxnidd,dxnidt,dxnida,dxnidz, &
                       dxniddd,dxniddt,dxnidda,dxniddz, &
                       dxnidtt,dxnidta,dxnidtz,dxnidaa, &
                       dxnidaz,dxnidzz

      double precision chi,chifac,dchifacdt,dchifacdz, &
                       dchifacdtt,dchifacdtz,dchifacdzz

      double precision saha, &
                       dsaha_dd,dsaha_dt,dsaha_da,dsaha_dz,dsaha_deta, &
                       dsaha_ddd,dsaha_ddt,dsaha_dda,dsaha_ddz, &
                       dsaha_dtt,dsaha_dta,dsaha_dtz,dsaha_daa, &
                       dsaha_daz,dsaha_dzz,dsaha_deta_dd,dsaha_deta_dt, &
                       dsaha_deta_da,dsaha_deta_dz,dsaha_deta2

      double precision sfac, &
                       dsfac_dd,dsfac_dt,dsfac_da,dsfac_dz,dsfac_deta, &
                       dsfac_ddd,dsfac_ddt,dsfac_dda,dsfac_ddz, &
                       dsfac_dtt,dsfac_dta,dsfac_dtz,dsfac_daa, &
                       dsfac_daz,dsfac_dzz,dsfac_deta_dd,dsfac_deta_dt, &
                       dsfac_deta_da,dsfac_deta_dz,dsfac_deta2

      double precision &
                       dzeff_dd,dzeff_dt,dzeff_da,dzeff_dz,dzeff_deta, &
                       dzeff_ddd,dzeff_ddt,dzeff_dda,dzeff_ddz, &
                       dzeff_dtt,dzeff_dta,dzeff_dtz,dzeff_daa, &
                       dzeff_daz,dzeff_dzz,dzeff_deta_dd,dzeff_deta_dt, &
                       dzeff_deta_da,dzeff_deta_dz,dzeff_deta2

      double precision &
                       dxne_dd,dxne_dt,dxne_da,dxne_dz,dxne_deta, &
                       dxne_ddd,dxne_ddt,dxne_dda,dxne_ddz, &
                       dxne_dtt,dxne_dta,dxne_dtz,dxne_daa, &
                       dxne_daz,dxne_dzz,dxne_deta_dd,dxne_deta_dt, &
                       dxne_deta_da,dxne_deta_dz,dxne_deta2


      double precision dxnefer_deta,dxnefer_dbeta, &
                       dxnefer_deta2,dxnefer_dbeta2, &
                       dxnefer_deta_dbeta, &
                       dxnpfer_detap, &
                       dxnpfer_detap2, &
                       dxnpfer_detap_dbeta, &
                       detap_deta,detap_dbeta, &
                       detap_deta2,detap_dbeta2, &
                       detap_deta_dbeta, &
                       dxnpfer_deta,dxnpfer_dbeta, &
                       dxnpfer_deta2,dxnpfer_dbeta2, &
                       dxnpfer_deta_dbeta

      double precision dxep_deta,dxep_dbeta, &
                       dxep_deta2,dxep_dbeta2, &
                       dxep_deta_dbeta

      double precision &
                       dzeffdd,dzeffdt,dzeffda,dzeffdz, &
                       dzeffddd,dzeffddt,dzeffdda,dzeffddz, &
                       dzeffdtt,dzeffdta,dzeffdtz,dzeffdaa, &
                       dzeffdaz,dzeffdzz

      double precision &
                       dpele_deta,dpele_dbeta, &
                       dpele_deta2,dpele_dbeta2,dpele_deta_dbeta, &
                       deele_deta,deele_dbeta, &
                       deele_deta2,deele_dbeta2,deele_deta_dbeta

      double precision &
                       dppos_deta,dppos_dbeta, &
                       dppos_deta2,dppos_dbeta2,dppos_deta_dbeta, &
                       dppos_detap,dppos_detap2,dppos_detap_dbeta, &
                       depos_deta,depos_dbeta, &
                       depos_deta2,depos_dbeta2,depos_deta_dbeta, &
                       depos_detap,depos_detap2,depos_detap_dbeta


      double precision &
                       deipddd,deipddt,deipdda,deipddz, &
                       deipdtt,deipdta,deipdtz,deipdaa, &
                       deipdaz,deipdzz, &
                       dsipddd,dsipddt,dsipdda,dsipddz, &
                       dsipdtt,dsipdta,dsipdtz,dsipdaa, &
                       dsipdaz,dsipdzz


      double precision xconst,pconst,econst,mecc,dbetadt,safe, &
                       positron_start
      parameter        (xconst  = 8.0d0 * pi * sqrt(2.0d0) * (me/h)**3 * clight**3, &
                        pconst  = xconst * 2.0d0/3.0d0 * me * clight**2, &
                        econst  = xconst * me * clight**2, &
                        mecc    = me * clight * clight, &
                        dbetadt = kerg/mecc, &
                        safe    = 0.005d0, &
                        positron_start = 0.02d0)

! some common factors
      ytot1   = 1.0d0/abar
      deni    = 1.0d0/den
      kt      = kerg * temp
      kti     = 1.0d0/kt
      beta    = kt/mecc
      beta12  = sqrt(beta)
      beta32  = beta * beta12
      beta52  = beta * beta32
      etaele  = aa


! ion number density in 1/cm**3
      xni     = avo * ytot1 * den
      dxnidd  = avo * ytot1
      dxnidt  = 0.0d0
      dxnida  = -xni * ytot1
      dxnidz  = 0.0d0

      dxniddd = 0.0d0
      dxniddt = 0.0d0
      dxnidda = -dxnidd*ytot1
      dxniddz = 0.0d0
      dxnidtt = 0.0d0
      dxnidta = 0.0d0
      dxnidtz = 0.0d0
      dxnidaa = -2.0d0 * dxnida * ytot1
      dxnidaz = 0.0d0
      dxnidzz = 0.0d0



! get the number density of free electrons
! saha is the ratio of the ground state to the ionized state
! this model is exact for a pure hydrogen composition
! denion is a crude pressure ionization model


! assume fully ionized
      chi        = 0.0d0
      chifac     = 0.0d0
      dchifacdt  = 0.0d0
      dchifacdz  = 0.0d0
      dchifacdtt = 0.0d0
      dchifacdtz = 0.0d0
      dchifacdzz = 0.0d0
      saha       = 0.0d0
      dsaha_dd   = 0.0d0
      dsaha_dt   = 0.0d0
      dsaha_da   = 0.0d0
      dsaha_dz   = 0.0d0
      dsaha_deta = 0.0d0
      dsaha_ddd  = 0.0d0
      dsaha_ddt  = 0.0d0
      dsaha_dda  = 0.0d0
      dsaha_ddz  = 0.0d0
      dsaha_dtt  = 0.0d0
      dsaha_dta  = 0.0d0
      dsaha_dtz  = 0.0d0
      dsaha_daa  = 0.0d0
      dsaha_daz  = 0.0d0
      dsaha_dzz  = 0.0d0
      dsaha_deta_dd = 0.0d0
      dsaha_deta_dt = 0.0d0
      dsaha_deta_da = 0.0d0
      dsaha_deta_dz = 0.0d0
      dsaha_deta2   = 0.0d0


! do a simple saha approach

      if (ionized .eq. 0) then

       denion     = 0.1d0
       chi        = hion * ev2erg * zbar
       chifac     = chi*kti
       dchifacdt  = -chifac/temp
       dchifacdz  = chifac/zbar
       dchifacdtt = -2.0d0*dchifacdt/temp
       dchifacdtz = -dchifacdz/temp
       dchifacdzz = 0.0d0

       yy         = chifac - den/denion

! completely neutral; set it for a fake convergence
       if (yy .gt. 200.0) then
        saha       = 1.0d90
        f          = 0.0d0
        df         = 1.0d0
        etaele     = -100.0d0
        if (mode .eq. 0) return

! ionization possible
       else if (yy .gt. -200.0) then
        ww   = min(200.0d0,chifac + etaele - den/denion)
        saha = 2.0d0 * exp(ww)
        if (ww .ne. 200.0d0) then

! first derivatives
         dsaha_dd   = -saha/denion
         dsaha_dt   = saha * dchifacdt
         dsaha_da   = 0.0d0
         dsaha_dz   = saha * dchifacdz
         dsaha_deta = saha

! second derivatives
         dsaha_ddd  = -dsaha_dd/denion
         dsaha_ddt  = -dsaha_dt/denion
         dsaha_dda  = -dsaha_da/denion
         dsaha_ddz  = -dsaha_dz/denion
         dsaha_dtt  = dsaha_dt*dchifacdt + saha*dchifacdtt
         dsaha_dta  = dsaha_da*dchifacdt
         dsaha_dtz  = dsaha_dz*dchifacdt + saha*dchifacdtz
         dsaha_daa  = 0.0d0
         dsaha_daz  = 0.0d0
         dsaha_dzz  = dsaha_dz*dchifacdz + saha * dchifacdzz
         dsaha_deta_dd = dsaha_dd
         dsaha_deta_dt = dsaha_dt
         dsaha_deta_da = dsaha_da
         dsaha_deta_dz = dsaha_dz
         dsaha_deta2   = dsaha_deta

        end if
       end if
      end if



! saha factor, the fraction ionized
      sfac       = 1.0d0/(1.0d0 + saha)

! first derivatives
      y          = -sfac*sfac
      dsfac_dd   = y*dsaha_dd
      dsfac_dt   = y*dsaha_dt
      dsfac_da   = y*dsaha_da
      dsfac_dz   = y*dsaha_dz
      dsfac_deta = y*dsaha_deta

! second derivatives
      ww         = -2.0d0*sfac
      dsfac_ddd  = ww*dsfac_dd*dsaha_dd + y*dsaha_ddd
      dsfac_ddt  = ww*dsfac_dt*dsaha_dd + y*dsaha_ddt
      dsfac_dda  = ww*dsfac_da*dsaha_dd + y*dsaha_dda
      dsfac_ddz  = ww*dsfac_dz*dsaha_dd + y*dsaha_ddz
      dsfac_dtt  = ww*dsfac_dt*dsaha_dt + y*dsaha_dtt
      dsfac_dta  = ww*dsfac_da*dsaha_dt + y*dsaha_dta
      dsfac_dtz  = ww*dsfac_dz*dsaha_dt + y*dsaha_dtz
      dsfac_daa  = ww*dsfac_da*dsaha_da + y*dsaha_daa
      dsfac_daz  = ww*dsfac_dz*dsaha_da + y*dsaha_daz
      dsfac_dzz  = ww*dsfac_dz*dsaha_dz + y*dsaha_dzz
      dsfac_deta_dd = ww*dsfac_dd*dsaha_deta + y*dsaha_deta_dd
      dsfac_deta_dt = ww*dsfac_dt*dsaha_deta + y*dsaha_deta_dt
      dsfac_deta_da = ww*dsfac_da*dsaha_deta + y*dsaha_deta_da
      dsfac_deta_dz = ww*dsfac_dz*dsaha_deta + y*dsaha_deta_dz
      dsfac_deta2   = ww*dsfac_deta*dsaha_deta + y*dsaha_deta2



! effective charge
      zeff       = zbar * sfac

! first derivatives
      dzeff_dd   = zbar * dsfac_dd
      dzeff_dt   = zbar * dsfac_dt
      dzeff_da   = zbar * dsfac_da
      dzeff_dz   = sfac
      dzeff_deta = zbar * dsfac_deta

! second derivatives
      dzeff_ddd  = zbar*dsfac_ddd
      dzeff_ddt  = zbar*dsfac_ddt
      dzeff_dda  = zbar*dsfac_dda
      dzeff_ddz  = zbar*dsfac_ddz
      dzeff_dtt  = zbar*dsfac_dtt
      dzeff_dta  = zbar*dsfac_dta
      dzeff_dtz  = zbar*dsfac_dtz
      dzeff_daa  = zbar*dsfac_daa
      dzeff_daz  = zbar*dsfac_daz
      dzeff_dzz  = dsfac_dz
      dzeff_deta_dd = zbar*dsfac_deta_dd
      dzeff_deta_dt = zbar*dsfac_deta_dt
      dzeff_deta_da = zbar*dsfac_deta_da
      dzeff_deta_dz = dsfac_deta + zbar*dsfac_deta_dz
      dzeff_deta2   = zbar * dsfac_deta2




! number density of free electrons
      xne        = xni * zeff

! first derivatives
      dxne_dd    = dxnidd * zeff + xni*dzeff_dd
      dxne_dt    = dxnidt * zeff + xni*dzeff_dt
      dxne_da    = dxnida * zeff + xni*dzeff_da
      dxne_dz    = dxnidz * zeff + xni*dzeff_dz
      dxne_deta  = xni*dzeff_deta

! second derivatives
      dxne_ddd = dxniddd*zeff + 2.0d0*dxnidd*dzeff_dd + xni*dzeff_ddd
      dxne_ddt = dxniddt*zeff + dxnidd*dzeff_dt &
                 + dxnidt*dzeff_dt + xni*dzeff_ddt
      dxne_dda = dxnidda*zeff + dxnidd*dzeff_da &
                 + dxnida*dzeff_dt + xni*dzeff_dda
      dxne_ddz = dxniddz*zeff + dxnidd*dzeff_dz &
                 + dxnidz*dzeff_dt + xni*dzeff_ddz
      dxne_dtt = dxnidtt*zeff + 2.0d0*dxnidt*dzeff_dt + xni*dzeff_dtt
      dxne_dta = dxnidta*zeff + dxnidt*dzeff_da &
                 + dxnida*dzeff_dt + xni*dzeff_dta
      dxne_dtz = dxnidtz*zeff + dxnidt*dzeff_dz &
                 + dxnidz*dzeff_dt + xni*dzeff_dtz
      dxne_daa = dxnidaa*zeff + 2.0d0*dxnida*dzeff_da + xni*dzeff_daa
      dxne_daz = dxnidaz*zeff + dxnida*dzeff_dz &
                 + dxnidz*dzeff_da + xni*dzeff_daz
      dxne_dzz = dxnidzz*zeff + 2.0d0*dxnidz*dzeff_dz + xni*dzeff_dzz
      dxne_deta_dd = dxnidd*dzeff_deta + xni*dzeff_deta_dd
      dxne_deta_dt = dxnidt*dzeff_deta + xni*dzeff_deta_dt
      dxne_deta_dz = dxnidz*dzeff_deta + xni*dzeff_deta_dz
      dxne_deta_da = dxnida*dzeff_deta + xni*dzeff_deta_da
      dxne_deta2   = xni*dzeff_deta2



! get the fermi-dirac integral electron contribution
      call dfermi(0.5d0, etaele, beta, f12, f12eta, f12beta, &
                                       f12eta2, f12beta2, f12etabeta)
      call dfermi(1.5d0, etaele, beta, f32, f32eta, f32beta, &
                                       f32eta2, f32beta2, f32etabeta)

      zz   = xconst * beta32
      ww   = xconst * 1.5d0 * beta12
      yy   = f12 + beta * f32
      dum1 = f12eta + beta*f32eta
      dum2 = f12beta + f32 + beta*f32beta

      xnefer             = zz * yy
      dxnefer_deta       = zz * dum1
      dxnefer_dbeta      = ww*yy + zz*dum2
      dxnefer_deta2      = zz * (f12eta2 + beta * f32eta2)
      dxnefer_dbeta2     = 0.5d0*ww/beta*yy + 2.0d0*ww*dum2 &
                         + zz*(f12beta2 + 2.0d0*f32beta + beta*f32beta2)
      dxnefer_deta_dbeta = ww*dum1 &
                         + zz*(f12etabeta + f32eta + beta*f32etabeta)



! if the temperature is not too low, get the positron contributions
! chemical equilibrium means etaele + etapos = eta_photon = 0.
      etapos             = 0.0d0
      detap_deta         = 0.0d0
      detap_dbeta        = 0.0d0
      detap_deta2        = 0.0d0
      detap_dbeta2       = 0.0d0
      detap_deta_dbeta   = 0.0d0
      xnpfer             = 0.0d0
      dxnpfer_deta       = 0.0d0
      dxnpfer_dbeta      = 0.0d0
      dxnpfer_deta2      = 0.0d0
      dxnpfer_dbeta2     = 0.0d0
      dxnpfer_deta_dbeta = 0.0d0

      if (beta .gt. positron_start) then
       etapos           = -aa - 2.0d0/beta
       detap_deta       = -1.0d0
       detap_deta2      = 0.0d0
       detap_dbeta      = 2.0d0/beta**2
       detap_dbeta2     = -4.0d0/beta**3
       detap_deta_dbeta = 0.0d0

       call dfermi(0.5d0, etapos, beta, f12, f12eta, f12beta, &
                                        f12eta2, f12beta2, f12etabeta)
       call dfermi(1.5d0, etapos, beta, f32, f32eta, f32beta, &
                                        f32eta2, f32beta2, f32etabeta)

       zz   = xconst * beta32
       ww   = xconst * 1.5d0 * beta12
       yy   = f12 + beta * f32
       dum1 = f12eta + beta*f32eta
       dum2 = f12beta + f32 + beta*f32beta

       xnpfer              = zz * yy
       dxnpfer_detap       = zz * dum1
       dxnpfer_dbeta      = ww*yy + zz*dum2
       dxnpfer_detap2      = zz * (f12eta2 + beta * f32eta2)
       dxnpfer_dbeta2     = 0.5d0*ww/beta*yy + 2.0d0*ww*dum2 &
                         + zz*(f12beta2 + 2.0d0*f32beta + beta*f32beta2)
       dxnpfer_detap_dbeta = ww*dum1 &
                         + zz*(f12etabeta + f32eta + beta*f32etabeta)


! convert the etap derivatives to eta derivatives
! all derived from the operator dxp = dxp/detap detap + dxp/dbeta dbeta

       dxnpfer_deta  = dxnpfer_detap * detap_deta
       dxnpfer_dbeta = dxnpfer_dbeta + dxnpfer_detap * detap_dbeta
       dxnpfer_deta2 = dxnpfer_detap2 * detap_deta**2 &
                       + dxnpfer_detap * detap_deta2
       dxnpfer_dbeta2 = dxnpfer_dbeta2 &
                       + 2.0d0 * dxnpfer_detap_dbeta * detap_dbeta &
                       + dxnpfer_detap2 * detap_dbeta**2 &
                       + dxnpfer_detap * detap_dbeta2
       dxnpfer_deta_dbeta = dxnpfer_detap2 * detap_dbeta * detap_deta &
                       + dxnpfer_detap_dbeta * detap_deta &
                       + dxnpfer_detap * detap_deta_dbeta

      end if



! charge neutrality means ne_ionization = ne_electrons - ne_positrons
      f  = xnefer - xnpfer - xne


! derivative of f with eta for newton-like root finders
      df = dxnefer_deta  - dxnpfer_deta  - dxne_deta



! if we are in root finder mode, return


      if (mode .eq. 0) return





! if we are not in root finder mode, polish off the calculation


! all the derivatives are in terms of eta and beta.
! we want to convert to temperature, density, abar and zbar derivatives.
! so, after the root find above on eta we have xne = xnefer - xnpfer.
! taking the derivative of this for property p
! dxne/deta deta + dxne/dp dp = dxnefer/deta deta + dxnefer/dbeta dbeta
! solving for the unknown eta derivative yields
! deta/dp = (dxne_dp - dxnefer/dbeta dbeta/dp) / (dxnefer/deta - dxne/deta)


      dxep_deta       = dxnefer_deta  - dxnpfer_deta
      dxep_dbeta      = dxnefer_dbeta - dxnpfer_dbeta
      dxep_deta2      = dxnefer_deta2  - dxnpfer_deta2
      dxep_dbeta2     = dxnefer_dbeta2 - dxnpfer_dbeta2
      dxep_deta_dbeta = dxnefer_deta_dbeta  - dxnpfer_deta_dbeta

      y          = 1.0d0/(dxep_deta - dxne_deta)

! the all important first derivatives of eta
      detadd = dxne_dd * y
      detadt = (dxne_dt - dxep_dbeta*dbetadt) * y
      detada = dxne_da * y
      detadz = dxne_dz * y

! second derivatives
      detaddd = dxne_ddd*y - detadd*y*(dxep_deta2 - dxne_deta2)*detadd
      detaddt = dxne_ddt*y - detadd*y*(dxep_deta2*detadt &
                + dxep_deta_dbeta*dbetadt - dxne_deta2*detadt)
      detadda = dxne_dda*y - detadd*y*(dxep_deta2 - dxne_deta2)*detada
      detaddz = dxne_ddz*y - detadd*y*(dxep_deta2 - dxne_deta2)*detadz

      detadtt = (dxne_dtt - (dxep_deta_dbeta*detadt &
                 + dxep_dbeta2*dbetadt)*dbetadt)*y &
                 - detadt*y*(dxep_deta2*detadt &
                         + dxep_deta_dbeta*dbetadt - dxne_deta2*detadt)

      detadta = (dxne_dta - dxep_deta_dbeta*detada*dbetadt)*y &
                 - detadt*y*(dxep_deta2*detada - dxne_deta2*detada)

      detadtz = (dxne_dtz - dxep_deta_dbeta*detadz*dbetadt)*y &
                 - detadt*y*(dxep_deta2*detadz - dxne_deta2*detadz)

      detadaa = dxne_daa*y - detada*y*(dxep_deta2 - dxne_deta2)*detada
      detadaz = dxne_daz*y - detada*y*(dxep_deta2 - dxne_deta2)*detadz
      detadzz = dxne_dzz*y - detadz*y*(dxep_deta2 - dxne_deta2)*detadz



! first derivatives of the effective charge
      dzeffdd = dzeff_deta*detadd + dzeff_dd
      dzeffdt = dzeff_deta*detadt + dzeff_dt
      dzeffda = dzeff_deta*detada + dzeff_da
      dzeffdz = dzeff_deta*detadz + dzeff_dz

! second derivatives
      dzeffddd = dzeff_deta_dd*detadd + dzeff_deta*detaddd + dzeff_ddd
      dzeffddt = dzeff_deta_dt*detadd + dzeff_deta*detaddt + dzeff_ddt
      dzeffdda = dzeff_deta_da*detadd + dzeff_deta*detadda + dzeff_dda
      dzeffddz = dzeff_deta_dz*detadd + dzeff_deta*detaddz + dzeff_ddz
      dzeffdtt = dzeff_deta_dt*detadt + dzeff_deta*detadtt + dzeff_dtt
      dzeffdta = dzeff_deta_da*detadt + dzeff_deta*detadta + dzeff_dta
      dzeffdtz = dzeff_deta_dz*detadt + dzeff_deta*detadtz + dzeff_dtz
      dzeffdaa = dzeff_deta_da*detada + dzeff_deta*detadaa + dzeff_daa
      dzeffdaz = dzeff_deta_dz*detada + dzeff_deta*detadaz + dzeff_daz
      dzeffdzz = dzeff_deta_dz*detadz + dzeff_deta*detadzz + dzeff_dzz


! first derivatives of the electron number density
      dxnedd = dxnidd * zeff + xni * dzeffdd
      dxnedt = dxnidt * zeff + xni * dzeffdt
      dxneda = dxnida * zeff + xni * dzeffda
      dxnedz = dxnidz * zeff + xni * dzeffdz

! second derivatives
      dxneddd = dxniddd*zeff+dxnidd*dzeffdd+dxnidd*dzeffdd+xni*dzeffddd
      dxneddt = dxniddt*zeff+dxnidd*dzeffdt+dxnidt*dzeffdd+xni*dzeffddt
      dxnedda = dxnidda*zeff+dxnidd*dzeffda+dxnida*dzeffdd+xni*dzeffdda
      dxneddz = dxniddz*zeff+dxnidd*dzeffdz+dxnidz*dzeffdd+xni*dzeffddz
      dxnedtt = dxnidtt*zeff+dxnidt*dzeffdt+dxnidt*dzeffdt+xni*dzeffdtt
      dxnedta = dxnidta*zeff+dxnidt*dzeffda+dxnida*dzeffdt+xni*dzeffdta
      dxnedtz = dxnidtz*zeff+dxnidt*dzeffdz+dxnidz*dzeffdt+xni*dzeffdtz
      dxnedaa = dxnidaa*zeff+dxnida*dzeffda+dxnida*dzeffda+xni*dzeffdaa
      dxnedaz = dxnidaz*zeff+dxnida*dzeffdz+dxnidz*dzeffda+xni*dzeffdaz
      dxnedzz = dxnidzz*zeff+dxnidz*dzeffdz+dxnidz*dzeffdz+xni*dzeffdzz



! first derivatives of the fermi integral electron number densities
      dxneferdd = dxnefer_deta * detadd
      dxneferdt = dxnefer_deta * detadt + dxnefer_dbeta * dbetadt
      dxneferda = dxnefer_deta * detada
      dxneferdz = dxnefer_deta * detadz

! second derivatives
      dxneferddd = dxnefer_deta2*detadd*detadd + dxnefer_deta*detaddd
      dxneferddt = dxnefer_deta2*detadt*detadd &
                   + dxnefer_deta_dbeta * dbetadt * detadd &
                   + dxnefer_deta *detaddt
      dxneferdda = dxnefer_deta2*detada*detadd + dxnefer_deta*detadda
      dxneferddz = dxnefer_deta2*detadz*detadd + dxnefer_deta*detaddz
      dxneferdtt = dxnefer_deta2*detadt*detadt + dxnefer_deta*detadtt &
                   + 2.0d0*dxnefer_deta_dbeta * detadt * dbetadt &
                   + dxnefer_dbeta2 * dbetadt * dbetadt
      dxneferdta = dxnefer_deta2 * detada * detadt &
                   + dxnefer_deta_dbeta * dbetadt * detada &
                   + dxnefer_deta * detadta
      dxneferdtz = dxnefer_deta2 * detadz * detadt &
                   + dxnefer_deta_dbeta * dbetadt * detadz &
                   + dxnefer_deta * detadtz
      dxneferdaa = dxnefer_deta2*detada*detada + dxnefer_deta*detadaa
      dxneferdaz = dxnefer_deta2*detadz*detada + dxnefer_deta*detadaz
      dxneferdzz = dxnefer_deta2*detadz*detadz + dxnefer_deta*detadzz



! first derivatives of the fermi integral positron number densities
      dxnpferdd = dxnpfer_deta * detadd
      dxnpferdt = dxnpfer_deta * detadt + dxnpfer_dbeta * dbetadt
      dxnpferda = dxnpfer_deta * detada
      dxnpferdz = dxnpfer_deta * detadz


! second derivatives
      dxnpferddd = dxnpfer_deta2*detadd*detadd + dxnpfer_deta*detaddd
      dxnpferddt = dxnpfer_deta2*detadt*detadd &
                   + dxnpfer_deta_dbeta * dbetadt * detadd &
                   + dxnpfer_deta *detaddt
      dxnpferdda = dxnpfer_deta2*detada*detadd + dxnpfer_deta*detadda
      dxnpferddz = dxnpfer_deta2*detadz*detadd + dxnpfer_deta*detaddz
      dxnpferdtt = dxnpfer_deta2*detadt*detadt + dxnpfer_deta*detadtt &
                   + 2.0d0*dxnpfer_deta_dbeta * detadt * dbetadt &
                   + dxnpfer_dbeta2 * dbetadt * dbetadt
      dxnpferdta = dxnpfer_deta2 * detada * detadt &
                   + dxnpfer_deta_dbeta * dbetadt * detada &
                   + dxnpfer_deta * detadta
      dxnpferdtz = dxnpfer_deta2 * detadz * detadt &
                   + dxnpfer_deta_dbeta * dbetadt * detadz &
                   + dxnpfer_deta * detadtz
      dxnpferdaa = dxnpfer_deta2*detada*detada + dxnpfer_deta*detadaa
      dxnpferdaz = dxnpfer_deta2*detadz*detada + dxnpfer_deta*detadaz
      dxnpferdzz = dxnpfer_deta2*detadz*detadz + dxnpfer_deta*detadzz





! now get the pressure and energy
! for the electrons

      call dfermi(1.5d0, etaele, beta, f32, f32eta, f32beta, &
                                       f32eta2, f32beta2, f32etabeta)
      call dfermi(2.5d0, etaele, beta, f52, f52eta, f52beta, &
                                       f52eta2, f52beta2, f52etabeta)

! pressure in erg/cm**3
      yy   = pconst * beta52
      ww   = pconst * 2.5d0 * beta32
      dum1 = f32 + 0.5d0*beta*f52
      dum2 = f32beta + 0.5d0*f52 + 0.5d0*beta*f52beta
      dum3 = f32eta + 0.5d0*beta*f52eta

      pele         = yy * dum1
      dpele_deta   = yy * dum3
      dpele_dbeta  = ww * dum1 + yy * dum2
      dpele_deta2  = yy * (f32eta2 + 0.5d0*beta*f52eta2)
      dpele_dbeta2 = pconst*beta12*3.75d0*dum1 + 2.0d0*ww*dum2 &
                     + yy * (f32beta2 + f52beta + 0.5d0*beta*f52beta2)
      dpele_deta_dbeta = ww*dum3 + yy*(f32etabeta +0.5d0*f52eta &
                     + 0.5d0*beta*f52etabeta)


! first derivatives of the electron pressure
      dpeledd = dpele_deta * detadd
      dpeledt = dpele_deta * detadt + dpele_dbeta * dbetadt
      dpeleda = dpele_deta * detada
      dpeledz = dpele_deta * detadz

! second derivatives
      dpeleddd = dpele_deta2*detadd*detadd + dpele_deta*detaddd
      dpeleddt = dpele_deta2*detadt*detadd &
                   + dpele_deta_dbeta * dbetadt * detadd &
                   + dpele_deta *detaddt
      dpeledda = dpele_deta2*detada*detadd + dpele_deta*detadda
      dpeleddz = dpele_deta2*detadz*detadd + dpele_deta*detaddz
      dpeledtt = dpele_deta2*detadt*detadt + dpele_deta*detadtt &
                   + 2.0d0*dpele_deta_dbeta * detadt * dbetadt &
                   + dpele_dbeta2 * dbetadt * dbetadt
      dpeledta = dpele_deta2 * detada * detadt &
                   + dpele_deta_dbeta * dbetadt * detada &
                   + dpele_deta * detadta
      dpeledtz = dpele_deta2 * detadz * detadt &
                   + dpele_deta_dbeta * dbetadt * detadz &
                   + dpele_deta * detadtz
      dpeledaa = dpele_deta2*detada*detada + dpele_deta*detadaa
      dpeledaz = dpele_deta2*detadz*detada + dpele_deta*detadaz
      dpeledzz = dpele_deta2*detadz*detadz + dpele_deta*detadzz


! energy in erg/cm**3
      yy   = econst * beta52
      ww   = econst * 2.5d0 * beta32
      dum1 = f32 + beta*f52
      dum2 = f32beta + f52 + beta*f52beta
      dum3 = f32eta + beta*f52eta

      eele        = yy * dum1
      deele_deta  = yy * dum3
      deele_dbeta = ww * dum1 + yy * dum2
      deele_deta2  = yy * (f32eta2 + beta*f52eta2)
      deele_dbeta2 = econst*beta12*3.75d0*dum1 + 2.0d0*ww*dum2 &
                      + yy * (f32beta2 + 2.0d0*f52beta + beta*f52beta2)
      deele_deta_dbeta = ww*dum3 &
                      + yy*(f32etabeta + f52eta + beta*f52etabeta)


! first derivatives of the electron energy
      deeledd = deele_deta * detadd
      deeledt = deele_deta * detadt + deele_dbeta * dbetadt
      deeleda = deele_deta * detada
      deeledz = deele_deta * detadz

! second derivatives
      deeleddd = deele_deta2*detadd*detadd + deele_deta*detaddd
      deeleddt = deele_deta2*detadt*detadd &
                   + deele_deta_dbeta * dbetadt * detadd &
                   + deele_deta *detaddt
      deeledda = deele_deta2*detada*detadd + deele_deta*detadda
      deeleddz = deele_deta2*detadz*detadd + deele_deta*detaddz
      deeledtt = deele_deta2*detadt*detadt + deele_deta*detadtt &
                   + 2.0d0*deele_deta_dbeta * detadt * dbetadt &
                   + deele_dbeta2 * dbetadt * dbetadt
      deeledta = deele_deta2 * detada * detadt &
                   + deele_deta_dbeta * dbetadt * detada &
                   + deele_deta * detadta
      deeledtz = deele_deta2 * detadz * detadt &
                   + deele_deta_dbeta * dbetadt * detadz &
                   + deele_deta * detadtz
      deeledaa = deele_deta2*detada*detada + deele_deta*detadaa
      deeledaz = deele_deta2*detadz*detada + deele_deta*detadaz
      deeledzz = deele_deta2*detadz*detadz + deele_deta*detadzz





! for the positrons
      ppos              = 0.0d0
      dppos_detap       = 0.0d0
      dppos_dbeta       = 0.0d0
      dppos_detap2      = 0.0d0
      dppos_dbeta2      = 0.0d0
      dppos_detap_dbeta = 0.0d0
      epos              = 0.0d0
      depos_detap       = 0.0d0
      depos_dbeta       = 0.0d0
      depos_detap2      = 0.0d0
      depos_dbeta2      = 0.0d0
      depos_detap_dbeta = 0.0d0


      if (beta .gt. positron_start) then
       call dfermi(1.5d0, etapos, beta, f32, f32eta, f32beta, &
                                        f32eta2, f32beta2, f32etabeta)
       call dfermi(2.5d0, etapos, beta, f52, f52eta, f52beta, &
                                        f52eta2, f52beta2, f52etabeta)

! pressure
       yy   = pconst * beta52
       ww   = pconst * 2.5d0 * beta32
       dum1 = f32 + 0.5d0*beta*f52
       dum2 = f32beta + 0.5d0*f52 + 0.5d0*beta*f52beta
       dum3 = f32eta + 0.5d0*beta*f52eta

       ppos         = yy * dum1
       dppos_detap  = yy * dum3
       dppos_dbeta  = ww * dum1 + yy * dum2
       dppos_detap2 = yy * (f32eta2 + 0.5d0*beta*f52eta2)
       dppos_dbeta2 = pconst*beta12*3.75d0*dum1 + 2.0d0*ww*dum2 &
                      + yy * (f32beta2 + f52beta + 0.5d0*beta*f52beta2)
       dppos_detap_dbeta = ww*dum3 + yy*(f32etabeta +0.5d0*f52eta &
                      + 0.5d0*beta*f52etabeta)

! energy
       yy   = econst * beta52
       ww   = econst * 2.5d0 * beta32
       dum1 = f32 + beta*f52
       dum2 = f32beta + f52 + beta*f52beta
       dum3 = f32eta + beta*f52eta

       epos        = yy * dum1
       depos_detap  = yy * dum3
       depos_dbeta = ww * dum1 + yy * dum2
       depos_detap2  = yy * (f32eta2 + beta*f52eta2)
       depos_dbeta2 = econst*beta12*3.75d0*dum1 + 2.0d0*ww*dum2 &
                       + yy * (f32beta2 + 2.0d0*f52beta + beta*f52beta2)
       depos_detap_dbeta = ww*dum3 &
                       + yy*(f32etabeta + f52eta + beta*f52etabeta)
      end if



! convert the etap derivatives to eta derivatives
! all derived from the operator dxp = dxp/detap detap + dxp/dbeta dbeta

       dppos_deta  = dppos_detap * detap_deta
       dppos_dbeta = dppos_dbeta + dppos_detap * detap_dbeta
       dppos_deta2 = dppos_detap2 * detap_deta**2 &
                       + dppos_detap * detap_deta2
       dppos_dbeta2 = dppos_dbeta2 &
                       + 2.0d0 * dppos_detap_dbeta * detap_dbeta &
                       + dppos_detap2 * detap_dbeta**2 &
                       + dppos_detap * detap_dbeta2
       dppos_deta_dbeta = dppos_detap2 * detap_dbeta * detap_deta &
                       + dppos_detap_dbeta * detap_deta &
                       + dppos_detap * detap_deta_dbeta

! first derivatives of the positron pressure
      dpposdd     = dppos_deta * detadd
      dpposdt     = dppos_deta * detadt + dppos_dbeta * dbetadt
      dpposda     = dppos_deta * detada
      dpposdz     = dppos_deta * detadz


! second derivatives
      dpposddd = dppos_deta2*detadd*detadd + dppos_deta*detaddd
      dpposddt = dppos_deta2*detadt*detadd &
                   + dppos_deta_dbeta * dbetadt * detadd &
                   + dppos_deta *detaddt
      dpposdda = dppos_deta2*detada*detadd + dppos_deta*detadda
      dpposddz = dppos_deta2*detadz*detadd + dppos_deta*detaddz
      dpposdtt = dppos_deta2*detadt*detadt + dppos_deta*detadtt &
                   + 2.0d0*dppos_deta_dbeta * detadt * dbetadt &
                   + dppos_dbeta2 * dbetadt * dbetadt
      dpposdta = dppos_deta2 * detada * detadt &
                   + dppos_deta_dbeta * dbetadt * detada &
                   + dppos_deta * detadta
      dpposdtz = dppos_deta2 * detadz * detadt &
                   + dppos_deta_dbeta * dbetadt * detadz &
                   + dppos_deta * detadtz
      dpposdaa = dppos_deta2*detada*detada + dppos_deta*detadaa
      dpposdaz = dppos_deta2*detadz*detada + dppos_deta*detadaz
      dpposdzz = dppos_deta2*detadz*detadz + dppos_deta*detadzz



! convert the etap derivatives to eta derivatives
! all derived from the operator dxp = dxp/detap detap + dxp/dbeta dbeta

       depos_deta  = depos_detap * detap_deta
       depos_dbeta = depos_dbeta + depos_detap * detap_dbeta
       depos_deta2 = depos_detap2 * detap_deta**2 &
                       + depos_detap * detap_deta2
       depos_dbeta2 = depos_dbeta2 &
                       + 2.0d0 * depos_detap_dbeta * detap_dbeta &
                       + depos_detap2 * detap_dbeta**2 &
                       + depos_detap * detap_dbeta2
       depos_deta_dbeta = depos_detap2 * detap_dbeta * detap_deta &
                       + depos_detap_dbeta * detap_deta &
                       + depos_detap * detap_deta_dbeta


! first derivatives of the positron energy
      deposdd     = depos_deta * detadd
      deposdt     = depos_deta * detadt + depos_dbeta * dbetadt
      deposda     = depos_deta * detada
      deposdz     = depos_deta * detadz

! second derivatives
      deposddd = depos_deta2*detadd*detadd + depos_deta*detaddd
      deposddt = depos_deta2*detadt*detadd &
                   + depos_deta_dbeta * dbetadt * detadd &
                   + depos_deta *detaddt
      deposdda = depos_deta2*detada*detadd + depos_deta*detadda
      deposddz = depos_deta2*detadz*detadd + depos_deta*detaddz
      deposdtt = depos_deta2*detadt*detadt + depos_deta*detadtt &
                   + 2.0d0*depos_deta_dbeta * detadt * dbetadt &
                   + depos_dbeta2 * dbetadt * dbetadt
      deposdta = depos_deta2 * detada * detadt &
                   + depos_deta_dbeta * dbetadt * detada &
                   + depos_deta * detadta
      deposdtz = depos_deta2 * detadz * detadt &
                   + depos_deta_dbeta * dbetadt * detadz &
                   + depos_deta * detadtz
      deposdaa = depos_deta2*detada*detada + depos_deta*detadaa
      deposdaz = depos_deta2*detadz*detada + depos_deta*detadaz
      deposdzz = depos_deta2*detadz*detadz + depos_deta*detadzz




! electron+positron pressure and its derivatives
! note: at high temperatures and low densities, dpepdd is very small
! and can go negative, so limit it to be positive definite

      pep     = pele    + ppos
      dpepdd  = max(dpeledd + dpposdd, 1.0d-30)
      dpepdt  = dpeledt + dpposdt
      dpepda  = dpeleda + dpposda
      dpepdz  = dpeledz + dpposdz
      dpepddd = dpeleddd + dpposddd
      dpepddt = dpeleddt + dpposddt
      dpepdda = dpeledda + dpposdda
      dpepddz = dpeleddz + dpposddz
      dpepdtt = dpeledtt + dpposdtt
      dpepdta = dpeledta + dpposdta
      dpepdtz = dpeledtz + dpposdtz
      dpepdaa = dpeledaa + dpposdaa
      dpepdaz = dpeledaz + dpposdaz
      dpepdzz = dpeledzz + dpposdzz


! electron+positron thermal energy and its derivatives
      eep     = eele    + epos
      deepdd  = deeledd + deposdd
      deepdt  = deeledt + deposdt
      deepda  = deeleda + deposda
      deepdz  = deeledz + deposdz
      deepddd = deeleddd + deposddd
      deepddt = deeleddt + deposddt
      deepdda = deeledda + deposdda
      deepddz = deeleddz + deposddz
      deepdtt = deeledtt + deposdtt
      deepdta = deeledta + deposdta
      deepdtz = deeledtz + deposdtz
      deepdaa = deeledaa + deposdaa
      deepdaz = deeledaz + deposdaz
      deepdzz = deeledzz + deposdzz



! electron entropy in erg/gr/kelvin and its derivatives
      y       = kerg*deni
      sele    = ((pele + eele)*kti - etaele*xnefer) * y


! first derivatives
      dseledd = ((dpeledd + deeledd)*kti - detadd*xnefer &
                  - etaele*dxneferdd)*y - sele*deni
      dseledt = ((dpeledt + deeledt)*kti - (pele + eele)/(kt*temp) &
                  - detadt*xnefer - etaele*dxneferdt)*y
      dseleda = ((dpeleda + deeleda)*kti - detada*xnefer &
                   - etaele*dxneferda)*y
      dseledz = ((dpeledz + deeledz)*kti - detadz*xnefer &
                   - etaele*dxneferdz)*y

! second derivatives
      dseleddd = ((dpeleddd + deeleddd)*kti - detaddd*xnefer &
                  - 2.0d0*detadd*dxneferdd - etaele*dxneferddd)*y &
                  - 2.0d0*((dpeledd + deeledd)*kti - detadd*xnefer &
                      - etaele*dxneferdd)*y*deni &
                  + 2.0d0*sele*deni**2
      dseleddt = ((dpeleddt + deeleddt)*kti &
                   - (dpeledd + deeledd)*kti/temp &
                   - detaddt*xnefer - detadd*dxneferdt &
                   - detadt*dxneferdd - etaele*dxneferddt)*y &
                   - dseledt*deni
      dseledda = ((dpeledda + deeledda)*kti &
                   - detadda*xnefer - detadd*dxneferda &
                   - detada*dxneferdd - etaele*dxneferdda)*y &
                   - dseleda*deni
      dseleddz = ((dpeleddz + deeleddz)*kti &
                   - detaddz*xnefer - detadd*dxneferdz &
                   - detadz*dxneferdd - etaele*dxneferddz)*y &
                   - dseledz*deni
      dseledtt = ((dpeledtt + deeledtt)*kti &
                 - 2.0d0*(dpeledt + deeledt)*kti/temp &
                 + 2.0d0*(pele + eele)*kti/temp**2 &
                 - detadtt*xnefer - 2.0d0*detadt*dxneferdt &
                 - etaele*dxneferdtt)*y
      dseledta = ((dpeledta + deeledta)*kti &
                  - (dpeleda + deeleda)*kti/temp &
                  - detadta*xnefer - detadt*dxneferda &
                  - detada*dxneferdt - etaele*dxneferdta)*y
      dseledtz = ((dpeledtz + deeledtz)*kti &
                  - (dpeledz + deeledz)*kti/temp &
                  - detadtz*xnefer - detadt*dxneferdz &
                  - detadz*dxneferdt - etaele*dxneferdtz)*y
      dseledaa = ((dpeledaa + deeledaa)*kti - detadaa*xnefer &
                 - 2.0d0*detada*dxneferda - etaele*dxneferdaa)*y
      dseledaz = ((dpeledaz + deeledaz)*kti - detadaz*xnefer &
                   - detada*dxneferdz &
                   - detadz*dxneferda  - etaele*dxneferdaz)*y
      dseledzz = ((dpeledzz + deeledzz)*kti - detadzz*xnefer &
                 - 2.0d0*detadz*dxneferdz - etaele*dxneferdzz)*y



! positron entropy in erg/gr/kelvin and its derivatives
      spos    = ((ppos + epos)/kt - etapos*xnpfer) * y

! first derivatives
      dsposdd = ((dpposdd + deposdd)*kti &
                 - detap_deta*detadd*xnpfer &
                 - etapos*dxnpferdd)*y - spos*deni
      dsposdt = ((dpposdt + deposdt)*kti - (ppos + epos)/(kt*temp) &
                 - (detap_deta*detadt + detap_dbeta*dbetadt)*xnpfer &
                 - etapos*dxnpferdt)*y
      dsposda = ((dpposda + deposda)*kti &
                 - detap_deta*detada*xnpfer &
                 - etapos*dxnpferda)*y
      dsposdz = ((dpposdz + deposdz)*kti - detap_deta*detadz*xnpfer &
                   - etapos*dxnpferdz)*y

! second derivatives
      dsposddd = ((dpposddd + deposddd)*kti &
                  - detap_deta*detaddd*xnpfer &
                  - 2.0d0*detap_deta*detadd*dxnpferdd &
                  - etapos*dxnpferddd)*y &
                  - 2.0d0*((dpposdd + deposdd)*kti &
                  - detap_deta*detadd*xnpfer &
                  - etapos*dxnpferdd)*y*deni &
                  + 2.0d0*spos*deni**2
      dsposddt = ((dpposddt + deposddt)*kti &
                   - (dpposdd + deposdd)*kti/temp &
                   - detap_deta*detaddt*xnpfer &
                   - detap_deta*detadd*dxnpferdt &
                   - (detap_deta*detadt + detap_dbeta*dbetadt)*dxnpferdd &
                   - etapos*dxnpferddt)*y &
                   - dsposdt*deni
      dsposdda = ((dpposdda + deposdda)*kti &
                   - detap_deta*detadda*xnpfer &
                   - detap_deta*detadd*dxnpferda &
                   - detap_deta*detada*dxnpferdd &
                   - etapos*dxnpferdda)*y &
                   - dsposda*deni
      dsposddz = ((dpposddz + deposddz)*kti &
                   - detap_deta*detaddz*xnpfer &
                   - detap_deta*detadd*dxnpferdz &
                   - detap_deta*detadz*dxnpferdd &
                   - etapos*dxnpferddz)*y &
                   - dsposdz*deni

      dsposdtt = ((dpposdtt + deposdtt)*kti &
                 - 2.0d0*(dpposdt + deposdt)*kti/temp &
                 + 2.0d0*(ppos + epos)*kti/temp**2 &
                 - (detap_deta2*detadt**2 &
                    + 2.0d0*detap_deta_dbeta*detadt*dbetadt &
                    + detap_deta*detadtt &
                    + detap_dbeta2*dbetadt**2)*xnpfer &
                 - 2.0d0*(detap_deta*detadt &
                            + detap_dbeta*dbetadt)*dxnpferdt &
                 - etapos*dxnpferdtt)*y

!      dsposdt = ((dpposdt + deposdt)*kti - (ppos + epos)/(kt*temp)
!     1           - (detap_deta*detadt + detap_dbeta*dbetadt)*xnpfer
!     2           - etapos*dxnpferdt)*y

      dsposdta = ((dpposdta + deposdta)*kti &
                  - (dpposda + deposda)*kti/temp &
                  - (detap_deta2*detadt*detada &
                      + detap_deta_dbeta*dbetadt*detada &
                      + detap_deta*detadta)*xnpfer &
                  - (detap_deta*detadt &
                            + detap_dbeta*dbetadt)*dxnpferda &
                  - detap_deta*detada*dxnpferdt &
                  - etapos*dxnpferdta)*y
      dsposdtz = ((dpposdtz + deposdtz)*kti &
                  - (dpposdz + deposdz)*kti/temp &
                  - (detap_deta2*detadt*detadz &
                      + detap_deta_dbeta*dbetadt*detadz &
                      + detap_deta*detadtz)*xnpfer &
                  - (detap_deta*detadt &
                            + detap_dbeta*dbetadt)*dxnpferdz &
                  - detap_deta*detadz*dxnpferdt &
                  - etapos*dxnpferdtz)*y
      dsposdaa = ((dpposdaa + deposdaa)*kti &
                  - detap_deta*detadaa*xnpfer &
                  - 2.0d0*detap_deta*detada*dxnpferda &
                  - etapos*dxnpferdaa)*y
      dsposdaz = ((dpposdaz + deposdaz)*kti &
                  - detap_deta*detadaz*xnpfer &
                  - detap_deta*detada*dxnpferdz &
                  - detap_deta*detadz*dxnpferda &
                  - etapos*dxnpferdaz)*y
      dsposdzz = ((dpposdzz + deposdzz)*kti &
                 - detap_deta*detadzz*xnpfer &
                 - 2.0d0*detap_deta*detadz*dxnpferdz &
                 - etapos*dxnpferdzz)*y


! and their sum
      sep      = sele + spos
      dsepdd   = dseledd + dsposdd
      dsepdt   = dseledt + dsposdt
      dsepda   = dseleda + dsposda
      dsepdz   = dseledz + dsposdz
      dsepddd  = dseleddd + dsposddd
      dsepddt  = dseleddt + dsposddt
      dsepdda  = dseledda + dsposdda
      dsepddz  = dseleddz + dsposddz
      dsepdtt  = dseledtt + dsposdtt
      dsepdta  = dseledta + dsposdta
      dsepdtz  = dseledtz + dsposdtz
      dsepdaa  = dseledaa + dsposdaa
      dsepdaz  = dseledaz + dsposdaz
      dsepdzz  = dseledzz + dsposdzz



! adjust for the rest mass energy of the positrons
      y        = 2.0d0 * mecc
      epos     = epos     + y * xnpfer
      deposdd  = deposdd  + y * dxnpferdd
      deposdt  = deposdt  + y * dxnpferdt
      deposda  = deposda  + y * dxnpferda
      deposdz  = deposdz  + y * dxnpferdz
      deposddd = deposddd + y * dxnpferddd
      deposddt = deposddt + y * dxnpferddt
      deposdda = deposdda + y * dxnpferdda
      deposddz = deposddz + y * dxnpferddz
      deposdtt = deposdtt + y * dxnpferdtt
      deposdta = deposdta + y * dxnpferdta
      deposdtz = deposdtz + y * dxnpferdtz
      deposdaa = deposdaa + y * dxnpferdaa
      deposdaz = deposdaz + y * dxnpferdaz
      deposdzz = deposdzz + y * dxnpferdzz


! and resum
      deepdd  = deeledd + deposdd
      deepdt  = deeledt + deposdt
      deepda  = deeleda + deposda
      deepdz  = deeledz + deposdz
      deepddd = deeleddd + deposddd
      deepddt = deeleddt + deposddt
      deepdda = deeledda + deposdda
      deepddz = deeleddz + deposddz
      deepdtt = deeledtt + deposdtt
      deepdta = deeledta + deposdta
      deepdtz = deeledtz + deposdtz
      deepdaa = deeledaa + deposdaa
      deepdaz = deeledaz + deposdaz
      deepdzz = deeledzz + deposdzz



! convert the electron-positron thermal energy in erg/cm**3
! to a specific thermal energy in erg/gr

      eele     = eele*deni
      deeledd  = (deeledd - eele)*deni
      deeledt  = deeledt*deni
      deeleda  = deeleda*deni
      deeledz  = deeledz*deni
      deeleddd = (deeleddd - 2.0d0*deeledd)*deni
      deeleddt = (deeleddt - deeledt)*deni
      deeledda = (deeledda - deeleda)*deni
      deeleddz = (deeleddz - deeledz)*deni
      deeledtt = deeledtt*deni
      deeledta = deeledta*deni
      deeledtz = deeledtz*deni
      deeledaa = deeledaa*deni
      deeledaz = deeledaz*deni
      deeledzz = deeledzz*deni

      epos     = epos*deni
      deposdd  = (deposdd - epos)*deni
      deposdt  = deposdt*deni
      deposda  = deposda*deni
      deposdz  = deposdz*deni
      deposddd = (deposddd - 2.0d0*deposdd)*deni
      deposddt = (deposddt - deposdt)*deni
      deposdda = (deposdda - deposda)*deni
      deposddz = (deposddz - deposdz)*deni
      deposdtt = deposdtt*deni
      deposdta = deposdta*deni
      deposdtz = deposdtz*deni
      deposdaa = deposdaa*deni
      deposdaz = deposdaz*deni
      deposdzz = deposdzz*deni

! and resum
      deepdd  = deeledd + deposdd
      deepdt  = deeledt + deposdt
      deepda  = deeleda + deposda
      deepdz  = deeledz + deposdz
      deepddd = deeleddd + deposddd
      deepddt = deeleddt + deposddt
      deepdda = deeledda + deposdda
      deepddz = deeleddz + deposddz
      deepdtt = deeledtt + deposdtt
      deepdta = deeledta + deposdta
      deepdtz = deeledtz + deposdtz
      deepdaa = deeledaa + deposdaa
      deepdaz = deeledaz + deposdaz
      deepdzz = deeledzz + deposdzz



! and take care of the ionization potential contributions
      if (potmult .eq. 0) then
       eip     = 0.0d0
       deipdd  = 0.0d0
       deipdt  = 0.0d0
       deipda  = 0.0d0
       deipdz  = 0.0d0
       deipddd = 0.0d0
       deipddt = 0.0d0
       deipdda = 0.0d0
       deipddz = 0.0d0
       deipdtt = 0.0d0
       deipdta = 0.0d0
       deipdtz = 0.0d0
       deipdaa = 0.0d0
       deipdaz = 0.0d0
       deipdzz = 0.0d0

       sip     = 0.0d0
       dsipdd  = 0.0d0
       dsipdt  = 0.0d0
       dsipda  = 0.0d0
       dsipdz  = 0.0d0
       dsipddd = 0.0d0
       dsipddt = 0.0d0
       dsipdda = 0.0d0
       dsipddz = 0.0d0
       dsipdtt = 0.0d0
       dsipdta = 0.0d0
       dsipdtz = 0.0d0
       dsipdaa = 0.0d0
       dsipdaz = 0.0d0
       dsipdzz = 0.0d0

      else
       eip     = chi * xne
       deipdd  = chi * dxnedd
       deipdt  = chi * dxnedt
       deipda  = chi * dxneda
       deipdz  = chi * dxnedz + hion*ev2erg*xne
       deipddd = chi * dxneddd
       deipddt = chi * dxneddt
       deipdda = chi * dxnedda
       deipddz = chi * dxneddz
       deipdtt = chi * dxnedtt
       deipdta = chi * dxnedta
       deipdtz = chi * dxnedtz
       deipdaa = chi * dxnedaa
       deipdaz = chi * dxnedaz
       deipdzz = chi * dxnedzz + 2.0d0*hion*ev2erg*dxnedz


! the ionization entropy in erg/gr/kelvin and its derivatives
       y       = kerg*deni
       sip     = eip*kti*y
       dsipdd  = deipdd*kti*y - sip*deni
       dsipdt  = (deipdt*kti - eip*kti/temp)*y
       dsipda  = deipda*kti*y
       dsipdz  = deipdz*kti*y
       dsipddd = deipddd*kti*y - dsipdd*deni + sip*deni*deni
       dsipddt = deipddt*kti*y - dsipdt*deni
       dsipddt = deipdda*kti*y - dsipda*deni
       dsipddt = deipddz*kti*y - dsipdz*deni
       dsipdtt = (deipdtt*kti - 2.0d0*deipdt*kti/temp &
                   + 2.0d0*eip*kti/temp**2)*y
       dsipdta = (deipdta*kti - deipda*kti/temp)*y
       dsipdtz = (deipdtz*kti - deipdz*kti/temp)*y
       dsipdaa = deipdaa*kti*y
       dsipdaz = deipdaz*kti*y
       dsipdzz = deipdzz*kti*y

!       sip    = (eip/kt - etaele*xne) * y
!       dsipdd = (deipdd/kt
!     1            - detadd*xne)*y
!     2            - etaele*dxnedd*y
!     3            - sip*deni
!       dsipdt = (deipdt/kt
!     1             - detadt*xne
!     2             - etaele*dxnedt
!     3             - eip/(kt*temp))*y

! convert the ionization energy from erg/cm**3 to  erg/gr
       eip    = eip*deni
       deipdd = (deipdd - eip)*deni
       deipdt = deipdt*deni
       deipda = deipda*deni
       deipdz = deipdz*deni
       deipddd = (deipddd - 2.0d0*deipdd)*deni
       deipddt = (deipddt - deipdt)*deni
       deipdda = (deipdda - deipda)*deni
       deipddz = (deipddz - deipdz)*deni
       deipdtt = deipdtt*deni
       deipdta = deipdta*deni
       deipdtz = deipdtz*deni
       deipdaa = deipdaa*deni
       deipdaz = deipdaz*deni
       deipdzz = deipdzz*deni

! end of ionization energy block
      end if

      return
      end



      subroutine etages(xni,zbar,temp,eta)
      include 'implno.dek'
      include 'const.dek'

! this routine makes a damn good guess for the electron degeneracy
! parameter eta.
! input is the ion number density xni, average charge zbar,
! and temperature temp.
! output is a guess at the electron chemical potential eta


! declare the pass
      double precision  xni,zbar,temp,eta

! declare
      double precision xne,x,y,z,kt,beta,tmkt,xnefac

      double precision rt2,rt3,rtpi,cpf0,cpf1,cpf2,cpf3, &
                       twoth,fa0,forpi,mecc
      parameter        (rt2     = 1.4142135623730951d0, &
                        rt3     = 1.7320508075688772d0, &
                        rtpi    = 1.7724538509055159d0, &
                        cpf0    = h/(me*clight), &
                        cpf1    = 3.0d0/(8.0d0*pi) * cpf0**3, &
                        cpf2    = 4.0d0/cpf1, &
                        cpf3    = 2.0d0*rt3*rtpi/(rt2*cpf1), &
                        twoth   = 2.0d0/3.0d0, &
                        fa0     = 64.0d0/(9.0d0*pi), &
                        forpi   = 4.0d0 * pi, &
                        mecc    = me * clight * clight)

! notes: rt2=sqrt(2)  rt3=sqrt(3)  rtpi=sqrt(pi)


! for the purposes of guessing eta, assume full ionization
      xne   = xni * zbar
      kt    = kerg * temp
      beta  = kt/mecc


! number density of ionized electrons (c&g 24.354k) and number density at
! turning point (c&g 24.354i). if either of these exceed the number density
! as given by a saha equation, then pairs are important. set eta = -1/2.

      if (beta .ge. 1.0) then
       x = cpf2 * beta * beta
      else
       x = cpf3 * beta * (1.0d0 + 0.75d0*beta) * exp(-1.0d0/beta)
      end if
      if (x .ge. xne) then
       eta = -0.5d0


! get the dimensionless number density (c&g 24.313), if it is large apply the
! formula (c&g 24.309) to get a possible alfa, if not large do a two term
! binomial expansion on (c&g 24.309) to estimate eta.

      else
       z = (xne*cpf1)**twoth
       if (z .ge. 1.0e-6) then
        y = (sqrt(z + 1.0d0) - 1.0d0)/beta
       else
        y = z * (1.0d0 - z * 0.25d0) * 0.5d0/beta
       end if


! isolate the constant in front of the number density integral. if it is
! small enough run the divine approximation backwards with c&g 24.43. then
! join it smoothly with the lower limit.

       x = log10(xne**0.6d0/temp)
       if (x .le. 9.5) then
        z = ((1.0d0 + fa0*beta)*sqrt(1.0d0 + fa0*beta*0.5) - 1.0d0)/fa0
        tmkt    = 2.0d0 * me/h * kt/h
        xnefac  = forpi * tmkt * sqrt(tmkt)
        eta = -log(xnefac*rtpi*(0.5d0+0.75d0*z)/xne)
        if (x .ge. 8.5) eta = eta*(9.5d0-x) + y * (1.0d0 - (9.5d0-x))
       else
        eta = y
       end if
      end if

      return
      end






      subroutine zero_eos_vector(jlo,jhi)
      include 'implno.dek'
      include 'vector_eos.dek'
!
! this routine sets the pipeline eos vector to zero
!
! declare
      integer   jlo,jhi,j

      do j=jlo,jhi

        ptot_row(j) = 0.0d0
        dpt_row(j)  = 0.0d0
        dpd_row(j)  = 0.0d0
        dpa_row(j)  = 0.0d0
        dpz_row(j)  = 0.0d0
        dpdd_row(j) = 0.0d0
        dpdt_row(j) = 0.0d0
        dpda_row(j) = 0.0d0
        dpdz_row(j) = 0.0d0
        dptt_row(j) = 0.0d0
        dpta_row(j) = 0.0d0
        dptz_row(j) = 0.0d0
        dpaa_row(j) = 0.0d0
        dpaz_row(j) = 0.0d0
        dpzz_row(j) = 0.0d0

        etot_row(j) = 0.0d0
        det_row(j)  = 0.0d0
        ded_row(j)  = 0.0d0
        dea_row(j)  = 0.0d0
        dez_row(j)  = 0.0d0
        dedd_row(j) = 0.0d0
        dedt_row(j) = 0.0d0
        deda_row(j) = 0.0d0
        dedz_row(j) = 0.0d0
        dett_row(j) = 0.0d0
        deta_row(j) = 0.0d0
        detz_row(j) = 0.0d0
        deaa_row(j) = 0.0d0
        deaz_row(j) = 0.0d0
        dezz_row(j) = 0.0d0

        stot_row(j) = 0.0d0
        dst_row(j)  = 0.0d0
        dsd_row(j)  = 0.0d0
        dsa_row(j)  = 0.0d0
        dsz_row(j)  = 0.0d0
        dsdd_row(j) = 0.0d0
        dsdt_row(j) = 0.0d0
        dsda_row(j) = 0.0d0
        dsdz_row(j) = 0.0d0
        dstt_row(j) = 0.0d0
        dsta_row(j) = 0.0d0
        dstz_row(j) = 0.0d0
        dsaa_row(j) = 0.0d0
        dsaz_row(j) = 0.0d0
        dszz_row(j) = 0.0d0

        prad_row(j)    = 0.0d0
        dpradt_row(j)  = 0.0d0
        dpradd_row(j)  = 0.0d0
        dprada_row(j)  = 0.0d0
        dpradz_row(j)  = 0.0d0
        dpraddd_row(j) = 0.0d0
        dpraddt_row(j) = 0.0d0
        dpradda_row(j) = 0.0d0
        dpraddz_row(j) = 0.0d0
        dpradtt_row(j) = 0.0d0
        dpradta_row(j) = 0.0d0
        dpradtz_row(j) = 0.0d0
        dpradaa_row(j) = 0.0d0
        dpradaz_row(j) = 0.0d0
        dpradzz_row(j) = 0.0d0

        erad_row(j)    = 0.0d0
        deradt_row(j)  = 0.0d0
        deradd_row(j)  = 0.0d0
        derada_row(j)  = 0.0d0
        deradz_row(j)  = 0.0d0
        deraddd_row(j) = 0.0d0
        deraddt_row(j) = 0.0d0
        deradda_row(j) = 0.0d0
        deraddz_row(j) = 0.0d0
        deradtt_row(j) = 0.0d0
        deradta_row(j) = 0.0d0
        deradtz_row(j) = 0.0d0
        deradaa_row(j) = 0.0d0
        deradaz_row(j) = 0.0d0
        deradzz_row(j) = 0.0d0

        srad_row(j)    = 0.0d0
        dsradt_row(j)  = 0.0d0
        dsradd_row(j)  = 0.0d0
        dsrada_row(j)  = 0.0d0
        dsradz_row(j)  = 0.0d0
        dsraddd_row(j) = 0.0d0
        dsraddt_row(j) = 0.0d0
        dsradda_row(j) = 0.0d0
        dsraddz_row(j) = 0.0d0
        dsradtt_row(j) = 0.0d0
        dsradta_row(j) = 0.0d0
        dsradtz_row(j) = 0.0d0
        dsradaa_row(j) = 0.0d0
        dsradaz_row(j) = 0.0d0
        dsradzz_row(j) = 0.0d0

        pgas_row(j)    = 0.0d0
        dpgast_row(j)  = 0.0d0
        dpgasd_row(j)  = 0.0d0
        dpgasa_row(j)  = 0.0d0
        dpgasz_row(j)  = 0.0d0
        dpgasdd_row(j) = 0.0d0
        dpgasdt_row(j) = 0.0d0
        dpgasda_row(j) = 0.0d0
        dpgasdz_row(j) = 0.0d0
        dpgastt_row(j) = 0.0d0
        dpgasta_row(j) = 0.0d0
        dpgastz_row(j) = 0.0d0
        dpgasaa_row(j) = 0.0d0
        dpgasaz_row(j) = 0.0d0
        dpgaszz_row(j) = 0.0d0

        egas_row(j)    = 0.0d0
        degast_row(j)  = 0.0d0
        degasd_row(j)  = 0.0d0
        degasa_row(j)  = 0.0d0
        degasz_row(j)  = 0.0d0
        degasdd_row(j) = 0.0d0
        degasdt_row(j) = 0.0d0
        degasda_row(j) = 0.0d0
        degasdz_row(j) = 0.0d0
        degastt_row(j) = 0.0d0
        degasta_row(j) = 0.0d0
        degastz_row(j) = 0.0d0
        degasaa_row(j) = 0.0d0
        degasaz_row(j) = 0.0d0
        degaszz_row(j) = 0.0d0

        sgas_row(j)    = 0.0d0
        dsgast_row(j)  = 0.0d0
        dsgasd_row(j)  = 0.0d0
        dsgasa_row(j)  = 0.0d0
        dsgasz_row(j)  = 0.0d0
        dsgasdd_row(j) = 0.0d0
        dsgasdt_row(j) = 0.0d0
        dsgasda_row(j) = 0.0d0
        dsgasdz_row(j) = 0.0d0
        dsgastt_row(j) = 0.0d0
        dsgasta_row(j) = 0.0d0
        dsgastz_row(j) = 0.0d0
        dsgasaa_row(j) = 0.0d0
        dsgasaz_row(j) = 0.0d0
        dsgaszz_row(j) = 0.0d0

        pion_row(j)    = 0.0d0
        dpiont_row(j)  = 0.0d0
        dpiond_row(j)  = 0.0d0
        dpiona_row(j)  = 0.0d0
        dpionz_row(j)  = 0.0d0
        dpiondd_row(j) = 0.0d0
        dpiondt_row(j) = 0.0d0
        dpionda_row(j) = 0.0d0
        dpiondz_row(j) = 0.0d0
        dpiontt_row(j) = 0.0d0
        dpionta_row(j) = 0.0d0
        dpiontz_row(j) = 0.0d0
        dpionaa_row(j) = 0.0d0
        dpionaz_row(j) = 0.0d0
        dpionzz_row(j) = 0.0d0

        eion_row(j)    = 0.0d0
        deiont_row(j)  = 0.0d0
        deiond_row(j)  = 0.0d0
        deiona_row(j)  = 0.0d0
        deionz_row(j)  = 0.0d0
        deiondd_row(j) = 0.0d0
        deiondt_row(j) = 0.0d0
        deionda_row(j) = 0.0d0
        deiondz_row(j) = 0.0d0
        deiontt_row(j) = 0.0d0
        deionta_row(j) = 0.0d0
        deiontz_row(j) = 0.0d0
        deionaa_row(j) = 0.0d0
        deionaz_row(j) = 0.0d0
        deionzz_row(j) = 0.0d0

        sion_row(j)    = 0.0d0
        dsiont_row(j)  = 0.0d0
        dsiond_row(j)  = 0.0d0
        dsiona_row(j)  = 0.0d0
        dsionz_row(j)  = 0.0d0
        dsiondd_row(j) = 0.0d0
        dsiondt_row(j) = 0.0d0
        dsionda_row(j) = 0.0d0
        dsiondz_row(j) = 0.0d0
        dsiontt_row(j) = 0.0d0
        dsionta_row(j) = 0.0d0
        dsiontz_row(j) = 0.0d0
        dsionaa_row(j) = 0.0d0
        dsionaz_row(j) = 0.0d0
        dsionzz_row(j) = 0.0d0

        xni_row(j)    = 0.0d0
        dxnit_row(j)  = 0.0d0
        dxnid_row(j)  = 0.0d0
        dxnia_row(j)  = 0.0d0
        dxniz_row(j)  = 0.0d0
        dxnidd_row(j) = 0.0d0
        dxnidt_row(j) = 0.0d0
        dxnida_row(j) = 0.0d0
        dxnidz_row(j) = 0.0d0
        dxnitt_row(j) = 0.0d0
        dxnita_row(j) = 0.0d0
        dxnitz_row(j) = 0.0d0
        dxniaa_row(j) = 0.0d0
        dxniaz_row(j) = 0.0d0
        dxnizz_row(j) = 0.0d0

        xnim_row(j)   = 0.0d0

        pele_row(j)   = 0.0d0
        ppos_row(j)   = 0.0d0
        dpept_row(j)  = 0.0d0
        dpepd_row(j)  = 0.0d0
        dpepa_row(j)  = 0.0d0
        dpepz_row(j)  = 0.0d0
        dpepdd_row(j) = 0.0d0
        dpepdt_row(j) = 0.0d0
        dpepda_row(j) = 0.0d0
        dpepdz_row(j) = 0.0d0
        dpeptt_row(j) = 0.0d0
        dpepta_row(j) = 0.0d0
        dpeptz_row(j) = 0.0d0
        dpepaa_row(j) = 0.0d0
        dpepaz_row(j) = 0.0d0
        dpepzz_row(j) = 0.0d0

        eele_row(j)   = 0.0d0
        epos_row(j)   = 0.0d0
        deept_row(j)  = 0.0d0
        deepd_row(j)  = 0.0d0
        deepa_row(j)  = 0.0d0
        deepz_row(j)  = 0.0d0
        deepdd_row(j) = 0.0d0
        deepdt_row(j) = 0.0d0
        deepda_row(j) = 0.0d0
        deepdz_row(j) = 0.0d0
        deeptt_row(j) = 0.0d0
        deepta_row(j) = 0.0d0
        deeptz_row(j) = 0.0d0
        deepaa_row(j) = 0.0d0
        deepaz_row(j) = 0.0d0
        deepzz_row(j) = 0.0d0

        sele_row(j)   = 0.0d0
        spos_row(j)   = 0.0d0
        dsept_row(j)  = 0.0d0
        dsepd_row(j)  = 0.0d0
        dsepa_row(j)  = 0.0d0
        dsepz_row(j)  = 0.0d0
        dsepdd_row(j) = 0.0d0
        dsepdt_row(j) = 0.0d0
        dsepda_row(j) = 0.0d0
        dsepdz_row(j) = 0.0d0
        dseptt_row(j) = 0.0d0
        dsepta_row(j) = 0.0d0
        dseptz_row(j) = 0.0d0
        dsepaa_row(j) = 0.0d0
        dsepaz_row(j) = 0.0d0
        dsepzz_row(j) = 0.0d0

        xnem_row(j)   = 0.0d0
        xne_row(j)    = 0.0d0
        xnp_row(j)    = 0.0d0
        dxnet_row(j)  = 0.0d0
        dxned_row(j)  = 0.0d0
        dxnea_row(j)  = 0.0d0
        dxnez_row(j)  = 0.0d0
        dxnedd_row(j) = 0.0d0
        dxnedt_row(j) = 0.0d0
        dxneda_row(j) = 0.0d0
        dxnedz_row(j) = 0.0d0
        dxnett_row(j) = 0.0d0
        dxneta_row(j) = 0.0d0
        dxnetz_row(j) = 0.0d0
        dxneaa_row(j) = 0.0d0
        dxneaz_row(j) = 0.0d0
        dxnezz_row(j) = 0.0d0


        zeff_row(j)   = 0.0d0

        etaele_row(j) = 0.0d0
        detat_row(j)  = 0.0d0
        detad_row(j)  = 0.0d0
        detaa_row(j)  = 0.0d0
        detaz_row(j)  = 0.0d0
        detadd_row(j) = 0.0d0
        detadt_row(j) = 0.0d0
        detada_row(j) = 0.0d0
        detadz_row(j) = 0.0d0
        detatt_row(j) = 0.0d0
        detata_row(j) = 0.0d0
        detatz_row(j) = 0.0d0
        detaaa_row(j) = 0.0d0
        detaaz_row(j) = 0.0d0
        detazz_row(j) = 0.0d0

        etapos_row(j) = 0.0d0

        etaion_row(j) = 0.0d0
        detait_row(j)  = 0.0d0
        detaid_row(j)  = 0.0d0
        detaia_row(j)  = 0.0d0
        detaiz_row(j)  = 0.0d0
        detaidd_row(j) = 0.0d0
        detaidt_row(j) = 0.0d0
        detaida_row(j) = 0.0d0
        detaidz_row(j) = 0.0d0
        detaitt_row(j) = 0.0d0
        detaita_row(j) = 0.0d0
        detaitz_row(j) = 0.0d0
        detaiaa_row(j) = 0.0d0
        detaiaz_row(j) = 0.0d0
        detaizz_row(j) = 0.0d0

        eip_row(j)    = 0.0d0
        pip_row(j)    = 0.0d0
        sip_row(j)    = 0.0d0

        pcou_row(j)   = 0.0d0
        dpcoud_row(j) = 0.0d0
        dpcout_row(j) = 0.0d0
        dpcoua_row(j) = 0.0d0
        dpcouz_row(j) = 0.0d0

        ecou_row(j)   = 0.0d0
        decoud_row(j) = 0.0d0
        decout_row(j) = 0.0d0
        decoua_row(j) = 0.0d0
        decouz_row(j) = 0.0d0

        scou_row(j)   = 0.0d0
        dscoud_row(j) = 0.0d0
        dscout_row(j) = 0.0d0
        dscoua_row(j) = 0.0d0
        dscouz_row(j) = 0.0d0

        plasg_row(j)  = 0.0d0

        cv_gas_row(j)    = 0.0d0
        dcv_gasdd_row(j) = 0.0d0
        dcv_gasdt_row(j) = 0.0d0
        dcv_gasda_row(j) = 0.0d0
        dcv_gasdz_row(j) = 0.0d0

        cp_gas_row(j)    = 0.0d0
        dcp_gasdd_row(j) = 0.0d0
        dcp_gasdt_row(j) = 0.0d0
        dcp_gasda_row(j) = 0.0d0
        dcp_gasdz_row(j) = 0.0d0

        gam1_gas_row(j)    = 0.0d0
        dgam1_gasdd_row(j) = 0.0d0
        dgam1_gasdt_row(j) = 0.0d0
        dgam1_gasda_row(j) = 0.0d0
        dgam1_gasdz_row(j) = 0.0d0

        gam2_gas_row(j)    = 0.0d0
        dgam2_gasdd_row(j) = 0.0d0
        dgam2_gasdt_row(j) = 0.0d0
        dgam2_gasda_row(j) = 0.0d0
        dgam2_gasdz_row(j) = 0.0d0

        gam3_gas_row(j)    = 0.0d0
        dgam3_gasdd_row(j) = 0.0d0
        dgam3_gasdt_row(j) = 0.0d0
        dgam3_gasda_row(j) = 0.0d0
        dgam3_gasdz_row(j) = 0.0d0

        nabad_gas_row(j)  = 0.0d0
        dnab_gasdd_row(j) = 0.0d0
        dnab_gasdt_row(j) = 0.0d0
        dnab_gasda_row(j) = 0.0d0
        dnab_gasdz_row(j) = 0.0d0

        cs_gas_row(j)    = 0.0d0
        dcs_gasdd_row(j) = 0.0d0
        dcs_gasdt_row(j) = 0.0d0
        dcs_gasda_row(j) = 0.0d0
        dcs_gasdz_row(j) = 0.0d0

        cv_row(j)      = 0.0d0
        dcvdd_row(j)   = 0.0d0
        dcvdt_row(j)   = 0.0d0
        dcvda_row(j)   = 0.0d0
        dcvdz_row(j)   = 0.0d0

        cp_row(j)      = 0.0d0
        dcpdd_row(j)   = 0.0d0
        dcpdt_row(j)   = 0.0d0
        dcpda_row(j)   = 0.0d0
        dcpdz_row(j)   = 0.0d0

        gam1_row(j)    = 0.0d0
        dgam1dd_row(j) = 0.0d0
        dgam1dt_row(j) = 0.0d0
        dgam1da_row(j) = 0.0d0
        dgam1dz_row(j) = 0.0d0

        gam2_row(j)    = 0.0d0
        dgam2dd_row(j) = 0.0d0
        dgam2dt_row(j) = 0.0d0
        dgam2da_row(j) = 0.0d0
        dgam2dz_row(j) = 0.0d0

        gam3_row(j)    = 0.0d0
        dgam3dd_row(j) = 0.0d0
        dgam3dt_row(j) = 0.0d0
        dgam3da_row(j) = 0.0d0
        dgam3dz_row(j) = 0.0d0

        nabad_row(j)  = 0.0d0
        dnabdd_row(j) = 0.0d0
        dnabdt_row(j) = 0.0d0
        dnabda_row(j) = 0.0d0
        dnabdz_row(j) = 0.0d0

        cs_row(j)     = 0.0d0
        dcsdd_row(j)  = 0.0d0
        dcsdt_row(j)  = 0.0d0
        dcsda_row(j)  = 0.0d0
        dcsdz_row(j)  = 0.0d0

        xn_row(j)     = 0.0d0
        xp_row(j)     = 0.0d0
        xa_row(j)     = 0.0d0
        xhv_row(j)    = 0.0d0
        xmuhat_row(j) = 0.0d0

        dse_row(j)    = 0.0d0
        dpe_row(j)    = 0.0d0
        dsp_row(j)    = 0.0d0

      enddo
      return
      end






! routine dfermi gets the fermi-dirac functions and their derivaties
! routine fdfunc1 forms the integrand of the fermi-dirac functions
! routine fdfunc2 same as fdfunc but with the change of variable z**2=x
! routine dqleg010 does 10 point gauss-legendre integration  9 fig accuracy
! routine dqleg020 does 20 point gauss-legendre integration 14 fig accuracy
! routine dqleg040 does 40 point gauss-legendre integration 18 fig accuracy
! routine dqleg080 does 80 point gauss-legendre integration 32 fig accuracy
! routine dqlag010 does 10 point gauss-laguerre integration  9 fig accuracy
! routine dqlag020 does 20 point gauss-laguerre integration 14 fig accuracy
! routine dqlag040 does 40 point gauss-laguerre integration 18 fig accuracy
! routine dqlag080 does 80 point gauss-laguerre integration 32 fig accuracy




      subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta, &
                        fdeta2,fdtheta2,fdetadtheta)
      include 'implno.dek'

! this routine computes the fermi-dirac integrals f(dk,eta,theta) of
! index dk, with degeneracy parameter eta and relativity parameter theta.

! input is dk the double precision index of the fermi-dirac function,
! eta the degeneracy parameter, and theta the relativity parameter.

!. output is fd is computed by applying three 10-point gauss-legendre
! and one 10-point gauss-laguerre rules over four appropriate subintervals.
! fdeta is the derivative of f with respect to eta,
! fdthetha is the derivative of f with theta.

! reference: j.m. aparicio, apjs 117, 632 1998


! declare
      external         fdfunc1,fdfunc2
      double precision dk,eta,theta,fd,fdeta,fdtheta, &
                       fdeta2,fdtheta2,fdetadtheta, &
                       d,sg,a1,b1,c1,a2,b2,c2,d2,e2,a3,b3,c3,d3,e3, &
                       eta1,xi,xi2,x1,x2,x3,s1,s2,s3,s12,par(3), &
                       res(4),drde(4),drdt(4),drde2(4),drdt2(4),drdet(4)


!   parameters defining the location of the breakpoints for the
!   subintervals of integration:
      data d   / 3.3609d0 /
      data sg  / 9.1186d-2 /
      data a1  / 6.7774d0 /
      data b1  / 1.1418d0 /
      data c1  / 2.9826d0 /
      data a2  / 3.7601d0 /
      data b2  / 9.3719d-2 /
      data c2  / 2.1063d-2 /
      data d2  / 3.1084d1 /
      data e2  / 1.0056d0 /
      data a3  / 7.5669d0 /
      data b3  / 1.1695d0 /
      data c3  / 7.5416d-1 /
      data d3  / 6.6558d0 /
      data e3  /-1.2819d-1 /


!   integrand parameters:
      par(1)=dk
      par(2)=eta
      par(3)=theta


!   definition of xi:
      eta1=sg*(eta-d)
      if (eta1.le.5.d1) then
        xi=log(1.d0+exp(eta1))/sg
      else
        xi=eta-d
      endif
      xi2=xi*xi

!   definition of the x_i:
      x1=(a1  +b1*xi+c1*   xi2) &
        /(1.d0+c1*xi)
      x2=(a2  +b2*xi+c2*d2*xi2) &
        /(1.d0+e2*xi+c2*   xi2)
      x3=(a3  +b3*xi+c3*d3*xi2) &
        /(1.d0+e3*xi+c3*   xi2)

!   breakpoints:
      s1=x1-x2
      s2=x1
      s3=x1+x3
      s12=sqrt(s1)

!   quadrature integrations:

! 14 significant figure accuracy
      call dqleg020(fdfunc2, 0.d0,  s12, res(1), drde(1), drdt(1), &
                    drde2(1), drdt2(1), drdet(1), par,3)
      call dqleg020(fdfunc1,   s1,   s2, res(2), drde(2), drdt(2), &
                    drde2(2), drdt2(2), drdet(2), par,3)
      call dqleg020(fdfunc1,   s2,   s3, res(3), drde(3), drdt(3), &
                    drde2(3), drdt2(3), drdet(3), par,3)
      call dqlag020(fdfunc1,   s3, 1.d0, res(4), drde(4), drdt(4), &
                    drde2(4), drdt2(4), drdet(4), par,3)

! sum the contributions
      fd          = res(1)   + res(2)   + res(3)   + res(4)
      fdeta       = drde(1)  + drde(2)  + drde(3)  + drde(4)
      fdtheta     = drdt(1)  + drdt(2)  + drdt(3)  + drdt(4)
      fdeta2      = drde2(1) + drde2(2) + drde2(3) + drde2(4)
      fdtheta2    = drdt2(1) + drdt2(2) + drdt2(3) + drdt2(4)
      fdetadtheta = drdet(1) + drdet(2) + drdet(3) + drdet(4)
      return
      end




      subroutine fdfunc1(x,par,n,fd,fdeta,fdtheta, &
                         fdeta2,fdtheta2,fdetadtheta)
      include 'implno.dek'

! forms the fermi-dirac integrand and its first and second
! derivatives with eta and theta.

! input:
! x is the integration variable
! par(1) is the double precision index
! par(2) is the degeneravy parameter
! par(3) is the relativity parameter

! output:
! fd is the integrand
! fdeta is the first derivative with eta
! fdeta2 is the second derivative with eta
! fdtheta is the first derivative with theta
! fdtheta2 is the second derivative with theta
! fdetadtheta is the mixed second derivative

! declare the pass
      integer          n
      double precision x,par(n),fd, &
                       fdeta,fdeta2,fdtheta,fdtheta2,fdetadtheta

! local variables
      double precision dk,eta,theta, &
                       factor,xst,dxst,denom,denomi,denom2,xdk


! initialize
      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xdk   = x**dk
      xst   = 1.0d0 + 0.5d0*x*theta
      dxst  = sqrt(xst)

!   avoid overflow in the exponentials at large x
      if ((x-eta) .lt. 1.0d2) then
       factor      = exp(x-eta)
       denom       = factor + 1.0d0
       denomi      = 1.0d0/denom
       fd          = xdk * dxst * denomi
       fdeta       = fd * factor * denomi
       fdeta2      = (2.0d0 * factor * denomi - 1.0d0)*fdeta
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * x * denom2
       fdtheta2    = -fdtheta * x * denom2
       fdetadtheta = fdtheta * factor * denomi
      else
       factor      = exp(eta-x)
       fd          = xdk * dxst * factor
       fdeta       = fd
       fdeta2      = fd
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * x * denom2
       fdtheta2    = -fdtheta * x * denom2
       fdetadtheta = fdtheta
      endif

      return
      end




      subroutine fdfunc2(x,par,n,fd,fdeta,fdtheta, &
                         fdeta2,fdtheta2,fdetadtheta)
      include 'implno.dek'

! forms the fermi-dirac integrand and its first and second
! derivatives with eta and theta when the variable change z**2=x
! has been made.

! input:
! x is the integration variable
! par(1) is the double precision index
! par(2) is the degeneravy parameter
! par(3) is the relativity parameter

! output:
! fd is the integrand
! fdeta is the first derivative with eta
! fdeta2 is the second derivative with eta
! fdtheta is the first derivative with theta
! fdtheta2 is the second derivative with theta
! fdetadtheta is the mixed second derivative

! declare the pass
      integer          n
      double precision x,par(n),fd, &
                       fdeta,fdeta2,fdtheta,fdtheta2,fdetadtheta

! local variables
      double precision dk,eta,theta, &
                       factor,xst,dxst,denom,denomi,denom2,xdk,xsq


      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xsq   = x * x
      xdk   = x**(2.0d0 * dk + 1.0d0)
      xst   = 1.0d0 + 0.5d0 * xsq * theta
      dxst  = sqrt(xst)

!   avoid an overflow in the denominator at large x:
      if ((xsq-eta) .lt. 1.d2) then
       factor      = exp(xsq - eta)
       denom       = factor + 1.0d0
       denomi      = 1.0d0/denom
       fd          = 2.0d0 * xdk * dxst * denomi
       fdeta       = fd * factor * denomi
       fdeta2      = (2.0d0 * factor * denomi - 1.0d0)*fdeta
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * xsq * denom2
       fdtheta2    = -fdtheta * xsq * denom2
       fdetadtheta = fdtheta * factor * denomi
      else
       factor      = exp(eta - xsq)
       fd          = 2.0d0 * xdk * dxst * factor
       fdeta       = fd
       fdeta2      = fd
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * xsq * denom2
       fdtheta2    = -fdtheta * xsq * denom2
       fdetadtheta = fdtheta
      endif

      return
      end




      subroutine dqleg020(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 20 point gauss-legendre rule for the fermi-dirac function and
! its derivatives with respect to eta and theta.
! on input f is the name of the subroutine containing the integrand,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to subroutine f,
! and n is the length of the par array. on output result is the
! approximation from applying the 20-point gauss-legendre rule,
! drdeta is the derivative with respect to eta, and drdtheta is the
! derivative with respect to theta.
!
! note: since the number of nodes is even, zero is not an abscissa.
!
! declare
      external         f
      integer          j,n
      double precision a,b,result,drdeta,drdtheta, &
                       drdeta2,drdtheta2,drdetadtheta,par(n), &
                       absc1,absc2,center,hlfrun,wg(10),xg(10), &
                       fval(2),deta(2),dtheta(2), &
                       deta2(2),dtheta2(2),detadtheta(2)


! the abscissae and weights are given for the interval (-1,1).
! xg     - abscissae of the 20-point gauss-legendre rule
!          for half of the usual run (-1,1), i.e.
!          the positive nodes of the 20-point rule
! wg     - weights of the 20-point gauss rule.
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   7.65265211334973337546404093988382110d-2 /
      data xg (  2) /   2.27785851141645078080496195368574624d-1 /
      data xg (  3) /   3.73706088715419560672548177024927237d-1 /
      data xg (  4) /   5.10867001950827098004364050955250998d-1 /
      data xg (  5) /   6.36053680726515025452836696226285936d-1 /
      data xg (  6) /   7.46331906460150792614305070355641590d-1 /
      data xg (  7) /   8.39116971822218823394529061701520685d-1 /
      data xg (  8) /   9.12234428251325905867752441203298113d-1 /
      data xg (  9) /   9.63971927277913791267666131197277221d-1 /
      data xg ( 10) /   9.93128599185094924786122388471320278d-1 /

      data wg (  1) /   1.52753387130725850698084331955097593d-1 /
      data wg (  2) /   1.49172986472603746787828737001969436d-1 /
      data wg (  3) /   1.42096109318382051329298325067164933d-1 /
      data wg (  4) /   1.31688638449176626898494499748163134d-1 /
      data wg (  5) /   1.18194531961518417312377377711382287d-1 /
      data wg (  6) /   1.01930119817240435036750135480349876d-1 /
      data wg (  7) /   8.32767415767047487247581432220462061d-2 /
      data wg (  8) /   6.26720483341090635695065351870416063d-2 /
      data wg (  9) /   4.06014298003869413310399522749321098d-2 /
      data wg ( 10) /   1.76140071391521183118619623518528163d-2 /


!           list of major variables
!           -----------------------
!
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 20-point gauss formula


      center       = 0.5d0 * (a+b)
      hlfrun       = 0.5d0 * (b-a)
      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0
      do j=1,10
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)

        call f(absc1, par, n, fval(1), deta(1), dtheta(1), &
               deta2(1), dtheta2(1), detadtheta(1))
        call f(absc2, par, n, fval(2), deta(2), dtheta(2), &
               deta2(2), dtheta2(2), detadtheta(2))

        result       = result    + (fval(1)    + fval(2))*wg(j)
        drdeta       = drdeta    + (deta(1)    + deta(2))*wg(j)
        drdtheta     = drdtheta  + (dtheta(1)  + dtheta(2))*wg(j)
        drdeta2      = drdeta2   + (deta2(1)   + deta2(2))*wg(j)
        drdtheta2    = drdtheta2 + (dtheta2(1) + dtheta2(2))*wg(j)
        drdetadtheta = drdetadtheta+(detadtheta(1)+detadtheta(2))*wg(j)
      enddo

      result       = result * hlfrun
      drdeta       = drdeta * hlfrun
      drdtheta     = drdtheta * hlfrun
      drdeta2      = drdeta2 * hlfrun
      drdtheta2    = drdtheta2 * hlfrun
      drdetadtheta = drdetadtheta * hlfrun
      return
      end







      subroutine dqlag020(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 20 point gauss-laguerre rule for the fermi-dirac function.
! on input f is the external function defining the integrand
! f(x)=g(x)*w(x), where w(x) is the gaussian weight
! w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to the function f,
! and n is the length of the par array. on output result is the
! approximation from applying the 20-point gauss-laguerre rule.
! since the number of nodes is even, zero is not an abscissa.
!
! declare
      external         f
      integer          j,n
      double precision a,b,result,drdeta,drdtheta, &
                       drdeta2,drdtheta2,drdetadtheta,par(n), &
                       absc,wg(20),xg(20),fval,deta,dtheta, &
                       deta2,dtheta2,detadtheta


! the abscissae and weights are given for the interval (0,+inf).
! xg     - abscissae of the 20-point gauss-laguerre rule
! wg     - weights of the 20-point gauss rule. since f yet
!          includes the weight function, the values in wg
!          are actually exp(xg) times the standard
!          gauss-laguerre weights
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   7.05398896919887533666890045842150958d-2 /
      data xg (  2) /   3.72126818001611443794241388761146636d-1 /
      data xg (  3) /   9.16582102483273564667716277074183187d-1 /
      data xg (  4) /   1.70730653102834388068768966741305070d0 /
      data xg (  5) /   2.74919925530943212964503046049481338d0 /
      data xg (  6) /   4.04892531385088692237495336913333219d0 /
      data xg (  7) /   5.61517497086161651410453988565189234d0 /
      data xg (  8) /   7.45901745367106330976886021837181759d0 /
      data xg (  9) /   9.59439286958109677247367273428279837d0 /
      data xg ( 10) /   1.20388025469643163096234092988655158d1 /
      data xg ( 11) /   1.48142934426307399785126797100479756d1 /
      data xg ( 12) /   1.79488955205193760173657909926125096d1 /
      data xg ( 13) /   2.14787882402850109757351703695946692d1 /
      data xg ( 14) /   2.54517027931869055035186774846415418d1 /
      data xg ( 15) /   2.99325546317006120067136561351658232d1 /
      data xg ( 16) /   3.50134342404790000062849359066881395d1 /
      data xg ( 17) /   4.08330570567285710620295677078075526d1 /
      data xg ( 18) /   4.76199940473465021399416271528511211d1 /
      data xg ( 19) /   5.58107957500638988907507734444972356d1 /
      data xg ( 20) /   6.65244165256157538186403187914606659d1 /

      data wg (  1) /   1.81080062418989255451675405913110644d-1 /
      data wg (  2) /   4.22556767878563974520344172566458197d-1 /
      data wg (  3) /   6.66909546701848150373482114992515927d-1 /
      data wg (  4) /   9.15352372783073672670604684771868067d-1 /
      data wg (  5) /   1.16953970719554597380147822239577476d0 /
      data wg (  6) /   1.43135498592820598636844994891514331d0 /
      data wg (  7) /   1.70298113798502272402533261633206720d0 /
      data wg (  8) /   1.98701589079274721410921839275129020d0 /
      data wg (  9) /   2.28663578125343078546222854681495651d0 /
      data wg ( 10) /   2.60583472755383333269498950954033323d0 /
      data wg ( 11) /   2.94978373421395086600235416827285951d0 /
      data wg ( 12) /   3.32539578200931955236951937421751118d0 /
      data wg ( 13) /   3.74225547058981092111707293265377811d0 /
      data wg ( 14) /   4.21423671025188041986808063782478746d0 /
      data wg ( 15) /   4.76251846149020929695292197839096371d0 /
      data wg ( 16) /   5.42172604424557430380308297989981779d0 /
      data wg ( 17) /   6.25401235693242129289518490300707542d0 /
      data wg ( 18) /   7.38731438905443455194030019196464791d0 /
      data wg ( 19) /   9.15132873098747960794348242552950528d0 /
      data wg ( 20) /   1.28933886459399966710262871287485278d1 /


!           list of major variables
!           -----------------------
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 20-point gauss formula

      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0

      do j=1,20
       absc = a+b*xg(j)

       call f(absc, par, n, fval, deta, dtheta, &
              deta2,dtheta2,detadtheta)

       result       = result    + fval*wg(j)
       drdeta       = drdeta    + deta*wg(j)
       drdtheta     = drdtheta  + dtheta*wg(j)
       drdeta2      = drdeta2   + deta2*wg(j)
       drdtheta2    = drdtheta2 + dtheta2*wg(j)
       drdetadtheta = drdetadtheta + detadtheta*wg(j)
      enddo

      result       = result*b
      drdeta       = drdeta*b
      drdtheta     = drdtheta*b
      drdeta2      = drdeta2*b
      drdtheta2    = drdtheta2*b
      drdetadtheta = drdetadtheta*b
      return
      end



      subroutine pretty_eos_out(whose)
      include 'implno.dek'
      include 'vector_eos.dek'

! writes a pretty output for the eos tester


! declare the pass
      character*(*) whose


! local variables
      integer          i,j
      double precision ye,xcess,avo,kerg,xka
      parameter        (avo     = 6.0221417930d23, &
                        kerg    = 1.380650424d-16, &
                        xka = kerg*avo)


! popular formats
01    format(1x,t2,a,t11,a,t27,a,t43,a,t59,a,t75,a,t91,a,t107,a)
02    format(1x,t2,a,1p7e16.8)
03    format(1x,t2,a7,1pe12.4,t22,a7,1pe12.4, &
               t42,a7,1pe12.4,t62,a7,1pe12.4)
04    format(1x,t2,a,t11,'total',t24,'ion',t34,'e- + e+', &
             t58,'radiation',t70,'coulomb')
05    format(1x,t2,a,1p3e12.4,t56,1p2e12.4)
06    format(1x,t2,a,a,1pe12.4, &
                t30,a,a,1pe12.4, &
                t58,a,a,1pe12.4)



! loop over the pipeline
      do j=jlo_eos,jhi_eos


! the input
      write(6,03) 'temp  =',temp_row(j),'den   =',den_row(j), &
                  'abar  =',abar_row(j),'zbar  =',zbar_row(j)

      ye = zbar_row(1)/abar_row(1)
      xcess = 1.0d0 - 2.0d0*ye
      write(6,03) 'ye    =',ye,'xcess =',xcess
      write(6,*) ' '


! and the output

       write(6,01)  whose,'value','d/dd','d/dt','d/da','d/dz'

       write(6,02) 'p tot=',ptot_row(j), &
                    dpd_row(j),dpt_row(j),dpa_row(j),dpz_row(j)
       write(6,02) 'p gas=',pgas_row(j), &
                 dpgasd_row(j),dpgast_row(j),dpgasa_row(j),dpgasz_row(j)
       write(6,02) 'p rad=',prad_row(j), &
                 dpradd_row(j),dpradt_row(j),dprada_row(j),dpradz_row(j)
       write(6,02) 'p ion=',pion_row(j), &
                dpiond_row(j),dpiont_row(j),dpiona_row(j),dpionz_row(j)
       write(6,02) 'p  e-=',pele_row(j), &
                dpepd_row(j),dpept_row(j),dpepa_row(j),dpepz_row(j)
       write(6,02) 'p  e+=',ppos_row(j)
       write(6,02) 'p cou=',pcou_row(j), &
                dpcoud_row(j),dpcout_row(j),dpcoua_row(j),dpcouz_row(j)


       write(6,*)  ' '
       write(6,02) 'e tot=',etot_row(j), &
                    ded_row(j),det_row(j),dea_row(j),dez_row(j)
       write(6,02) 'e gas=',egas_row(j), &
                 degasd_row(j),degast_row(j),degasa_row(j),degasz_row(j)
       write(6,02) 'e rad=',erad_row(j), &
                deradd_row(j),deradt_row(j),derada_row(j),deradz_row(j)
       write(6,02) 'e ion=',eion_row(j), &
                deiond_row(j),deiont_row(j),deiona_row(j),deionz_row(j)
       write(6,02) 'e  e-=',eele_row(j), &
                deepd_row(j),deept_row(j),deepa_row(j),deepz_row(j)
       write(6,02) 'e  e+=',epos_row(j)
       write(6,02) 'e cou=',ecou_row(j), &
                decoud_row(j),decout_row(j),decoua_row(j),decouz_row(j)

       write(6,*)  ' '
       write(6,02) 's tot=',stot_row(j), &
                    dsd_row(j),dst_row(j),dsa_row(j),dsz_row(j)
       write(6,02) 's/xka=',stot_row(j)/xka, &
             dsd_row(j)/xka,dst_row(j)/xka,dsa_row(j)/xka,dsz_row(j)/xka
       write(6,02) 's gas=',sgas_row(j), &
                 dsgasd_row(j),dsgast_row(j),dsgasa_row(j),dsgasz_row(j)
       write(6,02) 's rad=',srad_row(j), &
                dsradd_row(j),dsradt_row(j),dsrada_row(j),dsradz_row(j)
       write(6,02) 's ion=',sion_row(j), &
                dsiond_row(j),dsiont_row(j),dsiona_row(j),dsionz_row(j)
       write(6,02) 's  e-=',sele_row(j), &
                dsepd_row(j),dsept_row(j),dsepa_row(j),dsepz_row(j)
       write(6,02) 's  e+=',spos_row(j)
       write(6,02) 's cou=',scou_row(j), &
                dscoud_row(j),dscout_row(j),dscoua_row(j),dscouz_row(j)


! specific heats, and ratio of electostatic to thermal energy
! the 3 gammas and the sound speed for both the gas and the total
       write(6,*)  ' '
       write(6,02) 'cv  =',cv_row(j)/(kerg*avo)*abar_row(1), &
                    dcvdd_row(j),dcvdt_row(j), &
                    dcvda_row(j),dcvdz_row(j)
       write(6,02) 'cp  =',cp_row(j), &
                    dcpdd_row(j),dcpdt_row(j), &
                    dcpda_row(j),dcpdz_row(j)
       write(6,02) 'gam1=',gam1_row(j), &
                    dgam1dd_row(j),dgam1dt_row(j), &
                    dgam1da_row(j),dgam1dz_row(j)
       write(6,02) 'gam2=',gam2_row(j), &
                    dgam2dd_row(j),dgam2dt_row(j), &
                    dgam2da_row(j),dgam2dz_row(j)
       write(6,02) 'gam3=',gam3_row(j), &
                    dgam3dd_row(j),dgam3dt_row(j), &
                    dgam3da_row(j),dgam3dz_row(j)
       write(6,02) 'cs  =',cs_row(j), &
                    dcsdd_row(j),dcsdt_row(j), &
                    dcsda_row(j),dcsdz_row(j)

       write(6,*)  ' '
       write(6,02) 'cvgas=',cv_gas_row(j)/(kerg*avo)*abar_row(1), &
                    dcv_gasdd_row(j),dcv_gasdt_row(j), &
                    dcv_gasda_row(j),dcv_gasdz_row(j)
       write(6,02) 'cpgas=',cp_gas_row(j), &
                    dcp_gasdd_row(j),dcp_gasdt_row(j), &
                    dcp_gasda_row(j),dcp_gasdz_row(j)
       write(6,02) 'g1gas=',gam1_gas_row(j), &
                    dgam1_gasdd_row(j),dgam1_gasdt_row(j), &
                    dgam1_gasda_row(j),dgam1_gasdz_row(j)
       write(6,02) 'g2gas=',gam2_gas_row(j), &
                    dgam2_gasdd_row(j),dgam2_gasdt_row(j), &
                    dgam2_gasda_row(j),dgam2_gasdz_row(j)
       write(6,02) 'g3gas=',gam3_gas_row(j), &
                    dgam3_gasdd_row(j),dgam3_gasdt_row(j), &
                    dgam3_gasda_row(j),dgam3_gasdz_row(j)
       write(6,02) 'csgas=',cs_gas_row(j), &
                    dcs_gasdd_row(j),dcs_gasdt_row(j), &
                    dcs_gasda_row(j),dcs_gasdz_row(j)


! the thermodynamic consistency relations, these should all be
! at the floating point limit of zero
       write(6,*) ' '
       write(6,03) 'maxw1 =',dse_row(j),'maxw2 =',dpe_row(j), &
                   'maxw3 =',dsp_row(j)

! number density of ions and its derivatives
       write(6,03) 'xni   =',xni_row(j),  'xnim  =',xnim_row(j)
       write(6,03) 'dxnidd=',dxned_row(j),'dxnidt=',dxnet_row(j), &
                   'dxnida=',dxnea_row(j),'dxnidz=',dxnez_row(j)

! ion chemical potential and its derivatives
       write(6,03) 'etaion=',etaion_row(j)
       write(6,03) 'detaid=',detaid_row(j),'detait=',detait_row(j), &
                   'detaia=',detaia_row(j),'detaiz=',detaiz_row(j)


! number density of electrons+positrons and its derivatives
       write(6,03) 'xnele =',xne_row(j),'xnpos =',xnp_row(j), &
                   'xnem  =',xnem_row(j)
       write(6,03) 'dxnedd=',dxned_row(j),'dxnedt=',dxnet_row(j), &
                   'dxneda=',dxnea_row(j),'dxnedz=',dxnez_row(j)


! electron chemical potential, positron chemical potential and its derivatives
       write(6,03) 'etaele=',etaele_row(j),'etapos=',etapos_row(j)
       write(6,03) 'detadd=',detad_row(j),'detadt=',detat_row(j), &
                   'detada=',detaa_row(j),'detadz=',detaz_row(j)

       write(6,03) 'zeff  =',zeff_row(j), &
                   'ionzd =',zeff_row(j)/zbar_row(j), &
                   'plasg =',plasg_row(j)

! end of pipeline loop
      enddo

      return
      end




!
! Lorenzo Lodi, last modified 24 March 2016
!
! This is a companion test suite to the program 
! `DUO: a general program for calculating spectra of diatomic molecules'
! by Sergei N. Yurchenko, Lorenzo Lodi, Jonathan Tennyson and Andrey V. Stolyarov
! Computer Physics Communication (2016)
!
! The present test suite produces input files, feeds them to DUO and
! tries to determine if they ran successfully.
! For each test input file there is a certain number of subtests; typically we
! check that a few energies or other key values have their expected (correct) value.
! USAGE:
! Just run it in the same directory where the duo executable is located:
!
! ./duo_test.exe
!
! The DUO executable should be named 'duo.exe'. This can be changed by modifying
! the `executable' character variable at the very beginning of the program. 
! The test suite accepts numbers as command line parameters, indicating which 
! test files should be skipped. For example
! 
! ./duo_test.exe 1 2 4
!
! runs all tests EXCEPT tests number 1,2 and 4. Any other command line parameter
! apart from a number in the appropriate range is ignored. 
!
! INSTRUCTIONS TO ADD NEW SUBTESTS FOR AN EXISTING TEST FILE
! If you want to add further subtests for an existing test file,
! here's what you have to do.
! 1) go to the section immediately after the relative test input and update the line
!    `n_subtests(itest)='
! 2) Immediately below, add new `SubtestDescription' for the new subtests
! 3) In the do loop immediately below insert the new test, e.g. the new
!    output line you want to check
!
! INSTRUCTIONS TO ADD A NEW TEST FILE AND RELATIVE NEW SUBTESTS
! 1) Modify the parameter `n_test_files' at the beginning of the program
! 2) Use another file input section as a model: Copy, paste and modify.
! In particular:
! -update TestFileName(itest)
! -change nlines (number of lines expected in the output).
! -change n_subtests(itest) to the number of subtests for the new test file
! -change SubtestDescription(.) to some meaningful description of the subtest
! -modify the output lines to be checked. Leave out all spaces both on the left and
!  on the right.
!
! Note that the way results are checked is often rather primitive (comparison of whole lines)
! so that changes in the formatting of the output may cause a subtest to fail.
!
!*****************************************************************************
! I need a subroutine to run externally the DUO executable.
! There is no standard way of doing that in Fortran 2003 or earlier;
! Fortran 2008 added the intrinsic subroutine execute_command_line but not all
! compilers support it at this time. According to
! www.fortranplus.co.uk/.../fortran_2003_2008_compiler_support.pdf
! this feature is supported by: Absoft v14, Cray v8.4.0, gfortran 4.9,
! IBM 15.1, Intel 16.0, NAG 6.0. 
! It is NOT supported by: g95, HP, Pathscale, PGI, Oracle.
! When execute_command_line is not available as an intrisic, usually an
! equivalent routine is available under another name. For example,
! the Intel compiler before v15 has a routine called SYSTEM.
! The following 4-line subroutine is a wrapper to the version available
! in a specific compiler.
! Note that if the intrisic version is available, it will be used automatically.
subroutine execute_command_line( cmd )
  character(len=*) :: cmd
  call SYSTEM( cmd )
end subroutine execute_command_line
!*****************************************************************************

program test
implicit none
character(len=20), parameter :: executable='duo.exe'
character(len=20) :: FileInput, FileOutput
integer, parameter :: u1=100
character(len=500) :: cmd, line
integer :: itest
integer :: nlines
integer, parameter :: n_test_files=6
logical :: z_run_test_file(n_test_files) = .true.
integer :: n_subtests(n_test_files)=0
integer, parameter :: max_n_subtests=250 ! max number of sub-tests for each test input file.
character(len=200) :: TestFileName(n_test_files)
character(len=200) :: SubtestDescription(max_n_subtests)
integer :: iflag(max_n_subtests) 
integer :: i, ierr, itmp, n_test_files_to_run
logical :: zExists
integer, parameter :: dp=kind(1.d0)
real(kind=dp) :: value, thresh, ref
integer :: tot_n_subtests, n_failed_subtests=0
character(len=32) :: arg

! Reads command line parameters; the 0-th is the ./duo_test.exe command.
! This allows to run only specific tests.
i = 1
do
  call get_command_argument(i, arg)
  if (len_trim(arg) == 0) exit
  read(arg,*, iostat=ierr) itmp
  if(ierr /=0) cycle ! ignore argument if not an integer
  if(itmp <=0) cycle ! ignore if argument less than or equal to zero
  if(itmp > n_test_files) cycle ! ignore if argument greater than n_test_files
  z_run_test_file(itmp) = .false. ! disable specific test
!  write (*,'(A)') 'Command line argument: ' // trim(arg)
  i = i+1
enddo

n_test_files_to_run=0
do i=1, n_test_files
  if(z_run_test_file(i) ) n_test_files_to_run=n_test_files_to_run+1
enddo

write(*,'(A)') 'This program will run a battery of tests for the program DUO'
write(*,'(A)') 'Running all the tests should take less than about 30 seconds.'
write(*,'(A,I0,a,/)') 'A total number of ', n_test_files_to_run, ' test files will created and run.'


inquire(file=trim(executable), EXIST=zExists)
if( .not. zExists) then
  write(*,'(A)') 'Error: the DUO executable should be in the current directory!'
  write(*,'(A)') 'It should be called: ' // trim(executable)
  write(*,'(A)') 'Exiting now...'
  stop
endif

!***************************************************************
!**************** TEST 001 *************************************
!***************************************************************
1 continue
itest=1
if(z_run_test_file(itest) .eqv. .false. ) goto 2 ! skip test if flag is set
TestFileName(itest) = 'Simple harmonic oscillator'
write(FileInput,  '(a,i4.4,a)') 'duo_test_', itest, '.inp'
write(FileOutput, '(a,i4.4,a)') 'duo_test_', itest, '.out'
! FileOutput2='duo_test_001.log'

open(unit=u1, file=FileInput, status='unknown', action='write')

write(u1,'(a)') '( This is a test system with a harmonic potential                                   )'
write(u1,'(a)') '( The potential is V(r) = 1 cm^{-1}/ang^2 (r/ang-100)^2                             )'
write(u1,'(a)') '( The minimum is chosen to be at 100 angstroms so that we are very far from r=0     )'
write(u1,'(a)') '( The mass is set to give eigenvalues 0.5,1.5,2.5,...                               )'
write(u1,'(a)') '( Given (cgs system) h=6.62606957e-27, c=2.99792458e10, uma =1.660538921e-24 grams  )'
write(u1,'(a)') '(      then mass = h *1e16 / (4*pi^2*c) ; the 1e16 is conversion from cm^2 to angstroms^2 )'
write(u1,'(a)') ''
write(u1,'(a)') 'masses 33.71525838831536 33.71525838831536'
write(u1,'(a)') 'molecule fictitious'
write(u1,'(a)') ''
write(u1,'(a)') 'nstates 1'
write(u1,'(a)') 'jrot  0.0'
write(u1,'(a)') 'SYMMETRY Cs(M)'
write(u1,'(a)') ''
write(u1,'(a)') '( SolutionMethod 5PointDifferences)'
write(u1,'(a)') 'SolutionMethod Sinc'
write(u1,'(a)') ''
write(u1,'(a)') '( With SolutionMethod Sinc 110 points are enough to converge all levels to 1e-6 cm-1)'
write(u1,'(a)') '( With SolutionMethod 5PointDifferences 6900 points are enough to converge all levels to 1e-6 cm-1)'
write(u1,'(a)') ''
write(u1,'(a)') ''
write(u1,'(a)') 'grid'
write(u1,'(a)') '  npoints  120'
write(u1,'(a)') '  range  80.00, 120.00'
write(u1,'(a)') '  type 0'
write(u1,'(a)') 'end'
write(u1,'(a)') ''
write(u1,'(a)') 'DIAGONALIZER'
write(u1,'(a)') '  SYEV'
write(u1,'(a)') '  enermax  15000.0'
write(u1,'(a)') '  nroots 100'
write(u1,'(a)') 'end'
write(u1,'(a)') ''
write(u1,'(a)') 'CONTRACTION'
write(u1,'(a)') '  vib'
write(u1,'(a)') '  vmax  50'
write(u1,'(a)') 'END'
write(u1,'(a)') ''
write(u1,'(a)') 'poten 1'
write(u1,'(a)') 'name "Harmonic"'
write(u1,'(a)') 'type   polynomial'
write(u1,'(a)') 'lambda 0'
write(u1,'(a)') 'mult   1'
write(u1,'(a)') 'symmetry + '
write(u1,'(a)') 'units  cm-1'
write(u1,'(a)') 'values'
write(u1,'(a)') 'c0 0.0     0.0000000'
write(u1,'(a)') 're 0.0   100.0000000'
write(u1,'(a)') 'c1 0.0     0.0000000'
write(u1,'(a)') 'c2 0.0     0.2500000'
write(u1,'(a)') 'end'
write(u1,'(a)') ''
close(u1)

 cmd = "./" // trim(executable) // " < "  // trim(FileInput) // " > " // trim(FileOutput) !// " 2> " // trim(FileOutput2)
call execute_command_line(cmd)

open(unit=u1, file=FileOutput, status='old', action='read')

! The output should have 508 lines but we try to read more than that, to allow for
! future further output
nlines = 550
iflag = 0 !0 means the test has failed

n_subtests(itest)=5 
SubtestDescription(1) = 'Check mass'
SubtestDescription(2) = 'Check pot. r=80'
SubtestDescription(3) = 'Check v=25 energy'
SubtestDescription(4) = 'Check v=49 energy'
SubtestDescription(5) = 'Check ZPE'


do i=1, nlines
   read(u1, '(a500)', iostat=ierr) line
   if(ierr /=0) exit
   line = ADJUSTL(line) ! remove leading spaces if present
   if( trim(line) == 'Reduced mass is      16.8576291941577         atomic mass units (Daltons)')          iflag(1) = 1
   if( trim(line) == '80.00000000        1.00000000E+02')                                                  iflag(2) = 1
   if( trim(line) == '0.0   25         24.000000   1    24   0     0.0     0.0     0.0   +    ||Harmonic') iflag(3) = 1
   if( trim(line) == '0.0   50         49.000000   1    49   0     0.0     0.0     0.0   +    ||Harmonic') iflag(4) = 1
   if( trim(line) == 'Zero point energy (ZPE) =           0.500000')                                       iflag(5) = 1
enddo
close(u1)

do i=1, n_subtests(itest)
  write(line,'(A, I6,A)') 'Input file ' // trim(FileInput) //  " ( "  // trim(TestFileName(itest)) // " )" // &
                                          ' , subtest ', i, " (" // trim(SubtestDescription(i)) // ") : "
  write(*,'(a120)',advance='no') line
  if(iflag(i) == 1) then
     write(*,'(a)') 'Passed'
  else
      write(*,'(a)') 'FAILED! <---'
      n_failed_subtests=n_failed_subtests+1
  endif

enddo


!***************************************************************
!**************** END OF TEST 001 ******************************
!***************************************************************


!***************************************************************
!**************** TEST 002 *************************************
!***************************************************************
2 continue
itest=2
if(z_run_test_file(itest) .eqv. .false. ) goto 3 ! skip test if flag is set
TestFileName(itest) = 'Double well'
write(FileInput,  '(a,i4.4,a)') 'duo_test_', itest, '.inp'
write(FileOutput, '(a,i4.4,a)') 'duo_test_', itest, '.out'

open(unit=u1, file=FileInput, status='unknown', action='write')

write(u1,'(a)') 'masses 40000. 40000.'
write(u1,'(a)') ''
write(u1,'(a)') 'nstates 1'
write(u1,'(a)') 'jrot  0.0'
write(u1,'(a)') ''
write(u1,'(a)') 'SolutionMethod Sinc'
write(u1,'(a)') ''
write(u1,'(a)') 'grid'
write(u1,'(a)') '  npoints  350'
write(u1,'(a)') '  range  6.00, 14.00'
write(u1,'(a)') '  type 0'
write(u1,'(a)') 'end'
write(u1,'(a)') ''
write(u1,'(a)') 'CONTRACTION'
write(u1,'(a)') '  vib'
write(u1,'(a)') '  vmax  200'
write(u1,'(a)') 'END'
write(u1,'(a)') ''
write(u1,'(a)') 'DIAGONALIZER'
write(u1,'(a)') '  SYEV'
write(u1,'(a)') '  enermax  15000.0'
write(u1,'(a)') '  nroots 200'
write(u1,'(a)') 'end'
write(u1,'(a)') ''
write(u1,'(a)') ''
write(u1,'(a)') 'poten 1'
write(u1,'(a)') 'name "DWell"'
write(u1,'(a)') 'type   grid'
write(u1,'(a)') 'lambda 0'
write(u1,'(a)') 'mult   1'
write(u1,'(a)') 'symmetry +'
write(u1,'(a)') 'units  cm-1'
write(u1,'(a)') 'values'
write(u1,'(a)') '6.      12.380822009198916'
write(u1,'(a)') '6.02    12.25108063829405'
write(u1,'(a)') '6.04    12.12198960007542'
write(u1,'(a)') '6.06    11.993548894543036'
write(u1,'(a)') '6.08    11.865758521696883'
write(u1,'(a)') '6.1     11.738618481536978'
write(u1,'(a)') '6.12    11.612128774063308'
write(u1,'(a)') '6.14    11.486289399275883'
write(u1,'(a)') '6.16    11.361100357174692'
write(u1,'(a)') '6.18    11.236561647759745'
write(u1,'(a)') '6.2     11.112673271031035'
write(u1,'(a)') '6.22    10.989435226988569'
write(u1,'(a)') '6.24    10.866847515632339'
write(u1,'(a)') '6.26    10.744910136962353'
write(u1,'(a)') '6.28    10.623623090978601'
write(u1,'(a)') '6.3     10.502986377681095'
write(u1,'(a)') '6.32    10.382999997069824'
write(u1,'(a)') '6.34    10.263663949144798'
write(u1,'(a)') '6.36    10.144978233906006'
write(u1,'(a)') '6.38    10.026942851353459'
write(u1,'(a)') '6.4     9.909557801487148'
write(u1,'(a)') '6.42    9.792823084307082'
write(u1,'(a)') '6.44    9.676738699813251'
write(u1,'(a)') '6.46    9.561304648005665'
write(u1,'(a)') '6.48    9.446520928884313'
write(u1,'(a)') '6.5     9.332387542449206'
write(u1,'(a)') '6.52    9.218904488700339'
write(u1,'(a)') '6.54    9.106071767637708'
write(u1,'(a)') '6.5600000000000005      8.993889379261315'
write(u1,'(a)') '6.58    8.882357323571167'
write(u1,'(a)') '6.6     8.77147560056726'
write(u1,'(a)') '6.62    8.661244210249588'
write(u1,'(a)') '6.64    8.551663152618161'
write(u1,'(a)') '6.66    8.44273242767297'
write(u1,'(a)') '6.68    8.334452035414023'
write(u1,'(a)') '6.7     8.226821975841311'
write(u1,'(a)') '6.72    8.119842248954843'
write(u1,'(a)') '6.74    8.01351285475461'
write(u1,'(a)') '6.76    7.9078337932406235'
write(u1,'(a)') '6.78    7.80280506441287'
write(u1,'(a)') '6.8     7.698426668271363'
write(u1,'(a)') '6.82    7.59469860481609'
write(u1,'(a)') '6.84    7.4916208740470625'
write(u1,'(a)') '6.86    7.38919347596427'
write(u1,'(a)') '6.88    7.2874164105677215'
write(u1,'(a)') '6.9     7.186289677857409'
write(u1,'(a)') '6.92    7.08581327783334'
write(u1,'(a)') '6.94    6.985987210495508'
write(u1,'(a)') '6.96    6.886811475843919'
write(u1,'(a)') '6.98    6.788286073878567'
write(u1,'(a)') '7.      6.690411004599458'
write(u1,'(a)') '7.02    6.59318626800659'
write(u1,'(a)') '7.04    6.4966118640999575'
write(u1,'(a)') '7.0600000000000005      6.400687792879565'
write(u1,'(a)') '7.08    6.305414054345416'
write(u1,'(a)') '7.1     6.210790648497509'
write(u1,'(a)') '7.12    6.116817575335838'
write(u1,'(a)') '7.140000000000001       6.023494834860409'
write(u1,'(a)') '7.16    5.9308224270712255'
write(u1,'(a)') '7.18    5.838800351968285'
write(u1,'(a)') '7.2     5.747428609551585'
write(u1,'(a)') '7.22    5.656707199821139'
write(u1,'(a)') '7.24    5.566636122776939'
write(u1,'(a)') '7.26    5.477215378419001'
write(u1,'(a)') '7.28    5.388444966747329'
write(u1,'(a)') '7.3     5.300324887761944'
write(u1,'(a)') '7.32    5.212855141462859'
write(u1,'(a)') '7.34    5.126035727850119'
write(u1,'(a)') '7.36    5.0398666469237625'
write(u1,'(a)') '7.38    4.954347898683878'
write(u1,'(a)') '7.4     4.869479483130566'
write(u1,'(a)') '7.42    4.7852614002640035'
write(u1,'(a)') '7.4399999999999995      4.701693650084434'
write(u1,'(a)') '7.46    4.618776232592216'
write(u1,'(a)') '7.48    4.536509147787893'
write(u1,'(a)') '7.5     4.454892395672251'
write(u1,'(a)') '7.52    4.373925976246431'
write(u1,'(a)') '7.54    4.293609889512097'
write(u1,'(a)') '7.5600000000000005      4.2139441354716665'
write(u1,'(a)') '7.58    4.134928714128632'
write(u1,'(a)') '7.6     4.056563625488009'
write(u1,'(a)') '7.62    3.9788488695570017'
write(u1,'(a)') '7.640000000000001       3.9017844463459137'
write(u1,'(a)') '7.66    3.8253703558694268'
write(u1,'(a)') '7.68    3.7496065981483797'
write(u1,'(a)') '7.7     3.6744931732122645'
write(u1,'(a)') '7.72    3.6000300811026826'
write(u1,'(a)') '7.74    3.5262173218780806'
write(u1,'(a)') '7.76    3.4530548956203324'
write(u1,'(a)') '7.78    3.3805428024436863'
write(u1,'(a)') '7.8     3.308681042507057'
write(u1,'(a)') '7.82    3.2374696160306797'
write(u1,'(a)') '7.84    3.16690852331875'
write(u1,'(a)') '7.86    3.096997764789926'
write(u1,'(a)') '7.88    3.027737341018425'
write(u1,'(a)') '7.9     2.9591272527890324'
write(u1,'(a)') '7.92    2.8911675011705906'
write(u1,'(a)') '7.9399999999999995      2.823858087613669'
write(u1,'(a)') '7.96    2.757199014079982'
write(u1,'(a)') '7.98    2.6911902832131456'
write(u1,'(a)') '8.      2.6258318985631006'
write(u1,'(a)') '8.02    2.561123864879933'
write(u1,'(a)') '8.04    2.4970661884970355'
write(u1,'(a)') '8.06    2.433658877828778'
write(u1,'(a)') '8.08    2.3709019440143524'
write(u1,'(a)') '8.1     2.308795401747307'
write(u1,'(a)') '8.120000000000001       2.2473392703401207'
write(u1,'(a)') '8.14    2.186533575084777'
write(u1,'(a)') '8.16    2.126378348984446'
write(u1,'(a)') '8.18    2.0668736349485237'
write(u1,'(a)') '8.2     2.008019488563279'
write(u1,'(a)') '8.22    1.949815981574398'
write(u1,'(a)') '8.24    1.892263206245755'
write(u1,'(a)') '8.26    1.8353612807914201'
write(u1,'(a)') '8.280000000000001       1.7791103561159405'
write(u1,'(a)') '8.3     1.7235106241412999'
write(u1,'(a)') '8.32    1.6685623280484492'
write(u1,'(a)') '8.34    1.6142657748172542'
write(u1,'(a)') '8.36    1.5606213505108306'
write(u1,'(a)') '8.379999999999999       1.5076295388191203'
write(u1,'(a)') '8.4     1.4552909434515573'
write(u1,'(a)') '8.42    1.4036063150494356'
write(u1,'(a)') '8.44    1.3525765833739836'
write(u1,'(a)') '8.46    1.3022028956152505'
write(u1,'(a)') '8.48    1.2524866617571508'
write(u1,'(a)') '8.5     1.2034296080234808'
write(u1,'(a)') '8.52    1.155033839514993'
write(u1,'(a)') '8.54    1.107301913224531'
write(u1,'(a)') '8.56    1.0602369226813584'
write(u1,'(a)') '8.58    1.0138425955213168'
write(u1,'(a)') '8.6     0.9681234052997508'
write(u1,'(a)') '8.620000000000001       0.9230846988522471'
write(u1,'(a)') '8.64    0.8787328404555184'
write(u1,'(a)') '8.66    0.8350753739388327'
write(u1,'(a)') '8.68    0.7921212037359368'
write(u1,'(a)') '8.7     0.7498807956386083'
write(u1,'(a)') '8.72    0.7083663977067426'
write(u1,'(a)') '8.74    0.667592281396984'
write(u1,'(a)') '8.76    0.6275750024842248'
write(u1,'(a)') '8.780000000000001       0.5883336807614332'
write(u1,'(a)') '8.8     0.5498902968081878'
write(u1,'(a)') '8.82    0.5122700033154052'
write(u1,'(a)') '8.84    0.4755014475444629'
write(u1,'(a)') '8.86    0.4396171004884444'
write(u1,'(a)') '8.879999999999999       0.40465358720224703'
write(u1,'(a)') '8.9     0.3706520115921036'
write(u1,'(a)') '8.92    0.3376582677253648'
write(u1,'(a)') '8.94    0.30572332846584654'
write(u1,'(a)') '8.96    0.2749035009932244'
write(u1,'(a)') '8.98    0.2452606375673128'
write(u1,'(a)') '9.      0.21686228879685562'
write(u1,'(a)') '9.02    0.18978178572023555'
write(u1,'(a)') '9.04    0.16409823625920256'
write(u1,'(a)') '9.06    0.13989642112714581'
write(u1,'(a)') '9.08    0.1172665741224295'
write(u1,'(a)') '9.1     0.09630403197645122'
write(u1,'(a)') '9.120000000000001       0.07710873961419758'
write(u1,'(a)') '9.14    0.059784597874899176'
write(u1,'(a)') '9.16    0.04443864247661817'
write(u1,'(a)') '9.18    0.031180045324005357'
write(u1,'(a)') '9.2     0.020118932171036597'
write(u1,'(a)') '9.22    0.01136501416093202'
write(u1,'(a)') '9.24    0.005026034853676509'
write(u1,'(a)') '9.26    0.001206038975180806'
write(u1,'(a)') '9.280000000000001       3.4742145088537065e-6'
write(u1,'(a)') '9.3     0.0015091428648649349'
write(u1,'(a)') '9.32    0.005804025833349027'
write(u1,'(a)') '9.34    0.012957007392782134'
write(u1,'(a)') '9.36    0.023022534853237343'
write(u1,'(a)') '9.379999999999999       0.03603825290967934'
write(u1,'(a)') '9.4     0.052022657579507636'
write(u1,'(a)') '9.42    0.07097281917541098'
write(u1,'(a)') '9.44    0.0928622274583435'
write(u1,'(a)') '9.46    0.11763881478141017'
write(u1,'(a)') '9.48    0.1452232144803451'
write(u1,'(a)') '9.5     0.1755073118238141'
write(u1,'(a)') '9.52    0.20835314337059976'
write(u1,'(a)') '9.54    0.24359219749254543'
write(u1,'(a)') '9.56    0.28102516405808825'
write(u1,'(a)') '9.58    0.3204221748287472'
write(u1,'(a)') '9.6     0.3615235680528527'
write(u1,'(a)') '9.620000000000001       0.4040412011576378'
write(u1,'(a)') '9.64    0.4476603245124995'
write(u1,'(a)') '9.66    0.49204201718898816'
write(u1,'(a)') '9.68    0.5368261727563602'
write(u1,'(a)') '9.7     0.581635009751608'
write(u1,'(a)') '9.72    0.6260770679131223'
write(u1,'(a)') '9.74    0.6697516379585761'
write(u1,'(a)') '9.76    0.7122535600255839'
write(u1,'(a)') '9.780000000000001       0.7531783142840695'
write(u1,'(a)') '9.8     0.7921273170648412'
write(u1,'(a)') '9.82    0.828713327493262'
write(u1,'(a)') '9.84    0.8625658633911788'
write(u1,'(a)') '9.86    0.8933365213798322'
write(u1,'(a)') '9.879999999999999       0.9207040948781807'
write(u1,'(a)') '9.9     0.9443793851653894'
write(u1,'(a)') '9.92    0.9641096048992092'
write(u1,'(a)') '9.94    0.9796822804022559'
write(u1,'(a)') '9.96    0.990928568506014'
write(u1,'(a)') '9.98    0.9977259155532991'
write(u1,'(a)') '10.     1.'
write(u1,'(a)') '10.02   0.9977259155532991'
write(u1,'(a)') '10.04   0.990928568506014'
write(u1,'(a)') '10.06   0.9796822804022559'
write(u1,'(a)') '10.08   0.9641096048992092'
write(u1,'(a)') '10.1    0.9443793851653894'
write(u1,'(a)') '10.120000000000001      0.9207040948781807'
write(u1,'(a)') '10.14   0.8933365213798322'
write(u1,'(a)') '10.16   0.8625658633911788'
write(u1,'(a)') '10.18   0.828713327493262'
write(u1,'(a)') '10.2    0.7921273170648412'
write(u1,'(a)') '10.219999999999999      0.7531783142840695'
write(u1,'(a)') '10.24   0.7122535600255839'
write(u1,'(a)') '10.26   0.6697516379585761'
write(u1,'(a)') '10.280000000000001      0.6260770679131186'
write(u1,'(a)') '10.3    0.581635009751608'
write(u1,'(a)') '10.32   0.5368261727563602'
write(u1,'(a)') '10.34   0.49204201718898816'
write(u1,'(a)') '10.36   0.4476603245124995'
write(u1,'(a)') '10.379999999999999      0.4040412011576378'
write(u1,'(a)') '10.4    0.3615235680528527'
write(u1,'(a)') '10.42   0.3204221748287472'
write(u1,'(a)') '10.440000000000001      0.2810251640580848'
write(u1,'(a)') '10.46   0.24359219749254543'
write(u1,'(a)') '10.48   0.20835314337059976'
write(u1,'(a)') '10.5    0.1755073118238141'
write(u1,'(a)') '10.52   0.1452232144803451'
write(u1,'(a)') '10.54   0.11763881478141017'
write(u1,'(a)') '10.56   0.0928622274583435'
write(u1,'(a)') '10.58   0.07097281917541098'
write(u1,'(a)') '10.600000000000001      0.05202265757950601'
write(u1,'(a)') '10.620000000000001      0.03603825290967934'
write(u1,'(a)') '10.64   0.023022534853237343'
write(u1,'(a)') '10.66   0.012957007392782134'
write(u1,'(a)') '10.68   0.005804025833349027'
write(u1,'(a)') '10.7    0.0015091428648649349'
write(u1,'(a)') '10.719999999999999      3.4742145088537065e-6'
write(u1,'(a)') '10.74   0.001206038975180806'
write(u1,'(a)') '10.76   0.005026034853676509'
write(u1,'(a)') '10.780000000000001      0.011365014160932652'
write(u1,'(a)') '10.8    0.020118932171036597'
write(u1,'(a)') '10.82   0.031180045324005357'
write(u1,'(a)') '10.84   0.04443864247661817'
write(u1,'(a)') '10.86   0.059784597874899176'
write(u1,'(a)') '10.879999999999999      0.07710873961419758'
write(u1,'(a)') '10.9    0.09630403197645122'
write(u1,'(a)') '10.92   0.1172665741224295'
write(u1,'(a)') '10.940000000000001      0.1398964211271479'
write(u1,'(a)') '10.96   0.16409823625920256'
write(u1,'(a)') '10.98   0.18978178572023555'
write(u1,'(a)') '11.     0.21686228879685562'
write(u1,'(a)') '11.02   0.2452606375673128'
write(u1,'(a)') '11.04   0.2749035009932244'
write(u1,'(a)') '11.06   0.30572332846584654'
write(u1,'(a)') '11.08   0.3376582677253648'
write(u1,'(a)') '11.100000000000001      0.3706520115921065'
write(u1,'(a)') '11.120000000000001      0.40465358720224703'
write(u1,'(a)') '11.14   0.4396171004884444'
write(u1,'(a)') '11.16   0.4755014475444629'
write(u1,'(a)') '11.18   0.5122700033154052'
write(u1,'(a)') '11.2    0.5498902968081878'
write(u1,'(a)') '11.219999999999999      0.5883336807614332'
write(u1,'(a)') '11.24   0.6275750024842248'
write(u1,'(a)') '11.26   0.667592281396984'
write(u1,'(a)') '11.280000000000001      0.7083663977067463'
write(u1,'(a)') '11.3    0.7498807956386083'
write(u1,'(a)') '11.32   0.7921212037359368'
write(u1,'(a)') '11.34   0.8350753739388327'
write(u1,'(a)') '11.36   0.8787328404555184'
write(u1,'(a)') '11.379999999999999      0.9230846988522471'
write(u1,'(a)') '11.4    0.9681234052997508'
write(u1,'(a)') '11.42   1.0138425955213168'
write(u1,'(a)') '11.440000000000001      1.0602369226813626'
write(u1,'(a)') '11.46   1.107301913224531'
write(u1,'(a)') '11.48   1.155033839514993'
write(u1,'(a)') '11.5    1.2034296080234808'
write(u1,'(a)') '11.52   1.2524866617571508'
write(u1,'(a)') '11.54   1.3022028956152505'
write(u1,'(a)') '11.56   1.3525765833739836'
write(u1,'(a)') '11.58   1.4036063150494356'
write(u1,'(a)') '11.600000000000001      1.455290943451562'
write(u1,'(a)') '11.620000000000001      1.5076295388191203'
write(u1,'(a)') '11.64   1.5606213505108306'
write(u1,'(a)') '11.66   1.6142657748172542'
write(u1,'(a)') '11.68   1.6685623280484492'
write(u1,'(a)') '11.7    1.7235106241412999'
write(u1,'(a)') '11.719999999999999      1.7791103561159405'
write(u1,'(a)') '11.74   1.8353612807914201'
write(u1,'(a)') '11.76   1.892263206245755'
write(u1,'(a)') '11.780000000000001      1.9498159815744032'
write(u1,'(a)') '11.8    2.008019488563279'
write(u1,'(a)') '11.82   2.0668736349485237'
write(u1,'(a)') '11.84   2.126378348984446'
write(u1,'(a)') '11.86   2.186533575084777'
write(u1,'(a)') '11.879999999999999      2.2473392703401207'
write(u1,'(a)') '11.9    2.308795401747307'
write(u1,'(a)') '11.92   2.3709019440143524'
write(u1,'(a)') '11.940000000000001      2.4336588778287838'
write(u1,'(a)') '11.96   2.4970661884970355'
write(u1,'(a)') '11.98   2.561123864879933'
write(u1,'(a)') '12.     2.6258318985631006'
write(u1,'(a)') '12.02   2.6911902832131456'
write(u1,'(a)') '12.04   2.7571990140799794'
write(u1,'(a)') '12.06   2.823858087613669'
write(u1,'(a)') '12.08   2.8911675011705906'
write(u1,'(a)') '12.100000000000001      2.9591272527890387'
write(u1,'(a)') '12.120000000000001      3.027737341018428'
write(u1,'(a)') '12.14   3.0969977647899296'
write(u1,'(a)') '12.16   3.16690852331875'
write(u1,'(a)') '12.18   3.2374696160306797'
write(u1,'(a)') '12.2    3.3086810425070543'
write(u1,'(a)') '12.219999999999999      3.3805428024436828'
write(u1,'(a)') '12.24   3.4530548956203324'
write(u1,'(a)') '12.26   3.5262173218780806'
write(u1,'(a)') '12.280000000000001      3.6000300811026853'
write(u1,'(a)') '12.3    3.674493173212268'
write(u1,'(a)') '12.32   3.7496065981483797'
write(u1,'(a)') '12.34   3.8253703558694268'
write(u1,'(a)') '12.36   3.9017844463459137'
write(u1,'(a)') '12.379999999999999      3.978848869556998'
write(u1,'(a)') '12.4    4.056563625488009'
write(u1,'(a)') '12.42   4.134928714128632'
write(u1,'(a)') '12.440000000000001      4.213944135471674'
write(u1,'(a)') '12.46   4.2936098895121'
write(u1,'(a)') '12.48   4.373925976246431'
write(u1,'(a)') '12.5    4.454892395672251'
write(u1,'(a)') '12.52   4.536509147787893'
write(u1,'(a)') '12.54   4.618776232592213'
write(u1,'(a)') '12.56   4.701693650084434'
write(u1,'(a)') '12.58   4.7852614002640035'
write(u1,'(a)') '12.600000000000001      4.869479483130573'
write(u1,'(a)') '12.620000000000001      4.954347898683881'
write(u1,'(a)') '12.64   5.039866646923767'
write(u1,'(a)') '12.66   5.126035727850119'
write(u1,'(a)') '12.68   5.212855141462859'
write(u1,'(a)') '12.7    5.30032488776194'
write(u1,'(a)') '12.719999999999999      5.388444966747325'
write(u1,'(a)') '12.74   5.477215378419001'
write(u1,'(a)') '12.76   5.566636122776939'
write(u1,'(a)') '12.780000000000001      5.656707199821143'
write(u1,'(a)') '12.8    5.74742860955159'
write(u1,'(a)') '12.82   5.838800351968285'
write(u1,'(a)') '12.84   5.9308224270712255'
write(u1,'(a)') '12.86   6.023494834860409'
write(u1,'(a)') '12.879999999999999      6.116817575335833'
write(u1,'(a)') '12.9    6.210790648497509'
write(u1,'(a)') '12.92   6.305414054345416'
write(u1,'(a)') '12.940000000000001      6.400687792879573'
write(u1,'(a)') '12.96   6.496611864099962'
write(u1,'(a)') '12.98   6.59318626800659'
write(u1,'(a)') '13.     6.690411004599458'
write(u1,'(a)') '13.02   6.788286073878567'
write(u1,'(a)') '13.04   6.886811475843914'
write(u1,'(a)') '13.06   6.985987210495512'
write(u1,'(a)') '13.08   7.08581327783334'
write(u1,'(a)') '13.100000000000001      7.186289677857417'
write(u1,'(a)') '13.120000000000001      7.287416410567726'
write(u1,'(a)') '13.14   7.389193475964274'
write(u1,'(a)') '13.16   7.4916208740470625'
write(u1,'(a)') '13.18   7.59469860481609'
write(u1,'(a)') '13.2    7.698426668271358'
write(u1,'(a)') '13.219999999999999      7.8028050644128655'
write(u1,'(a)') '13.24   7.9078337932406235'
write(u1,'(a)') '13.26   8.01351285475461'
write(u1,'(a)') '13.280000000000001      8.119842248954846'
write(u1,'(a)') '13.3    8.226821975841315'
write(u1,'(a)') '13.32   8.334452035414023'
write(u1,'(a)') '13.34   8.44273242767297'
write(u1,'(a)') '13.36   8.551663152618158'
write(u1,'(a)') '13.379999999999999      8.661244210249585'
write(u1,'(a)') '13.4    8.77147560056726'
write(u1,'(a)') '13.42   8.882357323571167'
write(u1,'(a)') '13.440000000000001      8.993889379261326'
write(u1,'(a)') '13.46   9.106071767637712'
write(u1,'(a)') '13.48   9.218904488700339'
write(u1,'(a)') '13.5    9.332387542449206'
write(u1,'(a)') '13.52   9.446520928884313'
write(u1,'(a)') '13.54   9.56130464800566'
write(u1,'(a)') '13.56   9.676738699813257'
write(u1,'(a)') '13.58   9.792823084307082'
write(u1,'(a)') '13.600000000000001      9.909557801487159'
write(u1,'(a)') '13.620000000000001      10.026942851353466'
write(u1,'(a)') '13.64   10.144978233906013'
write(u1,'(a)') '13.66   10.263663949144798'
write(u1,'(a)') '13.68   10.382999997069824'
write(u1,'(a)') '13.7    10.50298637768109'
write(u1,'(a)') '13.719999999999999      10.623623090978596'
write(u1,'(a)') '13.74   10.744910136962353'
write(u1,'(a)') '13.76   10.866847515632339'
write(u1,'(a)') '13.780000000000001      10.989435226988574'
write(u1,'(a)') '13.8    11.11267327103104'
write(u1,'(a)') '13.82   11.236561647759745'
write(u1,'(a)') '13.84   11.361100357174692'
write(u1,'(a)') '13.86   11.486289399275877'
write(u1,'(a)') '13.879999999999999      11.612128774063303'
write(u1,'(a)') '13.9    11.738618481536978'
write(u1,'(a)') '13.92   11.865758521696883'
write(u1,'(a)') '13.940000000000001      11.993548894543041'
write(u1,'(a)') '13.96   12.121989600075425'
write(u1,'(a)') '13.98   12.25108063829405'
write(u1,'(a)') '14.     12.380822009198916'
write(u1,'(a)') 'end'
close(u1)

 cmd = "./" // trim(executable) // " < "  // trim(FileInput) // " > " // trim(FileOutput)
call execute_command_line(cmd)

open(unit=u1, file=FileOutput, status='old', action='read')

! The output should have 1566 lines but we try to read more than that, to allow for
! future further output
nlines = 1600
iflag = 0 !0 means the test has failed

n_subtests(itest)=6
SubtestDescription(1) = 'Check mass'
SubtestDescription(2) = 'Check pot. r=6'
SubtestDescription(3) = 'Check energy v=39'
SubtestDescription(4) = 'Check energy v=130'
SubtestDescription(5) = 'Check energy v=199'
SubtestDescription(6) = 'Check ZPE'

do i=1, nlines
   read(u1, '(a500)', iostat=ierr) line
   if(ierr /=0) exit
   line = ADJUSTL(line) ! remove leading spaces if present
   if( trim(line) == 'Reduced mass is      20000.0000000000         atomic mass units (Daltons)')        iflag(1) = 1
   if( trim(line) == '6.00000000        1.23808220E+01')                                                 iflag(2) = 1
   if( trim(line) == '0.0   40          1.712003   1    39   0     0.0     0.0     0.0   +    ||DWell')  iflag(3) = 1
   if( trim(line) == '0.0  131          6.318107   1   130   0     0.0     0.0     0.0   +    ||DWell')  iflag(4) = 1
   if( trim(line) == '0.0  200          9.897210   1   199   0     0.0     0.0     0.0   +    ||DWell')  iflag(5) = 1
   if( trim(line) == 'Zero point energy (ZPE) =           0.052749')                                     iflag(6) = 1


enddo
close(u1)

write(*,'(a)')
do i=1, n_subtests(itest)
  write(line,'(A, I6,A)') 'Input file ' // trim(FileInput) //  " ( "  // trim(TestFileName(itest)) // " )" // &
                                          ' , subtest ', i, " (" // trim(SubtestDescription(i)) // ") : "
  write(*,'(a120)',advance='no') line
  if(iflag(i) == 1) then
     write(*,'(a)') 'Passed'
  else
      write(*,'(a)') 'FAILED! <---'
      n_failed_subtests=n_failed_subtests+1
  endif

enddo



!***************************************************************
!**************** END OF TEST 002 ******************************
!***************************************************************


!***************************************************************
!**************** TEST 003 *************************************
!***************************************************************
3 continue
itest=3
if(z_run_test_file(itest) .eqv. .false. ) goto 4 ! skip test if flag is set
TestFileName(itest) = 'Morse potential (interpolated)'
write(FileInput,  '(a,i4.4,a)') 'duo_test_', itest, '.inp'
write(FileOutput, '(a,i4.4,a)') 'duo_test_', itest, '.out'

open(unit=u1, file=FileInput, status='unknown', action='write')


write(u1,'(a)') 'masses 1.00000 1.000000'
write(u1,'(a)') 'molecule fictitious'
write(u1,'(a)') ''
write(u1,'(a)') '(Total number of states taken into account)'
write(u1,'(a)') 'nstates 1'
write(u1,'(a)') ''
write(u1,'(a)') '(Total angular momentum quantum - a value or an interval)'
write(u1,'(a)') 'jrot 0.0'
write(u1,'(a)') ''
write(u1,'(a)') 'SYMMETRY Cs(M)'
write(u1,'(a)') ''
write(u1,'(a)') '(SolutionMethod 5PointDifferences )'
write(u1,'(a)') 'SolutionMethod Sinc'
write(u1,'(a)') ''
write(u1,'(a)') '(Defining the integration grid)'
write(u1,'(a)') 'grid'
write(u1,'(a)') ' npoints 50'
write(u1,'(a)') ' range 0.30, 4.50'
write(u1,'(a)') ' type 0'
write(u1,'(a)') 'end'
write(u1,'(a)') ''
write(u1,'(a)') 'DIAGONALIZER'
write(u1,'(a)') ' SYEV'
write(u1,'(a)') ' enermax 15000.0'
write(u1,'(a)') ' nroots 100'
write(u1,'(a)') 'end'
write(u1,'(a)') ''
write(u1,'(a)') 'CONTRACTION'
write(u1,'(a)') ' vib'
write(u1,'(a)') ' vmax 10'
write(u1,'(a)') 'END'
write(u1,'(a)') ''
write(u1,'(a)') '( sampling of 40000*(1-exp( -(r-2)) )^2 )'
write(u1,'(a)') 'poten 1'
write(u1,'(a)') 'name "Morse"'
write(u1,'(a)') 'type grid'
write(u1,'(a)') 'lambda 0'
write(u1,'(a)') 'mult 1'
write(u1,'(a)') 'symmetry +'
write(u1,'(a)') 'units cm-1'
write(u1,'(a)') 'units angstroms'
write(u1,'(a)') 'interpolationType CubicSplines'
write(u1,'(a)') 'values'
write(u1,'(a)') '1.45 21506.3995'
write(u1,'(a)') '1.50 16833.5715'
write(u1,'(a)') '1.55 12919.1496'
write(u1,'(a)') '1.60 9675.6613'
write(u1,'(a)') '1.65 7024.7044'
write(u1,'(a)') '1.70 4896.0474'
write(u1,'(a)') '1.75 3226.8175'
write(u1,'(a)') '1.80 1960.7673'
write(u1,'(a)') '1.85 1047.6129'
write(u1,'(a)') '1.90 442.4369'
write(u1,'(a)') '1.95 105.1490'
write(u1,'(a)') '2.00 0.0000'
write(u1,'(a)') '2.05 95.1428'
write(u1,'(a)') '2.10 362.2367'
write(u1,'(a)') '2.15 776.0907'
write(u1,'(a)') '2.20 1314.3416'
write(u1,'(a)') '2.25 1957.1637'
write(u1,'(a)') '2.30 2687.0078'
write(u1,'(a)') '2.35 3488.3650'
write(u1,'(a)') '2.40 4347.5549'
write(u1,'(a)') '2.45 5252.5343'
write(u1,'(a)') '2.50 6192.7249'
write(u1,'(a)') '2.55 7158.8585'
write(u1,'(a)') '2.60 8142.8376'
write(u1,'(a)') '2.65 9137.6096'
write(u1,'(a)') '2.70 10137.0543'
write(u1,'(a)') '2.75 11135.8822'
write(u1,'(a)') '2.80 12129.5436'
write(u1,'(a)') '2.85 13114.1464'
write(u1,'(a)') '2.90 14086.3827'
write(u1,'(a)') '2.95 15043.4629'
write(u1,'(a)') '3.00 15983.0560'
write(u1,'(a)') '3.05 16903.2372'
write(u1,'(a)') '3.10 17802.4396'
write(u1,'(a)') '3.15 18679.4122'
write(u1,'(a)') '3.20 19533.1812'
write(u1,'(a)') '3.25 20363.0162'
write(u1,'(a)') '3.30 21168.3997'
write(u1,'(a)') '3.35 21948.9997'
write(u1,'(a)') '3.40 22704.6454'
write(u1,'(a)') '3.45 23435.3058'
write(u1,'(a)') '3.50 24141.0699'
write(u1,'(a)') 'end'
write(u1,'(a)') ''
write(u1,'(a)') '( totaly fictitious dipoles)'
write(u1,'(a)') 'dipole 1 1'
write(u1,'(a)') 'name "Simple-dipole"'
write(u1,'(a)') 'type grid'
write(u1,'(a)') 'lambda 0 0'
write(u1,'(a)') 'mult 1 1'
write(u1,'(a)') 'values'
write(u1,'(a)') '1.0000 0.334512'
write(u1,'(a)') '1.2000 0.203446'
write(u1,'(a)') '1.4000 0.082607'
write(u1,'(a)') '1.6000 -0.011786'
write(u1,'(a)') '1.8000 -0.073148'
write(u1,'(a)') '2.0000 -0.102422'
write(u1,'(a)') '2.2000 -0.105441'
write(u1,'(a)') '2.4000 -0.090370'
write(u1,'(a)') '2.6000 -0.065617'
write(u1,'(a)') '2.8000 -0.038387'
write(u1,'(a)') '3.0000 -0.013911'
write(u1,'(a)') '3.2000 0.004751'
write(u1,'(a)') '3.4000 0.016490'
write(u1,'(a)') 'end'
write(u1,'(a)') ''
write(u1,'(a)') 'END'
write(u1,'(a)') ''
write(u1,'(a)') 'poten 1'
write(u1,'(a)') 'name "null"'
write(u1,'(a)') 'type grid'
write(u1,'(a)') 'lambda 0'
write(u1,'(a)') 'mult 1'
write(u1,'(a)') 'symmetry +'
write(u1,'(a)') 'units cm-1'
write(u1,'(a)') 'units angstroms'
write(u1,'(a)') 'values'
write(u1,'(a)') '1.1 0.00000'
write(u1,'(a)') '1.2 0.00000'
write(u1,'(a)') '1.3 0.00000'
write(u1,'(a)') '1.4 0.00000'
write(u1,'(a)') '1.5 0.00000'
write(u1,'(a)') 'end'
close(u1)

 cmd = "./" // trim(executable) // " < "  // trim(FileInput) // " > " // trim(FileOutput)
call execute_command_line(cmd)

open(unit=u1, file=FileOutput, status='old', action='read')

! The output should have 463 lines but we try to read more than that, to allow for
! future further output
nlines = 500
iflag = 0 !0 means the test has failed

n_subtests(itest)=6
SubtestDescription(1) = 'Check mass'
SubtestDescription(2) = 'Check pot. r=0.3'
SubtestDescription(3) = 'Check dip. r=0.3'
SubtestDescription(4) = 'Check energy v=1'
SubtestDescription(5) = 'Check energy v=6'
SubtestDescription(6) = 'Check ZPE'

do i=1, nlines
   read(u1, '(a500)', iostat=ierr) line
   if(ierr /=0) exit
   line = ADJUSTL(line) ! remove leading spaces if present
   if( trim(line) == 'Reduced mass is     0.500000000000000         atomic mass units (Daltons)')    iflag(1) = 1
   if( trim(line) == '0.30000000        5.58881619E+05')                                             iflag(2) = 1
   if( trim(line) == '0.30000000         0.27357595')                                                iflag(3) = 1
   if( trim(line) == '0.0    2       2255.168415   1     1   0     0.0     0.0     0.0   +    ||Morse') iflag(4) = 1
   if( trim(line) == '0.0    7      12519.314535   1     6   0     0.0     0.0     0.0   +    ||Morse') iflag(5) = 1
   if( trim(line) == 'Zero point energy (ZPE) =        1152.863941')                                 iflag(6) = 1
enddo
close(u1)

write(*,'(a)')
do i=1, n_subtests(itest)
  write(line,'(A, I6,A)') 'Input file ' // trim(FileInput) //  " ( "  // trim(TestFileName(itest)) // " )" // &
                                          ' , subtest ', i, " (" // trim(SubtestDescription(i)) // ") : "
  write(*,'(a120)',advance='no') line
  if(iflag(i) == 1) then
     write(*,'(a)') 'Passed'
  else
      write(*,'(a)') 'FAILED! <---'
      n_failed_subtests=n_failed_subtests+1
  endif

enddo


!***************************************************************
!**************** END OF TEST 003 ******************************
!***************************************************************


!***************************************************************
!**************** TEST 004 *************************************
!***************************************************************
4 continue
itest=4
if(z_run_test_file(itest) .eqv. .false. ) goto 5 ! skip test if flag is set
TestFileName(itest) = 'Single well'
write(FileInput,  '(a,i4.4,a)') 'duo_test_', itest, '.inp'
write(FileOutput, '(a,i4.4,a)') 'duo_test_', itest, '.out'

open(unit=u1, file=FileInput, status='unknown', action='write')

write(u1, '(a)') 'masses 40000. 40000.'
write(u1, '(a)') ''
write(u1, '(a)') 'nstates 1'
write(u1, '(a)') 'jrot 0.0'
write(u1, '(a)') ''
write(u1, '(a)') 'SolutionMethod Sinc'
write(u1, '(a)') ''
write(u1, '(a)') 'grid'
write(u1, '(a)') ' npoints 200'
write(u1, '(a)') ' range 7.00, 13.00'
write(u1, '(a)') ' type 0'
write(u1, '(a)') 'end'
write(u1, '(a)') ''
write(u1, '(a)') 'CONTRACTION'
write(u1, '(a)') ' vib'
write(u1, '(a)') ' vmax 50'
write(u1, '(a)') 'END'
write(u1, '(a)') ''
write(u1, '(a)') 'DIAGONALIZER'
write(u1, '(a)') ' SYEV'
write(u1, '(a)') ' enermax 15000.0'
write(u1, '(a)') ' nroots 100'
write(u1, '(a)') 'end'
write(u1, '(a)') ''
write(u1, '(a)') ''
write(u1, '(a)') 'poten 1'
write(u1, '(a)') 'name "SWell"'
write(u1, '(a)') 'type polynomial'
write(u1, '(a)') 'lambda 0'
write(u1, '(a)') 'mult 1'
write(u1, '(a)') 'symmetry +'
write(u1, '(a)') 'units cm-1'
write(u1, '(a)') 'values'
write(u1, '(a)') 'a0 0.00000'
write(u1, '(a)') 're 9.278986556699559'
write(u1, '(a)') 'a1 0.00000'
write(u1, '(a)') 'a2 3.380822009198914'
write(u1, '(a)') 'end '
close(u1)

 cmd = "./" // trim(executable) // " < "  // trim(FileInput) // " > " // trim(FileOutput)
call execute_command_line(cmd)

open(unit=u1, file=FileOutput, status='old', action='read')

! The output should have 573 lines but we try to read more than that, to allow for
! future further output
nlines = 600
iflag = 0 !0 means the test has failed

n_subtests(itest)=5
SubtestDescription(1) = 'Check mass'
SubtestDescription(2) = 'Check pot. r=7'
SubtestDescription(3) = 'Check energy v=9'
SubtestDescription(4) = 'Check energy v=49'
SubtestDescription(5) = 'Check ZPE'

do i=1, nlines
   read(u1, '(a500)', iostat=ierr) line
   if(ierr /=0) exit
   line = ADJUSTL(line) ! remove leading spaces if present
   if( trim(line) == 'Reduced mass is      20000.0000000000         atomic mass units (Daltons)')    iflag(1) = 1
   if( trim(line) == '7.00000000        1.75592448E+01')                                             iflag(2) = 1
   if( trim(line) == '0.0   10          0.960875   1     9   0     0.0     0.0     0.0   +    ||SWell') iflag(3) = 1
   if( trim(line) == '0.0   50          5.231431   1    49   0     0.0     0.0     0.0   +    ||SWell') iflag(4) = 1
   if( trim(line) == 'Zero point energy (ZPE) =           0.053382') iflag(5) = 1
enddo
close(u1)

write(*,'(a)')
do i=1, n_subtests(itest)
  write(line,'(A, I6,A)') 'Input file ' // trim(FileInput) //  " ( "  // trim(TestFileName(itest)) // " )" // &
                                          ' , subtest ', i, " (" // trim(SubtestDescription(i)) // ") : "
  write(*,'(a120)',advance='no') line
  if(iflag(i) == 1) then
     write(*,'(a)') 'Passed'
  else
      write(*,'(a)') 'FAILED! <---'
      n_failed_subtests=n_failed_subtests+1
  endif

enddo
!***************************************************************
!**************** END OF TEST 004 ******************************
!***************************************************************


!***************************************************************
!**************** TEST 005 *************************************
!***************************************************************
5 continue
itest=5
if(z_run_test_file(itest) .eqv. .false. ) goto 6 ! skip test if flag is set
TestFileName(itest) = 'Single well with intensities'
write(FileInput,  '(a,i4.4,a)') 'duo_test_', itest, '.inp'
write(FileOutput, '(a,i4.4,a)') 'duo_test_', itest, '.out'

open(unit=u1, file=FileInput, status='unknown', action='write')
write(u1,'(a)') 'masses 40000. 40000. '
write(u1,'(a)') ' '
write(u1,'(a)') 'nstates 1 '
write(u1,'(a)') 'jrot 0.0 '
write(u1,'(a)') ' '
write(u1,'(a)') 'SolutionMethod Sinc '
write(u1,'(a)') ' '
write(u1,'(a)') 'symmetry Cs(M) '
write(u1,'(a)') ' '
write(u1,'(a)') 'grid '
write(u1,'(a)') ' npoints 350 '
write(u1,'(a)') ' range 6.00, 14.00 '
write(u1,'(a)') ' type 0 '
write(u1,'(a)') 'end '
write(u1,'(a)') ' '
write(u1,'(a)') 'CONTRACTION '
write(u1,'(a)') ' vib '
write(u1,'(a)') ' vmax 100 '
write(u1,'(a)') 'END '
write(u1,'(a)') ' '
write(u1,'(a)') 'DIAGONALIZER '
write(u1,'(a)') ' SYEV '
write(u1,'(a)') ' enermax 15000.0 '
write(u1,'(a)') ' nroots 100 '
write(u1,'(a)') 'end '
write(u1,'(a)') ' '
write(u1,'(a)') ' '
write(u1,'(a)') 'poten 1 '
write(u1,'(a)') 'name "DWell" '
write(u1,'(a)') 'type polynomial '
write(u1,'(a)') 'lambda 0 '
write(u1,'(a)') 'mult 1 '
write(u1,'(a)') 'symmetry + '
write(u1,'(a)') 'units cm-1 '
write(u1,'(a)') 'values '
write(u1,'(a)') 'a0 -0.6258317155998449 '
write(u1,'(a)') 're 10.000000 '
write(u1,'(a)') 'a1 0.00000 '
write(u1,'(a)') 'a2 0.8129158577999225 '
write(u1,'(a)') 'end '
write(u1,'(a)') ' '
write(u1,'(a)') ' '
write(u1,'(a)') 'dipole 1 1 '
write(u1,'(a)') 'name "barrier" '
write(u1,'(a)') 'type grid '
write(u1,'(a)') 'lambda 0 0 '
write(u1,'(a)') 'mult 1 1 '
write(u1,'(a)') 'symmetry + + '
write(u1,'(a)') 'values '
write(u1,'(a)') '6. 2.607526611678407e-28 '
write(u1,'(a)') '6.02 4.93721848873338e-28 '
write(u1,'(a)') '6.04 9.318504119925838e-28 '
write(u1,'(a)') '6.06 1.7531549978712018e-27 '
write(u1,'(a)') '6.08 3.287794649354648e-27 '
write(u1,'(a)') '6.1 6.146095613816321e-27 '
write(u1,'(a)') '6.12 1.1452602773924518e-26 '
write(u1,'(a)') '6.14 2.1272539607860593e-26 '
write(u1,'(a)') '6.16 3.938625984832693e-26 '
write(u1,'(a)') '6.18 7.269095930870267e-26 '
write(u1,'(a)') '6.2 1.3372922153091908e-25 '
write(u1,'(a)') '6.22 2.452350211324485e-25 '
write(u1,'(a)') '6.24 4.482795443163941e-25 '
write(u1,'(a)') '6.26 8.168185855293816e-25 '
write(u1,'(a)') '6.28 1.4835854237791893e-24 '
write(u1,'(a)') '6.3 2.6860232167683993e-24 '
write(u1,'(a)') '6.32 4.847493376142648e-24 '
write(u1,'(a)') '6.34 8.720370640704862e-24 '
write(u1,'(a)') '6.36 1.5637341528698e-23 '
write(u1,'(a)') '6.38 2.7951244978987925e-23 '
write(u1,'(a)') '6.4 4.9802328585657665e-23 '
write(u1,'(a)') '6.42 8.845215122566763e-23 '
write(u1,'(a)') '6.44 1.5659482572846938e-22 '
write(u1,'(a)') '6.46 2.7634822125305835e-22 '
write(u1,'(a)') '6.48 4.861230308135457e-22 '
write(u1,'(a)') '6.5 8.52404979276005e-22 '
write(u1,'(a)') '6.52 1.4898962143063935e-21 '
write(u1,'(a)') '6.54 2.595830392925519e-21 '
write(u1,'(a)') '6.5600000000000005 4.508238328598644e-21 '
write(u1,'(a)') '6.58 7.804546610854037e-21 '
write(u1,'(a)') '6.6 1.3467864800875917e-20 '
write(u1,'(a)') '6.62 2.316648130629384e-20 '
write(u1,'(a)') '6.64 3.972205053939896e-20 '
write(u1,'(a)') '6.66 6.789120264058192e-20 '
write(u1,'(a)') '6.68 1.1566596936355828e-19 '
write(u1,'(a)') '6.7 1.9643006461708697e-19 '
write(u1,'(a)') '6.72 3.3252214610449074e-19 '
write(u1,'(a)') '6.74 5.611041072370278e-19 '
write(u1,'(a)') '6.76 9.437926227322794e-19 '
write(u1,'(a)') '6.78 1.5824134533609102e-18 '
write(u1,'(a)') '6.8 2.6446829322028795e-18 '
write(u1,'(a)') '6.82 4.405929219681732e-18 '
write(u1,'(a)') '6.84 7.316639835535116e-18 '
write(u1,'(a)') '6.86 1.2111448759406027e-17 '
write(u1,'(a)') '6.88 1.998438430277145e-17 '
write(u1,'(a)') '6.9 3.286969746739775e-17 '
write(u1,'(a)') '6.92 5.38903369092568e-17 '
write(u1,'(a)') '6.94 8.807169368199645e-17 '
write(u1,'(a)') '6.96 1.4347361973795417e-16 '
write(u1,'(a)') '6.98 2.329796670135002e-16 '
write(u1,'(a)') '7. 3.771153782467911e-16 '
write(u1,'(a)') '7.02 6.084721747184384e-16 '
write(u1,'(a)') '7.04 9.786276016714005e-16 '
write(u1,'(a)') '7.0600000000000005 1.5689332106422145e-15 '
write(u1,'(a)') '7.08 2.507273492423265e-15 '
write(u1,'(a)') '7.1 3.994010662977935e-15 '
write(u1,'(a)') '7.12 6.342010962936676e-15 '
write(u1,'(a)') '7.140000000000001 1.0038180821206298e-14 '
write(u1,'(a)') '7.16 1.5837743251486305e-14 '
write(u1,'(a)') '7.18 2.4908171061980403e-14 '
write(u1,'(a)') '7.2 3.9048165804186394e-14 '
write(u1,'(a)') '7.22 6.101964763801575e-14 '
write(u1,'(a)') '7.24 9.504931904568e-14 '
write(u1,'(a)') '7.26 1.4758376379681298e-13 '
write(u1,'(a)') '7.28 2.2842225340107596e-13 '
write(u1,'(a)') '7.3 3.524102235109213e-13 '
write(u1,'(a)') '7.32 5.419620089198879e-13 '
write(u1,'(a)') '7.34 8.308056600413588e-13 '
write(u1,'(a)') '7.36 1.269522221007627e-12 '
write(u1,'(a)') '7.38 1.9337103640460677e-12 '
write(u1,'(a)') '7.4 2.935978024136269e-12 '
write(u1,'(a)') '7.42 4.4434923196863535e-12 '
write(u1,'(a)') '7.4399999999999995 6.7035727179349165e-12 '
write(u1,'(a)') '7.46 1.008088065865524e-11 '
write(u1,'(a)') '7.48 1.511126701864662e-11 '
write(u1,'(a)') '7.5 2.2579459600128797e-11 '
write(u1,'(a)') '7.52 3.363074323137583e-11 '
write(u1,'(a)') '7.54 4.993092876877705e-11 '
write(u1,'(a)') '7.5600000000000005 7.389466445025507e-11 '
write(u1,'(a)') '7.58 1.0901010951332925e-10 '
write(u1,'(a)') '7.6 1.602989686558041e-10 '
write(u1,'(a)') '7.62 2.349659536955435e-10 '
write(u1,'(a)') '7.640000000000001 3.433123333747632e-10 '
write(u1,'(a)') '7.66 5.000162639319345e-10 '
write(u1,'(a)') '7.68 7.259205878802915e-10 '
write(u1,'(a)') '7.7 1.050520070490049e-9 '
write(u1,'(a)') '7.72 1.5154089655265611e-9 '
write(u1,'(a)') '7.74 2.1790420350857133e-9 '
write(u1,'(a)') '7.76 3.123285044451754e-9 '
write(u1,'(a)') '7.78 4.462393792287464e-9 '
write(u1,'(a)') '7.8 6.355276080754993e-9 '
write(u1,'(a)') '7.82 9.022174196551332e-9 '
write(u1,'(a)') '7.84 1.2767276152388232e-8 '
write(u1,'(a)') '7.86 1.800924697418566e-8 '
write(u1,'(a)') '7.88 2.532229788157695e-8 '
write(u1,'(a)') '7.9 3.549122017130936e-8 '
write(u1,'(a)') '7.92 4.958485016114015e-8 '
write(u1,'(a)') '7.9399999999999995 6.905376104322296e-8 '
write(u1,'(a)') '7.96 9.58596695753113e-8 '
write(u1,'(a)') '7.98 1.3264618851501676e-7 '
write(u1,'(a)') '8. 1.8296325617914135e-7 '
write(u1,'(a)') '8.02 2.515609602314059e-7 '
write(u1,'(a)') '8.04 3.4477269567450517e-7 '
write(u1,'(a)') '8.06 4.710128362241596e-7 '
write(u1,'(a)') '8.08 6.41420563303087e-7 '
write(u1,'(a)') '8.1 8.706894305848934e-7 '
write(u1,'(a)') '8.120000000000001 1.1781319228460575e-6 '
write(u1,'(a)') '8.14 1.589040011902242e-6 '
write(u1,'(a)') '8.16 2.136416874014602e-6 '
write(u1,'(a)') '8.18 2.8631719044812897e-6 '
write(u1,'(a)') '8.2 3.824891373509501e-6 '
write(u1,'(a)') '8.22 5.0933209704149755e-6 '
write(u1,'(a)') '8.24 6.760724560335833e-6 '
write(u1,'(a)') '8.26 8.945316218782198e-6 '
write(u1,'(a)') '8.280000000000001 0.000011798000498015955 '
write(u1,'(a)') '8.3 0.000015510699371040517 '
write(u1,'(a)') '8.32 0.000020326593793680243 '
write(u1,'(a)') '8.34 0.000026552663632388537 '
write(u1,'(a)') '8.36 0.000034574972002626226 '
write(u1,'(a)') '8.379999999999999 0.000044877208846031165 '
write(u1,'(a)') '8.4 0.00005806308360181939 '
write(u1,'(a)') '8.42 0.00007488323755347042 '
write(u1,'(a)') '8.44 0.00009626743193601047 '
write(u1,'(a)') '8.46 0.00012336285680121463 '
write(u1,'(a)') '8.48 0.0001575794960560621 '
write(u1,'(a)') '8.5 0.000200643573500087 '
write(u1,'(a)') '8.52 0.0002546601898866553 '
write(u1,'(a)') '8.54 0.0003221863380591267 '
write(u1,'(a)') '8.56 0.0004063155472853246 '
write(u1,'(a)') '8.58 0.0005107754533982561 '
write(u1,'(a)') '8.6 0.0006400396117467881 '
write(u1,'(a)') '8.620000000000001 0.0007994548579217274 '
write(u1,'(a)') '8.64 0.0009953854686280915 '
write(u1,'(a)') '8.66 0.0012353752731370166 '
write(u1,'(a)') '8.68 0.0015283287051959568 '
write(u1,'(a)') '8.7 0.0018847115565826776 '
write(u1,'(a)') '8.72 0.002316771887195828 '
write(u1,'(a)') '8.74 0.0028387811536722926 '
write(u1,'(a)') '8.76 0.003467295130908433 '
write(u1,'(a)') '8.780000000000001 0.004221433611875735 '
write(u1,'(a)') '8.8 0.0051231771761455274 '
write(u1,'(a)') '8.82 0.006197678514638444 '
write(u1,'(a)') '8.84 0.007473584888731782 '
write(u1,'(a)') '8.86 0.008983367291508983 '
write(u1,'(a)') '8.879999999999999 0.01076365077786742 '
write(u1,'(a)') '8.9 0.012855539254042904 '
write(u1,'(a)') '8.92 0.015304926787379906 '
write(u1,'(a)') '8.94 0.01816278624169761 '
write(u1,'(a)') '8.96 0.021485424796674712 '
write(u1,'(a)') '8.98 0.025334694712119032 '
write(u1,'(a)') '9. 0.029778146596777928 '
write(u1,'(a)') '9.02 0.03488911148903424 '
write(u1,'(a)') '9.04 0.040746697310637486 '
write(u1,'(a)') '9.06 0.04743568477497995 '
write(u1,'(a)') '9.08 0.05504630768042008 '
write(u1,'(a)') '9.1 0.06367390275835835 '
write(u1,'(a)') '9.120000000000001 0.07341841493378391 '
write(u1,'(a)') '9.14 0.08438374504592218 '
write(u1,'(a)') '9.16 0.09667692881283788 '
write(u1,'(a)') '9.18 0.11040713813918195 '
write(u1,'(a)') '9.2 0.12568449877893018 '
write(u1,'(a)') '9.22 0.1426187218753048 '
write(u1,'(a)') '9.24 0.1613175509882864 '
write(u1,'(a)') '9.26 0.18188503084378785 '
write(u1,'(a)') '9.280000000000001 0.20441960913087515 '
write(u1,'(a)') '9.3 0.2290120881427486 '
write(u1,'(a)') '9.32 0.2557434487865101 '
write(u1,'(a)') '9.34 0.2846825753349806 '
write(u1,'(a)') '9.36 0.3158839150982333 '
write(u1,'(a)') '9.379999999999999 0.3493851127712329 '
write(u1,'(a)') '9.4 0.3852046643713807 '
write(u1,'(a)') '9.42 0.4233396402113618 '
write(u1,'(a)') '9.44 0.4637635300521322 '
write(u1,'(a)') '9.46 0.5064242662467984 '
write(u1,'(a)') '9.48 0.5512424821310913 '
write(u1,'(a)') '9.5 0.5981100629736783 '
write(u1,'(a)') '9.52 0.646889045333342 '
write(u1,'(a)') '9.54 0.6974109175819261 '
write(u1,'(a)') '9.56 0.7494763695878685 '
write(u1,'(a)') '9.58 0.8028555331126858 '
write(u1,'(a)') '9.6 0.8572887464047096 '
write(u1,'(a)') '9.620000000000001 0.9124878668911746 '
write(u1,'(a)') '9.64 0.9681381449414747 '
write(u1,'(a)') '9.66 1.023900659627162 '
write(u1,'(a)') '9.68 1.0794153045174928 '
write(u1,'(a)') '9.7 1.1343042981494593 '
write(u1,'(a)') '9.72 1.1881761802614537 '
write(u1,'(a)') '9.74 1.2406302415711463 '
write(u1,'(a)') '9.76 1.2912613222161533 '
write(u1,'(a)') '9.780000000000001 1.3396649023663985 '
write(u1,'(a)') '9.8 1.3854423983526893 '
write(u1,'(a)') '9.82 1.4282065693003891 '
write(u1,'(a)') '9.84 1.4675869330313456 '
write(u1,'(a)') '9.86 1.5032350861667987 '
write(u1,'(a)') '9.879999999999999 1.5348298221257066 '
write(u1,'(a)') '9.9 1.562081942187235 '
write(u1,'(a)') '9.92 1.5847386590091346 '
write(u1,'(a)') '9.94 1.602587498914021 '
write(u1,'(a)') '9.96 1.615459618733379 '
write(u1,'(a)') '9.98 1.623232464810024 '
write(u1,'(a)') '10. 1.625831715599845 '
write(u1,'(a)') '10.02 1.623232464810024 '
write(u1,'(a)') '10.04 1.615459618733379 '
write(u1,'(a)') '10.06 1.602587498914021 '
write(u1,'(a)') '10.08 1.5847386590091346 '
write(u1,'(a)') '10.1 1.562081942187235 '
write(u1,'(a)') '10.120000000000001 1.5348298221257066 '
write(u1,'(a)') '10.14 1.5032350861667987 '
write(u1,'(a)') '10.16 1.4675869330313456 '
write(u1,'(a)') '10.18 1.4282065693003891 '
write(u1,'(a)') '10.2 1.3854423983526893 '
write(u1,'(a)') '10.219999999999999 1.3396649023663985 '
write(u1,'(a)') '10.24 1.2912613222161533 '
write(u1,'(a)') '10.26 1.2406302415711463 '
write(u1,'(a)') '10.280000000000001 1.188176180261449 '
write(u1,'(a)') '10.3 1.1343042981494593 '
write(u1,'(a)') '10.32 1.0794153045174928 '
write(u1,'(a)') '10.34 1.023900659627162 '
write(u1,'(a)') '10.36 0.9681381449414747 '
write(u1,'(a)') '10.379999999999999 0.9124878668911746 '
write(u1,'(a)') '10.4 0.8572887464047096 '
write(u1,'(a)') '10.42 0.8028555331126858 '
write(u1,'(a)') '10.440000000000001 0.7494763695878638 '
write(u1,'(a)') '10.46 0.6974109175819261 '
write(u1,'(a)') '10.48 0.646889045333342 '
write(u1,'(a)') '10.5 0.5981100629736783 '
write(u1,'(a)') '10.52 0.5512424821310913 '
write(u1,'(a)') '10.54 0.5064242662467984 '
write(u1,'(a)') '10.56 0.4637635300521322 '
write(u1,'(a)') '10.58 0.4233396402113618 '
write(u1,'(a)') '10.600000000000001 0.3852046643713774 '
write(u1,'(a)') '10.620000000000001 0.3493851127712329 '
write(u1,'(a)') '10.64 0.3158839150982333 '
write(u1,'(a)') '10.66 0.2846825753349806 '
write(u1,'(a)') '10.68 0.2557434487865101 '
write(u1,'(a)') '10.7 0.2290120881427486 '
write(u1,'(a)') '10.719999999999999 0.20441960913087515 '
write(u1,'(a)') '10.74 0.18188503084378785 '
write(u1,'(a)') '10.76 0.1613175509882864 '
write(u1,'(a)') '10.780000000000001 0.1426187218753032 '
write(u1,'(a)') '10.8 0.12568449877893018 '
write(u1,'(a)') '10.82 0.11040713813918195 '
write(u1,'(a)') '10.84 0.09667692881283788 '
write(u1,'(a)') '10.86 0.08438374504592218 '
write(u1,'(a)') '10.879999999999999 0.07341841493378391 '
write(u1,'(a)') '10.9 0.06367390275835835 '
write(u1,'(a)') '10.92 0.05504630768042008 '
write(u1,'(a)') '10.940000000000001 0.047435684774979316 '
write(u1,'(a)') '10.96 0.040746697310637486 '
write(u1,'(a)') '10.98 0.03488911148903424 '
write(u1,'(a)') '11. 0.029778146596777928 '
write(u1,'(a)') '11.02 0.025334694712119032 '
write(u1,'(a)') '11.04 0.021485424796674712 '
write(u1,'(a)') '11.06 0.01816278624169761 '
write(u1,'(a)') '11.08 0.015304926787379906 '
write(u1,'(a)') '11.100000000000001 0.01285553925404271 '
write(u1,'(a)') '11.120000000000001 0.01076365077786742 '
write(u1,'(a)') '11.14 0.008983367291508983 '
write(u1,'(a)') '11.16 0.007473584888731782 '
write(u1,'(a)') '11.18 0.006197678514638444 '
write(u1,'(a)') '11.2 0.0051231771761455274 '
write(u1,'(a)') '11.219999999999999 0.004221433611875735 '
write(u1,'(a)') '11.24 0.003467295130908433 '
write(u1,'(a)') '11.26 0.0028387811536722926 '
write(u1,'(a)') '11.280000000000001 0.0023167718871957847 '
write(u1,'(a)') '11.3 0.0018847115565826776 '
write(u1,'(a)') '11.32 0.0015283287051959568 '
write(u1,'(a)') '11.34 0.0012353752731370166 '
write(u1,'(a)') '11.36 0.0009953854686280915 '
write(u1,'(a)') '11.379999999999999 0.0007994548579217274 '
write(u1,'(a)') '11.4 0.0006400396117467881 '
write(u1,'(a)') '11.42 0.0005107754533982561 '
write(u1,'(a)') '11.440000000000001 0.0004063155472853159 '
write(u1,'(a)') '11.46 0.0003221863380591267 '
write(u1,'(a)') '11.48 0.0002546601898866553 '
write(u1,'(a)') '11.5 0.000200643573500087 '
write(u1,'(a)') '11.52 0.0001575794960560621 '
write(u1,'(a)') '11.54 0.00012336285680121463 '
write(u1,'(a)') '11.56 0.00009626743193601047 '
write(u1,'(a)') '11.58 0.00007488323755347042 '
write(u1,'(a)') '11.600000000000001 0.00005806308360181805 '
write(u1,'(a)') '11.620000000000001 0.000044877208846031165 '
write(u1,'(a)') '11.64 0.000034574972002626226 '
write(u1,'(a)') '11.66 0.000026552663632388537 '
write(u1,'(a)') '11.68 0.000020326593793680243 '
write(u1,'(a)') '11.7 0.000015510699371040517 '
write(u1,'(a)') '11.719999999999999 0.000011798000498015955 '
write(u1,'(a)') '11.74 8.945316218782198e-6 '
write(u1,'(a)') '11.76 6.760724560335833e-6 '
write(u1,'(a)') '11.780000000000001 5.093320970414848e-6 '
write(u1,'(a)') '11.8 3.824891373509501e-6 '
write(u1,'(a)') '11.82 2.8631719044812897e-6 '
write(u1,'(a)') '11.84 2.136416874014602e-6 '
write(u1,'(a)') '11.86 1.589040011902242e-6 '
write(u1,'(a)') '11.879999999999999 1.1781319228460575e-6 '
write(u1,'(a)') '11.9 8.706894305848934e-7 '
write(u1,'(a)') '11.92 6.41420563303087e-7 '
write(u1,'(a)') '11.940000000000001 4.7101283622414625e-7 '
write(u1,'(a)') '11.96 3.4477269567450517e-7 '
write(u1,'(a)') '11.98 2.515609602314059e-7 '
write(u1,'(a)') '12. 1.8296325617914135e-7 '
write(u1,'(a)') '12.02 1.3264618851501676e-7 '
write(u1,'(a)') '12.04 9.585966957531267e-8 '
write(u1,'(a)') '12.06 6.905376104322296e-8 '
write(u1,'(a)') '12.08 4.958485016114015e-8 '
write(u1,'(a)') '12.100000000000001 3.5491220171308225e-8 '
write(u1,'(a)') '12.120000000000001 2.532229788157659e-8 '
write(u1,'(a)') '12.14 1.8009246974185343e-8 '
write(u1,'(a)') '12.16 1.2767276152388232e-8 '
write(u1,'(a)') '12.18 9.022174196551332e-9 '
write(u1,'(a)') '12.2 6.355276080755083e-9 '
write(u1,'(a)') '12.219999999999999 4.462393792287544e-9 '
write(u1,'(a)') '12.24 3.123285044451754e-9 '
write(u1,'(a)') '12.26 2.1790420350857133e-9 '
write(u1,'(a)') '12.280000000000001 1.5154089655265396e-9 '
write(u1,'(a)') '12.3 1.0505200704900302e-9 '
write(u1,'(a)') '12.32 7.259205878802915e-10 '
write(u1,'(a)') '12.34 5.000162639319345e-10 '
write(u1,'(a)') '12.36 3.433123333747632e-10 '
write(u1,'(a)') '12.379999999999999 2.3496595369554764e-10 '
write(u1,'(a)') '12.4 1.602989686558041e-10 '
write(u1,'(a)') '12.42 1.0901010951332925e-10 '
write(u1,'(a)') '12.440000000000001 7.389466445025272e-11 '
write(u1,'(a)') '12.46 4.9930928768776163e-11 '
write(u1,'(a)') '12.48 3.363074323137583e-11 '
write(u1,'(a)') '12.5 2.2579459600128797e-11 '
write(u1,'(a)') '12.52 1.511126701864662e-11 '
write(u1,'(a)') '12.54 1.0080880658655418e-11 '
write(u1,'(a)') '12.56 6.7035727179349165e-12 '
write(u1,'(a)') '12.58 4.4434923196863535e-12 '
write(u1,'(a)') '12.600000000000001 2.9359780241361543e-12 '
write(u1,'(a)') '12.620000000000001 1.933710364046033e-12 '
write(u1,'(a)') '12.64 1.2695222210076e-12 '
write(u1,'(a)') '12.66 8.308056600413588e-13 '
write(u1,'(a)') '12.68 5.419620089198879e-13 '
write(u1,'(a)') '12.7 3.5241022351092756e-13 '
write(u1,'(a)') '12.719999999999999 2.284222534010808e-13 '
write(u1,'(a)') '12.74 1.4758376379681298e-13 '
write(u1,'(a)') '12.76 9.504931904568e-14 '
write(u1,'(a)') '12.780000000000001 6.101964763801467e-14 '
write(u1,'(a)') '12.8 3.904816580418557e-14 '
write(u1,'(a)') '12.82 2.4908171061980403e-14 '
write(u1,'(a)') '12.84 1.5837743251486305e-14 '
write(u1,'(a)') '12.86 1.0038180821206298e-14 '
write(u1,'(a)') '12.879999999999999 6.3420109629368105e-15 '
write(u1,'(a)') '12.9 3.994010662977935e-15 '
write(u1,'(a)') '12.92 2.507273492423265e-15 '
write(u1,'(a)') '12.940000000000001 1.5689332106421478e-15 '
write(u1,'(a)') '12.96 9.786276016713798e-16 '
write(u1,'(a)') '12.98 6.084721747184384e-16 '
write(u1,'(a)') '13. 3.771153782467911e-16 '
write(u1,'(a)') '13.02 2.329796670135002e-16 '
write(u1,'(a)') '13.04 1.4347361973795723e-16 '
write(u1,'(a)') '13.06 8.807169368199457e-17 '
write(u1,'(a)') '13.08 5.38903369092568e-17 '
write(u1,'(a)') '13.100000000000001 3.286969746739635e-17 '
write(u1,'(a)') '13.120000000000001 1.9984384302771026e-17 '
write(u1,'(a)') '13.14 1.211144875940577e-17 '
write(u1,'(a)') '13.16 7.316639835535116e-18 '
write(u1,'(a)') '13.18 4.405929219681732e-18 '
write(u1,'(a)') '13.2 2.6446829322029547e-18 '
write(u1,'(a)') '13.219999999999999 1.582413453360944e-18 '
write(u1,'(a)') '13.24 9.437926227322794e-19 '
write(u1,'(a)') '13.26 5.611041072370278e-19 '
write(u1,'(a)') '13.280000000000001 3.3252214610448366e-19 '
write(u1,'(a)') '13.3 1.9643006461708278e-19 '
write(u1,'(a)') '13.32 1.1566596936355828e-19 '
write(u1,'(a)') '13.34 6.789120264058192e-20 '
write(u1,'(a)') '13.36 3.972205053939981e-20 '
write(u1,'(a)') '13.379999999999999 2.3166481306294333e-20 '
write(u1,'(a)') '13.4 1.3467864800875917e-20 '
write(u1,'(a)') '13.42 7.804546610854037e-21 '
write(u1,'(a)') '13.440000000000001 4.50823832859842e-21 '
write(u1,'(a)') '13.46 2.5958303929254635e-21 '
write(u1,'(a)') '13.48 1.4898962143063935e-21 '
write(u1,'(a)') '13.5 8.52404979276005e-22 '
write(u1,'(a)') '13.52 4.861230308135457e-22 '
write(u1,'(a)') '13.54 2.763482212530662e-22 '
write(u1,'(a)') '13.56 1.5659482572846491e-22 '
write(u1,'(a)') '13.58 8.845215122566763e-23 '
write(u1,'(a)') '13.600000000000001 4.9802328585655185e-23 '
write(u1,'(a)') '13.620000000000001 2.795124497898713e-23 '
write(u1,'(a)') '13.64 1.5637341528697554e-23 '
write(u1,'(a)') '13.66 8.720370640704862e-24 '
write(u1,'(a)') '13.68 4.847493376142648e-24 '
write(u1,'(a)') '13.7 2.6860232167684757e-24 '
write(u1,'(a)') '13.719999999999999 1.4835854237792313e-24 '
write(u1,'(a)') '13.74 8.168185855293816e-25 '
write(u1,'(a)') '13.76 4.482795443163941e-25 '
write(u1,'(a)') '13.780000000000001 2.4523502113244156e-25 '
write(u1,'(a)') '13.8 1.3372922153091621e-25 '
write(u1,'(a)') '13.82 7.269095930870267e-26 '
write(u1,'(a)') '13.84 3.938625984832693e-26 '
write(u1,'(a)') '13.86 2.1272539607861196e-26 '
write(u1,'(a)') '13.879999999999999 1.1452602773924845e-26 '
write(u1,'(a)') '13.9 6.146095613816321e-27 '
write(u1,'(a)') '13.92 3.287794649354648e-27 '
write(u1,'(a)') '13.940000000000001 1.753154997871152e-27 '
write(u1,'(a)') '13.96 9.318504119925572e-28 '
write(u1,'(a)') '13.98 4.93721848873338e-28 '
write(u1,'(a)') '14. 2.607526611678407e-28 '
write(u1,'(a)') 'end '
write(u1,'(a)') ' '
write(u1,'(a)') 'INTENSITY '
write(u1,'(a)') 'absorption '
write(u1,'(a)') 'thresh_intes 1e-30 '
write(u1,'(a)') 'thresh_line 1e-30 '
write(u1,'(a)') 'temperature 300 '
write(u1,'(a)') 'qstat 1.0000000 '
write(u1,'(a)') 'gns 1.0 1.0 '
write(u1,'(a)') 'zpe 0.00 '
write(u1,'(a)') 'selection (rules) 1 1 '
write(u1,'(a)') 'linelist duo_test_0005'
write(u1,'(a)') 'J, 0, 2 '
write(u1,'(a)') 'freq-window 0.0 25.0 '
write(u1,'(a)') 'energy low -0.001 16.000, upper -0.00 16.0000 '
write(u1,'(a)') 'END '
close(u1)

 cmd = "./" // trim(executable) // " < "  // trim(FileInput) // " > " // trim(FileOutput)
call execute_command_line(cmd)

open(unit=u1, file=FileOutput, status='old', action='read')

! The output should have 15990 lines but we try to read more than that, to allow for
! future further output
nlines = 16000
iflag = 0 !0 means the test has failed

n_subtests(itest)=25
SubtestDescription(1) = 'Check mass'
SubtestDescription(2) = 'Check pot. r=6'
SubtestDescription(3) = 'Check dip. r=10'
SubtestDescription(4) = 'Check energy J=0 v=10'
SubtestDescription(5) = 'Check energy J=0 v=99'

! check a few transition dipoles. I chose a few reference values
! arbitrarily but so that a few are strong, a few weak and a few very weak.
SubtestDescription( 6) = 'Check <  0|mu|  0>'
SubtestDescription( 7) = 'Check < 50|mu| 50>'
SubtestDescription( 8) = 'Check < 99|mu| 99>'

SubtestDescription( 9) = 'Check <  2|mu| 20>'
SubtestDescription(10) = 'Check < 15|mu| 29>'
SubtestDescription(11) = 'Check < 72|mu| 31>'
SubtestDescription(12) = 'Check < 74|mu| 23>'

SubtestDescription(13) = 'Check energy J=1 v=10'
SubtestDescription(14) = 'Check energy J=1 v=99'

SubtestDescription(15) = 'Check energy J=2 v=10'
SubtestDescription(16) = 'Check energy J=2 v=99'

SubtestDescription(17) = 'Check ZPE'
SubtestDescription(18) = 'Check trans 1'
SubtestDescription(19) = 'Check trans 2'
SubtestDescription(20) = 'Check trans 3'
SubtestDescription(21) = 'Check states file, line=50'
SubtestDescription(22) = 'Check states file, line=300'

SubtestDescription(23) = 'Check states file, large A'
SubtestDescription(24) = 'Check states file, medium A'
SubtestDescription(25) = 'Check states file, small A'



do i=1, nlines
   read(u1, '(a500)', iostat=ierr) line
   if(ierr /=0) exit
   line = ADJUSTL(line) ! remove leading spaces if present
   if( trim(line) == 'Reduced mass is      20000.0000000000         atomic mass units (Daltons)')    iflag(1) = 1
   if( trim(line) == '6.00000000        1.23808220E+01')                                             iflag(2) = 1
   if( trim(line) == '10.01146132         1.62497765')                                               iflag(3) = 1
   if( trim(line) == '11          0.523523 [    1  10 ] DWell')                                      iflag(4) = 1
   if( trim(line) == '100          5.182881 [    1  99 ] DWell')                                     iflag(5) = 1
!    if( trim(line) == '< 1,   0|                                 barrier     | 1,   0> =   1.530265649165520E+00') iflag(6) = 1


   thresh = 1e-12_dp  ! maximum allowed relative error in transition dipoles to reference value
   if( line(1:65) == '< 1,   0|                                 barrier     | 1,   0> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(1.530265649165520E+00_dp-abs(value))/1.530265649165520E+00_dp  < thresh) iflag(6) = 1
   endif



   if( line(1:65) == '< 1,  50|                                 barrier     | 1,  50> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(2.597142764335410E-01_dp-abs(value))/ 2.597142764335410E-01_dp < thresh) iflag(7) = 1
   endif


   if( line(1:65) == '< 1,  99|                                 barrier     | 1,  99> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(1.830321915314403E-01_dp-abs(value))/ 1.830321915314403E-01_dp < thresh) iflag(8) = 1
   endif

   
   thresh = 7e-8_dp  ! maximum allowed relative error in transition dipoles to reference value
   if( line(1:65) == '< 1,   2|                                 barrier     | 1,  20> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(2.339516678018351E-08_dp-abs(value)) / 2.339516678018351E-08_dp < thresh) iflag(9) = 1
   endif


   thresh = 2e-11_dp  ! maximum allowed relative error in transition dipoles to reference value
   if( line(1:65) == '< 1,  15|                                 barrier     | 1,  29> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(2.808262465924235E-04_dp-abs(value)) / 2.808262465924235E-04_dp < thresh) iflag(10) = 1
   endif

   thresh = 3e-16_dp  ! maximum allowed value for transition dipole
   if( line(1:65) == '< 1,  72|                                 barrier     | 1,  31> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(value)  < thresh) iflag(11) = 1
   endif

   thresh = 4e-16_dp  ! maximum allowed value for transition dipole
   if( line(1:65) == '< 1,  74|                                 barrier     | 1,  23> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(value)  < thresh) iflag(12) = 1
   endif

!    if( trim(line) == '< 1,  99|                                 barrier     | 1,  99> =   1.830321915314403E-01') iflag(8) = 1
!    if( trim(line) == '< 1,   2|                                 barrier     | 1,  20> =  -2.339516678018351E-08') iflag(9) = 1
!    if( trim(line) == '< 1,  15|                                 barrier     | 1,  29> =   2.808262465924235E-04') iflag(10) = 1
!    if( trim(line) == '< 1,  72|                                 barrier     | 1,  31> =   6.542520951923322E-17') iflag(11) = 1
!    if( trim(line) == '< 1,  74|                                 barrier     | 1,  23> =  -4.904802364064624E-19') iflag(12) = 1

   if( trim(line) == '1.0   11          0.523540   1    10   0     0.0     0.0     0.0   -    ||DWell')  iflag(13) = 1
   if( trim(line) == '1.0  100          5.182900   1    99   0     0.0     0.0     0.0   -    ||DWell')  iflag(14) = 1

   if( trim(line) == '2.0   11          0.523574   1    10   0     0.0     0.0     0.0   +    ||DWell')  iflag(15) = 1
   if( trim(line) == '2.0  100          5.182937   1    99   0     0.0     0.0     0.0   +    ||DWell')  iflag(16) = 1

   if( trim(line) == 'Zero point energy (ZPE) =          -0.599656')  iflag(17) = 1

   if( trim(line) == '2.0  A''    <-  1.0  A"    R       0.0001 <-      0.0000      0.0000     4.68342625E+00' &
         // '   1.12746728E-20   1.06377437E-29 (  1   0  0     0.0     0.0 )<-(  1   0  0     0.0     0.0 )') iflag(18) = 1

   if( trim(line) == '1.0  A"    <-  2.0  A''    P       3.9264 <-      3.1935      0.7329     3.70190714E-04' &
         // '   1.52347558E-11   3.90221073E-25 (  1  75  0     0.0     0.0 )<-(  1  61  0     0.0     0.0 )') iflag(19) = 1

   if( trim(line) == '1.0  A"    <-  0.0  A''    R       1.8323 <-      0.9947      0.8377     2.10118040E-08' &
        // '   1.29103950E-15   2.92323356E-29 (  1  35  0     0.0     0.0 )<-(  1  19  0     0.0     0.0 )')  iflag(20) = 1

enddo
close(u1)


! check the states file
open(unit=u1, file='duo_test_0005.states', status='old', action='read')
nlines = 300
do i=1, nlines
   read(u1, '(a500)', iostat=ierr) line
   if(ierr /=0) exit
   line = ADJUSTL(line) ! remove leading spaces if present
   if( trim(line) == '50     2.565264      1       0 + e DWell       49  0       0       0')    iflag(21) = 1
   if( trim(line) == '300     5.182937      5       2 + e DWell       99  0       0       0')    iflag(22) = 1
enddo
close(u1)


! check the trans file
open(unit=u1, file='duo_test_0005.trans', status='old', action='read')
nlines = 3825
do i=1, nlines
   read(u1, '(a500)', iostat=ierr) line
   if(ierr /=0) exit
   line = ADJUSTL(line) ! remove leading spaces if present
   if( trim(line) == '100         192  2.0985E-10            0.418800')    iflag(23) = 1
   if( trim(line) == '166         252  7.2970E-12            0.732897')    iflag(24) = 1
   if( trim(line) == '102           2  9.2440E-22            0.000017')    iflag(25) = 1
enddo
close(u1)

write(*,'(a)')
do i=1, n_subtests(itest)
  write(line,'(A, I6,A)') 'Input file ' // trim(FileInput) //  " ( "  // trim(TestFileName(itest)) // " )" // &
                                          ' , subtest ', i, " (" // trim(SubtestDescription(i)) // ") : "
  write(*,'(a120)',advance='no') line
  if(iflag(i) == 1) then
     write(*,'(a)') 'Passed'
  else
      write(*,'(a)') 'FAILED! <---'
      n_failed_subtests=n_failed_subtests+1
  endif

enddo
!***************************************************************
!**************** END OF TEST 005 ******************************
!***************************************************************

!***************************************************************
!**************** TEST 006 *************************************
!***************************************************************
6 continue
itest=6
if(z_run_test_file(itest) .eqv. .false. ) goto 7 ! skip test if flag is set
TestFileName(itest) = 'AlO complex case'
write(FileInput,  '(a,i4.4,a)') 'duo_test_', itest, '.inp'
write(FileOutput, '(a,i4.4,a)') 'duo_test_', itest, '.out'

open(unit=u1, file=FileInput, status='unknown', action='write')

write(u1,'(a)') '(AlO molecule: Coupled rovibronic calculations) '
write(u1,'(a)') '  '
write(u1,'(a)') 'masses 26.981539  17.9991610 '
write(u1,'(a)') 'molecule AlO '
write(u1,'(a)') '  '
write(u1,'(a)') 'nstates 3  '
write(u1,'(a)') '  '
write(u1,'(a)') 'jrot 0.5 - 2.5  '
write(u1,'(a)') '  '
write(u1,'(a)') 'grid  '
write(u1,'(a)') '  npoints 201  '
write(u1,'(a)') '  range  0.02 2.70  '
write(u1,'(a)') '  type 0  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') 'DIAGONALIZER '
write(u1,'(a)') ' SYEVR  '
write(u1,'(a)') ' nroots  50  '
write(u1,'(a)') ' ZPE 0.0  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') 'CONTRACTION  '
write(u1,'(a)') '  vib '
write(u1,'(a)') '  vmax  30 30 30  '
write(u1,'(a)') 'END  '
write(u1,'(a)') '  '
write(u1,'(a)') 'CHECKPOINT  '
write(u1,'(a)') ' eigenvectors save  '
write(u1,'(a)') 'END  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') '(All potential functions must appear first followed by the couplings)  '
write(u1,'(a)') '  '
write(u1,'(a)') 'poten 1  '
write(u1,'(a)') 'name "X2Sigma+"  '
write(u1,'(a)') 'symmetry +  '
write(u1,'(a)') 'type  EMO  '
write(u1,'(a)') 'lambda 0  '
write(u1,'(a)') 'mult  2.00  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'V0  0.00000000000000E+00  '
write(u1,'(a)') 'RE  0.16189091308069E+01  '
write(u1,'(a)') 'DE  0.39435485441416E+05  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'PL  0.40000000000000E+01  '
write(u1,'(a)') 'PR  0.40000000000000E+01  '
write(u1,'(a)') 'NL  0.20000000000000E+01  '
write(u1,'(a)') 'NR  0.60000000000000E+01  '
write(u1,'(a)') 'B0  0.19035467302005E+01  '
write(u1,'(a)') 'B1  -0.10019383862110E+00  '
write(u1,'(a)') 'B2  -0.71188810338016E-01  '
write(u1,'(a)') 'B3  -0.23780463075872E+00  '
write(u1,'(a)') 'B4  0.27801043427257E+00  '
write(u1,'(a)') 'B5  0.00000000000000E+00  '
write(u1,'(a)') 'B6  0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') 'poten 2  '
write(u1,'(a)') 'name "A2PI"  '
write(u1,'(a)') 'type  EMO  '
write(u1,'(a)') 'lambda 1  '
write(u1,'(a)') 'mult  2  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'V0  0.54058510643430E+04  '
write(u1,'(a)') 'RE  0.17690556534168E+01  '
write(u1,'(a)') 'DE  0.39435485441416E+05 link  1  1  3  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'PL  0.40000000000000E+01  '
write(u1,'(a)') 'PR  0.40000000000000E+01  '
write(u1,'(a)') 'NL  0.20000000000000E+01  '
write(u1,'(a)') 'NR  0.50000000000000E+01  '
write(u1,'(a)') 'B0  0.15264720726250E+01  '
write(u1,'(a)') 'B1  -0.16794482850527E+00  '
write(u1,'(a)') 'B2  0.58551367092287E-01  '
write(u1,'(a)') 'B3  0.76372699231568E-01  '
write(u1,'(a)') 'B4  -0.14954775136109E+00  '
write(u1,'(a)') 'B5  0.32678665276613E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') 'poten 3  '
write(u1,'(a)') 'name "B2Sigma+"  '
write(u1,'(a)') 'type  EMO  '
write(u1,'(a)') 'lambda 0  '
write(u1,'(a)') 'symmetry +  '
write(u1,'(a)') 'mult  2  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'V0  0.11688310378409E+05  '
write(u1,'(a)') 'RE  0.16680733690045E+01  '
write(u1,'(a)') 'DE  0.70000000000000E+05  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'PL  0.40000000000000E+01  '
write(u1,'(a)') 'PR  0.20000000000000E+01  '
write(u1,'(a)') 'NL  0.20000000000000E+01  '
write(u1,'(a)') 'NR  0.50000000000000E+01  '
write(u1,'(a)') 'B0  0.15100604333583E+01  '
write(u1,'(a)') 'B1  -0.14541739107331E+00  '
write(u1,'(a)') 'B2  0.45276864376232E-01  '
write(u1,'(a)') 'B3  0.12841783832665E+01  '
write(u1,'(a)') 'B4  0.17207107773002E+01  '
write(u1,'(a)') 'B5  0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') '(The couplings ref. numb. must correspond to the ref. numb. used in the PEC above) '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-orbit  2 2  '
write(u1,'(a)') 'name "<Lambda=+1,S=1 (a2Pi)|LSZ|+1 (a2Pi),S=1>" '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'lambda 1 1  '
write(u1,'(a)') 'sigma 0.5 0.5  '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'factor  1.0  (0, 1 or i)  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE    0.17700000000000E+01  '
write(u1,'(a)') 'RREF -0.10000000000000E+01  '
write(u1,'(a)') 'P     0.20000000000000E+01  '
write(u1,'(a)') 'NT    0.30000000000000E+01  '
write(u1,'(a)') 'B0    0.10389912869853E+01  '
write(u1,'(a)') 'B1   -0.39342225873616E-01  '
write(u1,'(a)') 'B2    0.00000000000000E+00  '
write(u1,'(a)') 'B3    0.00000000000000E+00  '
write(u1,'(a)') 'BINF  0.10000000000000E+01  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-orbit  1 2  '
write(u1,'(a)') 'name "<Lamda=0, (A 2Sigma+)|LSX|+1 (a 2Pi)>"  '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'lambda  0 -1 '
write(u1,'(a)') 'sigma  -0.5  0.5  (<  0 0.5 |  1  -0.5 > =  1.0)  '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'factor  1.0  (0, 1 or i)  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.20000000000000E+01  '
write(u1,'(a)') 'NT     0.30000000000000E+01  '
write(u1,'(a)') 'B0     0.88236331519677E+00  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'B2     0.00000000000000E+00  '
write(u1,'(a)') 'B3     0.00000000000000E+00  '
write(u1,'(a)') 'BINF   0.10000000000000E+01  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-orbit  3 2  '
write(u1,'(a)') 'name "<Lamda=0, (A 2Sigma+)|LSX|+1 (a 2Pi)>"  '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'lambda 0 -1  '
write(u1,'(a)') 'sigma -0.5 0.5  (<  0 0.5 |  1  -0.5 > =  1.0) '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'factor  1.0  (0, 1 or i)  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE      0.17700000000000E+01  '
write(u1,'(a)') 'RREF   -0.10000000000000E+01  '
write(u1,'(a)') 'P       0.20000000000000E+01  '
write(u1,'(a)') 'NT      0.30000000000000E+01  '
write(u1,'(a)') 'B0      0.12911045632094E+01  '
write(u1,'(a)') 'B1      0.00000000000000E+00  '
write(u1,'(a)') 'B2      0.00000000000000E+00  '
write(u1,'(a)') 'B3      0.00000000000000E+00  '
write(u1,'(a)') 'BINF    0.10000000000000E+01  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'LxLy  3 2  '
write(u1,'(a)') 'name "<2.1|Lx+iLy|-1>"  '
write(u1,'(a)') 'lambda  0 -1 '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'factor  1.0  (0, 1 or i)  '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE      0.17700000000000E+01  '
write(u1,'(a)') 'RREF   -0.10000000000000E+01  '
write(u1,'(a)') 'P       0.20000000000000E+01  '
write(u1,'(a)') 'NT      0.30000000000000E+01  '
write(u1,'(a)') 'B0      0.11325633215740E+01  '
write(u1,'(a)') 'B1      0.00000000000000E+00  '
write(u1,'(a)') 'B2      0.00000000000000E+00  '
write(u1,'(a)') 'B3      0.00000000000000E+00  '
write(u1,'(a)') 'BINF    0.10000000000000E+01  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'LxLy  1 2  '
write(u1,'(a)') 'name "<1.1|Lx+iLy|-1>"  '
write(u1,'(a)') 'lambda  0 -1 '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'factor  1.0  (0, 1 or i)  '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.20000000000000E+01  '
write(u1,'(a)') 'NT     0.30000000000000E+01  '
write(u1,'(a)') 'B0     0.69977284594431E+00  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'B2     0.00000000000000E+00  '
write(u1,'(a)') 'B3     0.00000000000000E+00  '
write(u1,'(a)') 'BINF   0.10000000000000E+01  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-spin 1 1  '
write(u1,'(a)') 'name "<X 2Sigma+|SS-O|X 2Sigma+>" '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'lambda  0  0 '
write(u1,'(a)') 'factor  1.0  (0, 1 or i)  '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.10000000000000E+01  '
write(u1,'(a)') 'NT     0.20000000000000E+01  '
write(u1,'(a)') 'B0     1.00000000000000E+00  '
write(u1,'(a)') 'B1     0.10000000000000E+00  '
write(u1,'(a)') 'B1     0.01000000000000E+00  '
write(u1,'(a)') 'BINF   0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-spin 2 2  '
write(u1,'(a)') 'name "<A|BROT|A>"  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'lambda 1 1  '
write(u1,'(a)') 'mult  2.0 2.0  '
write(u1,'(a)') 'factor 1  '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.10000000000000E+01  '
write(u1,'(a)') 'NT     0.20000000000000E+01  '
write(u1,'(a)') 'B0     1.50000000000000E+00  '
write(u1,'(a)') 'B1     0.20000000000000E+00  '
write(u1,'(a)') 'B1     0.02000000000000E+00  '
write(u1,'(a)') 'BINF   0.00300000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-spin 3 3  '
write(u1,'(a)') 'name "<B 2Sigma+|SS|B 2Sigma+>"  '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'lambda 0 0  '
write(u1,'(a)') 'mult  2.0 2.0  '
write(u1,'(a)') 'factor 1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.10000000000000E+01  '
write(u1,'(a)') 'NT     0.20000000000000E+01  '
write(u1,'(a)') 'B0     3.00000000000000E+00  '
write(u1,'(a)') 'B1     0.20000000000000E+00  '
write(u1,'(a)') 'B1     0.03000000000000E+00  '
write(u1,'(a)') 'BINF   0.00400000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-spin-o 1 1  '
write(u1,'(a)') 'name "<X 2Sigma+|SS-O|X 2Sigma+>" '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'lambda  0  0 '
write(u1,'(a)') 'factor  1.0  (0, 1 or i)  '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.10000000000000E+01  '
write(u1,'(a)') 'NT     0.20000000000000E+01  '
write(u1,'(a)') 'B0     1.00000000000000E+00  '
write(u1,'(a)') 'B1     0.10000000000000E+00  '
write(u1,'(a)') 'B1     0.01000000000000E+00  '
write(u1,'(a)') 'BINF   0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-spin-o 2 2  '
write(u1,'(a)') 'name "<A|BROT|A>"  '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'lambda 1 1  '
write(u1,'(a)') 'mult  2.0 2.0  '
write(u1,'(a)') 'factor 1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.10000000000000E+01  '
write(u1,'(a)') 'NT     0.20000000000000E+01  '
write(u1,'(a)') 'B0     2.00000000000000E+00  '
write(u1,'(a)') 'B1     0.20000000000000E+00  '
write(u1,'(a)') 'B1     0.02000000000000E+00  '
write(u1,'(a)') 'BINF   0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-spin-o 3 3  '
write(u1,'(a)') 'name "<B 2Sigma+|SS-o|B 2Sigma+>" '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'lambda 0 0  '
write(u1,'(a)') 'mult  2.0 2.0  '
write(u1,'(a)') 'factor 1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.10000000000000E+01  '
write(u1,'(a)') 'NT     0.20000000000000E+01  '
write(u1,'(a)') 'B0     3.00000000000000E+00  '
write(u1,'(a)') 'B1     0.30000000000000E+00  '
write(u1,'(a)') 'B1     0.03000000000000E+00  '
write(u1,'(a)') 'BINF   0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'bobrot 1 1  '
write(u1,'(a)') 'name "<X 2Sigma+|BROT|X 2Sigma+>" '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'lambda 0 0  '
write(u1,'(a)') 'mult  2.0 2.0  '
write(u1,'(a)') 'factor 1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.10000000000000E+01  '
write(u1,'(a)') 'NT     0.20000000000000E+01  '
write(u1,'(a)') 'B0     0.13461788781347E-02  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'BINF   0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'bobrot 2 2  '
write(u1,'(a)') 'name "<A|BROT|A>"  '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'lambda 1 1  '
write(u1,'(a)') 'mult  2.0 2.0  '
write(u1,'(a)') 'factor 1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P       0.10000000000000E+01  '
write(u1,'(a)') 'NT     0.20000000000000E+01  '
write(u1,'(a)') 'B0     0.11835367656259E-02  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'BINF   0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') 'bobrot 3 3  '
write(u1,'(a)') 'name "<B 2Sigma+|Brot|B 2Sigma+>" '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'lambda 0 0  '
write(u1,'(a)') 'mult  2.0 2.0  '
write(u1,'(a)') 'factor 1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.40000000000000E+01  '
write(u1,'(a)') 'NT     0.20000000000000E+01  '
write(u1,'(a)') 'B0     0.10802616405392E-02  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'BINF   0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-rot 1 1 '
write(u1,'(a)') 'name "<X 2Sigma+|Brot|X 2Sigma+>" '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'lambda 0 0  '
write(u1,'(a)') 'mult  2.0 2.0  '
write(u1,'(a)') 'factor 1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.40000000000000E+01  '
write(u1,'(a)') 'NT     0.20000000000000E+01  '
write(u1,'(a)') 'B0     0.67490157590221E-02  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'BINF   0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-rot 2 2 '
write(u1,'(a)') 'name "<A 2Pi|Brot|X 2Pi>"  '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'lambda 1 1  '
write(u1,'(a)') 'mult  2.0 2.0  '
write(u1,'(a)') 'factor 1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE    0.17700000000000E+01  '
write(u1,'(a)') 'RREF -0.10000000000000E+01  '
write(u1,'(a)') 'P     0.40000000000000E+01  '
write(u1,'(a)') 'NT    0.20000000000000E+01  '
write(u1,'(a)') 'B0   -0.20607819057186E-02  '
write(u1,'(a)') 'B1    0.00000000000000E+00  '
write(u1,'(a)') 'B1    0.00000000000000E+00  '
write(u1,'(a)') 'BINF  0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'spin-rot 3 3 '
write(u1,'(a)') 'name "<A 2Sigma+|Brot|A 2Sigma+>" '
write(u1,'(a)') 'type  BOBLEROY  '
write(u1,'(a)') 'units  cm-1  '
write(u1,'(a)') 'lambda 0 0  '
write(u1,'(a)') 'mult  2.0 2.0  '
write(u1,'(a)') 'factor 1  '
write(u1,'(a)') 'values  '
write(u1,'(a)') 'RE     0.17700000000000E+01  '
write(u1,'(a)') 'RREF  -0.10000000000000E+01  '
write(u1,'(a)') 'P      0.40000000000000E+01  '
write(u1,'(a)') 'NT     0.20000000000000E+01  '
write(u1,'(a)') 'B0    -0.16889491142793E-02  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'B1     0.00000000000000E+00  '
write(u1,'(a)') 'BINF   0.00000000000000E+00  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'dipole  1 1  '
write(u1,'(a)') 'name "<X 2 Sigma|DMZ|X 2Sigma>"  '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'lambda  0  0 '
write(u1,'(a)') 'type  grid  '
write(u1,'(a)') 'factor  1  (0, 1 or i)  '
write(u1,'(a)') 'values  '
write(u1,'(a)') '1.24  -3.74297513037805'
write(u1,'(a)') '1.28  -3.90022088914213'
write(u1,'(a)') '1.32  -4.02418555419016'
write(u1,'(a)') '1.36  -4.12350836856094'
write(u1,'(a)') '1.40  -4.20338358814335'
write(u1,'(a)') '1.44  -4.26703937696391'
write(u1,'(a)') '1.48  -4.31663015439479'
write(u1,'(a)') '1.52  -4.35378831191569'
write(u1,'(a)') '1.56  -4.37998902183456'
write(u1,'(a)') '1.60  -4.39678851842240'
write(u1,'(a)') '1.64  -4.40599303222464'
write(u1,'(a)') '1.68  -4.40980500098331'
write(u1,'(a)') '1.72  -4.41078909529458'
write(u1,'(a)') '1.76  -4.41187262963450'
write(u1,'(a)') '1.80  -4.41626982265248'
write(u1,'(a)') '1.84  -4.42726683963877'
write(u1,'(a)') '1.88  -4.44798871706648'
write(u1,'(a)') '1.92  -4.48128794843708'
write(u1,'(a)') '1.96  -4.52934774077280'
write(u1,'(a)') '2.00  -4.59367819738656'
write(u1,'(a)') '2.04  -4.67499527438795'
write(u1,'(a)') '2.08  -4.77340428527536'
write(u1,'(a)') '2.12  -4.88824474318141'
write(u1,'(a)') '2.16  -5.01846067351136'
write(u1,'(a)') '2.20  -5.16260559163532'
write(u1,'(a)') '2.24  -5.31910261703098'
write(u1,'(a)') '2.28  -5.48628753119164'
write(u1,'(a)') '2.32  -5.66251306624174'
write(u1,'(a)') '2.36  -5.84619637635847'
write(u1,'(a)') '2.40  -6.03585222025867'
write(u1,'(a)') '2.44  -6.23010362609406'
write(u1,'(a)') '2.48  -6.42763000440126'
write(u1,'(a)') '2.52  -6.62723621627096'
write(u1,'(a)') '2.56  -6.82765014627356'
write(u1,'(a)') '2.60  -7.02754356357453'
write(u1,'(a)') '2.64  -7.22541110033244'
write(u1,'(a)') '2.68  -7.41905358164998'
write(u1,'(a)') '2.72  -7.60513394499626'
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') 'dipole  3 3  '
write(u1,'(a)') 'name "<B 2Sigma|DMZ|B 2Sigma>"  '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'lambda  0  0 '
write(u1,'(a)') 'type  grid  '
write(u1,'(a)') 'factor  1  (0, 1 or i)  '
write(u1,'(a)') 'values  '
write(u1,'(a)') '1.24   0.39932900084131  '
write(u1,'(a)') '1.28   0.33554205793302  '
write(u1,'(a)') '1.32   0.223060158659815  '
write(u1,'(a)') '1.36   0.0729790578004454  '
write(u1,'(a)') '1.40  -0.107223348932292  '
write(u1,'(a)') '1.44  -0.312148681799442  '
write(u1,'(a)') '1.48  -0.53758230628379  '
write(u1,'(a)') '1.52  -0.779948587178871  '
write(u1,'(a)') '1.56  -1.03599334756638  '
write(u1,'(a)') '1.60  -1.30251334886105  '
write(u1,'(a)') '1.64  -1.57630524511931  '
write(u1,'(a)') '1.68  -1.85400838346955  '
write(u1,'(a)') '1.72  -2.1322068411598  '
write(u1,'(a)') '1.76  -2.40734123563592  '
write(u1,'(a)') '1.80  -2.67596972088375  '
write(u1,'(a)') '1.84  -2.93478609167239  '
write(u1,'(a)') '1.88  -3.18108361663149  '
write(u1,'(a)') '1.92  -3.41212317596315  '
write(u1,'(a)') '1.96  -3.62637323061642  '
write(u1,'(a)') '2.00  -3.82297020156609  '
write(u1,'(a)') '2.04  -4.00212217679762  '
write(u1,'(a)') '2.08  -4.16386541769391  '
write(u1,'(a)') '2.12  -4.30925300437009  '
write(u1,'(a)') '2.16  -4.43980032092697  '
write(u1,'(a)') '2.20  -4.55656984421131  '
write(u1,'(a)') '2.24  -4.66088770504403  '
write(u1,'(a)') '2.28  -4.75359029957750  '
write(u1,'(a)') '2.32  -4.83497713023923  '
write(u1,'(a)') '2.36  -4.90463485048954  '
write(u1,'(a)') '2.40  -4.96122299662847  '
write(u1,'(a)') '2.44  -5.00212226987129  '
write(u1,'(a)') '2.48  -5.02297083420418  '
write(u1,'(a)') '2.52  -5.01793668502258  '
write(u1,'(a)') '2.56  -4.97775692180470  '
write(u1,'(a)') '2.60  -4.88980704693652  '
write(u1,'(a)') '2.64  -4.73788132553041  '
write(u1,'(a)') '2.68  -4.50393753889427  '
write(u1,'(a)') '2.72  -4.17282164080206  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') 'dipole  3 1  '
write(u1,'(a)') 'name "<B 2Sigma|DMZ|X 2Sigma>"  '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'lambda  0 0  '
write(u1,'(a)') 'type  grid  '
write(u1,'(a)') 'factor  1  (0, 1 or i)  '
write(u1,'(a)') 'values  '
write(u1,'(a)') '1.24  2.46101727976436  '
write(u1,'(a)') '1.28  2.38249299414292  '
write(u1,'(a)') '1.32  2.32645492956901  '
write(u1,'(a)') '1.36  2.28584259259417  '
write(u1,'(a)') '1.4  2.25493150165221  '
write(u1,'(a)') '1.44  2.22915193482794  '
write(u1,'(a)') '1.48  2.20482812024951  '
write(u1,'(a)') '1.52  2.1789459941707  '
write(u1,'(a)') '1.56  2.14900519945455  '
write(u1,'(a)') '1.6  2.11285278444499  '
write(u1,'(a)') '1.64  2.06866640827698  '
write(u1,'(a)') '1.68  2.01494411145352  '
write(u1,'(a)') '1.72  1.95049022434868  '
write(u1,'(a)') '1.76  1.87452601250416  '
write(u1,'(a)') '1.8  1.78672171699014  '
write(u1,'(a)') '1.84  1.68729073838483  '
write(u1,'(a)') '1.88  1.5769733627308  '
write(u1,'(a)') '1.92  1.45707209210245  '
write(u1,'(a)') '1.96  1.32919484351265  '
write(u1,'(a)') '2  1.19523040666005 '
write(u1,'(a)') '2.04  1.057155421631  '
write(u1,'(a)') '2.08  0.916717362432375  '
write(u1,'(a)') '2.12  0.775337479478925  '
write(u1,'(a)') '2.16  0.634221913831869  '
write(u1,'(a)') '2.2  0.493818952175943  '
write(u1,'(a)') '2.24  0.354336630169673  '
write(u1,'(a)') '2.28  0.215403463007744  '
write(u1,'(a)') '2.32  0.0761523820812309  '
write(u1,'(a)') '2.36  -0.0647894112438063  '
write(u1,'(a)') '2.4  -0.209305141217806  '
write(u1,'(a)') '2.44  -0.359844327684191  '
write(u1,'(a)') '2.48  -0.519499731519158  '
write(u1,'(a)') '2.52  -0.691766668516972  '
write(u1,'(a)') '2.56  -0.881090665910492  '
write(u1,'(a)') '2.6  -1.09256798137821  '
write(u1,'(a)') '2.64  -1.33180776277106  '
write(u1,'(a)') '2.68  -1.60329727638649  '
write(u1,'(a)') '2.72  -1.9103760369686  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') 'dipole  2 2  '
write(u1,'(a)') 'name "<A 2Pi|DMZ|A 2Pi>"  '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'lambda  1  1 '
write(u1,'(a)') 'type  grid  '
write(u1,'(a)') 'factor  1  (0, 1 or i)  '
write(u1,'(a)') 'values  '
write(u1,'(a)') '1.24  0.498966372178846  '
write(u1,'(a)') '1.28  0.41788664128511  '
write(u1,'(a)') '1.32  0.320580986372765  '
write(u1,'(a)') '1.36  0.210290783610119  '
write(u1,'(a)') '1.4  0.0885674157125876  '
write(u1,'(a)') '1.44  -0.0434546546130722  '
write(u1,'(a)') '1.48  -0.184971621860844  '
write(u1,'(a)') '1.52  -0.335310755004951  '
write(u1,'(a)') '1.56  -0.493840097992247  '
write(u1,'(a)') '1.6  -0.659972195728683  '
write(u1,'(a)') '1.64  -0.833122879542169  '
write(u1,'(a)') '1.68  -1.01271396562238  '
write(u1,'(a)') '1.72  -1.19818471828099  '
write(u1,'(a)') '1.76  -1.38892640282539  '
write(u1,'(a)') '1.8  -1.58441179955666  '
write(u1,'(a)') '1.84  -1.78407446599176  '
write(u1,'(a)') '1.88  -1.98736845904451  '
write(u1,'(a)') '1.92  -2.1937605708703  '
write(u1,'(a)') '1.96  -2.40273035698289  '
write(u1,'(a)') '2  -2.61384465469763  '
write(u1,'(a)') '2.04  -2.82648039224136  '
write(u1,'(a)') '2.08  -3.04025049234958  '
write(u1,'(a)') '2.12  -3.25472164161974  '
write(u1,'(a)') '2.16  -3.46948756515089  '
write(u1,'(a)') '2.2  -3.68416735238955  '
write(u1,'(a)') '2.24  -3.89841361151925  '
write(u1,'(a)') '2.28  -4.11190089210641  '
write(u1,'(a)') '2.32  -4.32433149203729  '
write(u1,'(a)') '2.36  -4.53543932469594  '
write(u1,'(a)') '2.4  -4.74497386967699  '
write(u1,'(a)') '2.44  -4.95270876247975  '
write(u1,'(a)') '2.48  -5.15809114271666  '
write(u1,'(a)') '2.52  -5.36168171610632  '
write(u1,'(a)') '2.56  -5.5627459057469  '
write(u1,'(a)') '2.6  -5.76035722859042  '
write(u1,'(a)') '2.64  -5.95573315410853  '
write(u1,'(a)') '2.68  -6.14767502546518  '
write(u1,'(a)') '2.72  -6.33580183472748  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') 'dipole  1 2  '
write(u1,'(a)') 'name "<X 2Sigma|DMY|A 2Pi''>"  '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'lambda  0  -1  '
write(u1,'(a)') 'type  grid  '
write(u1,'(a)') 'factor  1  (0, 1 or i)  '
write(u1,'(a)') 'values  '
write(u1,'(a)') '1.24  0.466652497254106  '
write(u1,'(a)') '1.28  0.519235784296098  '
write(u1,'(a)') '1.32  0.557000887026224  '
write(u1,'(a)') '1.36  0.585099776400804  '
write(u1,'(a)') '1.4  0.606206260668169  '
write(u1,'(a)') '1.44  0.621726139770924  '
write(u1,'(a)') '1.48  0.632560271972133  '
write(u1,'(a)') '1.52  0.639301409418364  '
write(u1,'(a)') '1.56  0.642346613246548  '
write(u1,'(a)') '1.6  0.641993225054214  '
write(u1,'(a)') '1.64  0.63845711960365  '
write(u1,'(a)') '1.68  0.63192303971895  '
write(u1,'(a)') '1.72  0.622555396163213  '
write(u1,'(a)') '1.76  0.610515638654141  '
write(u1,'(a)') '1.8  0.595989936436187  '
write(u1,'(a)') '1.84  0.579183817877558  '
write(u1,'(a)') '1.88  0.560333458911694  '
write(u1,'(a)') '1.92  0.539706651521383  '
write(u1,'(a)') '1.96  0.517589078944545  '
write(u1,'(a)') '2  0.494272094988858  '
write(u1,'(a)') '2.04  0.470062085905852  '
write(u1,'(a)') '2.08  0.445245842847587  '
write(u1,'(a)') '2.12  0.420078026929938  '
write(u1,'(a)') '2.16  0.394786739487255  '
write(u1,'(a)') '2.2  0.369568299498086  '
write(u1,'(a)') '2.24  0.344583275827807  '
write(u1,'(a)') '2.28  0.319961642123823  '
write(u1,'(a)') '2.32  0.295805967012094  '
write(u1,'(a)') '2.36  0.272194617243861  '
write(u1,'(a)') '2.4  0.249187394655897  '
write(u1,'(a)') '2.44  0.226831264049009  '
write(u1,'(a)') '2.48  0.205131424328731  '
write(u1,'(a)') '2.52  0.184189814660109  '
write(u1,'(a)') '2.56  0.164004001422798  '
write(u1,'(a)') '2.6  0.144544100304911  '
write(u1,'(a)') '2.64  0.126010926029143  '
write(u1,'(a)') '2.68  0.108424797708835  '
write(u1,'(a)') '2.72  0.0919812823367467  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') 'dipole  3 2  '
write(u1,'(a)') 'name "<B 2Sigma|DMY|A 2Pi''>"  '
write(u1,'(a)') 'spin  0.5 0.5  '
write(u1,'(a)') 'lambda  0  -1  '
write(u1,'(a)') 'type  grid  '
write(u1,'(a)') 'factor  1  (0, 1 or i)  '
write(u1,'(a)') 'values  '
write(u1,'(a)') '1.24  0.0621788893197797  '
write(u1,'(a)') '1.28  0.0604543543515064  '
write(u1,'(a)') '1.32  0.0579526573686149  '
write(u1,'(a)') '1.36  0.0541174338813709  '
write(u1,'(a)') '1.4  0.048864600954598  '
write(u1,'(a)') '1.44  0.0423003562076878  '
write(u1,'(a)') '1.48  0.0345296870533159  '
write(u1,'(a)') '1.52  0.0256455903662754  '
write(u1,'(a)') '1.56  0.0157235837705812  '
write(u1,'(a)') '1.6  0.00481805857422163  '
write(u1,'(a)') '1.64  -0.00694956555545681 '
write(u1,'(a)') '1.68  -0.0194743415204157  '
write(u1,'(a)') '1.72  -0.0325836972112064  '
write(u1,'(a)') '1.76  -0.0460646859139908  '
write(u1,'(a)') '1.8  -0.059659723687733  '
write(u1,'(a)') '1.84  -0.0730879475078161  '
write(u1,'(a)') '1.88  -0.0860710052898898  '
write(u1,'(a)') '1.92  -0.098345947481875  '
write(u1,'(a)') '1.96  -0.109698138532068  '
write(u1,'(a)') '2  -0.119986297036847  '
write(u1,'(a)') '2.04  -0.129126316333835  '
write(u1,'(a)') '2.08  -0.137184231755181  '
write(u1,'(a)') '2.12  -0.144231452349704  '
write(u1,'(a)') '2.16  -0.150420370022485  '
write(u1,'(a)') '2.2  -0.155957213982501  '
write(u1,'(a)') '2.24  -0.161083407707581  '
write(u1,'(a)') '2.28  -0.166061555837385  '
write(u1,'(a)') '2.32  -0.171177871792678  '
write(u1,'(a)') '2.36  -0.176705323517438  '
write(u1,'(a)') '2.4  -0.182970280838075  '
write(u1,'(a)') '2.44  -0.190259326877428  '
write(u1,'(a)') '2.48  -0.198880065431766  '
write(u1,'(a)') '2.52  -0.209085359006163  '
write(u1,'(a)') '2.56  -0.221219684359497  '
write(u1,'(a)') '2.6  -0.235492213412527  '
write(u1,'(a)') '2.64  -0.252053406765127  '
write(u1,'(a)') '2.68  -0.270683679241823  '
write(u1,'(a)') '2.72  -0.290911023423012  '
write(u1,'(a)') 'end  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') '  '
write(u1,'(a)') 'INTENSITY'
write(u1,'(a)') '  absorption'
write(u1,'(a)') '  THRESH_INTES  1e-20'
write(u1,'(a)') '  THRESH_LINE  1e-20'
write(u1,'(a)') '  TEMPERATURE  1700.0'
write(u1,'(a)') '  QSTAT  42075.33950 (partition function)  '
write(u1,'(a)') '  GNS  6.0 6.0  '
write(u1,'(a)') '  ZPE  0  (ZPE=0 means recompute)'
write(u1,'(a)') '  selection (rules) 1 1  '
write(u1,'(a)') '  linelist  duo_test_0006'
write(u1,'(a)') '  J,  0.5, 2.5  '
write(u1,'(a)') '  freq-window  0.0,  35000.0  '
write(u1,'(a)') '  energy low  -0.001, 15000.00, upper  -0.00, 50000.0  '
write(u1,'(a)') 'END  '
close(u1)

 cmd = "./" // trim(executable) // " < "  // trim(FileInput) // " > " // trim(FileOutput)
call execute_command_line(cmd)

open(unit=u1, file=FileOutput, status='old', action='read')

! The output should have 9780 lines but we try to read more than that, to allow for
! future further output
nlines = 10000
iflag = 0 !0 means the test has failed

n_subtests(itest)=52
SubtestDescription( 1) = 'Check mass'
SubtestDescription( 2) = 'Check pots, r=0.02'
SubtestDescription( 3) = 'Check pots, r=2.00'
SubtestDescription( 4) = 'Check spin-orbit, r=0.02'
SubtestDescription( 5) = 'Check spin-orbit, r=2.00'
SubtestDescription( 6) = 'Check L+, r=0.02'
SubtestDescription( 7) = 'Check L+, r=2.00'
SubtestDescription( 8) = 'Check spin-spin, r=0.02'
SubtestDescription( 9) = 'Check spin-spin, r=2.00'
SubtestDescription(10) = 'Check spin-spin-o, r=0.02'
SubtestDescription(11) = 'Check spin-spin-o, r=2.00'
SubtestDescription(12) = 'Check spin-rot, r=0.02'
SubtestDescription(13) = 'Check spin-rot, r=2.00'
SubtestDescription(14) = 'Check bob-rot, r=0.02'
SubtestDescription(15) = 'Check bob-rot, r=2.00'
SubtestDescription(16) = 'Check dipoles, r=0.02'
SubtestDescription(17) = 'Check dipoles, r=2.00'

SubtestDescription(18) = 'energy, J=0.0, n=50'
SubtestDescription(19) = 'energy, J=0.0, n=51'
SubtestDescription(20) = 'energy, J=0.0, n=52'
SubtestDescription(21) = 'energy, J=0.0, n=53'
SubtestDescription(22) = 'energy, J=0.0, n=54'
SubtestDescription(23) = 'energy, J=0.0, n=55'
SubtestDescription(24) = 'energy, J=0.0, n=56'

SubtestDescription(25) = 'energy, J=0.5, n=47, +'
SubtestDescription(26) = 'energy, J=0.5, n=48, +'
SubtestDescription(27) = 'energy, J=0.5, n=49, +'
SubtestDescription(28) = 'energy, J=0.5, n=50, +'

SubtestDescription(29) = 'energy, J=0.5, n=47, -'
SubtestDescription(30) = 'energy, J=0.5, n=48, -'
SubtestDescription(31) = 'energy, J=0.5, n=49, -'
SubtestDescription(32) = 'energy, J=0.5, n=50, -'
                                         
SubtestDescription(33) = 'Check ZPE'

SubtestDescription(34) = 'Check energy, J=2.5, n=46, +'
SubtestDescription(35) = 'Check energy, J=2.5, n=47, +'
SubtestDescription(36) = 'Check energy, J=2.5, n=48, +'
SubtestDescription(37) = 'Check energy, J=2.5, n=49, +'
SubtestDescription(38) = 'Check energy, J=2.5, n=50, +'

SubtestDescription(39) = 'Check energy, J=2.5, n=46, -'
SubtestDescription(40) = 'Check energy, J=2.5, n=47, -'
SubtestDescription(41) = 'Check energy, J=2.5, n=48, -'
SubtestDescription(42) = 'Check energy, J=2.5, n=49, -'
SubtestDescription(43) = 'Check energy, J=2.5, n=50, -'

SubtestDescription(44) = 'Check <n=1,v=29|mu|n=1,v=29> (strong)'
SubtestDescription(45) = 'Check <n=2,v=13|mu|n=2,v=14> (strong)'
SubtestDescription(46) = 'Check <n=1,v=13|mu|n=2,v= 8> (strong)'

SubtestDescription(47) = 'Check <n=1,v=20|mu|n=1,v=11> (medium)'
SubtestDescription(48) = 'Check <n=2,v=22|mu|n=2,v=14> (medium)'
SubtestDescription(49) = 'Check <n=3,v=26|mu|n=1,v= 7> (medium)'

SubtestDescription(50) = 'Check <n=3,v=26|mu|n=3,v= 0> (weak)'
SubtestDescription(51) = 'Check <n=1,v=29|mu|n=2,v= 0> (weak)'
SubtestDescription(52) = 'Check <n=3,v=25|mu|n=2,v= 0> (weak)'



! TO BE CONTINUED!!! Add tests of expectation values and line intensities

do i=1, nlines
   read(u1, '(a500)', iostat=ierr) line
   if(ierr /=0) exit
   line = ADJUSTL(line) ! remove leading spaces if present
   if( trim(line) == 'Reduced mass is      10.7967431473672         atomic mass units (Daltons)')     iflag( 1) = 1
   if( trim(line) == '0.02000000        1.73521828E+07        1.42460732E+07        1.40093879E+07')  iflag( 2) = 1
   if( trim(line) == '2.00320000        1.01610256E+04        8.35726996E+03        2.07413348E+04')  iflag( 3) = 1

                       
   if( trim(line) == '0.02000000         1.15662694         0.76475667         1.58213480')           iflag( 4) = 1
   if( trim(line) == '2.00320000         1.02994198         0.89684889         1.25525846')           iflag( 5) = 1

   if( trim(line) == '0.02000000         1.26509280         0.39962235')                              iflag( 6) = 1
   if( trim(line) == '2.00320000         1.11623971         0.73674229')                              iflag( 7) = 1

   if( trim(line) == '0.02000000         1.80321015         2.61466052         5.59906586')           iflag( 8) = 1
   if( trim(line) == '2.00320000         0.94402999         1.41914754         2.82653872')           iflag( 9) = 1

   if( trim(line) == '0.02000000         1.80321015         3.60642030         5.40963045')           iflag(10) = 1 
   if( trim(line) == '2.03000000         0.93799652         1.87599304         2.81398956')           iflag(11) = 1

   if( trim(line) == '0.02000000         0.01349803        -0.00412156        -0.00337790')           iflag(12) = 1
   if( trim(line) == '2.00320000         0.00511172        -0.00156084        -0.00127921')           iflag(13) = 1

   if( trim(line) == '0.02000000         0.00266228         0.00234063         0.00216052')           iflag(14) = 1
   if( trim(line) == '2.00320000         0.00126298         0.00111039         0.00081819')           iflag(15) = 1
   if( trim(line) == '0.02000000        -0.04297383         0.04297813         0.11494874         0.05435801' &
                                                     // '        -0.01035873         0.00278061')     iflag(16) = 1
   if( trim(line) == '2.00320000        -4.59955297        -3.83794226         1.18431009        -2.63080649' &
                                                     // '         0.49236413        -0.12075964')     iflag(17) = 1

   if( trim(line) == '50      17560.965168 [    1  22 ] X2Sigma+') iflag(18) = 1
   if( trim(line) == '51      17708.162811 [    2  20 ] A2PI')     iflag(19) = 1
   if( trim(line) == '52      17872.677641 [    3   7 ] B2Sigma+') iflag(20) = 1
   if( trim(line) == '53      18219.829807 [    1  23 ] X2Sigma+') iflag(21) = 1
   if( trim(line) == '54      18253.192321 [    2  21 ] A2PI')     iflag(22) = 1
   if( trim(line) == '55      18731.068664 [    3   8 ] B2Sigma+') iflag(23) = 1
   if( trim(line) == '56      18791.474893 [    2  22 ] A2PI')     iflag(24) = 1

   if( trim(line) == '0.5   47      16890.650882   1    21   0     0.5     0.5     0.5   +    ||X2Sigma+') iflag(25) = 1
   if( trim(line) == '0.5   48      17007.544019   3     6   0     0.5     0.5     0.5   +    ||B2Sigma+') iflag(26) = 1
   if( trim(line) == '0.5   49      17155.306040   2    19   1     0.5    -0.5     0.5   +    ||A2PI')     iflag(27) = 1
   if( trim(line) == '0.5   50      17560.965349   1    22   0     0.5     0.5     0.5   +    ||X2Sigma+') iflag(28) = 1

   if( trim(line) == '0.5   47      16891.609207   1    21   0     0.5    -0.5    -0.5   -    ||X2Sigma+')  iflag(29) = 1
   if( trim(line) == '0.5   48      17008.622805   3     6   0     0.5    -0.5    -0.5   -    ||B2Sigma+')  iflag(30) = 1
   if( trim(line) == '0.5   49      17155.305429   2    19  -1     0.5     0.5    -0.5   -    ||A2PI')      iflag(31) = 1
   if( trim(line) == '0.5   50      17561.913088   1    22   0     0.5    -0.5    -0.5   -    ||X2Sigma+')  iflag(32) = 1

   if( trim(line) == 'Zero point energy (ZPE) =         470.714781') iflag(33) = 1

   if( trim(line) == '2.5   46      13481.923484   3     2   0     0.5     0.5     0.5   +    ||B2Sigma+') iflag(34) = 1
   if( trim(line) == '2.5   47      13701.037026   2    13   1     0.5    -0.5     0.5   +    ||A2PI')     iflag(35) = 1
   if( trim(line) == '2.5   48      13703.803630   2    13   1     0.5     0.5     1.5   +    ||A2PI')     iflag(36) = 1
   if( trim(line) == '2.5   49      14095.950398   1    17   0     0.5     0.5     0.5   +    ||X2Sigma+') iflag(37) = 1
   if( trim(line) == '2.5   50      14294.905239   2    14   1     0.5    -0.5     0.5   +    ||A2PI')     iflag(38) = 1

   if( trim(line) == '2.5   46      13485.247141   3     2   0     0.5    -0.5    -0.5   -    ||B2Sigma+') iflag(39) = 1
   if( trim(line) == '2.5   47      13701.036443   2    13  -1     0.5     0.5    -0.5   -    ||A2PI')     iflag(40) = 1
   if( trim(line) == '2.5   48      13703.802283   2    13  -1     0.5    -0.5    -1.5   -    ||A2PI')     iflag(41) = 1
   if( trim(line) == '2.5   49      14098.952290   1    17   0     0.5    -0.5    -0.5   -    ||X2Sigma+') iflag(42) = 1
   if( trim(line) == '2.5   50      14294.904583   2    14  -1     0.5     0.5    -0.5   -    ||A2PI')     iflag(43) = 1


   thresh = 1.e-12_dp ; ref=5.004468733531663_dp
   if( line(1:65) == '< 1,  29|                <X 2 Sigma|DMZ|X 2Sigma>     | 1,  29> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(ref-abs(value))/ref < thresh)                                                                iflag(44) = 1
   endif

   thresh = 1.e-12_dp ; ref=9.216621586690346E-01_dp
   if( line(1:65) == '< 2,  13|                       <A 2Pi|DMZ|A 2Pi>     | 2,  14> =') then
     read(line(67:200),*,iostat=ierr) value
     write(*,*) 
     if( abs(ref-abs(value))/ref < thresh)                                                                iflag(45) = 1
   endif

   thresh = 2.e-12_dp ; ref=1.078494426033815E-01_dp
   if( line(1:65) == '< 1,  13|                   <X 2Sigma|DMY|A 2Pi''>     | 2,   8> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(ref-abs(value))/ref < thresh)                                                                iflag(46) = 1
   endif


   thresh = 2.e-11_dp ; ref=9.220145314074652E-05_dp
   if( line(1:65) == '< 1,  20|                <X 2 Sigma|DMZ|X 2Sigma>     | 1,  11> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(ref-abs(value))/ref < thresh)                                                                iflag(47) = 1
   endif


   thresh = 6.e-11_dp ; ref=9.165906260969210E-05_dp
   if( line(1:65) == '< 2,  22|                       <A 2Pi|DMZ|A 2Pi>     | 2,  14> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(ref-abs(value))/ref < thresh)                                                                iflag(48) = 1
   endif


   thresh = 5.e-11_dp ; ref=9.160887351997415E-05_dp
   if( line(1:65) == '< 3,  26|                 <B 2Sigma|DMZ|X 2Sigma>     | 1,   7> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(ref-abs(value))/ref < thresh)                                                                iflag(49) = 1
   endif


   thresh = 1.e-5_dp ; ref=1.257597592525911E-09_dp
   if( line(1:65) == '< 3,  26|                 <B 2Sigma|DMZ|B 2Sigma>     | 3,   0> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(ref-abs(value))/ref < thresh)                                                                iflag(50) = 1
   endif


   thresh = 1.e-5_dp ; ref=1.041117668951722E-09_dp
   if( line(1:65) == '< 1,  29|                   <X 2Sigma|DMY|A 2Pi''>     | 2,   0> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(ref-abs(value))/ref < thresh)                                                                iflag(51) = 1
   endif

   thresh = 1.e-5_dp ; ref=8.648291917628661E-10_dp
   if( line(1:65) == '< 3,  25|                   <B 2Sigma|DMY|A 2Pi''>     | 2,   0> =') then
     read(line(67:200),*,iostat=ierr) value
     if( abs(ref-abs(value))/ref < thresh)                                                                iflag(52) = 1
   endif


enddo
close(u1)

write(*,'(a)')
do i=1, n_subtests(itest)
  write(line,'(A, I6,A)') 'Input file ' // trim(FileInput) //  " ( "  // trim(TestFileName(itest)) // " )" // &
                                          ' , subtest ', i, " (" // trim(SubtestDescription(i)) // ") : "
  write(*,'(a120)',advance='no') line
  if(iflag(i) == 1) then
     write(*,'(a)') 'Passed'
  else
      write(*,'(a)') 'FAILED! <---'
      n_failed_subtests=n_failed_subtests+1
  endif

enddo

!***************************************************************
!**************** END OF TEST 006 ******************************
!***************************************************************
7 continue

!TO DO :
!-tests of emission spectra
!-tests of fitting (very important)
!-tests for all functional form
tot_n_subtests=0
do itest=1, n_test_files 
  tot_n_subtests=tot_n_subtests+n_subtests(itest)
enddo

write(*,*)
write(*,'(A, I6)') 'Total number of subtests      =', tot_n_subtests
write(*,'(A, I6, A, F6.1,A)') 'Number of failed subtests     =', n_failed_subtests, &
                              ' (', 100._dp*real(n_failed_subtests, dp)/real(tot_n_subtests,dp), '% )'
write(*,'(A, I6, A, F6.1,A)') 'Number of successful subtests =', tot_n_subtests-n_failed_subtests, &
                              ' (', 100._dp*real(tot_n_subtests-n_failed_subtests, dp)/real(tot_n_subtests,dp), '% )'
write(*,*)



end program test

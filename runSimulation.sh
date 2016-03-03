scp ./Fortran/ABP-DI_hysteresisTest.f90 marco@bubble.ucd.ie:
ssh -X marco@bubble.ucd.ie 'ifort -O3 -o /home/marco/vicsek /home/marco/ABP-DI_hysteresisTest.f90;nice -19 nohup /home/marco/vicsek > test.out '

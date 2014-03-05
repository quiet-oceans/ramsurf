      program ramdec
      real tl(1000,1000)
      CHARACTER(len=256) :: arg
      threshold = .3
      eps = 1e-10
      CALL get_command_argument(1, arg)
      open(unit=2,status='unknown',file=arg, form='unformatted')
      read(2)nry
      nx=1
    5 read(2,end=6)(tl(j,nx),j=1,nry)
      write(*,*)(tl(j,nx) ,j=1,nry)
      nx=nx+1
      go to 5
    6 end


      program ramdec
      real tl(1000,1000)
      real rtl(1000,1000)
      threshold = 3.5
      open(unit=2,status='unknown',file='ref.grid', form='unformatted')
      open(unit=3,status='unknown',file='tl.grid', form='unformatted')
      read(2)nry
      read(3)ny
      write(*,*) nry .ne. ny
      nx=1
    5 read(2,end=6)(rtl(j,nx),j=1,nry)
      read(3,end=6)(tl(j,nx),j=1,nry)
      write(*,*)(
     > abs(100. * (rtl(j,nx) - tl(j,nx)) / (rtl(j,nx))) .gt. threshold
     > ,j=1,nry)
      nx=nx+1
      go to 5
    6 end


      program ramsurfclr
c
c     Version 2.0
c
c     This code creates color TL plots from the output of RAM. The 
c     Postscript plot file is ramclr.ps. The subroutines emulate 
c     subroutines from the DISSPLA graphics package. The input flags 
c     provide options for color or grayscale, overlayed contours, 
c     a color bar, and the number of colors. In many cases, nice 
c     results can be obtained by using approximately ten colors.
c     This version contains improvements developed by Ed McDonald that
c     allows the output to be displayed with xv and reduces the
c     size of the output file.   
c
      real tl(1000,1000),xb(100),yb(100),vhue(1000),vten(1000),rb(1000),
     >   zb(1000),xc(100),yc(100),rsrf(1000),zsrf(1000)
      open(unit=1,status='old',file='ramsurf.in')
      open(unit=2,status='old',file='tl.grid',form='unformatted')
c
      write(*,*)' '
      write(*,*)'   enter flags:'
      write(*,*)' '
      write(*,*)'       icol (0=grays, 1=colors)'
      write(*,*)'       ictr (0=without contours, 1=with contours)'
      write(*,*)'       ibar (0=without color bar, 1=with color bar)'
      write(*,*)'       ncol (number of colors)'
      write(*,*)' '
      read(*,*)icol,ictr,ibar,ncol
c
      write(*,*)' '
      write(*,*)'   enter delr (km) and delz (m)'
      write(*,*)' '
      read(*,*)delx,dely
c
      write(*,*)' '
      write(*,*)'   enter tlmin, tlmax, dtlbar (dB)'
      write(*,*)' '
      read(*,*)tlmin,tlmax,deltl
      write(*,*)' '
      dtlctr=(tlmax-tlmin)/float(ncol)
c
      if(icol.eq.0)then
      sat=0.0
      ten=0.0
      dten=1.0/float(ncol-1)
      do 1 i=1,ncol
      vhue(i)=1.0
      vten(i)=ten
      ten=ten+dten
    1 continue
c
      else
      sat=1.0
      hue=0.0
      dhue=0.75/float(ncol-1)
      do 2 i=1,ncol
      vhue(i)=hue
      vten(i)=1.0
      hue=hue+dhue
    2 continue
      end if
c
      read(1,*)
      read(1,*)
      read(1,*)xmax
      read(1,*)zdum,dz,ndz,ymax
      read(1,*)
      xmax=xmax*0.001
c
      nsrf=1
   73 read(1,*)rsrf(nsrf),zsrf(nsrf)
      if(rsrf(nsrf).lt.0.0)go to 74
      rsrf(nsrf)=rsrf(nsrf)*0.001
      nsrf=nsrf+1
      go to 73
   74 rsrf(nsrf)=xmax
      zsrf(nsrf)=zsrf(nsrf-1)
c
      nb=1
    3 read(1,*)rb(nb),zb(nb)
      if(rb(nb).lt.0.0)go to 4
      rb(nb)=rb(nb)*0.001
      nb=nb+1
      go to 3
    4 rb(nb)=xmax
      zb(nb)=zb(nb-1)
c
      par1=1.0
      par2=float(ncol)/(tlmax-tlmin)
c
      read(2)ny
      nx=1
    5 read(2,end=6)(tl(j,nx),j=1,ny)
      nx=nx+1
      go to 5
c
    6 nx=nx-1
      dx=xmax/float(nx-1)
      dy=ymax/float(ny-1)
      call ngrid(tl,nx,ny,tlmin,tlmax,dtlctr)
c
      call comprs
      call area2d(6.0,4.0)
      call xname('Range (km)')
      call yname('Depth (m)')
      call frame
      call graf(0.0,delx,xmax,ymax,-dely,0.0)
c
c     Color between the contours in each pixel.
c
      fa=-1.0e15
      fb=tlmin+dtlctr
      ihue=1
      hue=vhue(ihue)
      ten=vten(ihue)
      call hwhsi(hue,sat,ten)
      call kolor(tl,fa,fb,xc,yc,nx,ny,xmin,ymin,dx,dy,ihue)
      do 7 ihue=2,ncol-1
      fa=fb
      fb=fb+dtlctr
      hue=vhue(ihue)
      ten=vten(ihue)
      call hwhsi(hue,sat,ten)
      call kolor(tl,fa,fb,xc,yc,nx,ny,xmin,ymin,dx,dy,ihue)
    7 continue
      fa=fb
      fb=1.0e15
      ihue=ncol
      hue=vhue(ihue)
      ten=vten(ihue)
      call hwhsi(hue,sat,ten)
      call kolor(tl,fa,fb,xc,yc,nx,ny,xmin,ymin,dx,dy,ihue)
c
      call hwhsi(0.0,0.0,0.0)
c
c     Draw the contours.
c
      if(ictr.eq.1)then
      tlc=tlmin+dtlctr
    8 call contr(tl,tlc,xc,yc,nc,nx,ny,xmin,ymin,dx,dy)
      tlc=tlc+dtlctr
      if(tlc.lt.tlmax-0.5*dtlctr)go to 8
      end if
c
c     Draw the surface and bathymetry.
c
      call hwhsi(0.0,0.0,0.0)
c
      call thkcrv(0.02)
      call curve(rb,zb,nb,0)
      call curve(rsrf,zsrf,nsrf,0)
      call endpl(0)
c
c     Draw the color bar.
c
      if(ibar.eq.0)go to 11
      call rearea(6.0,0.2)
      call xname2('Loss (dB re 1 m)')
      call frame
      call graf2(tlmin,deltl,tlmax,0.0,1.0,1.0)
c
      yb(1)=0.0
      yb(2)=0.0
      yb(3)=1.0
      yb(4)=1.0
c
      dtll=(tlmax-tlmin)/float(ncol)
      tll=tlmin+0.5*dtll
      do 9 j=1,ncol
      tl(1,j)=tll
      tll=tll+dtll
    9 continue
      i=1
c
      j=1
      xb(1)=tlmin
      xb(2)=tlmin
      xb(3)=tlmin
      xb(4)=tlmin
      jhue=ifix(par1+par2*(tl(1,j)-tlmin))
      if(jhue.gt.ncol)jhue=ncol
      if(jhue.lt.1)jhue=1
c
   10 j=j+1
      xb(2)=xb(2)+dtll
      xb(3)=xb(3)+dtll
c
      if(j.ge.ncol+1)then
      xb(2)=tlmax
      xb(3)=tlmax
      hue=vhue(jhue)
      ten=vten(jhue)
      call hwhsi(hue,sat,ten)
      call shade(xb,yb,4,90.0,0.0,1,0,0)
      yb(1)=yb(1)+1.0
      yb(2)=yb(2)+1.0
      yb(3)=yb(3)+1.0
      yb(4)=yb(4)+1.0
      go to 11
      end if
c
      ihue=ifix(par1+par2*(tl(1,j)-tlmin))
      if(ihue.gt.ncol)ihue=ncol
      if(ihue.lt.1)ihue=1
      if(ihue.eq.jhue)go to 10
c
      hue=vhue(jhue)
      ten=vten(jhue)
      call hwhsi(hue,sat,ten)
      call shade(xb,yb,4,90.0,0.0,1,0,0)
      xb(1)=xb(2)
      xb(4)=xb(3)
      jhue=ihue
      go to 10
c
   11 call hwhsi(0.0,0.0,0.0)
      call endpl(0)
      call donepl
c
      stop
      end
c
c     Add a perturbation when a contour passes through a grid point.  
c
      subroutine ngrid(tl,nx,ny,tlmin,tlmax,dtlctr)
      real tl(1000,1000)
      iseed=1999
      eps=sqrt(tlmin**2+tlmax**2)*0.0001
c
      tlc=tlmin+dtlctr
    1 do 4 i=1,ny
      do 3 j=1,nx
c
    2 diff=tl(i,j)-tlc
      if(diff.eq.0.0)then
c      tl(i,j)=tl(i,j)+eps*(1.0-2.0*ran(iseed))
      tl(i,j)=tl(i,j)+eps*0.57713
      go to 2
      end if
c
    3 continue
    4 continue
c
      tlc=tlc+dtlctr
      if(tlc.lt.tlmax-0.5*dtlctr)go to 1
c
      return
      end
c
      subroutine comprs
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      open(unit=40,status='unknown',file='ramclr.ps')
      write(40,'(a)')'%!'
      write(40,'(a)')'/center {'
      write(40,'(a)')'    dup stringwidth pop'
      write(40,'(a)')'    2 div neg 0 rmoveto'
      write(40,'(a)')'} bind def'
      write(40,'(a)')'/bdef {bind def} bind def'
      write(40,'(a)')'/ldef {load def} bdef'
      write(40,'(a)')'/d /stroke ldef /l /lineto ldef /m /moveto ldef'
      write(40,'(a)')'/c /sethsbcolor ldef /k /closepath ldef  '
      write(40,'(a)')'/f /fill ldef'
      ifrm=0
      idsh=0
c
      return
      end
c
      subroutine rearea(ww,hh)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      w=ww
      h=hh
      xx0=(612.0-72.0*w+36.0)/2.0
      yy0=(792.0-72.0*4.0+36.0)/2.0+72.0*4.2
      xx1=xx0+72.0*w
      yy1=yy0+72.0*h
c
      return
      end
c
      subroutine area2d(ww,hh)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      w=ww
      h=hh
      xx0=(612.0-72.0*w+36.0)/2.0
      yy0=(792.0-72.0*h+36.0)/2.0
      xx1=xx0+72.0*w
      yy1=yy0+72.0*h
c
      return
      end
c
      subroutine frame
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      ifrm=1
      return
      end
c
      subroutine graf2(xmin,delx,xmax,ymin,dely,ymax)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      eps=0.000001
c
      write(40,'(a)')'/Times-Roman findfont'
      write(40,'(a)')'14 scalefont'
      write(40,'(a)')'setfont'
c
      ax=xx0-w*72.0*xmin/(xmax-xmin)
      bx=w*72.0/(xmax-xmin)
      ay=yy0-h*72.0*ymin/(ymax-ymin)
      by=h*72.0/(ymax-ymin)
c
      mx=ifix(0.5+(xmax-xmin)/delx)
      xtick=xmin
      do 1 i=1,mx+1
      xx=ax+bx*xtick
      yy=yy0+0.2*72.0
      write(40,'(2f6.1,a)')xx,yy,' m'
      yy=yy0+0.2*72.0+9.0
      write(40,'(2f6.1,a)')xx,yy,' l'
      write(40,'(a)')'d'
      yy=yy0+0.2*72.0+12.0
      write(40,'(2f6.1,a)')xx,yy,' m'
      jtick=ifix(alog10(abs(xtick)+eps)+eps)
      jtick=max(0,jtick)
      if(xtick.gt.-eps)then
      itick=ifix(xtick+0.5)
      else
      itick=ifix(xtick-0.5)
      jtick=jtick+1
      end if
c
      if(jtick.eq.0)write(40,3)itick
      if(jtick.eq.1)write(40,4)itick
      if(jtick.eq.2)write(40,5)itick
      if(jtick.eq.3)write(40,6)itick
      if(jtick.eq.4)write(40,7)itick
c
      xtick=xtick+delx
    1 continue
c
    3 format('(',i1,')',' center show')
    4 format('(',i2,')',' center show')
    5 format('(',i3,')',' center show')
    6 format('(',i4,')',' center show')
    7 format('(',i5,')',' center show')
c
      return
      end
c
      subroutine graf(xmin,delx,xmax,ymin,dely,ymax)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      eps=0.000001
c
      write(40,'(a)')'/Times-Roman findfont'
      write(40,'(a)')'14 scalefont'
      write(40,'(a)')'setfont'
c
      ax=xx0-w*72.0*xmin/(xmax-xmin)
      bx=w*72.0/(xmax-xmin)
      ay=yy0-h*72.0*ymin/(ymax-ymin)
      by=h*72.0/(ymax-ymin)
c
      mx=ifix(0.5+(xmax-xmin)/delx)
      xtick=xmin
      do 1 i=1,mx+1
      xx=ax+bx*xtick
      yy=yy0
      write(40,'(2f6.1,a)')xx,yy,' m'
      yy=yy0-9.0
      write(40,'(2f6.1,a)')xx,yy,' l'
      write(40,'(a)')'d'
      yy=yy0-22.0
      write(40,'(2f6.1,a)')xx,yy,' m'
      jtick=ifix(alog10(abs(xtick)+eps)+eps)
      jtick=max(0,jtick)
      if(xtick.gt.-eps)then
      itick=ifix(xtick+0.5)
      else
      itick=ifix(xtick-0.5)
      jtick=jtick+1
      end if
      if(jtick.eq.0)write(40,3)itick
      if(jtick.eq.1)write(40,4)itick
      if(jtick.eq.2)write(40,5)itick
      if(jtick.eq.3)write(40,6)itick
      if(jtick.eq.4)write(40,7)itick
      xtick=xtick+delx
    1 continue
c
      my=ifix(0.5+(ymax-ymin)/dely)
      ytick=ymin
      do 2 i=1,my+1
      xx=xx0
      yy=ay+by*ytick
      write(40,'(2f6.1,a)')xx,yy,' m'
      xx=xx0-9.0
      write(40,*)xx,yy,' l'
      write(40,'(a)')'d'
      xx=xx0-12.0
      write(40,*)xx,yy,' m'
      jtick=ifix(alog10(abs(ytick)+eps)+eps)
      jtick=max(0,jtick)
      if(ytick.gt.-eps)then
      itick=ifix(ytick+0.5)
      else
      itick=ifix(ytick-0.5)
      jtick=jtick+1
      end if
      write(40,'(a)')'90 rotate'
      if(jtick.eq.0)write(40,3)itick
      if(jtick.eq.1)write(40,4)itick
      if(jtick.eq.2)write(40,5)itick
      if(jtick.eq.3)write(40,6)itick
      if(jtick.eq.4)write(40,7)itick
      write(40,'(a)')'270 rotate'
      ytick=ytick+dely
    2 continue
c
    3 format('(',i1,')',' center show')
    4 format('(',i2,')',' center show')
    5 format('(',i3,')',' center show')
    6 format('(',i4,')',' center show')
    7 format('(',i5,')',' center show')
c
      return
      end
c
      subroutine curve(x,y,n,icrv)
      real x(5000),y(5000)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
c     Solid curve.
c
      if(idsh.eq.0)then
      xx=ax+bx*x(1)
      yy=ay+by*y(1)
      write(40,'(2f6.1,a)')xx,yy,' m'
      do 1 i=2,n
      xx=ax+bx*x(i)
      yy=ay+by*y(i)
      write(40,'(2f6.1,a)')xx,yy,' l'
    1 continue
      write(40,'(a)')'d'
      return
      end if
c
c     Dashed curve.
c
      s=0.0
      dels=10.0
      iblck=1
      xx=ax+bx*x(1)
      yy=ay+by*y(1)
      write(40,'(2f6.1,a)')xx,yy,' m'
      xxold=xx
      yyold=yy
c
      do 3 i=2,n
      xx=ax+bx*x(i)
      yy=ay+by*y(i)
c
    2 dx=xx-xxold
      dy=yy-yyold
      ds=sqrt(dx**2+dy**2)
      s=s+ds
c
      if(s.le.dels)then
      xxold=xx
      yyold=yy
      if(iblck.eq.1)write(40,'(2f6.1,a)')xx,yy,' l'
      go to 3
      end if
c
      frac=1.0-(s-dels)/ds
      xxold=xxold+frac*dx
      yyold=yyold+frac*dy
c
      if(iblck.eq.1)then
      write(40,'(2f6.1,a)')xxold,yyold,' l'
      write(40,'(a)')'d'
      iblck=0
      dels=3.0
c
      else
      iblck=1
      dels=10.0
      end if
c
      write(40,'(2f6.1,a)')xxold,yyold,' m'
      s=0.0
      go to 2
c
    3 continue
c
      if(iblck.eq.1)write(40,'(a)')'d'
c
      return
      end
c
      subroutine xname2(xchar)
      character xout*80,xchar*(*)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      xx=0.5*(xx0+xx1)
      yy=yy0+0.2*72.0+32.0
      xout='('//xchar//') center show'
      write(40,'(a)')'/Times-Roman findfont'
      write(40,'(a)')'18 scalefont'
      write(40,'(a)')'setfont'
      write(40,'(2f6.1,a)')xx, yy,' m'
      write(40,'(a)')xout
c
      return
      end
c
      subroutine xname(xchar)
      character xout*80,xchar*(*)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      xx=0.5*(xx0+xx1)
      yy=yy0-44.0
      xout='('//xchar//') center show'
      write(40,'(a)')'/Times-Roman findfont'
      write(40,'(a)')'18 scalefont'
      write(40,'(a)')'setfont'
      write(40,'(2f6.1,a)')xx, yy,' m'
      write(40,'(a)')xout
c
      return
      end
c
      subroutine yname(ychar)
      character yout*80,ychar*(*)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      xx=xx0-32.0
      yy=0.5*(yy0+yy1)
      yout='('//ychar//') center show'
      write(40,'(a)')'/Times-Roman findfont'
      write(40,'(a)')'18 scalefont'
      write(40,'(a)')'setfont'
      write(40,'(2f6.1,a)')xx, yy,' m'
      write(40,'(a)')'90 rotate'
      write(40,'(a)')yout
      write(40,'(a)')'270 rotate'
c
      return
      end
c
      subroutine dash
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      idsh=1-idsh
      return
      end
c
      subroutine shade(xb,yb,nc,dum1,dum2,idum2,idum3,idum4)
      real xb(100),yb(100)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      xx=ax+bx*xb(1)
      yy=ay+by*yb(1)
      write(40,'(2f6.1,a)')xx,yy,' m'
      do 1 ic=2,nc
      xx=ax+bx*xb(ic)
      yy=ay+by*yb(ic)
      write(40,'(2f6.1,a)')xx,yy,' l'
    1 continue
      write(40,'(a)')'k f'
c
      return
      end
c
      subroutine thkcrv(thknss)
      write(40,*)thknss*72.0,' setlinewidth'
      return
      end
c
      subroutine hwhsi(hue,sat,brt)
      write(40,'(3f5.2,a)')hue,sat,brt,' c'
      return
      end
c
c     Construct contour curves.
c
      subroutine contr(f,f0,xc,yc,nc,nx,ny,xmin,ymin,dx,dy)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      real f(1000,1000),xc(4),yc(4)
c
c     Draw the contour curve segments inside each rectangle.
c
      do 2 i=1,ny-1
      yy=float(i-1)*dy
c
      do 1 j=1,nx-1
      xx=float(j-1)*dx
      nc=0
c
      f1=f(i,j)
      f2=f(i+1,j)
      f3=f(i,j+1)
      f4=f(i+1,j+1)
c
      det12=(f1-f0)*(f2-f0)
      det24=(f2-f0)*(f4-f0)
      det13=(f1-f0)*(f3-f0)
      det34=(f3-f0)*(f4-f0)
c
      if(det13.le.0.0)then
      nc=nc+1
      xc(nc)=xx+dx*(f0-f1)/(f3-f1)
      yc(nc)=yy
      end if
c
      if(det34.le.0.0)then
      nc=nc+1
      xc(nc)=xx+dx
      yc(nc)=yy+dy*(f0-f3)/(f4-f3)
      end if
c
      if(det24.le.0.0)then
      nc=nc+1
      xc(nc)=xx+dx*(f0-f2)/(f4-f2)
      yc(nc)=yy+dy
      end if
c
      if(det12.le.0.0)then
      nc=nc+1
      xc(nc)=xx
      yc(nc)=yy+dy*(f0-f1)/(f2-f1)
      end if
c
      if(nc.ge.2)call cnect(f,f1,f2,f3,f4,xc,yc,nc,dx,dy,xx,yy)
c
    1 continue
    2 continue
c
      return
      end
c
c     Connect the contour points.
c
      subroutine cnect(f,f1,f2,f3,f4,xc,yc,nc,dx,dy,xx,yy)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      real f(1000,1000),xc(4),yc(4)
c
      if(nc.lt.4)then
      xxc=ax+bx*xc(1)
      yyc=ay+by*yc(1)
      write(40,'(2f6.1,a)')xxc,yyc,' m'
      xxc=ax+bx*xc(2)
      yyc=ay+by*yc(2)
      write(40,'(2f6.1,a)')xxc,yyc,' l'
      write(40,'(a)')' d'
      return
      end if
c
c     The saddle point case. 
c     Approximate f(x,y) by a*x + b*y + c*x*y + d. 
c     The mappings x = x(f0,y) and y = y(f0,x) are
c     singular at y = -a/c and x = -b/c. 
c
      b=(xx+dx)*(f2-f1)-xx*(f4-f3)
      c=f1-f2-f3+f4
      xsng=-b/c
      det=(xc(1)-xsng)*(xc(2)-xsng)
c
      if(det.lt.0.0)then
      xxc=ax+bx*xc(1)
      yyc=ay+by*yc(1)
      write(40,'(2f6.1,a)')xxc,yyc,' m'
      xxc=ax+bx*xc(4)
      yyc=ay+by*yc(4)
      write(40,'(2f6.1,a)')xxc,yyc,' l'
      write(40,'(a)')' d'
      xxc=ax+bx*xc(2)
      yyc=ay+by*yc(2)
      write(40,'(2f6.1,a)')xxc,yyc,' m'
      xxc=ax+bx*xc(3)
      yyc=ay+by*yc(3)
      write(40,'(2f6.1,a)')xxc,yyc,' l'
      write(40,'(a)')' d'
c
      else
      xxc=ax+bx*xc(1)
      yyc=ay+by*yc(1)
      write(40,'(2f6.1,a)')xxc,yyc,' m'
      xxc=ax+bx*xc(2)
      yyc=ay+by*yc(2)
      write(40,'(2f6.1,a)')xxc,yyc,' l'
      write(40,'(a)')' d'
      xxc=ax+bx*xc(3)
      yyc=ay+by*yc(3)
      write(40,'(2f6.1,a)')xxc,yyc,' m'
      xxc=ax+bx*xc(4)
      yyc=ay+by*yc(4)
      write(40,'(2f6.1,a)')xxc,yyc,' l'
      write(40,'(a)')' d'
      end if
c
      return
      end
c
c     Color between the contours.
c
      subroutine kolor(f,fa,fb,xc,yc,nx,ny,xmin,ymin,dx,dy,ihue)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      real f(1000,1000),xc(100),yc(100)
      integer iab(100)
c
      do 6 i=1,ny-1
      yy=float(i-1)*dy
c
      do 5 j=1,nx-1
      xx=float(j-1)*dx
c
      call sides(f,f1,f2,f3,f4,fa,fb,i,j,xc,yc,xx,yy,dx,dy,
     >   nc,iab)
c
      if(nc.lt.3)go to 5
c
      if((nc.ge.3).and.(nc.le.5))then
      call shade(xc,yc,nc,90.0,0.0,1,0,0)
      go to 5
      end if
c
c     Saddle points are possible for nc > 5.
c
      do 1 ic=1,nc
      xc(8+ic)=xc(ic)
      yc(8+ic)=yc(ic)
    1 continue
c
      jc=0
      kc=1
      isadl=0
    2 jc=jc+1
      xc(jc)=xc(8+kc)
      yc(jc)=yc(8+kc)
      kc=kc+1
      if(kc.gt.nc)go to 3
      if((iab(kc).eq.0).or.(iab(kc).ne.iab(jc)))go to 2
c
      xc1=xc(jc)
      xc2=xc(kc)
      call ksadl(f1,f2,f3,f4,xx,yy,dx,dy,det,xc1,xc2)
      if(det.gt.0.0)go to 2
      isadl=kc
      if(iab(kc+1).eq.0)then
      kc=kc+3
      else
      kc=kc+4
      end if
      if(kc.le.nc)go to 2
    3 call shade(xc,yc,jc,90.0,0.0,1,0,0)
c
      if(isadl.gt.0)then
      jc=nc-jc
      kc=isadl
      do 4 ic=1,jc
      xc(ic)=xc(8+kc)
      yc(ic)=yc(8+kc)
      kc=kc+1
    4 continue
      call shade(xc,yc,jc,90.0,0.0,1,0,0)
      end if
c
    5 continue
    6 continue
c
      return
      end
c
c     Sort out the saddles for coloring.
c
      subroutine ksadl(f1,f2,f3,f4,xx,yy,dx,dy,det,xc1,xc2)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      b=(xx+dx)*(f2-f1)-xx*(f4-f3)
      c=f1-f2-f3+f4
      xsng=-b/c
      det=(xc1-xsng)*(xc2-xsng)
c
      return
      end
c
c     Obtain the nodal points on each side.
c
      subroutine sides(f,f1,f2,f3,f4,fa,fb,i,j,xc,yc,xx,yy,dx,dy,
     >   nc,iab)
      real f(1000,1000),xc(100),yc(100)
      integer iab(100)
c
      nc=0
      f1=f(i,j)
      f2=f(i+1,j)
      f3=f(i,j+1)
      f4=f(i+1,j+1)
c
      det=(f1-fa)*(f1-fb)
      if(det.le.0.0)then
      icrn1=1
      nc=nc+1
      iab(nc)=0
      xc(nc)=xx
      yc(nc)=yy
      end if
c
      deta=(f1-fa)*(f3-fa)
      detb=(f1-fb)*(f3-fb)
      if(deta.lt.0.0)then
      nc=nc+1
      iab(nc)=1
      xc(nc)=xx+dx*abs((fa-f1)/(f3-f1))
      yc(nc)=yy
      end if
      if(detb.lt.0.0)then
      nc=nc+1
      iab(nc)=2
      xc(nc)=xx+dx*abs((fb-f1)/(f3-f1))
      yc(nc)=yy
      end if
      if((deta.lt.0.0).and.(detb.lt.0.0))then
      if(xc(nc-1).gt.xc(nc))then
      temp=xc(nc-1)
      itemp=iab(nc-1)
      xc(nc-1)=xc(nc)
      iab(nc-1)=iab(nc)
      xc(nc)=temp
      iab(nc)=itemp
      end if
      end if
c
      det=(f3-fa)*(f3-fb)
      if(det.le.0.0)then
      icrn3=1
      nc=nc+1
      iab(nc)=0
      xc(nc)=xx+dx
      yc(nc)=yy
      end if
c
      deta=(f3-fa)*(f4-fa)
      detb=(f3-fb)*(f4-fb)
      if(deta.lt.0.0)then
      nc=nc+1
      iab(nc)=1
      xc(nc)=xx+dx
      yc(nc)=yy+dy*abs((fa-f3)/(f4-f3))
      end if
      if(detb.lt.0.0)then
      nc=nc+1
      iab(nc)=2
      xc(nc)=xx+dx
      yc(nc)=yy+dy*abs((fb-f3)/(f4-f3))
      end if
      if((deta.lt.0.0).and.(detb.lt.0.0))then
      if(yc(nc-1).gt.yc(nc))then
      temp=yc(nc-1)
      itemp=iab(nc-1)
      yc(nc-1)=yc(nc)
      iab(nc-1)=iab(nc)
      yc(nc)=temp
      iab(nc)=itemp
      end if
      end if
c
      det=(f4-fa)*(f4-fb)
      if(det.le.0.0)then
      icrn4=1
      nc=nc+1
      iab(nc)=0
      xc(nc)=xx+dx
      yc(nc)=yy+dy
      end if
c
      deta=(f4-fa)*(f2-fa)
      detb=(f4-fb)*(f2-fb)
      if(deta.lt.0.0)then
      nc=nc+1
      iab(nc)=1
      xc(nc)=xx+dx*abs((fa-f2)/(f4-f2))
      yc(nc)=yy+dy
      end if
      if(detb.lt.0.0)then
      nc=nc+1
      iab(nc)=2
      xc(nc)=xx+dx*abs((fb-f2)/(f4-f2))
      yc(nc)=yy+dy
      end if
      if((deta.lt.0.0).and.(detb.lt.0.0))then
      if(xc(nc-1).lt.xc(nc))then
      temp=xc(nc-1)
      itemp=iab(nc-1)
      xc(nc-1)=xc(nc)
      iab(nc-1)=iab(nc)
      xc(nc)=temp
      iab(nc)=itemp
      end if
      end if
c
      det=(f2-fa)*(f2-fb)
      if(det.le.0.0)then
      icrn2=1
      nc=nc+1
      iab(nc)=0
      xc(nc)=xx
      yc(nc)=yy+dy
      end if
c
      deta=(f1-fa)*(f2-fa)
      detb=(f1-fb)*(f2-fb)
      if(deta.lt.0.0)then
      nc=nc+1
      iab(nc)=1
      xc(nc)=xx
      yc(nc)=yy+dy*abs((fa-f1)/(f2-f1))
      end if
      if(detb.lt.0.0)then
      nc=nc+1
      iab(nc)=1
      xc(nc)=xx
      yc(nc)=yy+dy*abs((fb-f1)/(f2-f1))
      end if
      if((deta.lt.0.0).and.(detb.lt.0.0))then
      if(yc(nc-1).lt.yc(nc))then
      temp=yc(nc-1)
      itemp=iab(nc-1)
      yc(nc-1)=yc(nc)
      iab(nc-1)=iab(nc)
      yc(nc)=temp
      iab(nc)=itemp
      end if
      end if
c
      return
      end
c
      subroutine endpl(iplot)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      if(ifrm.eq.1)then
      write(40,'(a)')'2.5 setlinewidth'
      write(40,'(2f6.1,a)')xx0,yy0,' m'
      write(40,'(2f6.1,a)')xx1,yy0,' l'
      write(40,'(2f6.1,a)')xx1,yy1,' l'
      write(40,'(2f6.1,a)')xx0,yy1,' l'
      write(40,'(2f6.1,a)')xx0,yy0,' l'
      write(40,'(2f6.1,a)')xx1,yy0,' l'
      write(40,'(a)')'d'
      write(40,'(a)')'1 setlinewidth'
      end if
c
      return
      end
c
      subroutine donepl
      write(40,'(a)')'showpage'
      return
      end

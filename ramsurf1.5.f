      program ramsurf
c
c     This version of ram handles variable surface height, including
c     rough surfaces and propagation in beach environments [J. Acoust.
c     Soc. Am. 97, 2767-2770 (1995)]. The input file is similar to the 
c     input file for ram but contains a block of data for the location 
c     of the surface. The computational grid extends from z=0 to z=zmax. 
c     The pressure release surface is located at z=zsrf(r), which is 
c     linearly interpolated between input points just like the bathymetry. 
c     The inputs zsrf must be greater than or equal to zero. As in ramgeo, 
c     the layering in the sediment parallels the bathymetry. 
c
c     ******************************************************************
c     ***** Range-dependent Acoustic Model, Version 1.5b, 19-Oct-00 ****
c     ******************************************************************
c
      complex ci,ksq,ksqb,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2,ksqw
      real k0
      INTEGER argc,dirLen,inLen,outLen,resLen
      CHARACTER input*1024,output*1024,directory*1024,resultat*1024
c
c     mr=bathymetry points, mz=depth grid, mp=pade terms.
c
      parameter (mr=505,mz=20002,mp=10)
      dimension rb(mr),zb(mr),cw(mz),cb(mz),rhob(mz),attn(mz),attw(mz),
     >   alpb(mz),f1(mz),f2(mz),f3(mz),ksq(mz),ksqw(mz),ksqb(mz),u(mz),
     >   v(mz),tlg(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),s1(mz,mp),
     >   s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp),rsrf(mr),zsrf(mr),
     >   alpw(mz)
c
      open(unit=1,status='old',file='ram.in')
      open(unit=2,status='unknown',file='tl.line')
      open(unit=3,status='unknown',file='tl.grid',form='unformatted')
c
      call setup(mr,mz,nz,mp,np,ns,mdr,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,
     >   dz,pi,eta,eps,omega,rmax,c0,k0,ci,r,rp,rs,rb,zb,cw,cb,rhob,
     >   attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,
     >   pd1,pd2,tlg,rsrf,zsrf,izsrf,isrf,attw)
c
c     March the acoustic field out in range.
c
    1 call updat(mr,mz,nz,mp,np,iz,ib,dr,dz,eta,omega,rmax,c0,k0,ci,r,
     >   rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,
     >   r1,r2,r3,s1,s2,s3,pd1,pd2,rsrf,zsrf,izsrf,isrf,attw)
      call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
      r=r+dr
      call outpt(mz,mdr,ndr,ndz,iz,nzplt,lz,ir,dir,eps,r,f3,u,tlg)
      if(r.lt.rmax)go to 1
c
      close(1)
      close(2)
      close(3)
c
      stop
      end
c
c     Initialize the parameters, acoustic field, and matrices.
c
      subroutine setup(mr,mz,nz,mp,np,ns,mdr,ndr,ndz,iz,nzplt,lz,ib,ir,
     >   dir,dr,dz,pi,eta,eps,omega,rmax,c0,k0,ci,r,rp,rs,rb,zb,cw,cb,
     >   rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,
     >   s3,pd1,pd2,tlg,rsrf,zsrf,izsrf,isrf,attw)
      complex ci,u(mz),v(mz),ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),
     >   r3(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp),
     >   ksqw(mz)
      real k0,rb(mr),zb(mr),cw(mz),cb(mz),rhob(mz),attn(mz),alpw(mz),
     >   alpb(mz),f1(mz),f2(mz),f3(mz),tlg(mz),rsrf(mr),
     >   zsrf(mr),attw(mz)
c
      read(1,*)
      read(1,*)freq,zs,zr
      read(1,*)rmax,dr,ndr
      read(1,*)zmax,dz,ndz,zmplt
      read(1,*)c0,np,ns,rs
c
      i=1
    1 read(1,*)rsrf(i),zsrf(i)
      if(rsrf(i).lt.0.0)go to 2
      i=i+1
      go to 1
    2 rsrf(i)=2.0*rmax
      zsrf(i)=zsrf(i-1)
c
      i=1
    3 read(1,*)rb(i),zb(i)
      if(rb(i).lt.0.0)go to 4
      i=i+1
      go to 3
    4 rb(i)=2.0*rmax
      zb(i)=zb(i-1)
c
      pi=4.0*atan(1.0)
      ci=cmplx(0.0,1.0)
      eta=1.0/(40.0*pi*alog10(exp(1.0)))
      eps=1.0e-20
      ib=1
      isrf=1
      mdr=0
      r=dr
      omega=2.0*pi*freq
      ri=1.0+zr/dz
      ir=ifix(ri)
      dir=ri-float(ir)
      k0=omega/c0
      nz=zmax/dz-0.5
      nzplt=zmplt/dz-0.5
c
      z=zsrf(1)
      izsrf=1.0+z/dz
c
      z=zb(1)
      iz=1.0+z/dz
      iz=max(2,iz)
      iz=min(nz,iz)
      if(rs.lt.dr)rs=2.0*rmax
c
      if(nz+2.gt.mz)then
      write(*,*)'   Need to increase parameter mz to ',nz+2
      stop
      end if
      if(np.gt.mp)then
      write(*,*)'   Need to increase parameter mp to ',np
      stop
      end if
      if(i.gt.mr)then
      write(*,*)'   Need to increase parameter mr to ',i
      stop
      end if
c
      do 5 j=1,mp
      r3(1,j)=0.0
      r1(nz+2,j)=0.0
    5 continue
      do 6 i=1,nz+2
      u(i)=0.0
      v(i)=0.0
    6 continue
      lz=0
      do 7 i=ndz,nzplt,ndz
      lz=lz+1
    7 continue
      write(3)lz
c
c     The initial profiles and starting field.
c
      call profl(mz,nz,ci,dz,eta,omega,rmax,c0,k0,rp,cw,cb,rhob,attn,
     >   alpw,alpb,ksqw,ksqb,attw)
      call selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,pi,c0,k0,rhob,alpw,alpb,ksq,
     >   ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2,izsrf)
      call outpt(mz,mdr,ndr,ndz,iz,nzplt,lz,ir,dir,eps,r,f3,u,tlg)
c
c     The propagation matrices.
c
      call epade(mp,np,ns,1,k0,c0,dr,pd1,pd2)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
     >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2,izsrf)
c
      return
      end
c
c     Set up the profiles.
c
      subroutine profl(mz,nz,ci,dz,eta,omega,rmax,c0,k0,rp,cw,cb,rhob,
     >   attn,alpw,alpb,ksqw,ksqb,attw)
      complex ci,ksqb(mz),ksqw(mz)
      real k0,cw(mz),cb(mz),rhob(mz),attn(mz),alpw(mz),alpb(mz),
     >   attw(mz)
c
      call zread(mz,nz,dz,cw)
      call zread(mz,nz,dz,attw)
      call zread(mz,nz,dz,cb)
      call zread(mz,nz,dz,rhob)
      call zread(mz,nz,dz,attn)
      rp=2.0*rmax
      read(1,*,end=1)rp
c
    1 do 2 i=1,nz+2
c      ksqw(i)=(omega/cw(i))**2-k0**2
      ksqw(i)=((omega/cw(i))*(1.0+ci*eta*attw(i)))**2-k0**2
      ksqb(i)=((omega/cb(i))*(1.0+ci*eta*attn(i)))**2-k0**2
      alpw(i)=sqrt(cw(i)/c0)
      alpb(i)=sqrt(rhob(i)*cb(i)/c0)
    2 continue
c
      return
      end
c
c     Profile reader and interpolator.
c
      subroutine zread(mz,nz,dz,prof)
      real prof(mz)
c
      do 1 i=1,nz+2
      prof(i)=-1.0
    1 continue
      read(1,*)zi,profi
      prof(1)=profi
      i=1.5+zi/dz
      prof(i)=profi
      iold=i
    2 read(1,*)zi,profi
      if(zi.lt.0.0)go to 3
      i=1.5+zi/dz
      if(i.eq.iold)i=i+1
      prof(i)=profi
      iold=i
      go to 2
    3 prof(nz+2)=prof(i)
      i=1
      j=1
    4 i=i+1
      if(prof(i).lt.0.0)go to 4
      if(i-j.eq.1)go to 6
      do 5 k=j+1,i-1
      prof(k)=prof(j)+float(k-j)*(prof(i)-prof(j))/float(i-j)
    5 continue
    6 j=i
      if(j.lt.nz+2)go to 4
c
      return
      end
c
c     The tridiagonal matrices.
c
      subroutine matrc(mz,nz,mp,np,iz,jz,dz,k0,rhob,alpw,alpb,ksq,ksqw,
     >   ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2,izsrf)
      complex d1,d2,d3,rfact,ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),
     >   r3(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp),
     >   ksqw(mz)
      real k0,rhob(mz),f1(mz),f2(mz),f3(mz),alpw(mz),alpb(mz)
c
      a1=k0**2/6.0
      a2=2.0*k0**2/3.0
      a3=k0**2/6.0
      cfact=0.5/dz**2
      dfact=1.0/12.0
c
      do 1 i=1,iz
      f1(i)=1.0/alpw(i)
      f2(i)=1.0
      f3(i)=alpw(i)
      ksq(i)=ksqw(i)
    1 continue
c
      ii=1
      do 2 i=iz+1,nz+2
      f1(i)=rhob(ii)/alpb(ii)
      f2(i)=1.0/rhob(ii)
      f3(i)=alpb(ii)
      ksq(i)=ksqb(ii)
      ii=ii+1
    2 continue
c
      do 4 i=2,nz+1
c
c     Discretization by Galerkin's method.
c
      c1=cfact*f1(i)*(f2(i-1)+f2(i))*f3(i-1)
      c2=-cfact*f1(i)*(f2(i-1)+2.0*f2(i)+f2(i+1))*f3(i)
      c3=cfact*f1(i)*(f2(i)+f2(i+1))*f3(i+1)
      d1=c1+dfact*(ksq(i-1)+ksq(i))
      d2=c2+dfact*(ksq(i-1)+6.0*ksq(i)+ksq(i+1))
      d3=c3+dfact*(ksq(i)+ksq(i+1))
c
      do 3 j=1,np
      r1(i,j)=a1+pd2(j)*d1
      r2(i,j)=a2+pd2(j)*d2
      r3(i,j)=a3+pd2(j)*d3
      s1(i,j)=a1+pd1(j)*d1
      s2(i,j)=a2+pd1(j)*d2
      s3(i,j)=a3+pd1(j)*d3
    3 continue
    4 continue
c
c     The entries above the surface. 
c
      do 6 j=1,np
      do 5 i=1,izsrf
      r1(i,j)=0.0
      r2(i,j)=1.0
      r3(i,j)=0.0
      s1(i,j)=0.0
      s2(i,j)=0.0
      s3(i,j)=0.0
    5 continue
    6 continue
c
c     The matrix decomposition.
c
      do 8 j=1,np
      do 7 i=2,nz+1
      rfact=1.0/(r2(i,j)-r1(i,j)*r3(i-1,j))
      r1(i,j)=r1(i,j)*rfact
      r3(i,j)=r3(i,j)*rfact
      s1(i,j)=s1(i,j)*rfact
      s2(i,j)=s2(i,j)*rfact
      s3(i,j)=s3(i,j)*rfact
    7 continue
    8 continue
c
      return
      end
c
c     The tridiagonal solver.
c
      subroutine solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
      complex u(mz),v(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),s1(mz,mp),
     >   s2(mz,mp),s3(mz,mp)
c
      do 4 j=1,np
c
      do 1 i=2,nz+1
      v(i)=s1(i,j)*u(i-1)+s2(i,j)*u(i)+s3(i,j)*u(i+1)
    1 continue
c
      do 2 i=3,nz+1
      v(i)=v(i)-r1(i,j)*v(i-1)
    2 continue
c
      do 3 i=nz+1,2,-1
      u(i)=v(i)-r3(i,j)*u(i+1)
    3 continue
c
    4 continue
c
      return
      end
c
c     Matrix updates.
c
      subroutine updat(mr,mz,nz,mp,np,iz,ib,dr,dz,eta,omega,rmax,c0,k0,
     >   ci,r,rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,
     >   f3,r1,r2,r3,s1,s2,s3,pd1,pd2,rsrf,zsrf,izsrf,isrf,attw)
      complex ci,ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),
     >   s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp),ksqw(mz)
      real k0,rb(mr),zb(mr),attn(mz),cb(mz),rhob(mz),cw(mz),
     >   f1(mz),f2(mz),f3(mz),alpw(mz),alpb(mz),rsrf(mr),zsrf(mr),
     >   attw(mz)
c
c     Varying bathymetry.
c
      if(r.ge.rb(ib+1))ib=ib+1
      if(r.ge.rsrf(isrf+1))isrf=isrf+1
c
      jzsrf=izsrf
      z=zsrf(isrf)+(r+0.5*dr-rsrf(isrf))*
     >   (zsrf(isrf+1)-zsrf(isrf))/(rsrf(isrf+1)-rsrf(isrf))
      izsrf=1.0+z/dz
c
      jz=iz
      z=zb(ib)+(r+0.5*dr-rb(ib))*(zb(ib+1)-zb(ib))/(rb(ib+1)-rb(ib))
      iz=1.0+z/dz
      iz=max(2,iz)
      iz=min(nz,iz)
      if((iz.ne.jz).or.(izsrf.ne.jzsrf))call matrc(mz,nz,mp,np,iz,iz,dz,
     >   k0,rhob,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,
     >   pd2,izsrf)
c
c     Varying profiles.
c
      if(r.ge.rp)then
      call profl(mz,nz,ci,dz,eta,omega,rmax,c0,k0,rp,cw,cb,rhob,attn,
     >   alpw,alpb,ksqw,ksqb,attw)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
     >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2,izsrf)
      end if
c
c     Turn off the stability constraints.
c
      if(r.ge.rs)then
      ns=0
      rs=2.0*rmax
      call epade(mp,np,ns,1,k0,c0,dr,pd1,pd2)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
     >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2,izsrf)
      end if
c
      return
      end
c
c     The self-starter.
c
      subroutine selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,pi,c0,k0,rhob,alpw,
     >   alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2,
     >   izsrf)
      complex u(mz),v(mz),ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),
     >   r3(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp),
     >   ksqw(mz)
      real k0,rhob(mz),alpw(mz),alpb(mz),f1(mz),f2(mz),f3(mz)
c
c     Conditions for the delta function.
c
      si=1.0+zs/dz
      is=ifix(si)
      dis=si-float(is)
      u(is)=(1.0-dis)*sqrt(2.0*pi/k0)/(dz*alpw(is))
      u(is+1)=dis*sqrt(2.0*pi/k0)/(dz*alpw(is))
c
c     Divide the delta function by (1-X)**2 to get a smooth rhs.
c
      pd1(1)=0.0
      pd2(1)=-1.0
      call matrc(mz,nz,mp,1,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
     >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2,izsrf)
      call solve(mz,nz,mp,1,iz,u,v,r1,r2,r3,s1,s2,s3)
      call solve(mz,nz,mp,1,iz,u,v,r1,r2,r3,s1,s2,s3)
c
c     Apply the operator (1-X)**2*(1+X)**(-1/4)*exp(ci*k0*r*sqrt(1+X)).
c
      call epade(mp,np,ns,2,k0,c0,dr,pd1,pd2)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
     >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2,izsrf)
      call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
c
      return
      end
c
c     Output transmission loss.
c
      subroutine outpt(mz,mdr,ndr,ndz,iz,nzplt,lz,ir,dir,eps,r,f3,u,tlg)
      complex ur,u(mz)
      real f3(mz),tlg(mz)
c
      ur=(1.0-dir)*f3(ir)*u(ir)+dir*f3(ir+1)*u(ir+1)
      tl=-20.0*alog10(cabs(ur)+eps)+10.0*alog10(r+eps)
      write(2,*)r,tl
c
      mdr=mdr+1
      if(mdr.eq.ndr)then
      mdr=0
c
      j=0
      iflag=1
      do 1 i=ndz,nzplt,ndz
      ur=u(i)*f3(i)
      j=j+1
      tlg(j)=-20.0*alog10(cabs(ur)+eps)+10.0*alog10(r+eps)
c
    1 continue
      write(3)(tlg(j),j=1,lz)
      end if
c
      return
      end
c
c     The coefficients of the rational approximation.
c
      subroutine epade(mp,np,ns,ip,k0,c0,dr,pd1,pd2)
c
      implicit real*8 (a-h,o-z)
      complex*16 ci,z1,z2,g,dg,dh1,dh2,dh3,a,b
      complex*8 pd1(mp),pd2(mp)
      real*8 nu
      real*4 k0,c0,dr
      parameter (m=40)
      dimension bin(m,m),a(m,m),b(m),dg(m),dh1(m),dh2(m),dh3(m),fact(m)
      pi=4.0d0*datan(1.0d0)
      ci=dcmplx(0.0d0,1.0d0)
      sig=k0*dr
      n=2*np
c
      if(ip.eq.1)then
      nu=0.0d0
      alp=0.0d0
      else
      nu=1.0d0
      alp=-0.25d0
      end if
c
c     The factorials.
c
      fact(1)=1.0d0
      do 1 i=2,n
      fact(i)=dfloat(i)*fact(i-1)
    1 continue
c
c     The binomial coefficients.
c
      do 2 i=1,n+1
      bin(i,1)=1.0d0
      bin(i,i)=1.0d0
    2 continue
      do 4 i=3,n+1
      do 3 j=2,i-1
      bin(i,j)=bin(i-1,j-1)+bin(i-1,j)
    3 continue
    4 continue
c
      do 6 i=1,n
      do 5 j=1,n
      a(i,j)=0.0d0
    5 continue
    6 continue
c
c     The accuracy constraints.
c
      call deriv(m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
c
      do 7 i=1,n
      b(i)=dg(i+1)
    7 continue
      do 9 i=1,n
      if(2*i-1.le.n)a(i,2*i-1)=fact(i)
      do 8 j=1,i
      if(2*j.le.n)a(i,2*j)=-bin(i+1,j+1)*fact(j)*dg(i-j+1)
    8 continue
    9 continue
c
c     The stability constraints.
c
      if(ns.ge.1)then
      z1=-3.0d0
      b(n)=-1.0d0
      do 10 j=1,np
      a(n,2*j-1)=z1**j
      a(n,2*j)=0.0d0
   10 continue
      end if
c
      if(ns.ge.2)then
      z1=-1.5d0
      b(n-1)=-1.0d0
      do 11 j=1,np
      a(n-1,2*j-1)=z1**j
      a(n-1,2*j)=0.0d0
   11 continue
      end if
c
      call gauss(m,n,a,b)
c
      dh1(1)=1.0d0
      do 12 j=1,np
      dh1(j+1)=b(2*j-1)
   12 continue
      call fndrt(dh1,np,dh2,m)
      do 13 j=1,np
      pd1(j)=-1.0d0/dh2(j)
   13 continue
c
      dh1(1)=1.0d0
      do 14 j=1,np
      dh1(j+1)=b(2*j)
   14 continue
      call fndrt(dh1,np,dh2,m)
      do 15 j=1,np
      pd2(j)=-1.0d0/dh2(j)
   15 continue
c
      return
      end
c
c     The operator function.
c
      function g(ci,sig,x,alp,nu)
      complex*16 ci,g
      real*8 alp,sig,x,nu
      g=(1.0d0-nu*x)**2*cdexp(alp*dlog(1.0d0+x)+
     >   ci*sig*(-1.0d0+dsqrt(1.0d0+x)))
      return
      end
c
c     The derivatives of the operator function at x=0.
c
      subroutine deriv(m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
      implicit real*8 (a-h,o-z)
      complex*16 ci,dg(m),dh1(m),dh2(m),dh3(m)
      real*8 bin(m,m),nu
      ci=dcmplx(0.0d0,1.0d0)
c
      dh1(1)=0.5d0*ci*sig
      exp1=-0.5d0
      dh2(1)=alp
      exp2=-1.0d0
      dh3(1)=-2.0d0*nu
      exp3=-1.0d0
      do 1 i=2,n
      dh1(i)=dh1(i-1)*exp1
      exp1=exp1-1.0d0
      dh2(i)=dh2(i-1)*exp2
      exp2=exp2-1.0d0
      dh3(i)=-nu*dh3(i-1)*exp3
      exp3=exp3-1.0d0
    1 continue
c
      dg(1)=1.0d0
      dg(2)=dh1(1)+dh2(1)+dh3(1)
      do 3 i=2,n
      dg(i+1)=dh1(i)+dh2(i)+dh3(i)
      do 2 j=1,i-1
      dg(i+1)=dg(i+1)+bin(i,j)*(dh1(j)+dh2(j)+dh3(j))*dg(i-j+1)
    2 continue
    3 continue
c
      return
      end
c
c     Gaussian elimination.
c
      subroutine gauss(m,n,a,b)
      implicit real*8 (a-h,o-z)
      complex*16 a(m,m),b(m)
c
c     Downward elimination.
c
      do 4 i=1,n
      if(i.lt.n)call pivot(m,n,i,a,b)
      a(i,i)=1.0d0/a(i,i)
      b(i)=b(i)*a(i,i)
      if(i.lt.n)then
      do 1 j=i+1,n
      a(i,j)=a(i,j)*a(i,i)
    1 continue
      do 3 k=i+1,n
      b(k)=b(k)-a(k,i)*b(i)
      do 2 j=i+1,n
      a(k,j)=a(k,j)-a(k,i)*a(i,j)
    2 continue
    3 continue
      end if
    4 continue
c
c     Back substitution.
c
      do 6 i=n-1,1,-1
      do 5 j=i+1,n
      b(i)=b(i)-a(i,j)*b(j)
    5 continue
    6 continue
c
      return
      end
c
c     Rows are interchanged for stability.
c
      subroutine pivot(m,n,i,a,b)
      implicit real*8 (a-h,o-z)
      complex*16 temp,a(m,m),b(m)
c
      i0=i
      amp0=cdabs(a(i,i))
      do 1 j=i+1,n
      amp=cdabs(a(j,i))
      if(amp.gt.amp0)then
      i0=j
      amp0=amp
      end if
    1 continue
      if(i0.eq.i)return
c
      temp=b(i)
      b(i)=b(i0)
      b(i0)=temp
      do 2 j=i,n
      temp=a(i,j)
      a(i,j)=a(i0,j)
      a(i0,j)=temp
    2 continue
c
      return
      end
c
c     The root-finding subroutine. 
c
      subroutine fndrt(a,n,z,m)
      complex*16 a(m),z(m),root
      real*8 err
c
      if(n.eq.1)then
      z(1)=-a(1)/a(2)
      return
      end if
      if(n.eq.2)go to 4
c
      do 3 k=n,3,-1
c
c     Obtain an approximate root.
c
      root=0.0d0
      err=1.0d-12
      call guerre(a,k,m,root,err,1000)
c
c     Refine the root by iterating five more times.
c
      err=0.0d0
      call guerre(a,k,m,root,err,5)
      z(k)=root
c
c     Divide out the factor (z-root).
c
      do 1 i=k,1,-1
      a(i)=a(i)+root*a(i+1)
    1 continue
      do 2 i=1,k
      a(i)=a(i+1)
    2 continue
c
    3 continue
c
c     Solve the quadratic equation.
c
    4 z(2)=0.5*(-a(2)+sqrt(a(2)**2-4.0*a(1)*a(3)))/a(3)
      z(1)=0.5*(-a(2)-sqrt(a(2)**2-4.0*a(1)*a(3)))/a(3)
c
      return
      end
c
c     This subroutine finds a root of a polynomial of degree n > 2
c     by Laguerre's method.
c
      subroutine guerre(a,n,m,z,err,nter)
      complex*16 a(m),az(50),azz(50),z,dz,p,pz,pzz,f,g,h,ci
      real*8 amp1,amp2,rn,eps,err
      ci=cmplx(0.0d0,1.0d0)
      eps=1.0d-20
      rn=real(n)
c
c     The coefficients of p'(z) and p''(z).
c
      do 1 i=1,n
      az(i)=float(i)*a(i+1)
    1 continue
      do 2 i=1,n-1
      azz(i)=float(i)*az(i+1)
    2 continue
c
      iter=0
    3 p=a(n)+a(n+1)*z
      do 4 i=n-1,1,-1
      p=a(i)+z*p
    4 continue
      if(abs(p).lt.eps)return
c
      pz=az(n-1)+az(n)*z
      do 5 i=n-2,1,-1
      pz=az(i)+z*pz
    5 continue
c
      pzz=azz(n-2)+azz(n-1)*z
      do 6 i=n-3,1,-1
      pzz=azz(i)+z*pzz
    6 continue
c
c     The Laguerre perturbation.
c
      f=pz/p
      g=f**2-pzz/p
      h=sqrt((rn-1.0d0)*(rn*g-f**2))
      amp1=abs(f+h)
      amp2=abs(f-h)
      if(amp1.gt.amp2)then
      dz=-rn/(f+h)
      else
      dz=-rn/(f-h)
      end if
c
      iter=iter+1
c
c     Rotate by 90 degrees to avoid limit cycles. 
c
      jter=jter+1
      if(jter.eq.10)then
      jter=1
      dz=dz*ci
      end if
      z=z+dz
c
      if(iter.eq.100)then
      write(*,*)' '
      write(*,*)'   Laguerre method not converging.'
      write(*,*)'   Try a different combination of DR and NP.'
      write(*,*)' '
      stop
      end if
c
      if((abs(dz).gt.err).and.(iter.lt.nter))go to 3
c
      return
      end

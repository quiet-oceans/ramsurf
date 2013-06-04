      program rams
c
c     ******************************************************************
c     ***** Range-dependent Acoustic Model, Version 0.5s, 13-Sep-00 ****
c     ******************************************************************
c
c     This code was developed by Michael D. Collins at the Naval
c     Research Laboratory in Washington, DC. It solves range-dependent 
c     seismo-ocean acoustics problems with the split-step Pade 
c     algorithm [M. D. Collins, J. Acoust. Soc. Am. 93, 1736-1742 
c     (1993)]. 
c
c     Version 0.5s contains a correction to a bug in the dimension of
c     quantities passed to subroutines fndrt and guerre that Laurie
c     Fialkowski noticed. 
c
c     Version 0.4s contains a correction to a minor bug in subroutine
c     guerre that Dave King noticed (amp1 and amp2 were declared
c     twice) and a few other minor improvements. 
c
c     Version 0.3s contains a new root-finding subroutine.
c
c     Version 0.2s contains a correction of a bug in the call to solve,
c     which Peter Mignerey helped to locate. The bug only caused
c     problems with certain compilers. 
c
c     Version 0.1s contains two improvements:
c
c     (1) An improved self starter. Stability is improved by using the 
c     factor (1-i*X)**2 instead of (1+X)**2 to smooth the delta function. 
c     The factor (1+X)**2 is nearly singular for some problems involving
c     deep water and/or weak attenuation. Numerical problems associated 
c     with this singularity were detected by Eddie Scheer of Woods Hole 
c     Oceanographic Institute. 
c
c     (2) Elimination of underflow problems. A very small number is 
c     added to the solution in subroutine solve to prevent underflow,
c     which can adversely affect run time on some computers. This
c     improvement was suggested by Ed McDonald of the SACLANT Undersea
c     Research Centre. 
c
      complex ci,lamb,mub,u,v,r1,r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,
     >   s7,t1,t2,t3,t4,t5,t6,t7,pd1,pd2,g0
      real k0,lamw
c
c     mr=bathymetry points, mz=depth grid, mp=pade terms.
c
      parameter (mr=100,mz=10000,mp=10)
      dimension rb(mr),zb(mr),cw(mz),cp(mz),cs(mz),rhob(mz),attnp(mz),
     >   attns(mz),lamw(mz),lamb(mz),mub(mz),u(mz),v(mz),tlg(mz),
     >   r1(mz,mp),r2(mz,mp),r3(mz,mp),r4(mz,mp),r5(mz,mp),r6(mz,mp),
     >   r7(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),s4(mz,mp),s5(mz,mp),
     >   s6(mz,mp),s7(mz,mp),t1(6,mp),t2(6,mp),t3(6,mp),t4(6,mp),
     >   t5(6,mp),t6(6,mp),t7(6,mp),pd1(mp),pd2(mp)
c
      open(unit=1,status='old',file='rams.in')
      open(unit=2,status='unknown',file='tl.line')
      open(unit=3,status='unknown',file='tl.grid',form='unformatted')
c
      call setup(mr,mz,nz,mp,np,mdr,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,
     >   dz,pi,eta,eps,omega,rmax,c0,k0,ci,r,rp,rb,zb,cw,cp,cs,rhob,
     >   attnp,attns,lamw,lamb,mub,u,v,r1,r2,r3,r4,r5,r6,r7,s1,s2,s3,
     >   s4,s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,pd1,pd2,irot,theta,g0,tlg)
c
c     March the seismo-acoustic field out in range.
c
    1 call updat(mr,mz,nz,mp,np,iz,zs,ib,dr,dz,eta,omega,rmax,c0,k0,ci,
     >   r,rp,rb,zb,cw,cp,cs,rhob,attnp,attns,lamw,lamb,mub,u,r1,r2,r3,
     >   r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,pd1,
     >   pd2)
      call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,
     >   s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,g0)
      r=r+dr
      call outpt(mz,mdr,ndr,ndz,iz,nzplt,lz,ir,dir,eps,r,u,tlg)
      if(r.lt.rmax)go to 1
c
      close(1)
      close(2)
      close(3)
c
      stop
      end
c
c     Initialize the parameters, seismo-acoustic field, and matrices.
c
      subroutine setup(mr,mz,nz,mp,np,mdr,ndr,ndz,iz,nzplt,lz,ib,ir,
     >   dir,dr,dz,pi,eta,eps,omega,rmax,c0,k0,ci,r,rp,rb,zb,cw,cp,
     >   cs,rhob,attnp,attns,lamw,lamb,mub,u,v,r1,r2,r3,r4,r5,r6,r7,
     >   s1,s2,s3,s4,s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,pd1,pd2,irot,theta,
     >   g0,tlg)
      complex ci,g0,u(mz),v(mz),lamb(mz),mub(mz),r1(mz,mp),r2(mz,mp),
     >   r3(mz,mp),r4(mz,mp),r5(mz,mp),r6(mz,mp),r7(mz,mp),s1(mz,mp),
     >   s2(mz,mp),s3(mz,mp),s4(mz,mp),s5(mz,mp),s6(mz,mp),s7(mz,mp),
     >   t1(6,mp),t2(6,mp),t3(6,mp),t4(6,mp),t5(6,mp),t6(6,mp),
     >   t7(6,mp),pd1(mp),pd2(mp),nu
      real k0,rb(mr),zb(mr),cw(mz),cp(mz),cs(mz),rhob(mz),attnp(mz),
     >   attns(mz),lamw(mz),tlg(mz)
c
      read(1,*)
      read(1,*)freq,zs,zr
      read(1,*)rmax,dr,ndr
      read(1,*)zmax,dz,ndz,zmplt
      read(1,*)c0,np,irot,theta
c
      i=1
    1 read(1,*)rb(i),zb(i)
      if(rb(i).lt.0.0)go to 2
      i=i+1
      go to 1
    2 zb(i)=zb(i-1)
      rb(i)=rmax+2.0*dr
c
      pi=4.0*atan(1.0)
      ci=cmplx(0.0,1.0)
      nu=ci
      eta=1.0/(40.0*pi*alog10(exp(1.0)))
      eps=1.0e-20
      ib=1
      mdr=0
      r=dr
      omega=2.0*pi*freq
      ri=1.0+zr/dz
      ir=ifix(ri)
      dir=ri-float(ir)
      k0=omega/c0
      nz=zmax/dz-0.5
      nzplt=zmplt/dz-0.5
      z=zb(1)
      iz=z/dz
c
      do 3 i=1,2*nz+4
      u(i)=0.0
      v(i)=0.0
    3 continue
      lz=0
      do 4 i=1+ndz,nzplt,ndz
      lz=lz+1
    4 continue
      write(3)lz
c
c     The initial profiles and starting field.
c
      call profl(mz,nz,ci,dz,eta,omega,rmax,rp,cw,cp,cs,rhob,attnp,
     >   attns,lamw,lamb,mub)
      call selfs(mz,nz,ci,mp,np,iz,zs,dr,dz,pi,c0,k0,omega,rhob,lamw,
     >   lamb,mub,u,v,r1,r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7,t1,t2,
     >   t3,t4,t5,t6,t7,pd1,pd2,nu)
      call outpt(mz,mdr,ndr,ndz,iz,nzplt,lz,ir,dir,eps,r,u,tlg)
c
c     The propagation matrices.
c
      call epade(mp,np,1,1,k0,c0,dr,pd1,pd2,irot,theta,g0,nu)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,omega,rhob,lamw,lamb,mub,r1,
     >   r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,
     >   pd1,pd2)
c
      return
      end
c
c     Set up the profiles.
c
      subroutine profl(mz,nz,ci,dz,eta,omega,rmax,rp,cw,cp,cs,rhob,
     >   attnp,attns,lamw,lamb,mub)
      complex ci,mub(mz),lamb(mz)
      real cw(mz),cp(mz),cs(mz),rhob(mz),attnp(mz),attns(mz),lamw(mz)
c
      call zread(mz,nz,dz,cw)
      call zread(mz,nz,dz,cp)
      call zread(mz,nz,dz,cs)
      call zread(mz,nz,dz,rhob)
      call zread(mz,nz,dz,attnp)
      call zread(mz,nz,dz,attns)
      rp=2.0*rmax
      read(1,*,end=1)rp
c
    1 do 2 i=1,nz+2
      lamw(i)=cw(i)**2
      lamb(i)=rhob(i)*((cp(i)/(1.0+ci*eta*attnp(i)))**2-
     >   2.0*(cs(i)/(1.0+ci*eta*attns(i)))**2)
      mub(i)=rhob(i)*(cs(i)/(1.0+ci*eta*attns(i)))**2
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
c     Output transmission loss.
c
      subroutine outpt(mz,mdr,ndr,ndz,iz,nzplt,lz,ir,dir,eps,r,u,tlg)
      complex ur,u(mz)
      real tlg(mz)
c
      ur=(1.0-dir)*u(2*ir-1)+dir*u(2*ir+1)
      tl=-20.0*alog10(cabs(ur)+eps)+10.0*alog10(r+eps)
      write(2,*)r,tl
      write(*,*)r,tl
c
      mdr=mdr+1
      if(mdr.eq.ndr)then
      mdr=0
c
      j=0
      iflag=1
      do 1 i=1+ndz,nzplt,ndz
      ur=u(2*i-1)
      j=j+1
      tlg(j)=-20.0*alog10(cabs(ur)+eps)+10.0*alog10(r+eps)
c
c     Mark the ocean bottom.
c
c      if((i.gt.iz).and.(iflag.eq.1))then
c      tlg(j)=0.0
c      iflag=0
c      end if
c
    1 continue
      write(3)(tlg(j),j=1,lz)
      end if
c
      return
      end
c
c     Matrix updates and energy conservation.
c
      subroutine updat(mr,mz,nz,mp,np,iz,zs,ib,dr,dz,eta,omega,rmax,c0,
     >   k0,ci,r,rp,rb,zb,cw,cp,cs,rhob,attnp,attns,lamw,lamb,mub,u,r1,
     >   r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,
     >   pd1,pd2)
      complex ci,lamb(mz),mub(mz),u(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),
     >   r4(mz,mp),r5(mz,mp),r6(mz,mp),r7(mz,mp),s1(mz,mp),s2(mz,mp),
     >   s3(mz,mp),s4(mz,mp),s5(mz,mp),s6(mz,mp),s7(mz,mp),t1(6,mp),
     >   t2(6,mp),t3(6,mp),t4(6,mp),t5(6,mp),t6(6,mp),t7(6,mp),pd1(mp),
     >   pd2(mp)
      real k0,rb(mr),zb(mr),cw(mz),cp(mz),cs(mz),rhob(mz),attnp(mz),
     >   attns(mz),lamw(mz)
c
c     Varying bathymetry.
c
      jz=iz
      z=zb(ib)+(r+0.5*dr-rb(ib))*(zb(ib+1)-zb(ib))/(rb(ib+1)-rb(ib))
      iz=z/dz
      if(iz.ne.jz)then
      call matrc(mz,nz,mp,np,iz,jz,dz,k0,omega,rhob,lamw,lamb,mub,r1,
     >   r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,
     >   pd1,pd2)
c
c     An approximate energy-flux correction.
c
      if(iz.lt.jz)then
      fact=sqrt(cw(iz)**3)/sqrt(rhob(iz)*cp(iz)**3)
      do 1 i=iz+1,jz
      u(2*i+1)=u(2*i+1)*fact
    1 continue
c
      else
      fact=sqrt(rhob(iz)*cp(iz)**3)/sqrt(cw(iz)**3)
      do 2 i=jz+1,iz
      u(2*i+1)=u(2*i+1)*fact      
    2 continue
      end if
c
      end if
c
      if(r.ge.rb(ib+1))ib=ib+1
c
c     Varying profiles.
c
      if(r.ge.rp)then
      call profl(mz,nz,ci,dz,eta,omega,rmax,rp,cw,cp,cs,rhob,attnp,
     >   attns,lamw,lamb,mub)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,omega,rhob,lamw,lamb,mub,r1,
     >   r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,
     >   pd1,pd2)
      end if
c
      return
      end
c
c     The self starter.
c
      subroutine selfs(mz,nz,ci,mp,np,iz,zs,dr,dz,pi,c0,k0,omega,rhob,
     >   lamw,lamb,mub,u,v,r1,r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7,
     >   t1,t2,t3,t4,t5,t6,t7,pd1,pd2,nu)
      complex ci,g0,u(mz),v(mz),lamb(mz),mub(mz),r1(mz,mp),r2(mz,mp),
     >   r3(mz,mp),r4(mz,mp),r5(mz,mp),r6(mz,mp),r7(mz,mp),s1(mz,mp),
     >   s2(mz,mp),s3(mz,mp),s4(mz,mp),s5(mz,mp),s6(mz,mp),s7(mz,mp),
     >   t1(6,mp),t2(6,mp),t3(6,mp),t4(6,mp),t5(6,mp),t6(6,mp),
     >   t7(6,mp),pd1(mp),pd2(mp),nu
      real k0,rhob(mz),lamw(mz)
c
c     Conditions for the delta function.
c
      si=1.0+zs/dz
      is=ifix(si)
      dis=si-float(is)
      u(2*is-1)=(1.0-dis)*sqrt(2.0*pi/k0)/dz
      u(2*is+1)=dis*sqrt(2.0*pi/k0)/dz
c
c     Divide the delta function by (1-nu*X)**2 to get a smooth rhs.
c
      pd1(1)=0.0
      pd2(1)=-nu
      g0=(1.0,0.0)
      call matrc(mz,nz,mp,1,iz,iz,dz,k0,omega,rhob,lamw,lamb,mub,r1,
     >   r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,
     >   pd1,pd2)
      call solve(mz,nz,mp,1,iz,u,v,r1,r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,
     >   s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,g0)
      call solve(mz,nz,mp,1,iz,u,v,r1,r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,
     >   s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,g0)
c
c     Apply the operator (1-nu*X)**2*(1+X)**(-1/4)*exp(ci*k0*r*sqrt(1+X)).
c
      call epade(mp,np,1,0,k0,c0,dr,pd1,pd2,0,0.0,g0,nu)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,omega,rhob,lamw,lamb,mub,r1,
     >   r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,
     >   pd1,pd2)
      call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,
     >   s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,g0)
c
      return
      end
c
c     The heptadiagonal matrices.
c
      subroutine matrc(mz,nz,mp,np,iz,jz,dz,k0,omega,rhob,lamw,lamb,
     >   mub,r1,r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7,t1,t2,t3,t4,t5,
     >   t6,t7,pd1,pd2)
      complex l1,l2,l3,l4,l5,l6,l7,m1,m2,m3,m4,m5,m6,m7,a11,a12,a21,
     >   a22,a23,a31,a32,a33,a34,a35,lamb(mz),mub(mz),r1(mz,mp),
     >   r2(mz,mp),r3(mz,mp),r4(mz,mp),r5(mz,mp),r6(mz,mp),r7(mz,mp),
     >   s1(mz,mp),s2(mz,mp),s3(mz,mp),s4(mz,mp),s5(mz,mp),s6(mz,mp),
     >   s7(mz,mp),t1(6,mp),t2(6,mp),t3(6,mp),t4(6,mp),t5(6,mp),
     >   t6(6,mp),t7(6,mp),pd1(mp),pd2(mp)
      real k0,rhob(mz),lamw(mz)
c
c     New matrices when iz.eq.jz.
c
      if(iz.eq.jz)then
      ia=2
      ib=nz-2
c
c     Updated matrices when iz.ne.jz.
c
      else
      ia=min(iz,jz)
      ib=max(iz,jz)
      end if
c
c     Eq. (4) of jasa 86, 1459-1464 in fluid layer.
c
      do 2 i=ia-1,iz
      l1=0.0
      l2=(lamw(i)+lamw(i+1))/12.0
      l3=0.0
      l4=(lamw(i)+6.0*lamw(i+1)+lamw(i+2))/12.0
      l5=0.0
      l6=(lamw(i+1)+lamw(i+2))/12.0
      l7=0.0
      m1=0.0
      m2=lamw(i+1)/dz**2+omega**2/6.0+
     >   0.5*(lamw(i)-lamw(i+1))/dz**2+
     >   0.5*(lamw(i)-lamw(i+1))/dz**2
      m3=0.0
      m4=-2.0*lamw(i+1)/dz**2+2.0*omega**2/3.0+
     >   0.5*(2.0*lamw(i+1)-lamw(i)-lamw(i+2))/dz**2+
     >   0.5*(lamw(i)-2.0*lamw(i+1)+lamw(i+2))/dz**2
      m5=0.0
      m6=lamw(i+1)/dz**2+omega**2/6.0+
     >   0.5*(lamw(i+2)-lamw(i+1))/dz**2+
     >   0.5*(lamw(i+2)-lamw(i+1))/dz**2
      m7=0.0
c
      do 1 j=1,np
      r1(2*i-1,j)=k0**2*l1+pd2(j)*(m1-k0**2*l1)
      r2(2*i-1,j)=k0**2*l2+pd2(j)*(m2-k0**2*l2)
      r3(2*i-1,j)=k0**2*l3+pd2(j)*(m3-k0**2*l3)
      r4(2*i-1,j)=k0**2*l4+pd2(j)*(m4-k0**2*l4)
      r5(2*i-1,j)=k0**2*l5+pd2(j)*(m5-k0**2*l5)
      r6(2*i-1,j)=k0**2*l6+pd2(j)*(m6-k0**2*l6)
      r7(2*i-1,j)=k0**2*l7+pd2(j)*(m7-k0**2*l7)
c
      s1(2*i-1,j)=k0**2*l1+pd1(j)*(m1-k0**2*l1)
      s2(2*i-1,j)=k0**2*l2+pd1(j)*(m2-k0**2*l2)
      s3(2*i-1,j)=k0**2*l3+pd1(j)*(m3-k0**2*l3)
      s4(2*i-1,j)=k0**2*l4+pd1(j)*(m4-k0**2*l4)
      s5(2*i-1,j)=k0**2*l5+pd1(j)*(m5-k0**2*l5)
      s6(2*i-1,j)=k0**2*l6+pd1(j)*(m6-k0**2*l6)
      s7(2*i-1,j)=k0**2*l7+pd1(j)*(m7-k0**2*l7)
    1 continue
    2 continue
c
c     Eq. (2) of jasa 86, 1459-1464 in fluid layer.
c
      do 4 i=ia-1,iz
      m1=(-1.0/6.0)*(lamw(i)+2.0*lamw(i+1))/dz+
     >   (1.0/6.0)*(-lamw(i)+lamw(i+1))/dz
      m2=0.0
      m3=(1.0/6.0)*(lamw(i)-lamw(i+2))/dz+
     >   (1.0/3.0)*(-lamw(i)+lamw(i+2))/dz
      m4=omega**2
      m5=(1.0/6.0)*(lamw(i+2)+2.0*lamw(i+1))/dz+
     >   (1.0/6.0)*(lamw(i+2)-lamw(i+1))/dz
      m6=0.0
      m7=0.0
c
      do 3 j=1,np
      r1(2*i,j)=m1
      r2(2*i,j)=m2
      r3(2*i,j)=m3
      r4(2*i,j)=m4
      r5(2*i,j)=m5
      r6(2*i,j)=m6
      r7(2*i,j)=m7
      s1(2*i,j)=0.0
      s2(2*i,j)=0.0
      s3(2*i,j)=0.0
      s4(2*i,j)=0.0
      s5(2*i,j)=0.0
      s6(2*i,j)=0.0
      s7(2*i,j)=0.0
    3 continue    
    4 continue
c
c     Eq. (4) of jasa 86, 1459-1464 in solid layer.
c
      do 6 i=iz+1,ib+2
      l1=0.0
      l2=(lamb(i)+2.0*mub(i)+lamb(i+1)+2.0*mub(i+1))/12.0
      l3=(2.0/6.0)*(-mub(i)+mub(i+1))/dz
      l4=(lamb(i)+2.0*mub(i)+6.0*lamb(i+1)+12.0*mub(i+1)+
     >   lamb(i+2)+2.0*mub(i+2))/12.0
      l5=(2.0/3.0)*(-mub(i)+mub(i+2))/dz
      l6=(lamb(i+1)+2.0*mub(i+1)+lamb(i+2)+2.0*mub(i+2))/12.0
      l7=(2.0/6.0)*(mub(i+2)-mub(i+1))/dz
      m1=0.0
      m2=(lamb(i+1)+2.0*mub(i+1))/dz**2+
     >   omega**2*(rhob(i)+rhob(i+1))/12.0+
     >   0.5*(lamb(i)+2.0*mub(i)-lamb(i+1)-2.0*mub(i+1))/dz**2+
     >   0.5*(lamb(i)-lamb(i+1))/dz**2
      m3=(omega**2/6.0)*(-rhob(i)+rhob(i+1))/dz+
     >   2.0*(-mub(i)+mub(i+1))/dz**3
      m4=-2.0*(lamb(i+1)+2.0*mub(i+1))/dz**2+
     >   omega**2*(rhob(i)+6.0*rhob(i+1)+rhob(i+2))/12.0+
     >   0.5*(2.0*lamb(i+1)+4.0*mub(i+1)-lamb(i)-2.0*mub(i)-
     >   lamb(i+2)-2.0*mub(i+2))/dz**2+
     >   0.5*(lamb(i)-2.0*lamb(i+1)+lamb(i+2))/dz**2
      m5=(omega**2/3.0)*(-rhob(i)+rhob(i+2))/dz+
     >   2.0*(mub(i)-mub(i+2))/dz**3
      m6=(lamb(i+1)+2.0*mub(i+1))/dz**2+
     >   omega**2*(rhob(i+1)+rhob(i+2))/12.0+
     >   0.5*(lamb(i+2)+2.0*mub(i+2)-lamb(i+1)-2.0*mub(i+1))/dz**2+
     >   0.5*(lamb(i+2)-lamb(i+1))/dz**2
      m7=(omega**2/6.0)*(rhob(i+2)-rhob(i+1))/dz+
     >   2.0*(mub(i+2)-mub(i+1))/dz**3
c
      do 5 j=1,np
      r1(2*i-1,j)=k0**2*l1+pd2(j)*(m1-k0**2*l1)
      r2(2*i-1,j)=k0**2*l2+pd2(j)*(m2-k0**2*l2)
      r3(2*i-1,j)=k0**2*l3+pd2(j)*(m3-k0**2*l3)
      r4(2*i-1,j)=k0**2*l4+pd2(j)*(m4-k0**2*l4)
      r5(2*i-1,j)=k0**2*l5+pd2(j)*(m5-k0**2*l5)
      r6(2*i-1,j)=k0**2*l6+pd2(j)*(m6-k0**2*l6)
      r7(2*i-1,j)=k0**2*l7+pd2(j)*(m7-k0**2*l7)
c
      s1(2*i-1,j)=k0**2*l1+pd1(j)*(m1-k0**2*l1)
      s2(2*i-1,j)=k0**2*l2+pd1(j)*(m2-k0**2*l2)
      s3(2*i-1,j)=k0**2*l3+pd1(j)*(m3-k0**2*l3)
      s4(2*i-1,j)=k0**2*l4+pd1(j)*(m4-k0**2*l4)
      s5(2*i-1,j)=k0**2*l5+pd1(j)*(m5-k0**2*l5)
      s6(2*i-1,j)=k0**2*l6+pd1(j)*(m6-k0**2*l6)
      s7(2*i-1,j)=k0**2*l7+pd1(j)*(m7-k0**2*l7)
    5 continue
    6 continue
c
c     Eq. (2) of jasa 86, 1459-1464 in solid layer.
c
      do 8 i=iz+1,ib+2
      l1=0.0
      l2=(mub(i)+mub(i+1))/12.0
      l3=0.0
      l4=(mub(i)+6.0*mub(i+1)+mub(i+2))/12.0
      l5=0.0
      l6=(mub(i+1)+mub(i+2))/12.0
      l7=0.0
      m1=(-1.0/6.0)*(lamb(i)+mub(i)+2.0*lamb(i+1)+
     >   2.0*mub(i+1))/dz+
     >   (1.0/6.0)*(-lamb(i)+lamb(i+1))/dz
      m2=mub(i+1)/dz**2+
     >   omega**2*(rhob(i)+rhob(i+1))/12.0+
     >   (mub(i)-mub(i+1))/dz**2
      m3=(1.0/6.0)*(lamb(i)+mub(i)-lamb(i+2)-mub(i+2))/dz+
     >   (1.0/3.0)*(-lamb(i)+lamb(i+2))/dz
      m4=-2.0*mub(i+1)/dz**2+
     >   omega**2*(rhob(i)+6.0*rhob(i+1)+rhob(i+2))/12.0+
     >   (2.0*mub(i+1)-mub(i)-mub(i+2))/dz**2
      m5=(1.0/6.0)*(lamb(i+2)+mub(i+2)+2.0*lamb(i+1)+
     >   2.0*mub(i+1))/dz+
     >   (1.0/6.0)*(lamb(i+2)-lamb(i+1))/dz
      m6=mub(i+1)/dz**2+
     >   omega**2*(rhob(i+1)+rhob(i+2))/12.0+
     >   (mub(i+2)-mub(i+1))/dz**2
      m7=0.0
c
      do 7 j=1,np
      r1(2*i,j)=k0**2*l1+pd2(j)*(m1-k0**2*l1)
      r2(2*i,j)=k0**2*l2+pd2(j)*(m2-k0**2*l2)
      r3(2*i,j)=k0**2*l3+pd2(j)*(m3-k0**2*l3)
      r4(2*i,j)=k0**2*l4+pd2(j)*(m4-k0**2*l4)
      r5(2*i,j)=k0**2*l5+pd2(j)*(m5-k0**2*l5)
      r6(2*i,j)=k0**2*l6+pd2(j)*(m6-k0**2*l6)
      r7(2*i,j)=k0**2*l7+pd2(j)*(m7-k0**2*l7)
c
      s1(2*i,j)=k0**2*l1+pd1(j)*(m1-k0**2*l1)
      s2(2*i,j)=k0**2*l2+pd1(j)*(m2-k0**2*l2)
      s3(2*i,j)=k0**2*l3+pd1(j)*(m3-k0**2*l3)
      s4(2*i,j)=k0**2*l4+pd1(j)*(m4-k0**2*l4)
      s5(2*i,j)=k0**2*l5+pd1(j)*(m5-k0**2*l5)
      s6(2*i,j)=k0**2*l6+pd1(j)*(m6-k0**2*l6)
      s7(2*i,j)=k0**2*l7+pd1(j)*(m7-k0**2*l7)
    7 continue    
    8 continue
c
c     Continuity conditions for a fluid-solid interface.
c
      a11=lamw(iz)/lamw(iz+2)
      a12=-2.0*dz*omega**2/lamw(iz+2)
c
      a21=-dz*lamw(iz+1)/mub(iz+1)
      a22=dz*lamb(iz+1)/mub(iz+1)
      a23=1.0
c
      a31=-2.0*(mub(iz+1)+mub(iz+2))*
     >   lamw(iz+1)/(mub(iz+1)*lamb(iz))
      a32=2.0*(mub(iz+1)+mub(iz+2))*
     >   lamb(iz+1)/(mub(iz+1)*lamb(iz))
      a33=2.0*rhob(iz+1)*dz*omega**2/lamb(iz)-
     >   2.0*(mub(iz)+2.0*mub(iz+1)+mub(iz+2))/(dz*lamb(iz))
      a34=lamb(iz+2)/lamb(iz)
      a35=2.0*(mub(iz)+2.0*mub(iz+1)+mub(iz+2))/(dz*lamb(iz))
c
      do 9 j=1,np
c
      i=2*iz-1
      r2(i,j)=r2(i,j)+a11*r6(i,j)
      s2(i,j)=s2(i,j)+a11*s6(i,j)
      r7(i,j)=a12*r6(i,j)
      s7(i,j)=a12*s6(i,j)
      r6(i,j)=0.0
      s6(i,j)=0.0
c
      i=2*iz
      r1(i,j)=r1(i,j)+a11*r5(i,j)
      s1(i,j)=s1(i,j)+a11*s5(i,j)
      r6(i,j)=a12*r5(i,j)
      s6(i,j)=a12*s5(i,j)
      r5(i,j)=0.0
      s5(i,j)=0.0
c
      i=2*iz+1
      r5(i,j)=r5(i,j)+a33*r2(i,j)
      s5(i,j)=s5(i,j)+a33*s2(i,j)
      r6(i,j)=r6(i,j)+a34*r2(i,j)
      s6(i,j)=s6(i,j)+a34*s2(i,j)
      r7(i,j)=r7(i,j)+a35*r2(i,j)
      s7(i,j)=s7(i,j)+a35*s2(i,j)
      r4(i,j)=r4(i,j)+a32*r2(i,j)
      s4(i,j)=s4(i,j)+a32*s2(i,j)
      r2(i,j)=a31*r2(i,j)
      s2(i,j)=a31*s2(i,j)
c
      r2(i,j)=r2(i,j)+a21*r3(i,j)
      s2(i,j)=s2(i,j)+a21*s3(i,j)
      r4(i,j)=r4(i,j)+a22*r3(i,j)
      s4(i,j)=s4(i,j)+a22*s3(i,j)
      r7(i,j)=r7(i,j)+a23*r3(i,j)
      s7(i,j)=s7(i,j)+a23*s3(i,j)
      r3(i,j)=0.0
      s3(i,j)=0.0
c
      i=2*iz+2
      r3(i,j)=r3(i,j)+a32*r1(i,j)
      s3(i,j)=s3(i,j)+a32*s1(i,j)
      r4(i,j)=r4(i,j)+a33*r1(i,j)
      s4(i,j)=s4(i,j)+a33*s1(i,j)
      r5(i,j)=r5(i,j)+a34*r1(i,j)
      s5(i,j)=s5(i,j)+a34*s1(i,j)
      r6(i,j)=r6(i,j)+a35*r1(i,j)
      s6(i,j)=s6(i,j)+a35*s1(i,j)
      r1(i,j)=a31*r1(i,j)
      s1(i,j)=a31*s1(i,j)
c
      r1(i,j)=r1(i,j)+a21*r2(i,j)
      s1(i,j)=s1(i,j)+a21*s2(i,j)
      r3(i,j)=r3(i,j)+a22*r2(i,j)
      s3(i,j)=s3(i,j)+a22*s2(i,j)
      r6(i,j)=r6(i,j)+a23*r2(i,j)
      s6(i,j)=s6(i,j)+a23*s2(i,j)
      r2(i,j)=0.0
      s2(i,j)=0.0
    9 continue
c
c     The matrix decomposition.
c
      do 20 j=1,np
c
      do 13 i=2*ia-3,2*iz
      if(i.le.3)go to 10
      r2(i,j)=r2(i,j)-r1(i,j)*r5(i-3,j)
      r3(i,j)=r3(i,j)-r1(i,j)*r6(i-3,j)
      r4(i,j)=r4(i,j)-r1(i,j)*r7(i-3,j)
   10 if(i.le.2)go to 11
      r3(i,j)=r3(i,j)-r2(i,j)*r5(i-2,j)
      r4(i,j)=r4(i,j)-r2(i,j)*r6(i-2,j)
      r5(i,j)=r5(i,j)-r2(i,j)*r7(i-2,j)
   11 if(i.le.1)go to 12
      r4(i,j)=r4(i,j)-r3(i,j)*r5(i-1,j)
      r5(i,j)=r5(i,j)-r3(i,j)*r6(i-1,j)
      r6(i,j)=r6(i,j)-r3(i,j)*r7(i-1,j)
   12 r4(i,j)=1.0/r4(i,j)
      r5(i,j)=r5(i,j)*r4(i,j)
      r6(i,j)=r6(i,j)*r4(i,j)
      r7(i,j)=r7(i,j)*r4(i,j)
   13 continue
c
      i1=2*nz-2
      i2=2*nz-1
      i3=2*nz
      do 17 i=2*ib+4,2*iz+1,-1
      if(i.ge.i1)go to 14
      r6(i,j)=r6(i,j)-r7(i,j)*r3(i+3,j)
      r5(i,j)=r5(i,j)-r7(i,j)*r2(i+3,j)
      r4(i,j)=r4(i,j)-r7(i,j)*r1(i+3,j)
   14 if(i.ge.i2)go to 15
      r5(i,j)=r5(i,j)-r6(i,j)*r3(i+2,j)
      r4(i,j)=r4(i,j)-r6(i,j)*r2(i+2,j)
      r3(i,j)=r3(i,j)-r6(i,j)*r1(i+2,j)
   15 if(i.ge.i3)go to 16
      r4(i,j)=r4(i,j)-r5(i,j)*r3(i+1,j)
      r3(i,j)=r3(i,j)-r5(i,j)*r2(i+1,j)
      r2(i,j)=r2(i,j)-r5(i,j)*r1(i+1,j)
   16 r4(i,j)=1.0/r4(i,j)
      r3(i,j)=r3(i,j)*r4(i,j)
      r2(i,j)=r2(i,j)*r4(i,j)
      r1(i,j)=r1(i,j)*r4(i,j)
   17 continue
c
c     The six-by-six matrix near interface.
c
      i0=2*iz-3
      do 18 i=1,3
      t4(i,j)=1.0
      t5(i,j)=r5(i0+i,j)
      t6(i,j)=r6(i0+i,j)
      t7(i,j)=r7(i0+i,j)
   18 continue
c
      do 19 i=4,6
      t1(i,j)=r1(i0+i,j)
      t2(i,j)=r2(i0+i,j)
      t3(i,j)=r3(i0+i,j)
      t4(i,j)=1.0
      t5(i,j)=0.0
      t6(i,j)=0.0
      t7(i,j)=0.0
   19 continue
c
      t2(4,j)=t2(4,j)-t1(4,j)*t5(1,j)
      t3(4,j)=t3(4,j)-t1(4,j)*t6(1,j)
      t4(4,j)=t4(4,j)-t1(4,j)*t7(1,j)
c
      t3(4,j)=t3(4,j)-t2(4,j)*t5(2,j)
      t4(4,j)=t4(4,j)-t2(4,j)*t6(2,j)
      t5(4,j)=t5(4,j)-t2(4,j)*t7(2,j)
c
      t4(4,j)=t4(4,j)-t3(4,j)*t5(3,j)
      t5(4,j)=t5(4,j)-t3(4,j)*t6(3,j)
      t6(4,j)=t6(4,j)-t3(4,j)*t7(3,j)
c
      t5(4,j)=t5(4,j)/t4(4,j)
      t6(4,j)=t6(4,j)/t4(4,j)
      t4(4,j)=1.0/t4(4,j)
c
      t2(5,j)=t2(5,j)-t1(5,j)*t5(2,j)
      t3(5,j)=t3(5,j)-t1(5,j)*t6(2,j)
      t4(5,j)=t4(5,j)-t1(5,j)*t7(2,j)
c
      t3(5,j)=t3(5,j)-t2(5,j)*t5(3,j)
      t4(5,j)=t4(5,j)-t2(5,j)*t6(3,j)
      t5(5,j)=t5(5,j)-t2(5,j)*t7(3,j)
c
      t4(5,j)=t4(5,j)-t3(5,j)*t5(4,j)
      t5(5,j)=t5(5,j)-t3(5,j)*t6(4,j)
c
      t5(5,j)=t5(5,j)/t4(5,j)
      t4(5,j)=1.0/t4(5,j)
c
      t2(6,j)=t2(6,j)-t1(6,j)*t5(3,j)
      t3(6,j)=t3(6,j)-t1(6,j)*t6(3,j)
      t4(6,j)=t4(6,j)-t1(6,j)*t7(3,j)
c
      t3(6,j)=t3(6,j)-t2(6,j)*t5(4,j)
      t4(6,j)=t4(6,j)-t2(6,j)*t6(4,j)
c
      t4(6,j)=t4(6,j)-t3(6,j)*t5(5,j)
c
      t4(6,j)=1.0/t4(6,j)
c
   20 continue
c
      return
      end
c
c     The heptadiagonal solver.
c
      subroutine solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,r4,r5,r6,r7,s1,s2,
     >   s3,s4,s5,s6,s7,t1,t2,t3,t4,t5,t6,t7,g0)
      complex u(mz),v(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),r4(mz,mp),
     >   r5(mz,mp),r6(mz,mp),r7(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),
     >   s4(mz,mp),s5(mz,mp),s6(mz,mp),s7(mz,mp),t1(6,mp),t2(6,mp),
     >   t3(6,mp),t4(6,mp),t5(6,mp),t6(6,mp),t7(6,mp),g0
      eps=1.0e-30
c
      do 6 j=1,np
c
c     The right side.
c
      i=1
      v(i+2)=s2(i,j)*u(i)+s3(i,j)*u(i+1)+s4(i,j)*u(i+2)+
     >   s5(i,j)*u(i+3)+s6(i,j)*u(i+4)+s7(i,j)*u(i+5)+eps
      do 1 i=2,2*nz
      v(i+2)=s1(i,j)*u(i-1)+s2(i,j)*u(i)+s3(i,j)*u(i+1)+
     >   s4(i,j)*u(i+2)+s5(i,j)*u(i+3)+s6(i,j)*u(i+4)+
     >   s7(i,j)*u(i+5)+eps
    1 continue
c
      i1=2*iz-1
c
c     The elimination steps.
c
      i=1
      v(i+2)=(v(i+2)-r2(i,j)*v(i)-r3(i,j)*v(i+1))*r4(i,j)+eps
      do 2 i=2,2*iz
      v(i+2)=(v(i+2)-r1(i,j)*v(i-1)-r2(i,j)*v(i)-
     >   r3(i,j)*v(i+1))*r4(i,j)+eps
    2 continue
      i=2*nz
      v(i+2)=(v(i+2)-r6(i,j)*v(i+4)-r5(i,j)*v(i+3))*r4(i,j)+eps
      do 3 i=2*nz-1,2*iz+1,-1
      v(i+2)=(v(i+2)-r7(i,j)*v(i+5)-r6(i,j)*v(i+4)-
     >   r5(i,j)*v(i+3))*r4(i,j)+eps
    3 continue
c
c     The six-by-six system near the ocean bottom is solved.
c
      v(4+i1)=(v(4+i1)-t1(4,j)*v(1+i1)-t2(4,j)*v(2+i1)-
     >   t3(4,j)*v(3+i1))*t4(4,j)+eps
      v(5+i1)=(v(5+i1)-t1(5,j)*v(2+i1)-t2(5,j)*v(3+i1)-
     >   t3(5,j)*v(4+i1))*t4(5,j)+eps
      u(6+i1)=(v(6+i1)-t1(6,j)*v(3+i1)-t2(6,j)*v(4+i1)-
     >   t3(6,j)*v(5+i1))*t4(6,j)+eps
c
      u(5+i1)=v(5+i1)-t5(5,j)*u(6+i1)+eps
      u(4+i1)=v(4+i1)-t5(4,j)*u(5+i1)-t6(4,j)*u(6+i1)+eps
      u(3+i1)=v(3+i1)-t5(3,j)*u(4+i1)-t6(3,j)*u(5+i1)-
     >   t7(3,j)*u(6+i1)+eps
      u(2+i1)=v(2+i1)-t5(2,j)*u(3+i1)-t6(2,j)*u(4+i1)-
     >   t7(2,j)*u(5+i1)+eps
      u(1+i1)=v(1+i1)-t5(1,j)*u(2+i1)-t6(1,j)*u(3+i1)-
     >   t7(2,j)*u(4+i1)+eps
c
c     The back substitution steps.
c
      do 4 i=2*iz-3,1,-1
      u(i+2)=v(i+2)-r5(i,j)*u(i+3)-r6(i,j)*u(i+4)-
     >   r7(i,j)*u(i+5)+eps
    4 continue
      do 5 i=2*iz+4,2*nz
      u(i+2)=v(i+2)-r3(i,j)*u(i+1)-r2(i,j)*u(i)-
     >   r1(i,j)*u(i-1)+eps
    5 continue
    6 continue
c
      do 7 i=1,nz
      u(2*i+1)=g0*u(2*i+1)
      u(2*i+2)=g0*u(2*i+2)
    7 continue
c
      return
      end
c
c     The rotated Pade coefficients [J. Acoust. Soc. Am. 101, 760-766
c     (1997)].
c
      subroutine rpade(mp,np,k0,dr,pd1,pd2,theta,g0)
c
      complex ci,g0,tfact,rot0,rot1,rot2,pd1(mp),pd2(mp)
      real k0
      pi=4.0*atan(1.0)
      ci=cmplx(0.0,1.0)
      tfact=cexp(-ci*theta*pi/360.0)
      den=float(2*np+1)
      rot0=1.0
c
      do 1 j=1,np
c
c     The Pade coefficients.
c
      pade1=(2.0/den)*sin(float(j)*pi/den)**2
      pade2=cos(float(j)*pi/den)**2
c
c     The rotated Pade coefficients.
c
      rot1=tfact*pade1/(1.0+pade2*(tfact**2-1.0))**2
      rot2=tfact**2*pade2/(1.0+pade2*(tfact**2-1.0))
      rot0=rot0+pade1*(tfact**2-1.0)/(1.0+pade2*(tfact**2-1.0))
c
c     The Crank-Nicolson coefficients.
c
      pd1(j)=rot2+0.5*ci*k0*dr*rot1
      pd2(j)=rot2-0.5*ci*k0*dr*rot1
c
    1 continue
      rot0=rot0/tfact
      g0=cexp(ci*k0*dr*rot0)
c
      return
      end
c
c     The coefficients of the rational approximation.
c
      subroutine epade(mp,np,ns,ip,k0,c0,dr,pd1,pd2,irot,theta,g0,nu8)
c
      implicit real*8 (a-h,o-z)
      complex*16 ci,z1,z2,g,dg,dh1,dh2,dh3,a,b,nu
      complex*8 ci8,g0,pd1(mp),pd2(mp),nu8
      real*4 k0,c0,dr,theta
      parameter (m=40)
      dimension bin(m,m),a(m,m),b(m),dg(m),dh1(m),dh2(m),dh3(m),fact(m)
      pi=4.0d0*datan(1.0d0)
      ci=dcmplx(0.0d0,1.0d0)
      ci8=cmplx(0.0,1.0)
      sig=k0*dr
      n=2*np
      g0=cexp(ci8*k0*dr)
c
      if((ip.eq.1).and.(irot.eq.1))then
      call rpade(mp,np,k0,dr,pd1,pd2,theta,g0)
      return
      end if
c
      if(ip.eq.1)then
      nu=0.0d0
      alp=0.0d0
c
      else
      nu=nu8
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
      complex*16 ci,g,nu
      real*8 alp,sig,x
      g=(1.0d0-nu*x)**2*cdexp(alp*dlog(1.0d0+x)+
     >   ci*sig*(-1.0d0+dsqrt(1.0d0+x)))
      return
      end
c
c     The derivatives of the operator function at x=0.
c
      subroutine deriv(m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
      implicit real*8 (a-h,o-z)
      complex*16 ci,dg(m),dh1(m),dh2(m),dh3(m),nu
      real*8 bin(m,m)
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

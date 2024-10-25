       subroutine define_van_der_waals
     x  (idnode,safe,ntpatm,dlrpot,rvdw,maxvdw,ntpvdw,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining van der Waals potentials
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c     wl
c     2006/11/28 16:32:48
c     1.10
c     Exp
c     
c***********************************************************************

      use vdw_module 
      use parse_module

      implicit none

      logical safe
      character*8 atom1,atom2,keyword
      character*1 message(80)

      integer ntpvdw,ntpatm,fail,idum,ivdw,maxvdw
      integer itpvdw,keypot,numpar,katom1,katom2,jtpatm,keyvdw,i
      integer ntab,idnode,j
      real(8) dlrpot,engunit,rvdw
      real(8), allocatable :: parpot(:)


      allocate (parpot(mxpvdw))
c      ntpvdw=intstr(record,lenrec,idum)
      do ivdw=1,maxvdw
        
        lstvdw(ivdw)=0
        ltpvdw(ivdw)=-1
        
      enddo
      
      do itpvdw=1,ntpvdw
        
        do i=1,mxpvdw
          parpot(i)=0.d0
        enddo
        
        call getrec(safe,idnode,nfield)
c         if(.not.safe)return
        call copystring(record,message,80)
        call getword(atom1,record,8,lenrec)
        call getword(atom2,record,8,lenrec)
        call lowcase(record,lenrec-16)
        call getword(keyword,record,4,lenrec) 
        if(keyword(1:4).eq.'12-6') then
          keypot=1
          numpar=2
        elseif(keyword(1:4).eq.'lj  ') then
          keypot=2
          numpar=2
        elseif(keyword(1:4).eq.'nm  ') then
          keypot=3
          numpar=4
        elseif(keyword(1:4).eq.'buck') then
          keypot=4
          numpar=3
        elseif(keyword(1:4).eq.'bhm ') then
          keypot=5
          numpar=5
        elseif(keyword(1:4).eq.'hbnd') then
          keypot=6
          numpar=2
        elseif(keyword(1:4).eq.'snm ') then
          keypot=7
          numpar=5
        elseif(keyword(1:4).eq.'hcnm') then
          keypot=7
          numpar=5
        elseif(keyword(1:4).eq.'mors')then
          keypot=8
          numpar=3
        elseif(keyword(1:4).eq.'wca ')then
          keypot=9
          numpar=3
        elseif(keyword(1:4).eq.'tab ') then
          keypot=0
          numpar=0
c        else
c          if(idnode.eq.0) write(nrite,*) message
c          call error(idnode,452)
        endif

        parpot(1)=dblstr(record,lenrec,idum)
        parpot(2)=dblstr(record,lenrec,idum)
        parpot(3)=dblstr(record,lenrec,idum)
        parpot(4)=dblstr(record,lenrec,idum)
        parpot(5)=dblstr(record,lenrec,idum)
        
c        if(idnode.eq.0) 
c     x    write(nrite,"(16x,2a8,2x,a4,3x,1p,9e13.5)") 
c     x    atom1,atom2,keyword,(parpot(j),j=1,numpar)
        
        katom1=0
        katom2=0
        
        do jtpatm=1,ntpatm

          if(atom1.eq.unqatm(jtpatm))katom1=jtpatm
          if(atom2.eq.unqatm(jtpatm))katom2=jtpatm
        enddo

        keyvdw=loc8(katom1,katom2)

c     convert energies to internal unit

c        if(keyvdw.gt.maxvdw) call error(idnode,82)
        
        parpot(1)=parpot(1)*engunit
        
        if(keypot.eq.1) then
          
          parpot(2)=parpot(2)*engunit
          
        else if(keypot.eq.4) then
          
          parpot(3)=parpot(3)*engunit
          
        else if(keypot.eq.5) then
          
          parpot(4)=parpot(4)*engunit
          parpot(5)=parpot(5)*engunit
          
        else if(keypot.eq.6) then
          
          parpot(2)=parpot(2)*engunit
          
        endif

        lstvdw(keyvdw)=itpvdw
 
        ltpvdw(itpvdw)=keypot
        
        do i=1,mxpvdw
          
          prmvdw(itpvdw,i)=parpot(i)
          
        enddo
        
      enddo

c     generate nonbonded force arrays

      if(ntpvdw.gt.0)then
        call forgen(idnode,ntpvdw,dlrpot,rvdw)
      endif

c     check for unspecified atom-atom potentials
      
      ntab=(ntpatm*(ntpatm+1))/2
      
      if(ntpvdw.lt.ntab) then
        
        do i=1,ntab
          
          if(lstvdw(i).eq.0)then
            
            lstvdw(i)=ntpvdw+1
            
          endif
          
        enddo

c     define zero potential for undefined interactions
        
        do i=1,mxgrid
          
          vvv(i,ntpvdw+1)=0.d0
          
        enddo
        
      endif

      deallocate (parpot)

      return
      end

      subroutine forgen(idnode,ntpvdw,dlrpot,rcut)

c***********************************************************************
c     
c     dl_poly subroutine for generating potential energy and 
c     force arrays for van der waals forces only
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith may 1992.
c     
c     wl
c     2006/11/28 16:32:48
c     1.10
c     Exp
c     
c***********************************************************************
      
      use vdw_module
      
      implicit none

      integer i,ivdw,ntpvdw,idnode
      real(8) dlrpot,rcut,rrr,ann,amm,gam,bet,eps,rr0,aaa,bbb
      real(8) ccc,ddd,eee,sig,rho,rrc
CVAM
CVAM      call VTBEGIN(140, ierr)
CVAM

c     define grid resolution for potential arrays
      dlrpot=rcut/dble(mxgrid-4)

c     construct arrays for all types of short ranged  potential
      
      do ivdw=1,ntpvdw
        
        if(ltpvdw(ivdw).eq.1)then
          
c       12 - 6 potential
      
          aaa=prmvdw(ivdw,1)
          bbb=prmvdw(ivdw,2)
          
          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=(aaa/rrr**6-bbb)/rrr**6
            
          enddo
          
        else if(ltpvdw(ivdw).eq.2)then
          
c       lennard-jones potential
      
          eps=prmvdw(ivdw,1)
          sig=prmvdw(ivdw,2)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=4.d0*eps*(sig/rrr)**6*((sig/rrr)**6-1.d0)
            
          enddo
          
        else if(ltpvdw(ivdw).eq.3)then

c       n - m potential
      
          eps=prmvdw(ivdw,1)
          ann=max(prmvdw(ivdw,2),prmvdw(ivdw,3))
          amm=min(prmvdw(ivdw,2),prmvdw(ivdw,3))
          rr0=prmvdw(ivdw,4)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=eps/(ann-amm)*(amm*(rr0/rrr)**ann-
     x        ann*(rr0/rrr)**amm)
            
          enddo
          
        else if(ltpvdw(ivdw).eq.4)then
          
c       buckingham exp - 6 potential
      
          aaa=prmvdw(ivdw,1)
          rho=prmvdw(ivdw,2)
          ccc=prmvdw(ivdw,3)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=aaa*exp(-rrr/rho)-ccc/rrr**6
            
          enddo
          
        else if(ltpvdw(ivdw).eq.5)then
          
c       born-huggins-meyer exp - 6 - 8 potential
      
          aaa=prmvdw(ivdw,1)
          bbb=prmvdw(ivdw,2)
          ccc=prmvdw(ivdw,3)
          ddd=prmvdw(ivdw,4)
          eee=prmvdw(ivdw,5)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=aaa*exp(bbb*(ccc-rrr))-ddd/rrr**6-eee/rrr**8
            
          enddo
          
        else if(ltpvdw(ivdw).eq.6) then
          
c       Hydrogen-bond 12 - 10 potential
      
          aaa=prmvdw(ivdw,1)
          bbb=prmvdw(ivdw,2)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=aaa/rrr**12-bbb/rrr**10
            
          enddo
          
        else if(ltpvdw(ivdw).eq.7) then
          
c       shifted and force corrected n - m potential (w. smith)
      
          eps=prmvdw(ivdw,1)
          ann=prmvdw(ivdw,2)
          amm=prmvdw(ivdw,3)
          rr0=prmvdw(ivdw,4)
          rrc=prmvdw(ivdw,5)
          if(rrc.lt.1.d-6)rrc=rcut
 
c          if(ann.le.amm) call error(idnode,470)

          gam=rrc/rr0
c          if(gam.lt.1.d0) call error(idnode,468)
          bet=gam*((gam**(amm+1.d0)-1.d0)/(gam**(ann+1.d0)-1.d0))
     x      **(1.d0/(ann-amm))
          eps=-eps*(ann-amm)/(amm*(bet**ann)*(1.d0+(ann/gam-ann-1.d0)
     x      /gam**ann)-ann*(bet**amm)*(1.d0+(amm/gam-amm-1.d0)
     x      /gam**amm))

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            if(rrr.gt.rrc)then

              vvv(i,ivdw)=0.d0

            else

              vvv(i,ivdw)=eps/(ann-amm)*(amm*(bet**ann)*((rr0/rrr)**ann-
     x          (1.d0/gam)**ann)-ann*(bet**amm)*((rr0/rrr)**amm-
     x          (1.d0/gam)**amm)+ann*amm*((rrr/(gam*rr0)-1.d0)*
     x          ((bet/gam)**ann-(bet/gam)**amm)))

            endif

          enddo
          
        else if(ltpvdw(ivdw).eq.8) then
          
c       morse potential
          
          eps=prmvdw(ivdw,1)
          rr0=prmvdw(ivdw,2)
          sig=prmvdw(ivdw,3)
          
          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=eps*((1.d0-exp(-sig*(rrr-rr0)))**2-1.d0)
            
          enddo
          
        else if(ltpvdw(ivdw).eq.9) then
          
c       weeks-chandler-anderson potential
          
          eps=prmvdw(ivdw,1)
          sig=prmvdw(ivdw,2)
          rr0=prmvdw(ivdw,3)
          ddd=sig*2.d0**(1.d0/6.d0)
          
          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot-rr0
            if(rrr.gt.ddd)then
              
              vvv(i,ivdw)=0.d0

            else if(rrr.gt.0.d0)then
              
              vvv(i,ivdw)=4.d0*eps*(sig/rrr)**6*
     x          ((sig/rrr)**6-1.d0)+eps
            
            endif
              
          enddo
        endif 
      enddo
      return
      end

      subroutine srfrce(sittyp,ik,engsrp,rcut,dlrpot)

c***********************************************************************
c     
c     dl_poly subroutine for calculating short range force and
c     potential energy terms using verlet neighbour list
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       march 1992
c     
c     version 3
c     author    - t. forester    june  1993
c     stress tensor added t.forester may 1994
c     
c     wl
c     2006/11/28 16:32:48
c     1.10
c     Exp
c     
c***********************************************************************
      
c      use setup_module
c      use config_module
      use vdw_module
c      use property_module
c      use pair_module
      
      implicit none

      integer iatm,ik,m,jatm,k,l,sittyp
      real(8) engsrp,virsrp,rcut,dlrpot
      real(8) fi,rcsq,rdr,strs1,strs2,strs3,strs5,strs6,strs9,ai,aj
      real(8) ab,rrr,rsq,ppp,t1,t2,vk0,vk1,vk2,gk0,gk1,gk2,gamma,fx
      real(8) fy,fz,omega

      dimension fi(3)

c     set cutoff condition for pair forces

      rcsq=rcut**2
c     interpolation spacing
      
      rdr=1.d0/dlrpot

c     initialise potential energy
      
      engsrp=0.d0

c     store forces for iatm 
      
      ai=dble(sittyp)
   

c     start of primary loop for forces evaluation
      
      do m=1,ik
        
c     atomic and potential function indices
        jatm=ilist(m)
        aj=dble(ltype(jatm))
        if(ai.gt.aj) then
          ab=ai*(ai-1.d0)*0.5d0+aj+0.5d0
        else
          ab=aj*(aj-1.d0)*0.5d0+ai+0.5d0
        endif
        k=lstvdw(int(ab))
        
        if((ltpvdw(k).lt.100).and.(abs(vvv(1,k)).gt.1.d-10))then

c     apply truncation of potential
          rsq=rsqdf(m)

          if(rcsq.gt.rsq)then
            
            rrr=sqrt(rsq)
           
            l=int(rrr*rdr)
            ppp=rrr*rdr-dble(l)

c     calculate interaction energy using 3-point interpolation
            
            vk0=vvv(l,k)
            vk1=vvv(l+1,k)
            vk2=vvv(l+2,k)
            
            t1=vk0+(vk1-vk0)*ppp
            t2=vk1+(vk2-vk1)*(ppp-1.0d0)
            
            omega=t1+(t2-t1)*ppp*0.5d0
            engsrp=omega+engsrp


          endif
          
        endif
      enddo

      return
      end

      subroutine lrcorrect
     x  (idnode,imcon,keyfce,natms,ntpatm,ntpvdw,elrc,engunit,
     x  rcut,volm)
      
c*************************************************************************
c     
c     DL_POLY subroutine to evaluate long-range corrections to
c     pressure and energy in a periodic system.
c     
c     copyright daresbury laboratory 1993
c     
c     author -       t. forester may 1993
c     
c     wl
c     2006/11/28 16:32:48
c     1.10
c     Exp
c     
c***************************************************************************
      
c      use setup_module
      use vdw_module
c      use site_module
c      use config_module

      implicit none

      integer idnode,imcon,keyfce,natms,ntpatm,i,ka,ntpvdw
      integer ivdw,j,k
      real(8) elrc,engunit,virlrc,rcut,volm,twopi,plrc,eadd,padd
      real(8) denprd,aaa,bbb,ccc,ddd,eee,eps,sig,rr0,ann,amm
      twopi=2.0d0*pi
      
c     initalise counter arrays
      
      do i=1,ntpatm
        
        numtyp(i)=0
        numfrz(i)=0
        
      enddo

c     evaluate number density in system
      
      do i=1,natms
        
        ka=ltype(i)
        numtyp(ka)=numtyp(ka)+1
        if(lfreezesite(i).ne.0)numfrz(ka)=numfrz(ka)+1
        
      enddo

c     number densities
      
c      do i=1,ntpatm
c
c        dens(i)=dble(numtyp(i))/volm
c        
c      enddo
      
c     long range corrections to energy and pressure
      
      plrc=0.d0
      elrc=0.d0
      
      if(imcon.ne.0.and.imcon.ne.6.and.ntpvdw.gt.0) then 
        
        if(mod(keyfce,2).eq.0) then
          
          ivdw=0
          
          do i=1, ntpatm
            
            do j=1,i
              
              eadd=0.d0
              padd=0.d0
              
              ivdw=ivdw+1
              k=lstvdw(ivdw)
              
              if(ltpvdw(k).eq.0) then
                
c     tabulated potential
                
                eadd=prmvdw(k,1)
                padd=-prmvdw(k,2)
                
              else if(ltpvdw(k).eq.1) then
                
c     12-6 potential
                
                aaa=prmvdw(k,1)
                bbb=prmvdw(k,2)
                
                eadd=aaa/(9.d0*rcut**9)-bbb/(3.d0*rcut**3)
                padd=12.d0*aaa/(9.d0*rcut**9)-6.d0*bbb/(3.d0*rcut**3)
                
              else if(ltpvdw(k).eq.2) then
                
c     Lennard Jones potential
      
                eps=prmvdw(k,1)
                sig=prmvdw(k,2)
                eadd=4.d0*eps*(sig**12/(9.d0*rcut**9)-
     x            sig**6/(3.d0*rcut**3))
                padd=4.d0*eps*(12.d0*sig**12/(9.d0*rcut**9)-
     x            2.d0*sig**6/(rcut**3))
                
              else if(ltpvdw(k).eq.3) then
                
c     n - m potential
                
                eps=prmvdw(k,1)
                ann=prmvdw(k,2)
                amm=prmvdw(k,3)
                rr0=prmvdw(k,4)
                
                eadd=eps/(ann-amm)*(amm*rr0**ann/((ann-3.d0)*
     x            rcut**(ann-3.d0))-ann*rr0**amm/((amm-3.0d0)*
     x            rcut**(amm-3.d0)))
                padd=eps/(ann-amm)*ann*amm*(rr0**ann/((ann-3.d0)*
     x            rcut**(ann-3.d0))-rr0**amm/((amm-3.0d0)*
     x            rcut**(amm-3.d0)))
                
              else if(ltpvdw(k).eq.4) then
                
c     buckingham exp - 6 potential
      
                ccc=prmvdw(k,3)
      
                eadd=-ccc/(3.d0*rcut**3)
                padd=-2.d0*ccc/(rcut**3)
                
              else if(ltpvdw(k).eq.5) then
                
c     born huggins meyer exp -6 - 8  potential
      
                ddd=prmvdw(k,4)
                eee=prmvdw(k,5)
                
                eadd=-ddd/(3.d0*rcut**3)-eee/(5.d0*rcut**5)
                padd=-2.d0*ddd/(rcut**3)-8.d0*eee/(5.d0*rcut**5)
                
              else if(ltpvdw(k).eq.6) then
                
c     hydrogen bond  12 - 10 potential
      
                aaa=prmvdw(k,1)
                bbb=prmvdw(k,2)

                eadd=aaa/(9.d0*rcut**9)-bbb/(7.d0*rcut**7)
                padd=12.d0*aaa/(9.d0*rcut**9)-1.d1*bbb/(7.d0*rcut**7)
                
              endif
              
              if(i.ne.j) then
                eadd=eadd*2.d0
                padd=padd*2.d0
              endif
c             store eadd so this calculation does not need to be
c             done again during guest insertions and deletions.
              steadd(ivdw)=eadd

              denprd=twopi*(dble(numtyp(i))*dble(numtyp(j))-
     x          dble(numfrz(i))*dble(numfrz(j)))/volm**2
              elrc=elrc+volm*denprd*eadd
              plrc=plrc+denprd*padd/3.d0
              
            enddo
            
          enddo
          
        endif
        
      endif
      
c      if(idnode.eq.0) write(nrite,
c     x  "(/,/,'long range correction for: vdw energy  ',e15.6,/,
c     x  25x,': vdw pressure',e15.6)") elrc/engunit,plrc*prsunt
      
      return
      end

      subroutine gstlrcorrect
     x  (idnode,imcon,mol,keyfce,natms,ntpatm,ntpvdw,engunit,
     x  delrc,rcut,volm)
      
c*************************************************************************
c     
c     subroutine to evaluate long-range corrections to
c     pressure and energy in a periodic system for a single guest.
c    
c     
c***************************************************************************
      
      use vdw_module

      implicit none

      integer idnode,imcon,keyfce,natms,ntpatm,i,ka,ntpvdw
      integer ivdw,j,k,mol
      real(8) engunit,rcut,volm,twopi,eadd
      real(8) denprd,delrc
      twopi=2.0d0*pi


c     this part will be put before calling this subroutine in gcmc.f      
c      do i=1,natms
        
c        ka=ltpsit(mol,i)
c        numtyp(ka)=numtyp(ka)+1
c        if(lfzsite(mol,i).ne.0)numfrz(ka)=numfrz(ka)+1
c      enddo

c     long range corrections to energy and pressure
      
      delrc=0.d0
      if(imcon.ne.0.and.imcon.ne.6.and.ntpvdw.gt.0) then 
        
        if(mod(keyfce,2).eq.0) then
          ivdw=0
          
          do i=1, ntpatm
            
            do j=1,i
              
              eadd=0.d0              
              ivdw=ivdw+1
              
              eadd=steadd(ivdw)

              denprd=twopi*(dble(numtyp(i))*dble(numtyp(j))-
     x          dble(numfrz(i))*dble(numfrz(j)))/volm**2
              delrc=delrc+volm*denprd*eadd
              
            enddo
            
          enddo
          
        endif
        
      endif
    
      return
      end

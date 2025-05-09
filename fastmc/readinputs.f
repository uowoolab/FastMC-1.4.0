c Program to read field file ripped from define_system_module.f
      module readinputs

      use parse_module
      use utility_pack

      integer nummls,numsit 

      contains

      subroutine readfield
     &(idnode,ntpvdw,maxvdw,ntpatm,ntpmls,ntpguest,
     &ntpfram,totatm,rvdw,dlrpot,engunit,maxguest,maxatm)
c*************************************************************************
c
c     Subroutine to read the FIELD file
c
c*************************************************************************
      implicit none
      
      logical loop,loop2,safe,atmchk,blank
      logical lunits,lguest
      integer itmols,msite,jsite,ntpatm,totatm,maxvdw
      integer ifrz,irept,nrept,ksite,isite,ntpmls
      integer lsite,natms,maxguest,maxatm
      integer i,n,j,k,idum,idnode,ntpvdw,ntpguest,ntpfram
      character*8 junk
      character*8 atom1
      character*1 message(80)
      real(8) weight,charge,engunit,rvdw,dlrpot,sumchg
      real(8) comx,comy,comz,mass,gpress
      natms=0
      ntpguest=0
      ntpfram=0
      loop=.true.
      lguest=.false.
      safe=.true.
      blank=.true.
      lunits=.false.
      engunit=1.d0
      totatm=0

      if(idnode.eq.0) open(nfield,file='FIELD',status='old')
 
c     allocate guest pressures
      do i=1,ntpguest
        gstpress(i)=-1.d0
      enddo

      call getrec(safe,idnode,nfield)
      if(.not.safe)call abort_field_read(1,idnode,nfield)

      do while(loop)
         call getrec(safe,idnode,nfield)
         if(.not.safe)call abort_field_read(1,idnode,nfield)

c        convert to lowercase and strip leading blanks

         call lowcase(record, lenrec)
         call strip(record, lenrec) 

c         if (.not.safe)then
c           loop=.false.
c         endif

         if (findstring('units',record,idum))then
           lunits=.true.
           do i=6,lenrec
             if(record(i).ne.' ')blank=.false.
           enddo

           if (blank)then
              if(idnode.eq.0)write(nrite,'(a)')
     &           "energy units = internal"
            elseif(findstring('ev',record,idum))then
              engunit=9648.530821d0
              if(idnode.eq.0)write(nrite,'(a)')
     &           "energy units = electron volts"
            elseif(findstring('kcal',record,idum))then
              engunit=418.4d0
              if(idnode.eq.0)write(nrite,'(a)')
     &           "energy units = kcal"
            elseif(findstring('kj',record,idum))then
              engunit=1.d2
              if(idnode.eq.0)write(nrite,'(a)')
     &           "energy units = kj"
            elseif(findstring('k',record,idum))then
              engunit=boltz
              if(idnode.eq.0)write(nrite,'(a)')
     &           "energy units = boltzmann"
            elseif(findstring('internal',record,idum))then
              if(idnode.eq.0)write(nrite,'(a)')
     &           "energy units = internal"
            endif

c        van der waals specification
   
         elseif (findstring('vdw',record,idum))then
           call define_van_der_waals
     &     (idnode,safe,ntpatm,dlrpot,rvdw,maxvdw,ntpvdw,engunit)

c        molecular specification

         elseif (findstring('molecu',record,idum))then
           ntpmls=intstr(record,lenrec,idum)
           
           if(idnode.eq.0)
     &       write(nrite,"(/,'number of molecular types:',3x,i4)")
     &       ntpmls           

c        initialize total system charge
          
           sumchg=0.d0

c        initialize number of sites counter
    
           nsite=0

           do n=1,ntpmls

             if(idnode.eq.0)
     &         write(nrite,"(/,2x,'molecular species type:',3x,i4)")
     &         n

c        get name of molecular species
             
             call getrec(safe,idnode,nfield)
             if(.not.safe)call abort_field_read(1,idnode,nfield)
             call lowcase(record,lenrec)
c        break if commented out
             if(findstring('&guest',record,idum))then
               call getword(junk,record,6,lenrec)
               lguest=.true.
               ntpguest=ntpguest+1
               locguest(ntpguest)=n
             else
               ntpfram=ntpfram+1
               locfram(ntpfram)=n
             endif
             call strip(record,lenrec)
             call copystring(record,molnam(1,n),40)
             if(idnode.eq.0)
     &         write(nrite,"(2x,'name of species:',3x,40a1)")
     &         (molnam(i,n),i=1,40)

            
             loop2=.true.
             do while(loop2)
               call getrec(safe,idnode,nfield)
               if(.not.safe)call abort_field_read(1,idnode,nfield)

               call lowcase(record, lenrec)
               call strip(record, lenrec)

               ksite=0

               if (findstring('nummol',record,idum))then
                  nummls=intstr(record,lenrec,idum)
                  nummols(n)=nummls
                  if(idnode.eq.0)
     &              write(nrite,"(2x,'number of molecules:',
     &              3x,i4)") nummls
                
               elseif (findstring('atoms',record,idum))then
c        read in atom name, site number, mass, charge,  etc
                  numsit=intstr(record,lenrec,idum)
                  numatoms(n)=numsit

                  if(idnode.eq.0)then
                    write(nrite,"(2x,'number of atoms:',3x,i4)")
     &                 numsit
                    write(nrite,"(2x,'atomic characteristics:',
     &            /,21x,' site',5x,'name',10x,'mass',8x,'charge',
     &            4x,'repeat',4x,'freeze')")
                  endif

                  do isite=1,numsit
                    if(ksite.lt.numatoms(n))then

                      call getrec(safe,idnode,nfield)
                      if(.not.safe)call 
     &                abort_field_read(1,idnode,nfield)
                      call copystring(record,message,80)
                      call getword(atom1,record,8,lenrec)
                      weight=dblstr(record,lenrec,idum)
                      charge=dblstr(record,lenrec,idum)
                  
                      if(lguest)then
                       nrept=1
                       guestx(ntpguest,isite)=dblstr(record,
     &                        lenrec,idum)
                       guesty(ntpguest,isite)=dblstr(record,
     &                        lenrec,idum)
                       guestz(ntpguest,isite)=dblstr(record,
     &                        lenrec,idum)
                       ifrz=0
                      else
                       nrept=intstr(record,lenrec,idum)
                       ifrz=intstr(record,lenrec,idum)
                      endif
                      if(nrept.eq.0)nrept=1

                      do irept=1,nrept
c                       nsite iterates over all atom sites 
                        nsite=nsite+1
 
                        ksite=ksite+1
                        atmname(n,ksite)=atom1
                        atmwght(n,ksite)=weight
                        atmchg(n,ksite)=charge
                        lfzsite(n,ksite)=ifrz
                      enddo
                      if(idnode.eq.0)then
                        write(nrite,"(21x,i5,5x,a8,2f12.5,3i10)")
     &                   ksite,atom1,weight,charge,nrept,ifrz
                      endif
                  
c          establish a list of unique atom types (for vdw calc)                    
                      atmchk=.true.
                    
                      do jsite=1,ntpatm
                        if (atom1.eq.unqatm(jsite))then
                          atmchk=.false.
                          do irept=ksite,ksite-nrept+1,-1
                            ltpsit(n,irept)=jsite
                          enddo

                        endif
                      enddo
                      if (atmchk)then
                        ntpatm=ntpatm+1
                        unqatm(ntpatm)=atom1
                        do irept=ksite,ksite-nrept+1,-1
                          ltpsit(n,irept)=ntpatm
                        enddo
                      endif 
                    endif
                  enddo

               elseif (findstring('finish',record,idum))then
                  loop2=.false.

                  natms=natms+nummols(n)*numatoms(n)
                  if(natms.gt.maxatm+maxguest)call error(idnode,75)

                  if(lguest)then
c                 calculate com of the guest and subtract
c                 so the reference atoms are centred around
c                 the com.
                    mass=0.d0
                    comx=0.d0
                    comy=0.d0
                    comz=0.d0
                    lguest=.false.
                    do i=1,numsit
                      comx=comx+guestx(ntpguest,i)*atmwght(n,i)
                      comy=comy+guesty(ntpguest,i)*atmwght(n,i)
                      comz=comz+guestz(ntpguest,i)*atmwght(n,i)
                      mass=mass+atmwght(n,i)
                    enddo
                    comx=comx/mass
                    comy=comy/mass
                    comz=comz/mass
                    do i=1,numsit
                      guestx(ntpguest,i)=guestx(ntpguest,i)-comx
                      guesty(ntpguest,i)=guesty(ntpguest,i)-comy
                      guestz(ntpguest,i)=guestz(ntpguest,i)-comz
                    enddo
                    if(idnode.eq.0)then
                      write(nrite,"(/,2x,'initial guest pos:')")
                      do i=1,numsit
                        write(nrite,"(21x,'site:',2x,a4,5x,3f12.6)")
     &                   atmname(n,i),guestx(ntpguest,i),
     &                   guesty(ntpguest,i),guestz(ntpguest,i)
                      enddo
                    endif
                  endif
               endif
             enddo
             totatm=totatm+numsit*nummls 
           enddo 
         elseif(findstring('close',record,idum))then
           loop=.false.

         else
           if(idnode.eq.0)write(nrite,'(100a)')record
           call abort_field_read(2,idnode,nfield) 

         endif
      enddo

c     check system charge
      do itmols=1,ntpmls
        do lsite=1,numatoms(itmols)  
          sumchg=sumchg+dble(nummols(itmols)*atmchg(itmols,lsite)) 
        enddo
      enddo

      if(abs(sumchg).gt.1.0d-6)then
          if(idnode.eq.0)write(nrite,'(/,1x,a,f12.6,a)')
     &  ' *** warning - total system charge:',sumchg,' ***'
      endif

      if(idnode.eq.0)close(nfield)

      return
      end subroutine readfield

      subroutine readconfig
     &(idnode,mxnode,imcon,cfgname,levcfg,
     &ntpmls,maxmls,totatm,volm,rcut,celprp)
c*************************************************************************
c
c     Subroutine to read the CONFIG file
c
c*************************************************************************
      implicit none

      character*1 atname(8),cfgname(80)
      logical safe
      integer icfg,imcon,idnode,indatm,maxmls
      integer k,l,m,totatm,i,ntpmls
      integer idum,levcfg,mxnode
      real(8) xcoord,ycoord,zcoord,junk,volm,axx,test,rt3
      real(8) rcut,width
      real(8), dimension(10) :: celprp
     
      safe=.true.
      if(idnode.eq.0)open(nconfig,file='CONFIG',status='old')
c     read header info
      call getrec(safe,idnode,nconfig)
      call copystring(record,cfgname,80)

      call getrec(safe,idnode,nconfig)

      levcfg=intstr(record,lenrec,idum) 

      if(imcon.eq.0)then
        if(idnode.eq.0)write(nrite,'(3x, "WARNING - no periodic 
     &boundaries requested!")')
        do i=1,9
          cell(i)=0.d0
        enddo
        volum=0.d0
      else
        call getrec(safe,idnode,nconfig)
        if(.not.safe)call abort_config_read(1,idnode,nconfig)
        cell(1)=dblstr(record,lenrec,idum)
        cell(2)=dblstr(record,lenrec,idum)
        cell(3)=dblstr(record,lenrec,idum)
        call getrec(safe,idnode,nconfig)
        if(.not.safe)call abort_config_read(1,idnode,nconfig)
        cell(4)=dblstr(record,lenrec,idum)
        cell(5)=dblstr(record,lenrec,idum)
        cell(6)=dblstr(record,lenrec,idum)
        call getrec(safe,idnode,nconfig)
        if(.not.safe)call abort_config_read(1,idnode,nconfig)
        cell(7)=dblstr(record,lenrec,idum)
        cell(8)=dblstr(record,lenrec,idum)
        cell(9)=dblstr(record,lenrec,idum)
      endif
c     read atomic coordinates
      indatm=0
      safe=.true.

      do k=1,maxmls
        indatm=0
        do l=1,nummols(k)
          do m=1,numatoms(k)

            indatm=indatm+1

            framwkxxx(k,indatm)=0.d0
            framwkyyy(k,indatm)=0.d0
            framwkzzz(k,indatm)=0.d0

            if(idnode.eq.0)then
              if(levcfg.eq.0)then
                read(nconfig,'(8a1)',end=100)atname
                read(nconfig,'(3f20.0)',end=100)xcoord,ycoord,zcoord
  
              elseif(levcfg.eq.1)then
                read(nconfig,'(8a1)',end=100)atname
                read(nconfig,'(3f20.0)',end=100)xcoord,ycoord,zcoord
                read(nconfig,'(3f20.0)',end=100)junk,junk,junk
            
              elseif(levcfg.eq.2)then
                read(nconfig,'(8a1)',end=100)atname
                read(nconfig,'(3f20.0)',end=100)xcoord,ycoord,zcoord
                read(nconfig,'(3f20.0)',end=100)junk,junk,junk
                read(nconfig,'(3f20.0)',end=100)junk,junk,junk
              endif


              call strip(atname,8)

              if (atmname(k,m).eq.mkwd8(atname))then
                 framwkxxx(k,indatm)=xcoord
                 framwkyyy(k,indatm)=ycoord
                 framwkzzz(k,indatm)=zcoord 

              else
                 write(nrite,"(/,/,'unidentified atom label :',8a1,
     &             ': atom number ',i5)")atname,indatm
                 safe=.false.
              endif

            endif

            call gstate(safe)
            if(.not.safe)call error(idnode,25)
           
          enddo
        enddo
      enddo

c     sum coordinate arrays across all nodes
      if(mxnode.gt.1)then
        call gd2dsum(framwkxxx,ntpmls,nummols,numatoms,totatm)
        call gd2dsum(framwkyyy,ntpmls,nummols,numatoms,totatm)
        call gd2dsum(framwkzzz,ntpmls,nummols,numatoms,totatm)
      endif

c     check integrity of cell vectors
      if((imcon.eq.1).or.(imcon.eq.4).or.(imcon.eq.5))then
         axx=(abs(cell(1))+abs(cell(5)))/2.d0
         test=1.d-8*axx
         if(abs(cell(1)-axx).gt.test)call error(idnode,410)
         if(abs(cell(5)-axx).gt.test)call error(idnode,410)
         if(imcon.eq.5)then
           if(abs(cell(9)-axx*sqrt(2.d0)).gt.test)
     &       call error(idnode,410)
         else
           if(abs(cell(9)-axx).gt.test)call error(idnode,410)
         endif
      endif

      if(imcon.eq.7)then
        rt3=sqrt(3.d0)
        if(abs(cell(1)-rt3*cell(5)).ge.1.d-6)
     &    call error(idnode,410)
      endif

      if(imcon.eq.6)then
        if(abs(cell(3)).gt.1.d-10) call error(idnode,410)
        if(abs(cell(6)).gt.1.d-10) call error(idnode,410)
        if(abs(cell(7)).gt.1.d-10) call error(idnode,410)
        if(abs(cell(8)).gt.1.d-10) call error(idnode,410)
      endif

      if((imcon.eq.1).or.(imcon.eq.2).or.(imcon.eq.4).or.
     & (imcon.eq.5).or.(imcon.eq.7))then

        if(abs(cell(2)).gt.1.d-10) call error(idnode,410)
        if(abs(cell(3)).gt.1.d-10) call error(idnode,410)
        if(abs(cell(4)).gt.1.d-10) call error(idnode,410)
        if(abs(cell(6)).gt.1.d-10) call error(idnode,410)
        if(abs(cell(7)).gt.1.d-10) call error(idnode,410)
        if(abs(cell(8)).gt.1.d-10) call error(idnode,410)

      endif

      call dcell(cell,celprp)
      if(imcon.eq.0)then

        volm=0.d0

      elseif(imcon.eq.4)then

        volm=0.5d0*celprp(10)

      elseif(imcon.eq.5)then

        volm=0.5d0*celprp(10)

      elseif(imcon.eq.7)then

        volm=0.5d0*celprp(10)

      else

        volm=celprp(10)

      endif

      if(idnode.eq.0)then
        write(nrite,"(/,'simulation cell vectors:')")
        write(nrite,"(21x,3f12.6)")cell
        write(nrite,
     &     "(/,'system volume:',/,21x,1p,g22.12)")volm
      endif

c     check value of cutoff and reset if necessary

      if(imcon.gt.0)then
        
        width=min(celprp(7),celprp(8),celprp(9))/2.d0
        if(imcon.eq.4)width=sqrt(3.d0)*cell(1)/4.d0
        if(imcon.eq.5)width=cell(1)/2.d0
        if(imcon.eq.6)width=min(celprp(7),celprp(8))/2.d0

c     halt program if potential cutoff exceeds cell width
   
        if(rcut.gt.width) call error(idnode,95)
      endif

      if(idnode.eq.0)close(nconfig)

      return

100   call abort_config_read(2,idnode,nconfig)

      end subroutine readconfig

      subroutine readcontrol(idnode,lspe,temp,
     &ljob,mcsteps,eqsteps,celprp,ntpguest,lrestart,laccsample,lnumg,
     &nnumg,nhis,nwind,mcinsf, mcdelf, mcdisf, mcjmpf, mcflxf, mcswpf, 
     &mctraf, mcrotf, mcmjpf, disp_ratio, tran_ratio, rota_ratio, lfuga,
     &maxguest,maxatm,desorb)
c*************************************************************************
c
c     Subroutine to read the CONTROL file
c
c*************************************************************************
      implicit none
      
      character*1 directive(lenrec)
      character*1 sysname(80)

      logical safe,loop,loop2,ltemp,lewald,lspe,ljob,lnumg,lfuga,lwind
      logical lmcsteps,leqsteps,lprob,loop3,rprob,lrestart,laccsample
      logical mvspefid(ntpguest), delrspefid(ntpguest),lmaxg,lmaxf
      logical disp_ratiospefid(ntpguest), tran_ratiospefid(ntpguest)
      real(8) drdf,dzdn,zlen,temp,rcut
      integer idnode,idum,keyres,eqsteps,mcsteps,idguest,nhis,nnumg
      integer n,iprob,i,j,ntpsite,ntpguest,ngst,cprob,nwind,maxguest
      integer maxatm,desorb
      real(8) disp_ratio(ntpguest), tran_ratio(ntpguest), rota_ratio
      real(8), dimension(10) ::  celprp
      real(8) mcinsf(ntpguest), mcdelf(ntpguest), mcdisf(ntpguest)
      real(8) mcswpf(ntpguest), mcflxf(ntpguest), mcjmpf(ntpguest)
      real(8) mctraf(ntpguest), mcrotf(ntpguest), delr(ntpguest)
      real(8) mcmjpf(ntpguest)
      real(8) tempmcinsf, tempmcdelf, tempmcdisf, tempmcswpf
      real(8) tempmcflxf, tempmcjmpf, tempmctraf, tempmcrotf
      real(8) tempmcmjpf
      real(8) tempdelr,tempdisp_ratio,temptran_ratio

      data ltemp/.false./,lprob/.false./,loop3/.false./,lmaxg/.false./
      data lmcsteps/.false./,leqsteps/.false./,lwind/.false./

      iprob=0
      ngst=0 
      nhis=0
      desorb=0
      temp=0.d0
      drdf=0.05d0
      dzdn=0.05d0
      zlen=0.d0
      loop=.true.
      loop2=.false.
      lewald=.false.
c     allocate guest pressures
c     allocate guest mole fractions
      do i=1,ntpguest
        gstpress(i)=0.d0
        gstmolfract(i)=0.d0
        adspres(i)=0.d0
        adsguen(i)=0.d0
      enddo
      total_pressure = 0.d0
      packf=0.d0
      despre=0.d0
      if(idnode.eq.0)open (ncontrol,file='CONTROL',status='old')

c     allocate guest pressures
      mcinsf(:)=0.d0
      mcdelf(:)=0.d0
      mcdisf(:)=0.d0
      mcjmpf(:)=0.d0
      mcflxf(:)=0.d0
      mcswpf(:)=0.d0
      mctraf(:)=0.d0
      mcrotf(:)=0.d0
      mcmjpf(:)=0.d0
      mvspefid(:)=.false.
      delrspefid(:)=.false.
      disp_ratiospefid(:)=.false.
      tran_ratiospefid(:)=.false.
 
      call getrec(safe,idnode,ncontrol)
      if(.not.safe)call abort_control_read(1,idnode,ncontrol)
      call copystring(record,sysname,80)
      do while(loop)
        call getrec(safe,idnode,ncontrol)
        if(.not.safe)call abort_control_read(1,idnode,ncontrol)
        call lowcase(record,lenrec)
        call strip(record,lenrec)
        call copystring(record,directive,lenrec) 

        if(record(1).eq.'#'.or.record(1).eq.' ')then
c       record is commented out

        else if(findstring('restart',directive,idum))then
           keyres=1
           if(idnode.eq.0)write(nrite,'(/,3x,a)')
     &'***restart requested***'
           lrestart=.true.
        elseif (findstring('temp',directive,idum))then
           temp=dblstr(directive,lenrec,idum)
           ltemp=.true.
        elseif (findstring('sampling', directive, idum))then
           laccsample=.true.
           write(nrite,"(/,3x,'warning! Sampling accepted',
     &' steps only, this is not true Metropolis importance ',
     &'sampling. This simulation is CURSED!!!!')")
c       guest related stuff
        elseif (findstring('&guest',directive,idum))then
           ngst=ngst+1
           idguest=intstr(directive,lenrec,idum)
           loop2=.true.
           rprob=.false.
           do while(loop2)
             call getrec(safe,idnode,ncontrol)
             
             if(.not.safe)call abort_control_read(1,idnode,ncontrol)
             call lowcase(record,lenrec)
             call strip(record,lenrec)
             if(record(1).eq.'#'.or.record(1).eq.' ')then
c            record is commented out
             elseif(findstring('fract',record,idum))then
               gstmolfract(idguest)=dblstr(record,lenrec,idum)
             elseif(findstring('press',record,idum))then
               gstpress(idguest)=dblstr(record,lenrec,idum)*1.d5
c               if(gstpress(idguest).le.1.d-12)call error(idnode,387)
               if(idnode.eq.0)write(nrite,"(/,'Guest ',i4,
     &' pressure: ',f14.4,' bar')")idguest,(gstpress(idguest)/100000)
             elseif(findstring('probability',record,idum))then
               lprob=.true.
               rprob=.true.
               n=intstr(record,lenrec,idum)
               nprob(idguest)=n
               loop3=.true.
               cprob=0
               do while(loop3)
                 call getrec(safe,idnode,ncontrol)
                 if(.not.safe)call abort_control_read(1,idnode,ncontrol)
                 call strip(record,lenrec)
                 if(record(1).eq.'#'.or.record(1).eq.' ')then
                 else
                   iprob=iprob+1
                   cprob=cprob+1 
                   ntpsite=intstr(record,lenrec,idum)
                   nprobsites(iprob)=ntpsite
                   do j=1,ntpsite
                     lprobsites(iprob,j)=intstr(record,lenrec,idum)
                   enddo
                 endif
                 if(cprob.eq.n)loop3=.false.
               enddo
             ! For individual guests (priority)
             elseif(findstring('move',record,idum))then
               mvspefid(idguest) = .true.
               if (findstring('ins',record,idum))then
                 mcinsf(idguest) = dblstr(record,lenrec,idum)
               elseif (findstring('del',record,idum))then
                 mcdelf(idguest) = dblstr(record,lenrec,idum)
               elseif (findstring('dis',record,idum))then
                 mcdisf(idguest) = dblstr(record,lenrec,idum)
               elseif (findstring('jum',record,idum))then
                 mcjmpf(idguest) = dblstr(record,lenrec,idum)
               elseif (findstring('fle',record,idum))then
                 mcflxf(idguest) = dblstr(record,lenrec,idum)
               elseif (findstring('swa',record,idum))then
                 mcswpf(idguest) = dblstr(record,lenrec,idum)
               elseif (findstring('tra',record,idum))then
                 mctraf(idguest) = dblstr(record,lenrec,idum)
               elseif (findstring('rot',record,idum))then
                 mcrotf(idguest) = dblstr(record,lenrec,idum)
               elseif (findstring('mul',record,idum))then
                 mcmjpf(idguest) = dblstr(record,lenrec,idum)
               else
                 if(idnode.eq.0)
     &write(nrite,"(/, 'ignoring unknown move for guest', i4)") idguest
               endif
             elseif(findstring('delr',directive,idum))then
               delrspefid(idguest) = .true.
               delr(idguest) = dblstr(record,lenrec,idum)
             elseif(findstring('accep',record,idum))then
               if (findstring('dis',record,idum))then
                 disp_ratiospefid(idguest) = .true.
                 disp_ratio(idguest) = dblstr(record,lenrec,idum)
               elseif (findstring('tra',record,idum))then
                 tran_ratiospefid(idguest) = .true.
                 tran_ratio(idguest) = dblstr(record,lenrec,idum)
               endif
             elseif(findstring('&end',record,idum))then
               loop2=.false.
             endif    
           enddo
           if(.not.rprob)nprob(idguest)=0
        elseif (findstring('steps',directive,idum))then
           mcsteps=intstr(directive,lenrec,idum)
           lmcsteps=.true.
        elseif (findstring('pressure',directive,idum).or.
     &findstring('total press',directive,idum))then
           total_pressure=dblstr(directive,lenrec,idum)*1.d5
        elseif (findstring('numg',directive,idum))then
           lnumg=.true.
           nnumg=intstr(directive,lenrec,idum)
c          set default writing numguests.out to 1000 steps
           if(nnumg.le.0)nnumg=1000
        elseif (findstring('fuga',directive,idum))then
           lfuga=.true.
        elseif (findstring('hist',directive,idum))then
           nhis=intstr(directive,lenrec,idum)
        elseif (findstring('equilibr',directive,idum))then
           eqsteps=intstr(directive,lenrec,idum)
           if(eqsteps.ne.0)leqsteps=.true.
        elseif (findstring('uvt',directive,idum))then
          if(idnode.eq.0)write(nrite, 
     &"(/,'gcmc requested')")
        elseif (findstring('single point',directive,idum))then
          if(idnode.eq.0)write(nrite,
     &"(/,'single point energy calculation requested ')")
          lspe=.true.
        elseif (findstring('jobcontrol',directive,idum))then
c           ljob=.true.
c           Set ljob to false, debugging option
           ljob=.false.
        elseif (findstring('cutoff',directive,idum))then
           rcut=dblstr(directive,lenrec,idum)
        elseif (findstring('delr',directive,idum))then
          tempdelr = dblstr(directive,lenrec,idum)
          do i = 1, ntpguest
            if (delrspefid(i).eqv..false.) delr(i) = tempdelr
          enddo
        elseif (findstring('ewald',directive,idum))then
           
        elseif (findstring('averaging window',directive,idum))then
           nwind=intstr(directive,lenrec,idum)
           lwind=.true.
        elseif (findstring('max guest atoms',directive,idum))then
           maxguest=intstr(directive,lenrec,idum)
           lmaxg=.true.
        elseif (findstring('max framework atoms',directive,idum))then
           maxatm=intstr(directive,lenrec,idum)
           lmaxf=.true.
        elseif (findstring('desorb',directive,idum))then
           desorb=intstr(directive,lenrec,idum)
        elseif (findstring('adsorb press',directive,idum))then
           do i=1,ntpguest
              adspres(i)=dblstr(directive,lenrec,idum)
           enddo
        elseif (findstring('adsorb guests',directive,idum))then
           do i=1,ntpguest
              adsguen(i)=dblstr(directive,lenrec,idum)
           enddo
        elseif (findstring('packing fraction',directive,idum))then
           packf = dblstr(directive,lenrec,idum)
        elseif (findstring('finishment',directive,idum))then
           despre = dblstr(directive,lenrec,idum)
        ! For all guests not specified
        elseif (findstring('move',directive,idum))then
          if (findstring('ins',directive,idum))then
            tempmcinsf = dblstr(directive,lenrec,idum)
            do i = 1, ntpguest
              if (mvspefid(i).eqv..false.) mcinsf(i) = tempmcinsf
            enddo
          elseif (findstring('del',directive,idum))then
            tempmcdelf = dblstr(directive,lenrec,idum)
            do i = 1, ntpguest
              if (mvspefid(i).eqv..false.) mcdelf(i) = tempmcdelf
            enddo
          elseif (findstring('dis',directive,idum))then
            tempmcdisf = dblstr(directive,lenrec,idum)
            do i = 1, ntpguest
              if (mvspefid(i).eqv..false.) mcdisf(i) = tempmcdisf
            enddo
          elseif (findstring('jum',directive,idum))then
            tempmcjmpf = dblstr(directive,lenrec,idum)
            do i = 1, ntpguest
              if (mvspefid(i).eqv..false.) mcjmpf(i) = tempmcjmpf
            enddo
          elseif (findstring('fle',directive,idum))then
            tempmcflxf = dblstr(directive,lenrec,idum)
            do i = 1, ntpguest
              if (mvspefid(i).eqv..false.) mcflxf(i) = tempmcflxf
            enddo
          elseif (findstring('swa',directive,idum))then
            tempmcswpf = dblstr(directive,lenrec,idum)
            do i = 1, ntpguest
              if (mvspefid(i).eqv..false.) mcswpf(i) = tempmcswpf
            enddo
          elseif (findstring('tra',directive,idum))then
            tempmctraf = dblstr(directive,lenrec,idum)
            do i = 1, ntpguest
              if (mvspefid(i).eqv..false.) mctraf(i) = tempmctraf
            enddo
          elseif (findstring('rot',directive,idum))then
            tempmcrotf = dblstr(directive,lenrec,idum)
            do i = 1, ntpguest
              if (mvspefid(i).eqv..false.) mcrotf(i) = tempmcrotf
            enddo
          elseif (findstring('mul',directive,idum))then
            tempmcmjpf = dblstr(directive,lenrec,idum)
            do i = 1, ntpguest
              if (mvspefid(i).eqv..false.) mcmjpf(i) = tempmcmjpf
            enddo
          else
            if(idnode.eq.0)then
              write(nrite,"(/, 'ignoring unknown general move')")
            endif
          endif
        elseif (findstring('accep',directive,idum))then
          if (findstring('dis',directive,idum))then
            tempdisp_ratio = dblstr(directive,lenrec,idum)
            do i = 1, ntpguest
              if (disp_ratiospefid(i).eqv..false.)then
                disp_ratio(i) = tempdisp_ratio
              endif
            enddo
          elseif (findstring('tra',directive,idum))then
            temptran_ratio = dblstr(directive,lenrec,idum)
            do i = 1, ntpguest
              if (tran_ratiospefid(i).eqv..false.)then
                tran_ratio(i) = temptran_ratio
              endif
            enddo
          elseif (findstring('rot',directive,idum))then
            rota_ratio = dblstr(directive,lenrec,idum)
          endif
        elseif (findstring('grid',directive,idum))then
c 'grid' is already read in initscan
        elseif (findstring('finish',directive,idum))then
           loop=.false.
        else
           if(idnode.eq.0)write(nrite,"(/,/,100a1)")record
           call error(idnode,3)
        endif

      enddo
      
      if(ljob)then

        if(idnode.eq.0)then
           write(nrite,'(/,a)')
     &"Equilibrium steps will be determined by entering
     & 'START AVERAGING' in the file jobcontrol.in"
        endif
        eqsteps=0
      else
        if(.not.leqsteps)then
          eqsteps=0
          if(idnode.eq.0)write(nrite,'(/,a)')
     &'*** warning no equilibrium steps specified in CONTROL file,
     & starting production averaging ***'
        endif
      endif
      if(.not.lnumg)then
        lnumg=.true.
        if (eqsteps.gt.1000)then
            nnumg=eqsteps/1000
        else
            nnumg=1
        endif
        if(idnode.eq.0)then
          write(nrite,'(/,a,i10,a)')
     &"*** warning no 'numguests' found in CONTROL, writing to 
     &numguests.out every ",nnumg," gcmc steps***"
        endif
      else
        if(idnode.eq.0)write(nrite,'(/,a31,i6,a12)')
     &"writing to numguests.out every ",nnumg," gcmc 
     &steps."
      endif
      if(nhis.ne.0)then
        if(idnode.eq.0)write(nrite,'(/,a24,i5,a12)')
     &"writing to his.xyz every",abs(nhis)," gcmc steps."
      endif
      if(.not.lmcsteps)then
        mcsteps=1
        lspe=.true.
        if(idnode.eq.0)write(nrite,'(/,a)')
     &'*** warning no gcmc production steps 
     &specified in CONTROL file, assuming single point calculation ***'
      endif
      if(.not.lwind)then
          nwind=100000
          if(idnode.eq.0)write(nrite, '(/a,i10,a)')
     &'*** warning the averaging window was not specified in the
     & CONTROL file, defaulting to ',nwind,' steps per window.'
      endif
      if(.not.lmaxg)then
          lmaxg=.true.
          if(idnode.eq.0)write(nrite, '(/a,i10,a)')
     &'*** warning the max amount of guest atoms  not specified in the
     & CONTROL file, defaulting to ',maxguest,' atoms in GCMC.'
      endif
      if(.not.lmaxg)then
          lmaxf=.true.
          if(idnode.eq.0)write(nrite, '(/a,i10,a)')
     &'*** warning max amount of framework atoms  not specified in the
     & CONTROL file, defaulting to ',maxatm,' atoms in GCMC.'
      endif
      if((ngst.ne.ntpguest))call error(idnode,2312)      
      if(.not.ltemp)then
        call error(idnode,380)
      endif

      if(idnode.eq.0)close(ncontrol)
 
      return
      end subroutine

      end module readinputs


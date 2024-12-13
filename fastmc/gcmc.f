************************************************************
c      main gcmc program                                    *
c************************************************************
c      personal reference                                   *
c      imcon                                                *
c      0  -  no periodic boundaries                         *
c      1  -  cubic boundary conditions  (b.c.)              *
c      2  -  orthorhombic b.c.                              *
c      3  -  parallelpiped b.c.                             *
c      4  -  truncated octahedral b.c.                      *
c      5  -  rhombic dodecahedral b.c                       *
c      6  -  x-y parallelogram b.c. no periodicity in z     *
c      7  -  hexagonal prism b.c.                           *
c************************************************************


      use utility_pack
      use readinputs
      use vdw_module


      implicit none

      character*1 cfgname(80)      
      character*8 outdir,localdir
      character*2 mnth
      character*7 debug
      character(8):: date
      character(10) :: time
      character(5) :: zone
      character*25 outfile
      character*80 command
      character*25 outfile2
      character*21 outfile3
      character*18 outfile4
      character*70 string
      logical lgchk,lspe,ljob,lprob,fram,newmol
      logical lfuga,file_exists
      logical insert,delete,displace,lrestart,laccsample
      logical lnorm_bin
      logical jump, flex, swap
      logical tran, rota, multijump
      logical accepted,production,jobsafe,lnumg
      logical tick_tock_cycles(5)
      integer, dimension(1) :: myseed
      integer nfo(8)
      integer accept_ins,accept_del,totaccept,nnumg,nhis
      integer accept_jump, accept_flex, accept_swap
      integer accept_rota, x, y, z
      integer ins_count,del_count,buffer,idum,levcfg,nstat
      integer jump_count, flex_count, swap_count
      integer rota_count, ichain, inodes
      integer totatm,jatm,isite,imol,at,atmadd,iatm,ntprob,ntpsite
      integer itatm,indx,newld,gcmccount,prodcount,globalprod,ii
      integer ibuff,iprob,cprob,totfram,prevnodes,nwind,rollstat,windn
      integer n,k,p,i,j,ik,ka,jj,kk,l,mm,ierr,ksite,ntpatm,nang,gridsize
      integer kmax1,kmax2,kmax3,ntpvdw,maxvdw,maxmls,mxnode,idnode
      integer ntpfram,randchoice,nmols,ngsts,molecules,totsteps,nummove
      integer mxatm,imcon,keyfce,mxatyp,sittyp,mcsteps,eqsteps,vstat
      integer iguest,ntpguest,ntpmls,natms,mol,maxanglegrid,rollcount
      integer ins_rollcount,maxguest,framatm,maxatm,desorb
      integer ngrida,ngridb,ngridc,nguests,mxewld,mxebuf,istat,avcount
      integer totalguests,globalnguests,globalsteps,mvloop
      integer, allocatable :: prevmultijump(:)
      real(8), allocatable :: guestmol(:),guesteng(:)
      real(8), allocatable :: gridbuff(:),dqp(:)
      real(8), dimension(10) :: celprp
      real(8), dimension(9) :: rcell
      real(8) vdweng,ewld1test,ewld1,engsicold,ang,hyp,norm,origenergy
      real(8) det,engcpe,engacc,engac1,drewd,epsq,statvolm,volm
      real(8) engsrp,rand,randmov,chg,rande,delta
      real(8) stdQst,stdCv,rotangle,elrc,delrc,junk,stdevcv,stdevq
      real(8) tzero,timelp,engunit,rvdw,press,temp,beta
      real(8) dlrpot,rcut,eps,alpha,init,gpress,initdelr
      real(8) ewld1eng,ewld2eng,alen,blen,clen,ecoul,evdw
      real(8) ecoulg,evdwg,engcomm,dummy,selectcount,changerate
      real(8) ewld2sum,vdwsum,ewld3sum,comx,comy,comz
      real(8) dmolecules,molecules2,energy2,Q_st,eng,C_v,cv_old
      real(8) weight,tw,spenergy,dlpeng,estep,thrd,twothrd,totpres
      real(8) E,aN,EN,E2,N2,H_const,avgH,stdH,othermol,stdS,temppres
      real(8) avgN,stdN,avgE,stdE,avgEN,stdEN,avgN2,stdN2,avgE2,stdE2
      real(8) griddim,grvol,randa,randb,randc,rand1,rand2,rand3,rand4
      real(8), dimension(9) :: iabc
c Cumulative move probabilities, remeber to zero and include in normalising
      real(8) rota_ratio,ij_select,presrate,presfrac,presleft
      real(8) rota_rotangle, dis_delr, dis_rotangle
      integer, dimension(3) :: gridfactor

      data lgchk/.true./,insert/.false./,delete/.false./,
     &displace/.false./,accepted/.false./,production/.false./
      data multijump/.false./
      data jump/.false./,flex/.false./,swap/.false./
      data tran/.false./,rota/.false./,laccsample/.false./
      data lspe/.false./,ljob/.false./,jobsafe/.true./,lrestart/.false./
     &,lnumg/.false./
      data lfuga/.true./ ! change to false
      data accept_ins,accept_del,totaccept/0,0,0/
      data accept_jump, accept_flex, accept_swap/0,0,0/
      data ins_count,del_count,gcmccount,prevnodes/0,0,0,0/
      data jump_count, flex_count, swap_count/0,0,0/
      data rota_count, accept_rota/0,0/

      integer, parameter, dimension(3) :: revision = (/1, 4, 0 /)
c     TODO(pboyd): include error checking for the number of guests
c     coinciding between the CONTROL file and the FIELD file.
      tw=0.d0
      ecoul=0.d0
      evdw=0.d0
      ntpatm=0
      mxatm=0
      mxgrid=0
      gridsize=0
      ibuff=0
      nnumg=1
      rollcount=0
      ins_rollcount=0
      avcount = 0
c     averaging window to calculate errors
      nwind=100000
      windn=5
      selectcount=0.d0
      othermol=0.d0

c Default target acceptance ratios of 0.5      
c      disp_ratio(:) = 0.5d0  ! Guest specific ratios can only be
c      tran_ratio(:) = 0.5d0  ! defaulted after allocation
      rota_ratio = 0.5d0

c scoping issues
      elrc = 0
      delrc = 0

      thrd=1.d0/3.d0
      twothrd=2.d0/3.d0

c default length of a side of a grid box in angstroms (in CONTROL)
      griddim=0.1d0
c default supecell folding for grid points (in CONTROL)
      gridfactor = (/ 1, 1, 1 /)
      
c     global number of production steps over all nodes
      globalprod=0
c     local production steps (after equilibrium)
      prodcount=0
c     local insertion steps (after equilibrium) - to keep track 
c     newld is the number of ewald points in reciprocal space
      newld=0
c     mcsteps = number of production steps, energies and molecules
c counted towards final averages
      mcsteps=1
c     eqsteps = number of equilibrium steps, energies and molecules
c ignored
      eqsteps=0
c when running mc cycles keep track of if history to help averaging
      tick_tock_cycles = .false.
      
c     initialize communications
      call initcomms()
      call gsync()
      call timchk(0,tzero)

c     determine processor identities
      call machine(idnode,mxnode)

c     open main output file.
      if(idnode.eq.0)then

        open(nrite,file='OUTPUT')
        write(nrite,
     &"(/,20x,'FastMC version ',i1,'.',i1,'.',i1,/,/)")revision
        call date_and_time(date,time,zone,nfo)
        mnth=date(5:6)
        
        write(nrite,
     &"('Started : ',9x,a2,3x,a9,1x,a4,3x,a2,':',a2,a6,' GMT',/)")
     &date(7:8),month(mnth),date(1:4),time(1:2),time(3:4),zone
        write(nrite,"('Running on ',i4,' nodes',/,/)")mxnode
      endif
      call initscan
     &(idnode,imcon,volm,keyfce,rcut,eps,alpha,kmax1,kmax2,kmax3,lprob,
     & initdelr,rvdw,ntpguest,ntprob,ntpsite,ntpvdw,maxmls,mxatm,mxatyp,
     & griddim,gridfactor,maxguest,maxatm,despre)

c     mxebuf is used to allocate the arrays ckcsum,ckssum
      mxebuf=(2*kmax1+1)*(2*kmax2+1)*(2*kmax3+1)-1
c     mxewld is used to allocate other ewald arrays.
      mxewld=mxatm

      maxvdw=max(ntpvdw,(mxatyp*(mxatyp+1))/2)
c      write(*,"('maxmls:',i6,' mxatm: ',i6,' mxatyp: ',i6,' volm: ',
c     & f9.3,' kmax2: ', i6, ' kmax3: ',i6,' mxebuf: ',i6,' mxewld: ',
c     & i6,' ntpguest: ',i6,' rcut: ',f6.3,' rvdw: ',f6.3,' delr: ',
c     & f6.3)")maxmls,mxatm,mxatyp,volm,kmax2,kmax3,mxebuf,mxewld,
c     & ntpguest,rcut,rvdw,delr


      call alloc_config_arrays
     & (idnode,mxnode,maxmls,mxatm,mxatyp,volm,kmax2,kmax3,
     &mxebuf,mxewld,ntpguest,rcut,rvdw,initdelr,maxguest)

      delr(:) = initdelr
      accept_disp(:)=0.d0
      accept_tran(:)=0.d0
      disp_count(:)=0.d0
      tran_count(:)=0.d0

c Default target acceptance ratios of 0.5      
      disp_ratio(:) = 0.5d0
      tran_ratio(:) = 0.5d0
      rota_ratio = 0.5d0

c     avgwindow will reset to 0's once the average has been carried
c     over to varwindow.
c     varwindow will be a rolling variance calculation of the windowed
c     averages 
c     currently the order is:
c     chainstats(1) = production gcmc step count
c     chainstats(2) = rolling average number of guests <N>
c     chainstats(3) = rolling average energy <E>
c     chainstats(4) = rolling average for energy*guests <EN> 
c     chainstats(5) = rolling average for (number of guests)^2 <N2> 
c     chainstats(6) = rolling average for (energy)^2 <E2>
c     chainstats(7) = rolling average for <exp(-E/kb/T)>
c     chainstats(8) = rolling stdev value for <N>
c     chainstats(9) = rolling stdev value for <E> 
c     chainstats(10) = rolling stdev value for <EN> 
c     chainstats(11) = rolling stdev value for <N2> 
c     chainstats(12) = rolling stdev value for <E2> 
c     chainstats(13) = rolling stdev value for <Q_st>
c     chainstats(14) = rolling stdev value for <C_v>
c     chainstats(15) = rolling stdev value for <exp(-E/kb/T)>

      call alloc_vdw_arrays(idnode,maxvdw)
      call readfield
     &(idnode,ntpvdw,maxvdw,ntpatm,ntpmls,ntpguest,
     &ntpfram,totatm,rvdw,dlrpot,engunit,maxguest,maxatm)

      call readconfig(idnode,mxnode,imcon,cfgname,levcfg,
     &ntpmls,maxmls,totatm,volm,rcut,celprp)
      framatm=totatm

c     volume reported in m^3 for statistical calculations      
      statvolm=volm*1d-30

      call invert(cell,rcell,det)

      if(lprob)then
c      calculate number of grid points in a,b,and c directions
c      calculate the volume of a grid point (differs from griddim^3)
c NB these are allocated before assigning to guests in readconfig
        ngrida=gridfactor(1)*ceiling(celprp(1)/(griddim*gridfactor(1)))
        ngridb=gridfactor(2)*ceiling(celprp(2)/(griddim*gridfactor(2)))
        ngridc=gridfactor(3)*ceiling(celprp(3)/(griddim*gridfactor(3)))
        gridsize=ngrida*ngridb*ngridc
        if(idnode.eq.0)
     &write(nrite,"(/,' Probability grid size:',
     &i8,i8,i8)")ngrida,ngridb,ngridc
        allocate(gridbuff(gridsize))
        call alloc_prob_arrays(idnode,ntpguest,ntpsite,ntprob,gridsize)
        grvol=celprp(1)/dble(ngrida)*celprp(2)/dble(ngridb)*
     &  celprp(3)/dble(ngridc)
      else
c       this is in case we run into allocation problems later on

        ntprob=0
        gridsize=1
        call alloc_prob_arrays(idnode,ntpguest,ntpsite,ntprob,gridsize)
      endif
      call readcontrol(idnode,lspe,temp,ljob,mcsteps,eqsteps,
     &celprp,ntpguest,lrestart,laccsample,lnumg,nnumg,nhis,nwind,
     &mcinsf, mcdelf, mcdisf, mcjmpf, mcflxf, mcswpf, mctraf, mcrotf,
     &mcmjpf,disp_ratio,tran_ratio,rota_ratio,lfuga,maxguest,maxatm,
     &desorb,lnorm_bin)

c Normalise the move frequencies (individually)
      if (desorb.eq.1)then
        allocate(dqp(ntpguest))
        presleft=0.d0
        do i=1,ntpguest
           dqp(i)=0.d0
        enddo
      endif
      do i = 1, ntpguest
        mcmvnorm(i) = mcinsf(i)+mcdelf(i)+mcdisf(i)+mcjmpf(i)+mcflxf(i)+
     &mcswpf(i)+mctraf(i)+mcrotf(i)+mcmjpf(i)
        if(mcmvnorm(i).eq.0)then
          if(idnode.eq.0)
     &  write(nrite,"(/,'No move frequencies specified, defaulting to '
     &  'Grand Canonical')")
          mcmvnorm(i) = 1
          mcinsf(i) = 1.d0/3.d0
          mcdelf(i) = 1.d0/3.d0
          mcdisf(i) = 1.d0/3.d0
          mcjmpf(i) = 0.d0
          mcflxf(i) = 0.d0
          mcswpf(i) = 0.d0
          mctraf(i) = 0.d0
          mcrotf(i) = 0.d0
          mcmjpf(i) = 0.d0
        else
          if(idnode.eq.0)
     &write(nrite,"(/,'Move frequencies specified')")
        endif
c Now we normalise
        mcinsf(i) = mcinsf(i)/mcmvnorm(i)
        mcdelf(i) = mcinsf(i) + (mcdelf(i)/mcmvnorm(i))
        mcdisf(i) = mcdelf(i) + (mcdisf(i)/mcmvnorm(i))
        mcjmpf(i) = mcdisf(i) + (mcjmpf(i)/mcmvnorm(i))
        mcflxf(i) = mcjmpf(i) + (mcflxf(i)/mcmvnorm(i))
        mcswpf(i) = mcflxf(i) + (mcswpf(i)/mcmvnorm(i))
        mctraf(i) = mcswpf(i) + (mctraf(i)/mcmvnorm(i))
        mcrotf(i) = mctraf(i) + (mcrotf(i)/mcmvnorm(i))
        mcmjpf(i) = mcrotf(i) + (mcmjpf(i)/mcmvnorm(i))
        if(idnode.eq.0)
     &  write(nrite,"('Normalised move frequencies:',/, a18, i4,/,
     &  a18,f7.3,a18,f7.3,a18,f7.3,/,a18,f7.3,a18,f7.3,a18,f7.3,/,
     &  a18,f7.3,a18,f7.3,a18,f7.3)")
     &  'guest:', i,
     &  'insertion:',mcinsf(i), 'deletion:',mcdelf(i)-mcinsf(i),
     &  'displacement:',mcdisf(i)-mcdelf(i),
     &  'jumping:',mcjmpf(i)-mcdisf(i),
     &  'flexing:',mcflxf(i)-mcjmpf(i),'swapping:',mcswpf(i)-mcflxf(i),
     &  'translation:',mctraf(i)-mcswpf(i),
     &  'rotation:',mcrotf(i)-mctraf(i),
     &  'multi jumps:',mcmjpf(i)-mcrotf(i)

        if(idnode.eq.0)
     &  write(nrite,"('Target acceptance ratios:',/,
     &  a18,f7.3,a18,f7.3,a18,f7.3,/)")
     &  'displacement:',disp_ratio(i),'translation:',tran_ratio(i),
     &  'rotation:',rota_ratio
      enddo

      if(.not.lspe)then
        call sleep(idnode+1)
        init=duni(idnode)
        init=duni(idnode)

c     initialize jobcontrol file
c        if(idnode.eq.0)then
c          open(205,file='jobcontrol.in')
c          close(205)
c        endif


c       initialize rotangle 
        rotangle=pi/3.d0
   
c==========================================================================        
c       if restart requested then descend into the branch
c       and read the REVIVE and REVCON for the appropriate
c       arrays

        if(lrestart)then
c         do a check to see if the number of branch directories
c         matches the number of nodes.  adjust accordingly depending
c         on the situation.  
c         This is a rather shitty way of scanning branch directories
c         if other files are called "branch" in the working directory
c         problems happen.

          call revscan(idnode,prevnodes)
          if(prevnodes.eq.0)then
            if(idnode.eq.0)write(nrite,"(/,a85,/)")
     &"No branches found in the working directory, starting from
     &the CONFIG file provided"
            write(outdir,"('branch',i2.2)")idnode+1
            write(command,"('mkdir 'a8)")outdir
            call system(command)
c         if the restart calculation has more nodes than the previous
c         calculation
          elseif(idnode+1.gt.prevnodes)then
c           create the new directory
            write(outdir,"('branch',i2.2)")idnode+1
            write(command,"('mkdir 'a8)")outdir
            call system(command)

c           apply the values from other branches to the new node.
            write(localdir,"('branch',i2.2)")(mod(idnode+1,prevnodes)+1)
            
            call revread(localdir,production,ntpmls,totatm,ntpguest)

          elseif(idnode+1.le.prevnodes)then 
c         read the values from the first mxnodes
            write(localdir,"('branch',i2.2)")idnode+1
            outdir=localdir
            call revread(localdir,production,ntpmls,totatm,ntpguest)
          endif
c     write warning if the number of nodes does not correspond with the 
c     previous calculation
          if(idnode.eq.0.and.prevnodes.ne.mxnode)then
           write(nrite,"(/,3x,a35,i2,a36,i2,/)")"WARNING - previous 
     &calculation had",prevnodes," branches while the current job has ",
     &mxnode
           if(mxnode.lt.prevnodes)write(nrite,"(3x,a9,i2,a4,i2,a17,/)")
     &"Branches ",mxnode+1,"  - ", prevnodes,"  will be ignored"
           if(mxnode.gt.prevnodes)
     & write(nrite,"(3x,a9,i2,a4,i2,a30,i2,a3,i2,/)")
     &"Branches ", prevnodes+1,"  - ", mxnode,"  will take data from
     & branches ",1," - ",prevnodes
          endif
          if(production)then
            if(eqsteps.gt.0)then
              production=.false.
              if(idnode.eq.0)write(nrite,"(3x,a41,i7,a6,/)")
     &"Production averaging will continue after ",eqsteps," steps"
            else
              if(idnode.eq.0)write(nrite,"(3x,a59,/)")
     &"Production averaging will start at the beginning of the run"
            endif
          endif

c     local output for each node
c     This may only work on system specific machines
        else 
          write(outdir,"('branch',i2.2)")idnode+1
          write(command,"('mkdir ',a8)")outdir
          call system(command)
        endif
c=========================================================================
c      open necessary archiving files

        if(ntpguest.gt.1)then
          do i=1,ntpguest
            if(lnumg)then
              write(outfile,"(a8,'/numguests',i2.2,'.out')")outdir,i
              open(400+i,file=outfile)
              call FLUSH(400+i)
c              inquire(file=outfile,exist=file_exists)
c              if(file_exists.eq..false.)then
c                call sleep(5)  ! Wait for the file to appear
c              endif
c              open(400+i,file=outfile)
            endif
            if(abs(nhis).gt.0)then
              write(outfile4,"(a8,'/his',i2.2,'.xyz')")outdir,i
              open(500+i,file=outfile4)
            endif
          enddo
        else
          if(lnumg)then
            outfile=outdir // '/numguests.out'
            open(401,file=outfile)
          endif
          if(abs(nhis).gt.0)then
            outfile4=outdir // '/his.xyz'
            open(501,file=outfile4)
          endif
        endif
c        outfile2=outdir // '/runningstats.out'
c        open(202,file=outfile2)
c        outfile3=outdir // '/energies.out'
c        open(203,file=outfile3)
      endif

c     ins,del,dis store the last move made for each guest type
 
c     debugging.. need to see if all information gets to each
c     node
c      write(debug,"('debug',i2.2)")idnode
c      open(999,file=debug)
c      write(999,'(a20,i2,/)')'data for node ',idnode
c      write(999,'(a20,/,3f16.5,/,3f16.5,/,3f16.5)')'cell vectors',
c     & (cell(i),i=1,3),(cell(j),j=4,6),(cell(k),k=7,9)
 
      call guest_exclude(ntpguest)
      call condense(imcon,totatm,ntpmls,ntpfram,ntpguest)

c      write(999,"('atomic information',/)")
c      do i=1,totatm
c        write(999,'(3x,a4,3f16.6)')atomname(i),xxx(i),yyy(i),zzz(i)
c      enddo

c      write(999,"('some other info ',/)")
c      write(999,"('production?',3x,l)")production
c      write(999,"('number of guests',3x,i3)")ntpguest
c      write(999,"('number of equilibrium steps',3x,i9)")eqsteps
c      write(999,"('number of production steps ',3x,i9)")mcsteps
c      write(999,"('number of probability plots',3x,i9)")ntprob
c      close(999)
      engsrp=0.d0
      engcpe=0.d0
      engacc=0.d0
      engac1=0.d0

c     beta is a constant used in the acceptance criteria for the gcmc
c     moves
      beta=1.d0/(kboltz*temp)
c     this is the relative dielectric constant. default is 1
c     we don't need this....

      epsq=1.d0

c    create ewald interpolation arrays 
      call erfcgen(keyfce,alpha,rcut,drewd)

      spenergy=0.d0
      evdwg=0.d0
      ecoulg=0.d0

      call single_point
     &(imcon,idnode,keyfce,alpha,rcut,delr,totatm,totfram,ntpfram,
     &ntpguest,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,spenergy,evdw,ecoul,evdwg,ecoulg,dlpeng)
    
      if(idnode.eq.0)then
c        write(nrite,'(/,/,a35,f22.6)')'Configurational energy:
c     & ',spenergy
c        write(nrite,'(/,a35,f22.6)')'Initial framework energy :',
c     &(evdw+ecoul)/engunit
        write(nrite,'(a35,f22.6)')'Initial guest energy :',
     &spenergy
        write(nrite,'(a35,f22.6)')'van der Waals energy :',
     &evdwg/engunit
        write(nrite,'(a35,f22.6)')'Electrostatic energy :',
     &ecoulg/engunit
        write(nrite,'(a35,f22.6)')'Energy reported by DL_POLY:',
     &dlpeng

      endif

c     ewald1 arrays populated in single_point subroutine
      call ewald1
     &(imcon,engacc,engsicold,mxatm,volm,
     &alpha,kmax1,kmax2,kmax3,epsq,newld)

c     populate ewald3 array with values for each guest
c     this assumes the guests are rigid so the correction
c     to the reciprocal space ewald arrays is constant.
      fram=.false.
      do k=1,ntpguest
        ewld3sum=0.d0
        mol=locguest(k)
        natms=numatoms(mol)
        do i=1,natms-1
          ik=0
          do j=i+1,natms
            ik=ik+1
            jlist(ik)=j
            xdf(ik)=guestx(k,i)-guestx(k,j)
            ydf(ik)=guesty(k,i)-guesty(k,j)
            zdf(ik)=guestz(k,i)-guestz(k,j)
          enddo

          call images(imcon,ik,cell,xdf,ydf,zdf)

          chg=atmchg(mol,i)
          call ewald3(chg,mol,ik,alpha,engcpe,epsq,fram)
          ewld3sum=ewld3sum+engcpe
          
        enddo
        ewald3en(k)=ewld3sum
      enddo 

      nguests=0
      do i=1,ntpguest
c       this is incorrect because spenergy includes energies
c       calculated from all guests in the framework at the begining
c       where energy separates these values
        energy(i)=spenergy
        mol=locguest(i)
        nmols=nummols(mol)
        nguests=nguests+nmols
      enddo
     
      if(lspe)then
        if(lspe)lgchk=.false.
        call error(idnode,0)
        
      endif

   
      call timchk(1,tzero)
      if((eqsteps.eq.0).and.(.not.ljob))production=.true.


c******************************************************************
c
c       gcmc begins
c
c******************************************************************

c     start delrdisp as delr (read from the CONTROL file)
      do i = 1, ntpguest
        delrdisp(i) = delr(i)
        tran_delr(i) = delr(i)
      enddo
      rota_rotangle = rotangle
      call timchk(0,timelp)
      do while(lgchk)
        gcmccount=gcmccount+1

c     every so often, update the total number of prodcounts
c     across all nodes
        if(mod(gcmccount,400).eq.0)then

c         check jobcontrol
          if(ljob)then
            call jobcheck(idnode,jobsafe,ljob,production)
          endif

          if(.not.jobsafe)then
            call gstate(jobsafe)
            lgchk=.false. 
            call error(idnode,-2312)
          endif
          if (production)then
            call gisum2(prodcount,1,globalprod)
            if(mcsteps.lt.0)then
c             Check for cycles here
              tick_tock_cycles = cshift(tick_tock_cycles, -1)
              tick_tock_cycles(1) = .true.
c             sum up for all guests
              totalguests = 0
              do i=1,ntpguest
                mol=locguest(i)
                totalguests = totalguests + nummols(mol)
              enddo
              call gisum2(totalguests,1,globalnguests)
c             If the total number of production steps over all nodes
c             is less than the total number of guests (plus one)
c             multiplied by desired cycles flag this check as not done
              if((globalprod).lt.(-mcsteps*(globalnguests+1)))then
                tick_tock_cycles(1) = .false.
              endif
c             safe to end if every check passes
              if(all(tick_tock_cycles))then
                if(idnode.eq.0)write(nrite, "('Completed at least ',
     &i10,' production cycles for each guest')")
     &-mcsteps
                lgchk=.false.
              endif
            else
              if(globalprod.ge.mcsteps)lgchk=.false.
            endif
          elseif(.not.ljob)then
c           not production yet; test if we are equilibratied
            if(despre.gt.0)then
              call gisum2(prodcount,1,globalsteps)
              tick_tock_cycles = cshift(tick_tock_cycles, -1)
              tick_tock_cycles(1) = .true.
              totalguests = 0
              do i=1,ntpguest
                mol=locguest(i)
                totalguests = totalguests + nummols(mol)
              enddo
              call gisum2(totalguests,1,globalnguests)
              if((globalsteps).lt.(-eqsteps*(globalnguests+1)/4))then
                tick_tock_cycles(1) = .false.
              endif
              if(tick_tock_cycles(1))then
                totpres = sum(gstpress)
                do i=1,ntpguest
                  istat=1+14*(i-1)
                  mol=locguest(i)
                  temppres=(gstpress(mol)* (1 - packf) * volm / packf) /
     & (temp*13806874.7924278)
                  dqp(mol) = adsguen(mol)-chainstats(istat+1) + 
     & (adspres(mol)-temppres)
                  chainstats(istat+1) = 0
                enddo
                do i=1,ntpguest-1
                  mol=locguest(i)
                  presrate = dqp(mol) - sum(dqp) * gstpress(mol) /
     & totpres
                  gstpress(mol) = dqp(mol) * totpres / sum(dqp)
                  if(gstpress(mol).lt.0)then
                    gstpress(mol) = 0
                  endif
                  presleft = presleft + gstpress(mol)
                enddo
                gstpress(locguest(ntpguest)) = totpres - presleft
                if(gstpress(locguest(ntpguest)).lt.0)then
                  gstpress(locguest(ntpguest)) = 0
                endif
                changerate = totpres / sum(gstpress)
                do i=1,ntpguest
                  mol = locguest(i)
                  gstpress(mol) = gstpress(mol) * changerate
                enddo
                presleft=0.d0
                prodcount = 0
                if((abs(presrate).lt.despre) .or. ((gcmccount).gt.
     &(-50*eqsteps)))then
                  if(idnode.eq.0)write(nrite, "('Equalibrated the
     & presssure to use after ', i10,' steps')")gcmccount
                  do i=1,ntpguest
                    mol=locguest(i)
                    write(nrite,"(32x, '<P', i1, '>: ',f15.9)")
     &              i, gstpress(i) / 100000
                  enddo
                  gcmccount=0
                  despre=0
                endif
              endif
            endif
            if(eqsteps.gt.0)then
              if(gcmccount.ge.eqsteps.and.despre.eq.0)then
                if(idnode.eq.0)write(nrite, "('Completed at least ',
     &i10,' equlibration steps',/,'Starting production at ',i10)")
     &eqsteps, gcmccount
                production=.true.
              endif
            elseif(eqsteps.lt.0)then
              call gisum2(gcmccount,1,globalsteps)
c             Check for cycles here
              tick_tock_cycles = cshift(tick_tock_cycles, -1)
              tick_tock_cycles(1) = .true.
c             sum up for all guests
              totalguests = 0
              do i=1,ntpguest
                mol=locguest(i)
                totalguests = totalguests + nummols(mol)
              enddo
              call gisum2(totalguests,1,globalnguests)
c             If the total number of production steps over all nodes
c             is less than the total number of guests (plus one)
c             multiplied by desired cycles flag this check as not done
              if((globalsteps).lt.(-eqsteps*(globalnguests+1)))then
                tick_tock_cycles(1) = .false.
              endif
              if(despre.gt.0)then
                tick_tock_cycles(1) = .false.
              endif
c             safe to end if every check passes
              if(all(tick_tock_cycles))then
                if(idnode.eq.0)write(nrite, "('Completed at least ',
     &i10,' equilibration cycles for each guest',/,
     &'Starting production at',i10,' over all nodes')")
     &-eqsteps, globalsteps
                production=.true.
c               reset the cycles for production
                tick_tock_cycles=.false.
              endif
            endif
          endif
        endif

        if(mod(gcmccount,1000).eq.0)then
          call revive
     &(idnode,totatm,0,production,ntpguest,ntpmls,imcon,cfgname,
     &   delE(iguest),outdir)
        endif

c      randomly choose a guest type to move 
        if(ntpguest.gt.1)then
          iguest=floor(duni(idnode)*ntpguest)+1
        elseif(ntpguest.eq.1)then
          iguest=1
        endif

c Randomly decide which MC move to do
c Each guest has its specific set of moves
        randmov=duni(idnode)

        if(randmov.lt.mcinsf(iguest))then
          insert = .true.
        elseif(randmov.lt.mcdelf(iguest))then
          delete = .true.
        elseif(randmov.lt.mcdisf(iguest))then
          displace = .true.
        elseif(randmov.lt.mcjmpf(iguest))then
          jump = .true.
        elseif(randmov.lt.mcflxf(iguest))then
          flex = .true.
        elseif(randmov.lt.mcswpf(iguest))then
          swap = .true.
        elseif(randmov.lt.mctraf(iguest))then
          displace = .true.
          tran = .true.
        elseif(randmov.lt.mcrotf(iguest))then
          displace = .true.
          rota = .true.
        elseif(randmov.lt.mcmjpf(iguest))then
          multijump = .true.
        else
c Failover displace -- shouldn't reach here
          displace=.true.
        endif
        
        mol=locguest(iguest)

        natms=numatoms(mol)
        nmols=nummols(mol)

        if((nhis.ne.0).and.(mod(gcmccount,abs(nhis)).eq.0))then
c          write(202,'(a35,f20.15,a15,f15.10,a15,f15.10)')
c     &'displacement acceptance ratio: ',
c     &(dble(accept_disp(iguest))/dble(disp_count(iguest))),
c     &'delr: ',sum(delrdisp),'angle: ',rotangle
          if(nhis.lt.0)call hisarchive(ntpguest,gcmccount)
          if((nhis.gt.0).and.(prodcount.gt.0))call hisarchive
     &      (ntpguest,gcmccount) 
        endif
        if(nmols.ge.1.and.disp_count(iguest).ge.1)then
          if(mod(disp_count(iguest),100).eq.0)then
c update distance part of displacement on n00s
            if((accept_disp(iguest)/dble(disp_count(iguest))).gt.
     &         disp_ratio(iguest))then
              delrdisp(iguest)=delrdisp(iguest)*1.05d0
            else
              if(delrdisp(iguest).gt.delrmin)delrdisp(iguest)=
     &delrdisp(iguest)*0.95d0
            endif
          elseif(mod(disp_count(iguest),50).eq.0)then
c update rotation part of displacement on n50s
            if((accept_disp(iguest)/dble(disp_count(iguest))).gt.
     &         disp_ratio(iguest))then
              rotangle=rotangle*1.05d0
            else
              if(rotangle.gt.minangle)rotangle=rotangle*0.95d0
            endif
          endif
        endif
        if((nmols.ge.1).and.(tran_count(iguest).ge.1).and.
     &mod(tran_count(iguest),100).eq.0)then
c update distance moves on n00s
          if((accept_tran(iguest)/dble(tran_count(iguest))).gt.
     &       tran_ratio(iguest))then 
            tran_delr(iguest)=tran_delr(iguest)*1.05d0
          else
            if(tran_delr(iguest).gt.delrmin)tran_delr(iguest)=
     &tran_delr(iguest)*0.95d0
          endif
        endif
        if((nmols.ge.1).and.(rota_count.ge.1).and.
     &mod(rota_count,100).eq.0)then
c update rotation moves on n00s
          if((accept_rota/dble(rota_count)).gt.rota_ratio)then
            rota_rotangle=rota_rotangle*1.05d0
          else
            if(rota_rotangle.gt.minangle)
     &rota_rotangle=rota_rotangle*0.95d0
          endif
        endif

c       the following is added in case of initial conditions
        if((nmols.eq.0).or.((nmols.eq.1).and.(multijump)))then
          insert=.true.
          delete=.false.
          displace=.false.
          jump = .false.
          flex = .false.
          swap = .false.
          tran = .false.
          rota = .false.
          multijump = .false.
        endif


c***********************************************************************
c    
c             Insertion
c
c***********************************************************************
        if(insert)then
c c       calculate the ewald sums ewald1 and ewald2(only for new mol)
c c       store ewald1 and 2 separately for the next step.
c c       calculate vdw interaction (only for new mol)
          ins(iguest)=1
          ins_count=ins_count+1
          delE(iguest)=0.d0
          ewld1eng=0.d0
          ewld2eng=0.d0
          ewld2sum=0.d0
          ewld3sum=0.d0
          vdweng=0.d0
          vdwsum=0.d0


          if((totatm-framatm)+natms.ge.maxguest) call error(idnode,75)
          randa=duni(idnode)
          randb=duni(idnode)
          randc=duni(idnode)
          rand1=duni(idnode)
          rand2=duni(idnode)
          rand3=duni(idnode)
          rand4=duni(idnode)
          call random_ins(imcon,natms,totatm,iguest,rcut,delr(iguest),
     & randa,randb,randc,rand1,rand2,rand3,rand4)
          do i=1,natms
            ind(i)=0
          enddo
          call guestlistgen
     &(imcon,iguest,totatm,rcut,delr(iguest),
     &natms,newx,newy,newz)

          call ewald1_guest
     &(imcon,ewld1eng,natms,iguest,volm,alpha,
     &kmax1,kmax2,kmax3,epsq,insert,delete,displace,1)
 
c         calculate long range correction to vdw for the insertion
c         of an additional guest
          do i=1,natms
            ka=ltpsit(mol,i)
            numtyp(ka)=numtyp(ka)+1
            if(lfzsite(mol,i).ne.0)numfrz(ka)=numfrz(ka)+1
          enddo
          call gstlrcorrect(idnode,imcon,mol,keyfce,natms,ntpatm,maxvdw,
     &engunit,delrc,rcut,volm)

          delrc=delrc-elrc
c         temp atmadd=0
          atmadd=0
c         do vdw and ewald2 energy calculations for the new atoms

          do i=1,natms
            ik=0
            do j=1,gstlentry(i)
              ik=ik+1
              jatm=gstlist(i,j)
              ilist(j)=jatm
              xdf(ik)=newx(i)-xxx(jatm)
              ydf(ik)=newy(i)-yyy(jatm)
              zdf(ik)=newz(i)-zzz(jatm)
            enddo
            call images(imcon,ik,cell,xdf,ydf,zdf)
            do l=1,gstlentry(i)
              rsqdf(l)=xdf(l)**2+ydf(l)**2+zdf(l)**2
            enddo
c        figure out which index contains charges and ltype arrays
c        that match the guest...

            chg=atmchg(mol,i)
            sittyp=ltpsit(mol,i)

c         calc ewald2 interactions
            call ewald2
     & (chg,gstlentry(i),ewld2eng,drewd,rcut,epsq)
            ewld2sum=ewld2sum+ewld2eng
c         calc vdw interactions 

            call srfrce
     & (sittyp,gstlentry(i),vdweng,rcut,dlrpot)
            vdwsum=vdwsum+vdweng 

          enddo
c        the pairwise intramolecular coulombic correction
c        calculated for the guest at the begining of the 
c        program run.  Assumes constant bond distances. 

          ewld3sum=ewald3en(iguest) 

          gpress=gstpress(iguest)
          ngsts=nummols(mol)

          estep=(ewld1eng+ewld2sum+ewld3sum+vdwsum+delrc)/engunit
          accepted=.false.

          rande=duni(idnode)
          call energy_eval
     &(estep,rande,statvolm,gpress,ngsts,temp,beta,
     &displace,insert,delete,accepted)
          
          if(accepted)then
            delE(iguest)=estep
            accept_ins=accept_ins+1
            mm=natms*nummols(mol)
c           debugging time
c            if(production)then 
c              call debugging(idnode,totatm,levcfg,imcon,cfgname,
c     &eng,outdir,insert,delete,displace,1,delE,ewld1eng,ewld2sum,
c     &ewld3sum,vdwsum,delrc,elrc,engunit,iguest)
c            endif
c           update atomic coordinates
            do i=1,natms
              framwkxxx(mol,mm+i)=newx(i)
              framwkyyy(mol,mm+i)=newy(i)
              framwkzzz(mol,mm+i)=newz(i)
            enddo
c           update ewald1 sums
            do i=1,newld
              ckcsum(i)=ckcsnew(i)
              ckssum(i)=ckssnew(i)
            enddo
c           update long range correction
            elrc=elrc+delrc
c           update nummols,totatm, then condense everything to 1d arrays
            nummols(mol)=nummols(mol)+1
            totatm=totatm+natms
            call condense(imcon,totatm,ntpmls,ntpfram,ntpguest)
            energy(iguest)=energy(iguest)+delE(iguest)

c           debugging time
c            if(production)then
c              call debugging(idnode,totatm,levcfg,imcon,cfgname,
c     &eng,outdir,insert,delete,displace,2,delE,ewld1eng,ewld2sum,
c     &ewld3sum,vdwsum,delrc,elrc,engunit,iguest)
c            endif
          else
             delE(iguest)=0.d0

c            if not accepted, must return numtyp to original value,
c            this is to calculate the long range correction to the vdw
c            sum
             do i=1,natms
               ka=ltpsit(mol,i)
               numtyp(ka)=numtyp(ka)-1
               if(lfzsite(mol,i).ne.0)numfrz(ka)=numfrz(ka)-1
             enddo 
            
          endif
          insert=.false.


c***********************************************************************
c    
c             Deletion 
c
c***********************************************************************
        elseif(delete)then

c       calculate the ewald sum of the mol you wish to delete
c       ewald1,ewald2,vdw of the mol you wish to delete
          del(iguest)=1
          del_count=del_count+1
          randchoice=floor(duni(idnode)*nmols)+1
          delE(iguest)=0.d0        
          ewld1eng=0.d0
          ewld2eng=0.d0
          ewld2sum=0.d0
          ewld3sum=0.d0
          vdweng=0.d0
          vdwsum=0.d0

c         find which index the molecule "randchoice" is
c         assuming the list of atoms is build with the guests
c         first(which should be the case)
          atmadd=0
          iatm=0

c         adds indices in case of multiple guest type
c         NOT TESTED (pretty simple - should work)
          if(iguest.gt.1)then
            do n=1,iguest-1
              imol=locguest(n)
              atmadd=atmadd+numatoms(imol)*nummols(imol)
            enddo
          endif


          at=(randchoice-1)*natms+1
          indx=atmadd+at
          do k=at,at+natms-1
            iatm=iatm+1
            ind(iatm)=atmadd+k
            newx(iatm)=framwkxxx(mol,k)
            newy(iatm)=framwkyyy(mol,k)
            newz(iatm)=framwkzzz(mol,k)
        
          enddo

          call guestlistgen
     &(imcon,iguest,totatm,rcut,delr(iguest),
     &natms,newx,newy,newz)
 
c         calculate long range correction of the system less one
c         molecule of the guest.  the delta value will be calculated
c         by subtracting this value by the current energy (elrc)

          do i=1,natms
            ka=ltpsit(mol,i)
            numtyp(ka)=numtyp(ka)-1
            if(lfzsite(mol,i).ne.0)numfrz(ka)=numfrz(ka)-1
          enddo

          call gstlrcorrect(idnode,imcon,mol,keyfce,natms,ntpatm,maxvdw,
     &engunit,delrc,rcut,volm) 

          delrc=delrc-elrc
          
          call ewald1_guest
     &(imcon,ewld1eng,natms,iguest,volm,alpha,
     &kmax1,kmax2,kmax3,epsq,insert,delete,displace,1)
c         do vdw and ewald2 energy calculations for the new atoms
          
          do i=1,natms

            itatm=ind(i)
            ik=0 
            do j=1,gstlentry(i)
              ik=ik+1
              jatm=gstlist(i,j)
              ilist(j)=jatm
              xdf(ik)=newx(i)-xxx(jatm)
              ydf(ik)=newy(i)-yyy(jatm)
              zdf(ik)=newz(i)-zzz(jatm)
            enddo
            
            call images(imcon,ik,cell,xdf,ydf,zdf)
            do l=1,gstlentry(i)
              rsqdf(l)=xdf(l)**2+ydf(l)**2+zdf(l)**2
            enddo
            chg=atmcharge(itatm)
            sittyp=ltype(itatm)
            
            call ewald2
     & (chg,gstlentry(i),ewld2eng,drewd,rcut,epsq)
c         calc vdw interactions
            
            call srfrce
     & (sittyp,gstlentry(i),vdweng,rcut,dlrpot)
            vdwsum=vdwsum+vdweng
            ewld2sum=ewld2sum+ewld2eng
          enddo
c      calculate the pairwise intramolecular coulombic correction
c      (calculated at the begining - assumes constant bond distance)
          ewld3sum=ewald3en(iguest)

          gpress=gstpress(iguest)
          ngsts=nummols(mol)

          estep=(ewld1eng-ewld2sum-ewld3sum-vdwsum+delrc)/engunit

          accepted=.false.

          rande=duni(idnode)
          call energy_eval
     &(estep,rande,statvolm,gpress,ngsts,temp,beta,
     &displace,insert,delete,accepted)
         
c         the following occurs if the move is accepted.

          if(accepted)then
            delE(iguest)=estep
            accept_del=accept_del+1
            mm=natms*nummols(mol)
c           debugging time
c            if(production)then
c              call debugging(idnode,totatm,levcfg,imcon,cfgname,
c     &eng,outdir,insert,delete,displace,1,delE,ewld1eng,ewld2sum,
c     &ewld3sum,vdwsum,delrc,elrc,engunit,iguest)
c            endif
c           update atomic coordinates
            do i=at+natms,mm
              framwkxxx(mol,i-natms)=framwkxxx(mol,i)
              framwkyyy(mol,i-natms)=framwkyyy(mol,i)
              framwkzzz(mol,i-natms)=framwkzzz(mol,i)
            enddo
c           update ewald1 sums
            do i=1,newld
              ckcsum(i)=ckcsnew(i)
              ckssum(i)=ckssnew(i)
            enddo

c           update nummols,totatm, then condense everything to 1d arrays
            elrc=elrc+delrc
            nummols(mol)=nummols(mol)-1
            totatm=totatm-natms
            call condense(imcon,totatm,ntpmls,ntpfram,ntpguest)
c            call images(imcon,totatm,cell,xxx,yyy,zzz)
            energy(iguest)=energy(iguest)+delE(iguest)
c           debugging time
c            if(production)then
c              call debugging(idnode,totatm,levcfg,imcon,cfgname,
c     &eng,outdir,insert,delete,displace,2,delE,ewld1eng,ewld2sum,
c     &ewld3sum,vdwsum,delrc,elrc,engunit,iguest)
c            endif
          else
            delE(iguest)=0.d0
            do i=1,natms
              ka=ltpsit(mol,i)
              numtyp(ka)=numtyp(ka)+1
              if(lfzsite(mol,i).ne.0)numfrz(ka)=numfrz(ka)+1
            enddo
          endif

          delete=.false.

            
c***********************************************************************
c    
c           Displacement 
c
c***********************************************************************

        elseif(displace)then
          if(tran)then
            dis(iguest)=2
            tran_count(iguest) = tran_count(iguest) + 1
            randa=duni(idnode)
            randb=duni(idnode)
            randc=duni(idnode)
            rand1=0.d0
            rand2=0.d0
            rand3=0.d0
            rand4=0.d0
            dis_delr = tran_delr(iguest)
            dis_rotangle = 0.d0
          elseif(rota)then
            dis(iguest)=3
            rota_count = rota_count + 1
            randa=0.5d0
            randb=0.5d0
            randc=0.5d0
            rand1=duni(idnode)
            rand2=duni(idnode)
            rand3=duni(idnode)
            rand4=duni(idnode)
            dis_delr = 0.d0
            dis_rotangle = rota_rotangle
          else
            dis(iguest)=1
            disp_count(iguest)=disp_count(iguest)+1
            randa=duni(idnode)
            randb=duni(idnode)
            randc=duni(idnode)
            rand1=duni(idnode)
            rand2=duni(idnode)
            rand3=duni(idnode)
            rand4=duni(idnode)
            dis_delr = delrdisp(iguest)
            dis_rotangle = rotangle
          endif
          delE(iguest)=0.d0
          ewld1eng=0.d0
          ewld2eng=0.d0
          ewld2sum=0.d0
          ewld3sum=0.d0
          vdweng=0.d0
          vdwsum=0.d0
c         choose a molecule from the list

          randchoice=floor(duni(idnode)*nmols)+1
c         find which index the molecule "randchoice" is

          atmadd=0
          iatm=0
          
          if(iguest.gt.1)then
            do n=1,iguest-1
              imol=locguest(n)
              atmadd=atmadd+numatoms(imol)*nummols(imol)
            enddo
          endif
          at=(randchoice-1)*natms+1
          do k=at,at-1+natms
            iatm=iatm+1 
            ind(iatm)=atmadd+k
            newx(iatm)=framwkxxx(mol,k)
            newy(iatm)=framwkyyy(mol,k)
            newz(iatm)=framwkzzz(mol,k)
          enddo

          call guestlistgen
     &(imcon,iguest,totatm,rcut,delr(iguest),
     &natms,newx,newy,newz)
 
c    LOOP 1 - calculate the energy for the chosen guest in its orginal
c             place
          call ewald1_guest
     &(imcon,ewld1eng,natms,iguest,volm,alpha,
     &kmax1,kmax2,kmax3,epsq,insert,delete,displace,1)
c         do vdw and ewald2 energy calculations for the atoms
          do i=1,natms
            ik=0
            do j=1,gstlentry(i)
              ik=ik+1
              jatm=gstlist(i,j)
              ilist(j)=jatm
              xdf(ik)=newx(i)-xxx(jatm)
              ydf(ik)=newy(i)-yyy(jatm)
              zdf(ik)=newz(i)-zzz(jatm)
            enddo
            call images(imcon,ik,cell,xdf,ydf,zdf)
            do l=1,gstlentry(i)
              rsqdf(l)=xdf(l)**2+ydf(l)**2+zdf(l)**2
            enddo
            iatm=ind(i)
            chg=atmcharge(iatm)
            sittyp=ltype(iatm)
 
            call ewald2
     & (chg,gstlentry(i),ewld2eng,drewd,rcut,epsq)
            ewld2sum=ewld2sum-ewld2eng 
c         calc vdw interactions 
            call srfrce
     & (sittyp,gstlentry(i),vdweng,rcut,dlrpot)
            vdwsum=vdwsum-vdweng
          enddo

c         LOOP 2 - shift the newx,newy,newz coordinates and re-calculate
c         the new energy and subtract from the above energy.
c dis_delr and dis_rotangle depend on the type of move!
          call random_disp
     &(imcon,natms,mol,newx,newy,newz,dis_delr,cell,rcell,dis_rotangle,
     &randa,randb,randc,rand1,rand2,rand3,rand4)

          call guestlistgen
     &(imcon,iguest,totatm,rcut,delr(iguest),
     &natms,newx,newy,newz)
          call ewald1_guest
     &(imcon,ewld1eng,natms,iguest,volm,alpha,
     &kmax1,kmax2,kmax3,epsq,insert,delete,displace,2) 

c         do vdw and ewald2 energy calculations for the new atoms
          do i=1,natms

            ik=0
            do j=1,gstlentry(i)
              ik=ik+1
              jatm=gstlist(i,j)
              ilist(j)=jatm
              xdf(ik)=newx(i)-xxx(jatm)
              ydf(ik)=newy(i)-yyy(jatm)
              zdf(ik)=newz(i)-zzz(jatm)
            enddo
            call images(imcon,ik,cell,xdf,ydf,zdf)
            do l=1,gstlentry(i)
              rsqdf(l)=xdf(l)**2+ydf(l)**2+zdf(l)**2
            enddo
            itatm=ind(i)
            chg=atmcharge(itatm)
            sittyp=ltype(itatm)
 
            call ewald2
     & (chg,gstlentry(i),ewld2eng,drewd,rcut,epsq)
            ewld2sum=ewld2sum+ewld2eng
c         calc vdw interactions 

            call srfrce
     & (sittyp,gstlentry(i),vdweng,rcut,dlrpot)
            vdwsum=vdwsum+vdweng

          enddo
          gpress=gstpress(iguest)
          ngsts=nummols(mol)

          estep=(ewld1eng+ewld2sum+vdwsum)/engunit

          if(estep.lt.0.d0)then
            accepted=.true.
          else
            accepted=.false.
            rande=duni(idnode)
            call energy_eval
     &(estep,rande,statvolm,gpress,ngsts,temp,beta,
     &displace,insert,delete,accepted)
          endif

          if(accepted)then
            delE(iguest)=estep
c           debugging time
c            if(production)then
c              call debugging(idnode,totatm,levcfg,imcon,cfgname,
c     &eng,outdir,insert,delete,displace,1,delE,ewld1eng,ewld2sum,
c     &ewld3sum,vdwsum,delrc,elrc,engunit,iguest)
c            endif
            if(tran)then
              accept_tran(iguest) = accept_tran(iguest) + 1
            elseif(rota)then
              accept_rota = accept_rota + 1
            else
              accept_disp(iguest)=accept_disp(iguest)+1
            endif
            jj=0
            do j=at,at-1+natms
              jj=jj+1
              framwkxxx(mol,j)=newx(jj)
              framwkyyy(mol,j)=newy(jj)
              framwkzzz(mol,j)=newz(jj)
            enddo
            do i=1,newld
              ckcsum(i)=ckcsnew(i)
              ckssum(i)=ckssnew(i)
            enddo
            call condense(imcon,totatm,ntpmls,ntpfram,ntpguest)
c            call images(imcon,totatm,cell,xxx,yyy,zzz)
            energy(iguest)=energy(iguest)+delE(iguest)
c           debugging time
c            if(production)then
c              call debugging(idnode,totatm,levcfg,imcon,cfgname,
c     &eng,outdir,insert,delete,displace,2,delE,ewld1eng,ewld2sum,
c     &ewld3sum,vdwsum,delrc,elrc,engunit,iguest)
c            endif
          else
            delE(iguest)=0.d0
          endif
          displace=.false.
          tran = .false.
          rota = .false.


c***********************************************************************
c    
c          Long jumps (eg insertion+deletion)
c
c***********************************************************************

        elseif(jump)then
          jmp(iguest)=1
          jump_count=jump_count+1
          delE(iguest)=0.d0
          ewld1eng=0.d0
          ewld2eng=0.d0
          ewld2sum=0.d0
          ewld3sum=0.d0
          vdweng=0.d0
          vdwsum=0.d0
c         choose a molecule from the list

          randchoice=floor(duni(idnode)*nmols)+1
c         find which index the molecule "randchoice" is

          atmadd=0
          iatm=0
          
          if(iguest.gt.1)then
            do n=1,iguest-1
              imol=locguest(n)
              atmadd=atmadd+numatoms(imol)*nummols(imol)
            enddo
          endif
          at=(randchoice-1)*natms+1
          do k=at,at-1+natms
            iatm=iatm+1 
            ind(iatm)=atmadd+k
            newx(iatm)=framwkxxx(mol,k)
            newy(iatm)=framwkyyy(mol,k)
            newz(iatm)=framwkzzz(mol,k)
          enddo

          call guestlistgen
     &(imcon,iguest,totatm,rcut,delr(iguest),
     &natms,newx,newy,newz)
 
          call ewald1_guest
     &(imcon,ewld1eng,natms,iguest,volm,alpha,
     &kmax1,kmax2,kmax3,epsq,insert,delete,jump,1)

c         do vdw and ewald2 energy calculations for the new atoms
          do i=1,natms
            ik=0
            do j=1,gstlentry(i)
              ik=ik+1
              jatm=gstlist(i,j)
              ilist(j)=jatm
              xdf(ik)=newx(i)-xxx(jatm)
              ydf(ik)=newy(i)-yyy(jatm)
              zdf(ik)=newz(i)-zzz(jatm)
            enddo
            call images(imcon,ik,cell,xdf,ydf,zdf)
            do l=1,gstlentry(i)
              rsqdf(l)=xdf(l)**2+ydf(l)**2+zdf(l)**2
            enddo
            iatm=ind(i)
            chg=atmcharge(iatm)
            sittyp=ltype(iatm)
 
            call ewald2
     & (chg,gstlentry(i),ewld2eng,drewd,rcut,epsq)
            ewld2sum=ewld2sum-ewld2eng
c         calc vdw interactions 
            call srfrce
     & (sittyp,gstlentry(i),vdweng,rcut,dlrpot)
            vdwsum=vdwsum-vdweng
          enddo

c         shift the newx,newy,newz coordinates and re-calculate
c         the new energy and subtract from the above energy.

          randa=duni(idnode)
          randb=duni(idnode)
          randc=duni(idnode)
          rand1=duni(idnode)
          rand2=duni(idnode)
          rand3=duni(idnode)
          rand4=duni(idnode)
          
          call random_jump
     &(natms,mol,newx,newy,newz,randa,randb,randc,
     &rand1,rand2,rand3,rand4)
          
          call guestlistgen
     &(imcon,iguest,totatm,rcut,delr(iguest),
     &natms,newx,newy,newz)

          call ewald1_guest
     &(imcon,ewld1eng,natms,iguest,volm,alpha,
     &kmax1,kmax2,kmax3,epsq,insert,delete,jump,2) 

c         do vdw and ewald2 energy calculations for the new atoms
          do i=1,natms

            ik=0
            do j=1,gstlentry(i)
              ik=ik+1
              jatm=gstlist(i,j)
              ilist(j)=jatm
              xdf(ik)=newx(i)-xxx(jatm)
              ydf(ik)=newy(i)-yyy(jatm)
              zdf(ik)=newz(i)-zzz(jatm)
            enddo
            call images(imcon,ik,cell,xdf,ydf,zdf)
            do l=1,gstlentry(i)
              rsqdf(l)=xdf(l)**2+ydf(l)**2+zdf(l)**2
            enddo
            itatm=ind(i)
            chg=atmcharge(itatm)
            sittyp=ltype(itatm)
 
            call ewald2
     & (chg,gstlentry(i),ewld2eng,drewd,rcut,epsq)
            ewld2sum=ewld2sum+ewld2eng
c         calc vdw interactions 

            call srfrce
     & (sittyp,gstlentry(i),vdweng,rcut,dlrpot)
            vdwsum=vdwsum+vdweng

          enddo
          gpress=gstpress(iguest)
          ngsts=nummols(mol)

          estep=(ewld1eng+ewld2sum+vdwsum)/engunit

          if(estep.lt.0.d0)then
            accepted=.true.
          else
            accepted=.false.
            rande=duni(idnode)
            call energy_eval
     &(estep,rande,statvolm,gpress,ngsts,temp,beta,
     &jump,insert,delete,accepted)
          endif

          if(accepted)then
            delE(iguest)=estep
            accept_jump=accept_jump+1
            jj=0
            do j=at,at-1+natms
              jj=jj+1
              framwkxxx(mol,j)=newx(jj)
              framwkyyy(mol,j)=newy(jj)
              framwkzzz(mol,j)=newz(jj)
            enddo
            do i=1,newld
              ckcsum(i)=ckcsnew(i)
              ckssum(i)=ckssnew(i)
            enddo
            call condense(imcon,totatm,ntpmls,ntpfram,ntpguest)
            energy(iguest)=energy(iguest)+delE(iguest)
          else
            delE(iguest)=0.d0
          endif
 

          jump=.false.


c***********************************************************************
c    
c          Multiple long jumps (eg insertion+deletion n molecules before 
c                               energy evaluation)
c
c***********************************************************************

        elseif(multijump)then
          jump=.true.
          jump_count=jump_count+1
          delE(iguest)=0.d0
          estep=0.d0
          allocate(prevmultijump(nmols))
          prevmultijump(:)=0

c         store the original framework, if the move is rejected, the
c         original framework is restored
          origframwkxxx = framwkxxx
          origframwkyyy = framwkyyy
          origframwkzzz = framwkzzz

c         store original ewald1 sums, if the move is rejected, the
c         original value is restored
          ckcsorig = ckcsum
          ckssorig = ckssum

c         store original energy
          origenergy = energy(iguest)

c         choose the number of guest molecules to jump in one move
          !nummove = floor(duni(idnode)*nmols)+1
          nummove = 2
          jmp(iguest)=nummove
          do mvloop = 1, nummove
            ewld1eng=0.d0
            ewld2eng=0.d0
            ewld2sum=0.d0
            ewld3sum=0.d0
            vdweng=0.d0
            vdwsum=0.d0

c           choose a molecule that hasn't jumped in this move
            newmol=.true.
            do while(newmol)
              randchoice=floor(duni(idnode)*nmols)+1
              if(all(prevmultijump.ne.randchoice))then
                prevmultijump(mvloop)=randchoice
                newmol=.false.
              endif
            enddo

c           find which index the molecule "randchoice" is
            atmadd=0
            iatm=0
            if(iguest.gt.1)then
              do n=1,iguest-1
                imol=locguest(n)
                atmadd=atmadd+numatoms(imol)*nummols(imol)
              enddo
            endif
            at=(randchoice-1)*natms+1

            do k=at,at-1+natms
              iatm=iatm+1 
              ind(iatm)=atmadd+k
              newx(iatm)=framwkxxx(mol,k)
              newy(iatm)=framwkyyy(mol,k)
              newz(iatm)=framwkzzz(mol,k)
            enddo

            call guestlistgen
     &(imcon,iguest,totatm,rcut,delr(iguest),
     &natms,newx,newy,newz)

            call ewald1_guest
     &(imcon,ewld1eng,natms,iguest,volm,alpha,
     &kmax1,kmax2,kmax3,epsq,insert,delete,jump,1)

c           do vdw and ewald2 energy calculations for the new atoms
            do i=1,natms
              ik=0
              do j=1,gstlentry(i)
                ik=ik+1
                jatm=gstlist(i,j)
                ilist(j)=jatm
                xdf(ik)=newx(i)-xxx(jatm)
                ydf(ik)=newy(i)-yyy(jatm)
                zdf(ik)=newz(i)-zzz(jatm)
              enddo
              call images(imcon,ik,cell,xdf,ydf,zdf)
              do l=1,gstlentry(i)
                rsqdf(l)=xdf(l)**2+ydf(l)**2+zdf(l)**2
              enddo
              iatm=ind(i)
              chg=atmcharge(iatm)
              sittyp=ltype(iatm)
  
              call ewald2
     & (chg,gstlentry(i),ewld2eng,drewd,rcut,epsq)
              ewld2sum=ewld2sum-ewld2eng 
    
c           calc vdw interactions 
              call srfrce
     & (sittyp,gstlentry(i),vdweng,rcut,dlrpot)
              vdwsum=vdwsum-vdweng
            enddo
    
c           shift the newx,newy,newz coordinates and re-calculate
c           the new energy and subtract from the above energy.
            randa=duni(idnode)
            randb=duni(idnode)
            randc=duni(idnode)
            rand1=duni(idnode)
            rand2=duni(idnode)
            rand3=duni(idnode)
            rand4=duni(idnode)

            call random_jump
     &(natms,mol,newx,newy,newz,randa,randb,randc,
     &rand1,rand2,rand3,rand4)
            
            call guestlistgen
     &(imcon,iguest,totatm,rcut,delr(iguest),
     &natms,newx,newy,newz)

            call ewald1_guest
     &(imcon,ewld1eng,natms,iguest,volm,alpha,
     &kmax1,kmax2,kmax3,epsq,insert,delete,jump,2) 

c           do vdw and ewald2 energy calculations for the new atoms
            do i=1,natms
    
              ik=0
              do j=1,gstlentry(i)
                ik=ik+1
                jatm=gstlist(i,j)
                ilist(j)=jatm
                xdf(ik)=newx(i)-xxx(jatm)
                ydf(ik)=newy(i)-yyy(jatm)
                zdf(ik)=newz(i)-zzz(jatm)
              enddo
              call images(imcon,ik,cell,xdf,ydf,zdf)
              do l=1,gstlentry(i)
                rsqdf(l)=xdf(l)**2+ydf(l)**2+zdf(l)**2
              enddo
              itatm=ind(i)
              chg=atmcharge(itatm)
              sittyp=ltype(itatm)
  
              call ewald2
     & (chg,gstlentry(i),ewld2eng,drewd,rcut,epsq)
              ewld2sum=ewld2sum+ewld2eng
    
c           calc vdw interactions 
              call srfrce
     & (sittyp,gstlentry(i),vdweng,rcut,dlrpot)
              vdwsum=vdwsum+vdweng
    
            enddo
            gpress=gstpress(iguest)
            ngsts=nummols(mol)

c           sum the energy of multiple jumps for evaluation later
            estep=estep+(ewld1eng+ewld2sum+vdwsum)/engunit

c           store conditions after each jump for the next jump
            if (mvloop.lt.nummove)then
              delE(iguest)=estep
c             store the framework after each jump
              jj=0
              do j=at,at-1+natms
                jj=jj+1
                framwkxxx(mol,j)=newx(jj)
                framwkyyy(mol,j)=newy(jj)
                framwkzzz(mol,j)=newz(jj)
              enddo
c             update ewald1 sums
              do i=1,newld
                ckcsum(i)=ckcsnew(i)
                ckssum(i)=ckssnew(i)
              enddo
              call condense(imcon,totatm,ntpmls,ntpfram,ntpguest)
              energy(iguest)=origenergy+delE(iguest)
            endif

          enddo

c         energy evaluation after multiple jumps
          if(estep.lt.0.d0)then
            accepted=.true.
          else
            accepted=.false.
            rande=duni(idnode)
            call energy_eval
     &(estep,rande,statvolm,gpress,ngsts,temp,beta,
     &jump,insert,delete,accepted)
          endif

          if(accepted)then
c            print *, 'MULTIJUMP ACCEPTED, guest ',iguest 
            delE(iguest)=estep
            accept_jump=accept_jump+1
c            print *, 'energy+delE='
            jj=0
            do j=at,at-1+natms
              jj=jj+1
              framwkxxx(mol,j)=newx(jj)
              framwkyyy(mol,j)=newy(jj)
              framwkzzz(mol,j)=newz(jj)
            enddo
            do i=1,newld
              ckcsum(i)=ckcsnew(i)
              ckssum(i)=ckssnew(i)
            enddo
            call condense(imcon,totatm,ntpmls,ntpfram,ntpguest)
c            print *, origenergy, '+', delE(iguest)
            energy(iguest)=origenergy+delE(iguest)
c            print *, energy(iguest)

          else
            delE(iguest)=0.d0
c           restore original framework if move is rejected
            framwkxxx=origframwkxxx
            framwkyyy=origframwkyyy
            framwkzzz=origframwkzzz
c           restore original ewald1 sums if step is rejected
            ckcsum=ckcsorig
            ckssum=ckssorig
            call condense(imcon,totatm,ntpmls,ntpfram,ntpguest)
            energy(iguest)=origenergy
          endif


          jump = .false.
          multijump=.false.
          deallocate(prevmultijump)


c*************************************************************************
        elseif(swap)then
c not implemented :(
          swap = .false.
        elseif(flex)then
c not implemented :(
          flex = .false.       
        endif
        if (.not.production.and.despre.gt.0)then
          if((laccsample).and.(accepted).or.(.not.laccsample))then
            prodcount=prodcount+1
            do i=1,ntpguest
              mol = locguest(i)
              dmolecules=real(nummols(mol))
              istat=1+14*(i-1)
              delta = dmolecules - chainstats(istat+1)
              aN = chainstats(istat+1) + delta/prodcount
              chainstats(istat+1) = aN
            enddo
          endif
        endif
c=========================================================================
c       once GCMC move is done, check if production is requested
c       if so, store averages, probability plots etc..
c=========================================================================

        if(production)then 
          if((laccsample).and.(accepted).or.(.not.laccsample))then
            prodcount=prodcount+1
            rollcount=rollcount+1
            chainstats(1) = dble(prodcount)
            do i=1,ntpguest
              mol=locguest(i) 
              dmolecules=real(nummols(mol))
              molecules2=dmolecules*dmolecules
              energy2=energy(i)*energy(i)
              if(prodcount.le.500001)then
                  gasmol((i-1)*500000+prodcount-1) = dmolecules
                  gaseng((i-1)*500000+prodcount-1) = energy(i)
              endif
              istat=1+14*(i-1)
c           Rolling <N>
              delta = dmolecules - chainstats(istat+1)
              aN = chainstats(istat+1) + delta/prodcount
              chainstats(istat+1) = aN 
c           Rolling <E>
              delta = energy(i) - chainstats(istat+2)
              E = chainstats(istat+2) + delta/prodcount
              chainstats(istat+2) = E 
c           Rolling <EN>
              delta = energy(i)*dmolecules - chainstats(istat+3)
              EN = chainstats(istat+3) + delta/prodcount
              chainstats(istat+3) = EN 
c           Rolling <N^2>
              delta = molecules2 - chainstats(istat+4)
              N2 = chainstats(istat+4) + delta/prodcount
              chainstats(istat+4) = N2 
c           Rolling <E^2>
              delta = energy2 - chainstats(istat+5)
              E2 = chainstats(istat+5) + delta/prodcount
              chainstats(istat+5) = E2
c           Sampling the Henry's coefficient (This requires an
c              energy calculation between the guest and framework
c              only. This can be done, but will require a major
c              overhaul of the energy calculations at each step
c              namely, will need to split the energy between
c              framework - guest and guest - guest.
c              tricky with the ewald sums.
              if((ins(i).eq.1).and.(accepted))then
c             Rolling <exp(-E/kT)>
                H_const=dexp(-1.d0*delE(i)/kboltz/temp)
                ins_rollcount = ins_rollcount + 1
                delta = H_const - chainstats(istat+6)
                H_const = chainstats(istat+6) + delta/accept_ins
                chainstats(istat+6) = H_const
c             Rolling <exp(E/k/T)> for window
                avgwindow(rollstat+8)=
     &((ins_rollcount-1)*avgwindow(rollstat+8)+H_const)/
     &dble(ins_rollcount)
              endif
              rollstat=8*(i-1)
c           Rolling <N> for window
              avgwindow(rollstat+1)=
     &((rollcount-1)*avgwindow(rollstat+1)+dmolecules)/dble(rollcount)
c           Rolling <E> for window
              avgwindow(rollstat+2)=
     &((rollcount-1)*avgwindow(rollstat+2)+energy(i))/dble(rollcount)
c           Rolling <EN> for window
              avgwindow(rollstat+3)=
     &((rollcount-1)*avgwindow(rollstat+3)+
     &energy(i)*dmolecules)/dble(rollcount)
c           Rolling <N^2> for window
              avgwindow(rollstat+4)=
     &((rollcount-1)*avgwindow(rollstat+4)+molecules2)/dble(rollcount)
c           Rolling <E^2> for window
              avgwindow(rollstat+5)=
     &((rollcount-1)*avgwindow(rollstat+5)+energy2)/dble(rollcount)
            enddo
          endif
          if((mod(prodcount,nwind).eq.0).or.(.not.lgchk))then
c         store averages for standard deviations
c         reset windows to zero
            avcount = avcount + 1
            weight = dble(rollcount) / dble(nwind)
            do i=1,ntpguest
              istat = 1+(i-1)*14
              rollstat = (i-1)*8
c           get the rolled averages for all the variables
              aN = chainstats(istat+1) 
              E = chainstats(istat+2)
              EN = chainstats(istat+3)
              N2 = chainstats(istat+4)
              E2 = chainstats(istat+5)
c           compute C_v and Q_st for the windowed averages
              avgwindow(rollstat+6) = calc_Qst(E2, E, aN, N2, EN, temp)
              avgwindow(rollstat+7) = calc_Cv(E2, E, aN, N2, EN, temp)
              if(i>1)then
                aN = avgwindow(rollstat+1)
                do j=1,i-1
                  if(aN.eq.0)then
                    cycle
                  endif
                  ij_select = avgwindow(8*(j-1)+1) / aN
                  delta = weight * (ij_select - seltot(i,j))
                  seltot(i,j) = seltot(i,j) + delta/avcount
                  selwin(i,j) = selwin(i,j)+weight*(ij_select*ij_select-
     & selwin(i,j)) / (avcount-1+weight)
                  seltot(j,i) = ((avcount-1+weight) * selwin(i,j) - 
     &((avcount-1+weight)*seltot(i,j)*seltot(i,j))) / (avcount-1+weight)
                enddo
              endif
            enddo
            do i=1,ntpguest*8
              delta = weight * (avgwindow(i) - sumwindowav(i))
              sumwindowav(i) = sumwindowav(i) + delta/avcount
              varwindow(i) = varwindow(i) + 
     &delta*weight*(avgwindow(i) - sumwindowav(i))
              avgwindow(i)=0.d0
            enddo
            rollcount = 0
        endif
          if(lprob)then
            call storeprob(ntpguest,rcell,ngrida,ngridb,ngridc,
     &lnorm_bin)
          endif
        endif

c       increment the numguests.out storage buffer

c        ibuff=ibuff+1
        if(lnumg.and.(mod(gcmccount,nnumg).eq.0))then
c          write(*,'(4e18.5)')ewld1eng/engunit,ewld2sum/engunit,
c     &ewld3sum/engunit,vdwsum/engunit
c          do i=1,natms
c            write(*,'(a3,3f15.5)')atmname(mol,i),newx(i),newy(i),newz(i)

c          enddo
          do i=1,ntpguest
            mol=locguest(i)
            write(400+i,"(i9,i7,2f20.6,6i4)")
     &        gcmccount,nummols(mol),energy(i),
     &        delE(i),ins(i),del(i),dis(i),jmp(i),flx(i),swp(i)

          enddo
        endif
      
        delE(iguest)=0.d0
        ins(iguest)=0
        del(iguest)=0
        dis(iguest)=0
        jmp(iguest)=0
        flx(iguest)=0
        swp(iguest)=0
      enddo
c*************************************************************************
c     END OF GCMC RUN
c*************************************************************************
      call timchk(0,timelp)

c     run statistics on uptake, energies locally
c     then add the sums globally for a global weighted
c     average
      if(prodcount.gt.0)then
        weight = chainstats(1)
        do i=1,ntpguest
          istat = 1+(i-1)*14
          vstat = (i-1)*8
          if(weight.le.500000)then
              nwind = weight / windn
              avcount = windn
              do x=(1-ntpguest)*8+1,ntpguest*8
                 varwindow(x) = 0.d0
                 sumwindowav(x) = 0.d0
                 avgwindow(x) = 0.d0
              enddo
              do x=1,windn
                y=(i-1)*500000+(x-1)*nwind
                aN = sum(gasmol(y:y+nwind)) / nwind
                E = sum(gaseng(y:y+nwind)) / nwind
                EN = dot_product(gasmol(y:y+nwind),gaseng(y:y+nwind))/ 
     &                           nwind
                N2 = dot_product(gasmol(y:y+nwind),gasmol(y:y+nwind))/
     &                           nwind
                E2 = dot_product(gaseng(y:y+nwind),gaseng(y:y+nwind))/
     &                           nwind
                avgwindow(1) = aN
                avgwindow(2) = E
                avgwindow(3) = EN
                avgwindow(4) = N2
                avgwindow(5) = E2
                avgwindow(6) = calc_Qst(E2,E,aN,N2,EN,temp)
                avgwindow(7) = calc_Cv(E2,E,aN,N2,EN,temp)
                do z=1,7
                  delta = (avgwindow(z) - sumwindowav(z))
                  sumwindowav(z) = sumwindowav(z) + delta/x
                  varwindow(vstat+z) = varwindow(vstat+z) + delta*
     &                                 (avgwindow(z) - sumwindowav(z))
                  avgwindow(z) = 0.d0
                enddo
              enddo
              do x=1,6
                chainstats(istat+x) = sumwindowav(x)
              enddo
          endif
          avgN = chainstats(istat+1)
          avgE = chainstats(istat+2)
          avgEN = chainstats(istat+3)
          avgN2 = chainstats(istat+4)
          avgE2 = chainstats(istat+5)
          avgH = chainstats(istat+6)
          stdN = sqrt(varwindow(vstat+1)/avcount)
          stdE = sqrt(varwindow(vstat+2)/avcount)
          stdEN = sqrt(varwindow(vstat+3)/avcount)
          stdN2 = sqrt(varwindow(vstat+4)/avcount)
          stdE2 = sqrt(varwindow(vstat+5)/avcount)
          stdQst = sqrt(varwindow(vstat+6)/avcount)
          stdCv = sqrt(varwindow(vstat+7)/avcount)
          stdH = sqrt(varwindow(vstat+8)/avcount)
          chainstats(istat+7) = stdN
          chainstats(istat+8) = stdE
          chainstats(istat+9) = stdEN
          chainstats(istat+10) = stdN2
          chainstats(istat+11) = stdE2
          chainstats(istat+12) = stdQst
          chainstats(istat+13) = stdCv
          chainstats(istat+14) = stdH
          do j=1,i-1
            seltot(j,i) = sqrt(seltot(j,i))
          enddo
        enddo
      endif

c     sum all variables for final presentation

c     first send all local statistics to the master node.
c     via csend([messagetag],[variable to send],[length],[destination],
c     [dummy var])
c     message tag for chainstats is even
      if(idnode.gt.0)then
        call csend(idnode*2+1,chainstats,1+ntpguest*14,0,1)
      endif


c     this is final file writing stuff.. probably should put this in a 
c     subroutine so it looks less messy in the main program. 
      if(idnode.eq.0)then

         nodeweight(1)=weight
         write(nrite,"(/,a75,/,3x,'Data reported from node ',
     &i3,/,a75,/)")repeat('*',75),0,repeat('*',75)
         do i=1,ntpguest
           istat=1+(i-1)*14
           mol=locguest(i)
           avgN = chainstats(istat+1)
           avgE = chainstats(istat+2)
           avgEN = chainstats(istat+3)
           avgN2 = chainstats(istat+4)
           avgE2 = chainstats(istat+5)
           avgH = chainstats(istat+6)
           stdN = chainstats(istat+7) 
           stdE = chainstats(istat+8) 
           stdEN = chainstats(istat+9) 
           stdN2 = chainstats(istat+10) 
           stdE2 = chainstats(istat+11) 
           stdQst = chainstats(istat+12)
           stdCv = chainstats(istat+13)
           stdH = chainstats(istat+14)
c     isosteric heat calculation 
           Q_st = calc_Qst(avgE2, avgE, avgN, avgN2, avgEN, temp)
c     constant volume heat capacity
           C_v = calc_Cv(avgE2, avgE, avgN, avgN2, avgEN, temp)


           write(nrite,"(5x,'guest ',i2,': ',40a,/)")i,
     &(molnam(p,mol),p=1,40)
           write(nrite,"(5x,a60,f20.9,/,5x,a60,f20.9,/,
     &     5x,a60,f20.9,/,5x,a60,f20.9,/,5x,a60,f20.9,/,
     &     5x,a60,i20,/,5x,a60,f20.9,/,5x,a60,f20.9,/,
     &     5x,a60,f20.9,/,5x,a60,f20.9,/)")
c     &     ,5x,a60,f20.9,/,5x,a60,f20.9,/)")
     &     '<N>: ', avgN, 
     &     '<E>: ', avgE,
     &     '<E*N>: ', avgEN, 
     &     '<N*N>: ', avgN2,
     &     '<E*E>: ', avgE2,
     &     'Multiplier: ',prodcount,
     &     'Isosteric heat of adsorption (kcal/mol): ',Q_st,
     &     'Isosteric heat error: ', stdQst,
     &     'Heat capacity, Cv (kcal/mol/K): ', C_v,
     &     'Heat capacity error: ', stdCv
c     &     "Henry's Constant (mol/kcal): ", avgH/kboltz/temp,
c     &     "Henry's Constant error: ", stdH/kboltz/temp

           nstat = (i-1)*8
           node_avg(1,(nstat+1):(nstat+5)) =
     &                    chainstats((istat+1):(istat+5))
           node_avg(1, nstat+6) = Q_st
           node_avg(1, nstat+7) = C_v
           node_avg(1, nstat+8) = chainstats(istat+6)
           node_std(1,nstat+1:nstat+8) =
     &                   chainstats(istat+7:istat+14)
         enddo
         if(ntpguest>1)then
           do i=1,ntpguest
             do j=1,i-1
               write(nrite,"(5x, 'guest ',i2,' over guest ',i2,/)")
     &               j,i
               write(nrite,"(5x,a60,f20.9,/)")
     &         '<S>: ', seltot(i,j)
             enddo
           enddo
         endif
         tw=tw+weight
         do i=1,mxnode-1
c        recieve all data from other nodes (via crecv)
c        prodcount used for weighting the mean and stdev
           statbuff(1:1+ntpguest*14) = 0.d0
           call crecv(i*2+1,statbuff,1+ntpguest*14,1)
           weight=statbuff(1)
           nodeweight(i+1)=weight
           tw=tw+weight
           prodcount=prodcount+int(weight)
           write(nrite,"(/,a75,/,3x,'Data reported from node ',i3,
     &/,a75,/)")repeat('*',75),i,repeat('*',75)
           do j=1,ntpguest
             mol=locguest(j)
             write(nrite,"(5x,'guest ',i2,': ',40a,/)")j,
     &       (molnam(p,mol),p=1,40)
             istat=1+(j-1)*14

             avgN = statbuff(istat+1)
             avgE = statbuff(istat+2)
             avgEN = statbuff(istat+3)
             avgN2 = statbuff(istat+4)
             avgE2 = statbuff(istat+5)
             avgH = statbuff(istat+6)
             stdN = statbuff(istat+7)
             stdE = statbuff(istat+8)
             stdEN = statbuff(istat+9)
             stdN2 = statbuff(istat+10)
             stdE2 = statbuff(istat+11)
             stdQst = statbuff(istat+12)
             stdCv = statbuff(istat+13)
             stdH = statbuff(istat+14)

             Q_st = calc_Qst(avgE2, avgE, avgN, avgN2, avgEN, temp)
             C_v = calc_cv(avgE2, avgE, avgN, avgN2, avgEN, temp)
             write(nrite,"(5x,a60,f20.9,/,5x,a60,f20.9,/,
     &       5x,a60,f20.9,/,5x,a60,f20.9,/,5x,a60,f20.9,/,
     &       5x,a60,i20,/,5x,a60,f20.9,/,5x,a60,f20.9,/,
     &       5x,a60,f20.9,/,5x,a60,f20.9,/)")
c     &       ,5x,a60,f20.9,/,5x,a60,f20.9,/)")
     &       '<N>: ', avgN, 
     &       '<E>: ', avgE,
     &       '<E*N>: ', avgEN, 
     &       '<N*N>: ', avgN2,
     &       '<E*E>: ', avgE2,
     &       'Multiplier: ',prodcount,
     &       'Isosteric heat of adsorption (kcal/mol): ',Q_st,
     &       'Isosteric heat error: ', stdQst,
     &       'Heat capacity, Cv (kcal/mol/K): ', C_v,
     &       'Heat capacity error: ', stdCv
c     &       'Henrys Constant (mol/kcal): ', avgH/kboltz/temp,
c     &       'Henrys Constant error: ', stdH/kboltz/temp

             nstat = (j-1)*8
             node_avg(i+1,nstat+1:nstat+5) = statbuff(istat+1:istat+5)
             node_avg(i+1,nstat+6) = Q_st
             node_avg(i+1,nstat+7) = C_v
             node_avg(i+1,nstat+8) = statbuff(istat+6)
             node_std(i+1,nstat+1:nstat+8) = statbuff(istat+7:istat+14)
           enddo
         enddo
      endif

c     sum over all parallel nodes
      call gisum(accept_ins,1,buffer)
      call gisum(ins_count,1,buffer)

      call gisum(accept_del,1,buffer)
      call gisum(del_count,1,buffer)

      call gisum(accept_disp,ntpguest,buffer)
      call gisum(disp_count,ntpguest,buffer)

      call gisum(accept_jump,1,buffer)
      call gisum(jump_count,1,buffer)

      call gisum(accept_flex,1,buffer)
      call gisum(flex_count,1,buffer)

      call gisum(accept_swap,1,buffer)
      call gisum(swap_count,1,buffer)

      call gisum(accept_tran,ntpguest,buffer)
      call gisum(tran_count,ntpguest,buffer)

      call gisum(accept_rota,1,buffer)
      call gisum(rota_count,1,buffer)

      call gisum(gcmccount,1,buffer)

c     write final probability cube files
      
      cprob=0
      cell=cell*angs2bohr
      if(lprob.and.prodcount.gt.0)then
        do i=1,ntpguest
          iprob=0
          do j=1,nprob(i)
            cprob=cprob+1
            iprob=iprob+1
            call gdsum3(grid,cprob,ntprob,gridsize,gridbuff)
            if(idnode.eq.0)call writeprob
     &(i,cprob,iprob,cell,ntpguest,ntpfram,gridsize,
     &ngrida,ngridb,ngridc,prodcount)
            if(lnorm_bin.and.idnode.eq.0)call writeprobnormalbin
     &(i,cprob,iprob,cell,ntpguest,ntpfram,gridsize,
     &ngrida,ngridb,ngridc,prodcount)
          enddo
        enddo
      endif

      if(idnode.eq.0)then
        write(nrite,"(/,a17,i9,a15,f13.3,a8)")
     &'time elapsed for ',gcmccount,' gcmc steps : ',timelp,' seconds'
        write(nrite,"(/,a30,i9)")
     &'total accepted steps : ',accept_ins+accept_del+
     &sum(accept_disp)+accept_jump+accept_flex+accept_swap+
     &sum(accept_tran)+accept_rota
        if(ins_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'insertion ratio: ',dble(accept_ins)/dble(ins_count)
        if(del_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'deletion ratio: ',dble(accept_del)/dble(del_count)
        if(sum(disp_count).gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'displacement ratio: ',dble(sum(accept_disp))/
     &dble(sum(disp_count))
        if(jump_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'jump ratio: ',dble(accept_jump)/dble(jump_count)
        if(flex_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'flex ratio: ',dble(accept_flex)/dble(flex_count)
        if(swap_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'swap ratio: ',dble(accept_swap)/dble(swap_count)
        if(sum(tran_count).gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'translation ratio: ',dble(sum(accept_tran))/
     &dble(sum(tran_count))
        if(rota_count.gt.0)
     &write(nrite,"(/,3x,a21,f15.9)")
     &'rotation ratio: ',dble(accept_rota)/dble(rota_count)
        do i=1,ntpguest
          mol=locguest(i)

          avgN = 0.d0
          avgE = 0.d0
          avgEN = 0.d0
          avgN2 = 0.d0
          avgE2 = 0.d0
          avgH = 0.d0
          stdN = 0.d0
          stdE = 0.d0
          stdEN = 0.d0
          stdN2 = 0.d0
          stdE2 = 0.d0
          stdQst = 0.d0
          stdCv = 0.d0
          stdH = 0.d0
c         compute unions of averages and standard deviations
c          do kk=1,mxnode
c            write(nrite,"(/,3x,f15.6)")nodeweight(kk)
c          enddo
          call avunion(i,mxnode,avgN,avgE,avgEN,avgN2,avgE2, avgH)
c       isosteric heat of adsorption 
          Q_st = calc_Qst(avgE2, avgE, avgN, avgN2, avgEN,temp)
c       heat capacity
          C_v = calc_Cv(avgE2, avgE, avgN, avgN2, avgEN,temp)
          call stdunion
     &(i,mxnode,stdN,stdE,stdEN,stdN2,stdE2,stdQst,stdCv,stdH,
     &avgN,avgE,avgEN,avgN2,avgE2,Q_st,C_v,avgH)

          write(nrite,"(/,'final stats for guest ',i2,3x,40a,/)")
     &       i,(molnam(p,mol),p=1,40)
          write(nrite,"(/,a36,f15.6,a5,f12.3,/,
     &    a36,f15.6,a5,f12.3,/,a36,f15.6,a5,f12.3,/,
     &    a36,f15.6,a5,f12.3,/,a36,f15.6,a5,f12.3,/)")
     &    '<N>: ',avgN, ' +/- ', stdN,
     &    '<E>: ',avgE, ' +/- ', stdE,
     &    '<E*N>: ',avgEN, ' +/- ', stdEN,
     &    '<N*N>: ',avgN2, ' +/- ', stdN2,
     &    '<E*E>: ',avgE2, ' +/- ', stdE2
          write(nrite,"(/,a60,f15.6,/,a60,f15.6)")
     &    'average number of guests: ',avgN,
     &    'standard error: ',stdN 
          write(nrite,"(a60,f15.6,/,a60,f15.6,/,
     &                  a60,f15.6,/,a60,f15.6, 
     &                  /,a60,i15,/)")
c     &                  /,a60,f15.6,/,a60,f15.6,/,a60,i15,/)")
     &    'Isosteric heat of adsorption (kcal/mol): ',Q_st,
     &    'Isosteric heat error: ', stdQst,
     &    'Heat capacity, Cv (kcal/mol/K): ', C_v,
     &    'Heat capacity error: ', stdCv,
c     &    'Henrys Constant (mol/kcal): ', avgH/kboltz/temp,
c     &    'Henrys Constant error: ', stdH/kboltz/temp,
     &    'Total steps counted: ',int(tw)
        enddo
        if(ntpguest>1)then
          do i=1,ntpguest
            do j=1,i-1
            write(nrite,"(/'selectivity stats for guest ',i2,'/'i2)")
     &         j,i
              write(nrite,"(a36,f15.6,a5,f12.3)")
     &        '<S>:',seltot(i,j), '+/-', seltot(j,i)
            enddo
          enddo
          if(desorb.eq.1)then
            write(nrite,"(/'Equalibrated pressures')")
            do i=1,ntpguest
              write(nrite,"(32x, '<P', i1, '>: ',f15.9)")
     &              i, gstpress(i) / 100000
            enddo
          endif
        endif
      endif
      close(202)
      close(ncontrol)
      close(nconfig)
      close(nfield)
      do i=1,ntpguest
        if(lnumg)close(400+i)
        if(abs(nhis).gt.0)close(500+i)
      enddo
      if(idnode.eq.0)then
       close(nrite)
c       close(nang)
      endif
      call exitcomms()
      contains
      character*9 function month(date)
      implicit none
      character*2 date

      if(date.eq.'01')month='January'
      if(date.eq.'02')month='February'
      if(date.eq.'03')month='March'   
      if(date.eq.'04')month='April'   
      if(date.eq.'05')month='May'     
      if(date.eq.'06')month='June'
      if(date.eq.'07')month='July'
      if(date.eq.'08')month='August'
      if(date.eq.'09')month='September'
      if(date.eq.'10')month='October'
      if(date.eq.'11')month='November'
      if(date.eq.'12')month='December'

      return
      end function month
      subroutine avunion(iguest,mxnode,avgN,avgE,avgEN,avgN2,avgE2,avgH)

      implicit none
      real(8) avgE,avgN,avgEN,avgN2,avgE2,avgH,sumweight,weight
      integer iguest,i,node,mxnode,istat
      istat=(iguest-1)*8
      sumweight=0.d0
      do node=1,mxnode
        weight = nodeweight(node)
        sumweight = weight + sumweight
        avgN = avgN + weight*node_avg(node,istat+1)
        avgE = avgE + weight*node_avg(node,istat+2)
        avgEN = avgEN + weight*node_avg(node,istat+3)
        avgN2 = avgN2 + weight*node_avg(node,istat+4)
        avgE2 = avgE2 + weight*node_avg(node,istat+5)
        avgH = avgH + weight*node_avg(node,istat+8)
      enddo

      avgE = avgE/sumweight
      avgN = avgN/sumweight
      avgEN = avgEN/sumweight
      avgN2 = avgN2/sumweight
      avgE2 = avgE2/sumweight
      avgH = avgH/sumweight
      end subroutine avunion
       
      subroutine stdunion(iguest,mxnode,stdN,stdE,stdEN,stdN2,stdE2,
     &stdQst,stdCv,stdH,avgN,avgE,avgEN,avgN2,avgE2,Q_st,C_v,avgH)

      implicit none
      real(8) stdN,stdE,stdEN,stdN2,stdE2,stdQst,stdCv,stdH
      real(8) weight,sumweight
      real(8) avgN,avgE,avgEN,avgN2,avgE2,Q_st,C_v,avgH
      integer mxnode,iguest,istat,node
      
      istat = (iguest-1)*8
      sumweight=0.d0
      do node = 1,mxnode 
        weight = nodeweight(node)
        sumweight = nodeweight(node)+sumweight

        stdN = stdN + weight*
     &(node_std(node,istat+1)**2+node_avg(node,istat+1)**2)
        stdE = stdE + weight*
     &(node_std(node,istat+2)**2+node_avg(node,istat+2)**2)
        stdEN = stdEN + weight*
     &(node_std(node,istat+3)**2+node_avg(node,istat+3)**2)
        stdN2 = stdN2 + weight*
     &(node_std(node,istat+4)**2+node_avg(node,istat+4)**2)
        stdE2 = stdE2 + weight*
     &(node_std(node,istat+5)**2+node_avg(node,istat+5)**2)
        stdQst = stdQst + weight*
     &(node_std(node,istat+6)**2+node_avg(node,istat+6)**2)
        stdCv = stdCv + weight*
     &(node_std(node,istat+7)**2+node_avg(node,istat+7)**2)
        stdH = stdH + weight*
     &(node_std(node,istat+8)**2+node_avg(node,istat+8)**2)
      enddo
      stdN = sqrt((stdN/sumweight) - avgN**2)
      stdE = sqrt((stdE/sumweight) - avgE**2)
      stdEN = sqrt((stdEN/sumweight) - avgEN**2)
      stdN2 = sqrt((stdN2/sumweight) - avgN2**2)
      stdE2 = sqrt((stdE2/sumweight) - avgE2**2)
      stdQst = sqrt((stdQst/sumweight) - Q_st**2)
      stdCv = sqrt((stdCv/sumweight) - C_v**2)
      stdH = sqrt((stdH/sumweight) - avgH**2)
      end subroutine stdunion
      subroutine timchk(ktim,time)

c***********************************************************************
c     
c     dlpoly timing routine for time elapsed in seconds
c     copyright daresbury laboratory
c     author w.smith nov 2003
c
c     wl
c     2006/11/28 16:32:48
c     1.2
c     Exp
c     
c***********************************************************************
      implicit none

      logical init
      character*12 dat,tim,zon
      integer idnode,mynode,ktim,day
      real(8) time,told,tsum,tnow
      integer info(8)

      save init,idnode,told,tsum,day

      data init/.true./

   10 format(/,'time elapsed since job start = ',f15.3,' seconds',/)

      call date_and_time(dat,tim,zon,info)
      
      if(init)then

         tsum=0.d0
         time=0.d0
         day=info(3)
         idnode=mynode()
         told=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))
         init=.false.

      else 

         tnow=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))
         if(day.ne.info(3))then
           told=told-86400.d0
           day=info(3)
         endif
         tsum=tsum+tnow-told
         told=tnow
         time=tsum

      endif

      if(ktim.gt.0.and.idnode.eq.0) write(nrite,10)time

      return
      end subroutine timchk
      subroutine hisarchive(ntpguest,gcmccount)
c*****************************************************************************
c
c     writes an xyz file of the current geometries of the guest
c
c*****************************************************************************
      implicit none

      integer i,imol,nmols,natms,gcmccount,iatm,ntpguest

      do k=1,ntpguest
        imol=locguest(k)
        nmols=nummols(imol)
        natms=numatoms(imol)
        write(500+k,"(i6,/,'step ',i7)")nmols*natms,gcmccount

        iatm=0
        do i=1,nmols
           do j=1,natms
            iatm=iatm+1
            write(500+k,'(2x,a1,3f16.6)')atmname(imol,j),
     &  framwkxxx(imol,iatm),framwkyyy(imol,iatm),framwkzzz(imol,iatm)
          enddo
        enddo
      enddo

      end subroutine hisarchive
      subroutine single_point
     &(imcon,idnode,keyfce,alpha,rcut,delr,totatm,totfram,ntpfram,
     &ntpguest,volm,kmax1,kmax2,kmax3,epsq,dlrpot,ntpatm,maxvdw,
     &engunit,spenergy,evdw,ecoul,evdwg,ecoulg,dlpeng)
c*****************************************************************************
c
c     does a single point energy calculation over all atoms in the
c      system
c
c*****************************************************************************
      implicit none
      logical insert,delete,displace,fram
      integer i,ii,ik,j,jj,p,kmax1,kmax2,kmax3,imcon,keyfce,indx
      integer totatm,newld,sittyp,idnode,ntpatm,maxvdw,totfram
      integer k,mol,ntpguest,natms,nmols,ntpfram
      real(8) drewd,dlrpot,volm,epsq,alpha,rcut,delr(ntpguest),ecoul
      real(8) tempdelr
      real(8) engcpe,engsic,chg,engsrp,elrc,engunit,ecoulg,evdwg,evdw
      real(8) ewald1sum,ewald2sum,ewald3sum,vdwsum
      real(8) ewald1eng,ewald2eng,ewald3eng,vdweng
      real(8) spenergy,delrc,dlpeng

      data ewald1sum,ewald2sum,ewald3eng,vdwsum/0.d0,0.d0,0.d0,0.d0/ 
      data ewald1eng,ewald2eng,vdweng/0.d0,0.d0,0.d0/ 
 

      spenergy=0.d0
c     long range correction to short range forces.
c     this initializes arrays, a different subroutine
c     is called during the gcmc simulation
c     "gstlrcorrect"

      call lrcorrect(idnode,imcon,keyfce,totatm,ntpatm,maxvdw,
     &elrc,engunit,rcut,volm)

      call erfcgen(keyfce,alpha,rcut,drewd)
   
      call condensefram(totfram,ntpfram)
c      call images(imcon,totfram,cell,xxx,yyy,zzz)
      
c     reciprocal space ewald calculation
      call ewald1(imcon,engcpe,engsic,totfram,volm,alpha,
     &kmax1,kmax2,kmax3,epsq,newld)

      ewald1sum=engcpe
c      fram=.true.
c     generate neighbour list
c      call parlst(imcon,totfram,fram,rcut,delr)

c     start loop over atoms

c      do i=1,totfram
c        do k=1,lentry(i)
c
c          j=list(i,k)
c
c          ilist(k)=j
c
c          xdf(k)=xxx(i)-xxx(j)
c          ydf(k)=yyy(i)-yyy(j)
c          zdf(k)=zzz(i)-zzz(j)
c
c        enddo

c        call images(imcon,lentry(i),cell,xdf,ydf,zdf)
c
c        do k=1,lentry(i)
          
c          j=list(i,k)
c          rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
c
c        enddo

c       calculate pairwise interactions
        
c        chg=atmcharge(i)
c        call ewald2(chg,lentry(i),engcpe,drewd,rcut,epsq)
c
c        ewald2sum=ewald2sum+engcpe
c
c        sittyp=ltype(i)
c        call srfrce(sittyp,lentry(i),engsrp,rcut,dlrpot)
c
c        vdwsum=vdwsum+engsrp
c
c      enddo

c     ewald 3 calculation NOT done on framework.  this is 
c     because the framework is considered fixed, all 
c     sites are excluded from the energy calculation so
c     this is redundant
      fram=.false.
c     done framework energy calcs
c      ecoul=ewald1eng+ewald2sum
c      evdw=vdwsum
c      spenergy=(ecoul+evdw)/engunit

c     second time 'round do guest energy calcs
      
c     generate neighbour list
      call condense(imcon,totatm,ntpmls,ntpfram,ntpguest)
c      call images(imcon,totatm,cell,xxx,yyy,zzz)

      tempdelr=delr(1)
      call parlst(imcon,totatm,fram,rcut,tempdelr)
c     reciprocal space ewald calculation
      call ewald1(imcon,engcpe,engsic,totatm,volm,alpha,
     &kmax1,kmax2,kmax3,epsq,newld)

      ewald1eng=engcpe-ewald1sum
c      ewald1eng=engcpe

c     start loop over atoms
c     NEED to separate energies calculated for different
c     guests if multiple guests are requested.
c     these initial values will skew the enthalpy 
c     of adsorption otherwise.

c     the if statement checks if any guests are present in the
c     system prior to the gcmc sim.  If there aren't any
c     the following loops are ignored.
      if(totatm.gt.totfram)then
        do i=1,totatm
          do k=1,lentry(i)

            j=list(i,k)
 
            ilist(k)=j

            xdf(k)=xxx(i)-xxx(j)
            ydf(k)=yyy(i)-yyy(j)
            zdf(k)=zzz(i)-zzz(j)

          enddo

          call images(imcon,lentry(i),cell,xdf,ydf,zdf)

          do k=1,lentry(i)
          
            j=list(i,k)
            rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2

          enddo

c       calculate pairwise interactions
        
          chg=atmcharge(i)
          call ewald2(chg,lentry(i),engcpe,drewd,rcut,epsq)

          ewald2eng=ewald2eng+engcpe

          sittyp=ltype(i)
          call srfrce(sittyp,lentry(i),engsrp,rcut,dlrpot)
 
          vdweng=vdweng+engsrp

        enddo

        fram=.true.
        do i=1,totatm
 
          do k=1,nexatm(i)
 
            j=lexatm(i,k)
            jlist(k)=j
  
            xdf(k)=xxx(i)-xxx(j)
            ydf(k)=yyy(i)-yyy(j)
            zdf(k)=zzz(i)-zzz(j)

          enddo

          call images(imcon,nexatm(i),cell,xdf,ydf,zdf)
 
          chg=atmcharge(i)
          call ewald3(chg,mol,nexatm(i),alpha,engcpe,epsq,fram)
          ewald3eng=ewald3eng+engcpe  

        enddo
        fram=.false.
        ecoulg=ecoulg+ewald1eng+ewald2eng+ewald3eng

        evdwg=evdwg+vdweng+elrc

      endif
      spenergy=(ecoulg+evdwg)/engunit
c      write(*,'(f25.10)') ewald1eng+ewald1sum,ewald2eng,ewald3eng,vdweng
c     &, elrc/engunit

      dlpeng=(ewald1eng+ewald1sum+ewald2eng+vdweng+elrc+ewald3eng)
     & /engunit
      return
      end subroutine single_point

      function duni(idnode)
 
c*********************************************************************
c     
c     dl_poly random number generator based on the universal
c     random number generator of marsaglia, zaman and tsang
c     (stats and prob. lett. 8 (1990) 35-39.) it must be
c     called once to initialise parameters u,c,cd,cm
c     
c     copyright daresbury laboratory 1992
c     author -  w.smith         july 1992
c      
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c      
c*********************************************************************

      implicit none

    
      logical new
      integer i,j,k,l,m,ii,jj,ir,jr,idnode
      integer,dimension(8) :: values
      real(4) s,t,u,c,cd,cm,uni
      real(8) duni
      dimension u(97)

      data new/.true./

      save u,c,cd,cm,uni,ir,jr,new
CVAM
CVAM      call VTBEGIN(134, ierr)
CVAM

      if(new)then

c      initial values of i,j,k must be in range 1 to 178 (not all 1)
c      initial value of l must be in range 0 to 168.

        call date_and_time(values=values)
c       note, these date_and_time values can be 0
        i=mod(values(8)*3, 177) + 1
        j=mod(values(7)*23, 177) + 1
        k=mod(values(8)*5, 177) + 1
        l=mod(values(7)*3, 168)

c       This is in case we find the same problem arises
c        write(nrite,"(/,'values for date and time', 9i6,/)")values
        if(idnode.eq.0)write(nrite,"(/,'i,j,k,l values for duni()', 4i5,
     &/)")i,j,k,l
        ir=97
        jr=33
        new=.false.

        do 200 ii=1,97
          s=0.0
          t=0.5
          do 100 jj=1,24
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32)s=s+t
            t=0.5*t
  100     continue
          u(ii)=s
  200   continue
        c =  362436.0/16777216.0
        cd= 7654321.0/16777216.0
        cm=16777213.0/16777216.0
      else

c      calculate random number
       uni=u(ir)-u(jr)
       if(uni.lt.0.0)uni=uni+1.0
       u(ir)=uni
       ir=ir-1
       if(ir.eq.0)ir=97
       jr=jr-1
       if(jr.eq.0)jr=97
       c=c-cd
       if(c.lt.0.0)c=c+cm
       uni=uni-c
       if(uni.lt.0.0)uni=uni+1.0
       duni=dble(uni)
      endif
CVAM
CVAM     call VTEND(134, ierr)
CVAM
      return
      end function duni
      end

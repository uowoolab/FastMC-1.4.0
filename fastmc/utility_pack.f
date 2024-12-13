      module utility_pack
      use parse_module

      implicit none

      integer nsite
      integer mxlist,mxgrid,mxegrd
      real(8) volum, total_pressure, despre, packf
      real(8), allocatable :: xxx(:),framwkxxx(:,:),origframwkxxx(:,:)
      real(8), allocatable :: yyy(:),framwkyyy(:,:),origframwkyyy(:,:)
      real(8), allocatable :: zzz(:),framwkzzz(:,:),origframwkzzz(:,:)
      real(8), allocatable :: newx(:),newy(:),newz(:),selsum(:,:)
      real(8), allocatable :: guestx(:,:),guesty(:,:),guestz(:,:)
      real(8), allocatable :: atmcharge(:),atmchg(:,:),selvar(:,:)
      real(8), allocatable :: atmweight(:),atmwght(:,:),selwin(:,:)
      real(8), allocatable :: frambuff(:,:), selectivity(:,:)
      real(8), allocatable :: xdf(:),ydf(:),zdf(:),rsqdf(:)
      real(8), allocatable :: selwin2(:,:), seltot(:,:)
      character(8), allocatable :: unqatm(:),atomname(:),atmname(:,:)
      character*70, allocatable :: numgbuff(:,:)
      character*1, allocatable :: molnam(:,:)
      integer, allocatable :: lentry(:),list(:,:)
      integer, allocatable :: gstlentry(:),gstlist(:,:),locguest(:)
      integer, allocatable :: locfram(:)
      integer, allocatable :: lfreezesite(:),lfzsite(:,:)
      integer, allocatable :: nummols(:),numatoms(:)
      real(8), dimension(9) :: cell
      integer, allocatable :: lexatm(:,:),ltpsit(:,:)
      integer, allocatable :: ltype(:)
      integer, allocatable :: nexatm(:),noxatm(:),lexsit(:,:,:)
      integer, allocatable :: nexsit(:,:)
      integer, allocatable :: numfrz(:),numtyp(:),dens(:)
      integer, allocatable :: ins(:),del(:),dis(:),jmp(:),flx(:)
      integer, allocatable :: swp(:)
      integer, allocatable :: ind(:),ilist(:),jlist(:),angdist(:)
      integer, allocatable :: nprob(:),nprobsites(:),lprobsites(:,:)
      integer, allocatable :: accept_disp(:)
      integer, allocatable :: disp_count(:)
      integer, allocatable :: accept_tran(:)
      integer, allocatable :: tran_count(:)
      real(8), allocatable :: grid(:,:)
      real(8), allocatable :: grid_norm(:,:)
      real(8), allocatable :: dbuff(:),delE(:)
      real(8), allocatable :: statbuff(:),chainstats(:)
      real(8), allocatable :: gasmol(:), gaseng(:)
      real(8), allocatable :: ewald3en(:)
      real(8), allocatable :: energy(:),node_avg(:,:),nodeweight(:)
      real(8), allocatable :: node_std(:,:)
      real(8), allocatable :: avgwindow(:), varwindow(:), sumwindowav(:)
c     bunch of ewald arrays and parameters
      real(8), allocatable :: ckc(:),cks(:),clm(:),slm(:)
      real(8), allocatable :: ckcsum(:),ckssum(:)
      real(8), allocatable :: ckcsorig(:),ckssorig(:)
      real(8), allocatable :: ckcsnew(:),ckssnew(:)
      real(8), allocatable :: elc(:,:),els(:,:)
      real(8), allocatable :: emc(:,:),ems(:,:)
      real(8), allocatable :: erc(:),fer(:)
      real(8), allocatable :: enc(:,:),ens(:,:)
      real(8), allocatable :: mcinsf(:), mcdelf(:), mcdisf(:), mcjmpf(:)
      real(8), allocatable :: mcflxf(:), mcswpf(:), mctraf(:), mcrotf(:)
      real(8), allocatable :: mcmvnorm(:), mcmjpf(:)
      real(8), allocatable :: delrdisp(:), delr(:), disp_ratio(:)
      real(8), allocatable :: tran_ratio(:), tran_delr(:)
      real(8), allocatable :: adspres(:), adsguen(:)
c     Fugacity stuff
      real(8), allocatable :: gstpress(:)
      real(8), allocatable :: gstfuga(:), gstmolfract(:)
      real(8), allocatable :: Acc_factor(:), P_crit(:), T_crit(:)
      real(8), allocatable :: K_fug(:,:)

c     FIXED PARAMETERS

      integer, parameter :: mxexcl=50

c     max number of sites (atoms) for guest molecule
      integer, parameter :: mxguestsite=50
c     standard pi values

      real(8), parameter :: pi=3.141592653589793d0
      real(8), parameter :: sqrpi=1.7724538509055159d0

c     angstroms to bohr converter for .cube files

      real(8), parameter :: angs2bohr=1.889725989d0


c     min angle for gcmc rotation
      real(8), parameter :: minangle=pi/18.d0
c      real(8), parameter :: maxangle=pi/3.d0
c     min translation for gcmc displacement
      real(8), parameter :: delrmin=0.1d0
      real(8), parameter :: delrmax=3.d0

c     conversion factor for coulombic terms in internal units
c     i.e. (unit(charge)**2/(4 pi eps0 unit(length))/unit(energy)

      real(8), parameter :: r4pie0=138935.4835d0

c     boltzmann constant in internal units

      real(8), parameter :: iboltz=8.31451115d-1

c     boltzmann constant in kcal/mol/K units 

      real(8), parameter :: kboltz=1.9872041d-3

c     boltzmann constant in kg m^2 / (s^2 K)
      real(8), parameter :: boltz=1.3806503d-23 

c     Gas constant in terms of kcal / (mol K)
      real(8), parameter :: Rgas=8.314472d0/4184.d0

c     planck's constant in internal units

      real(8), parameter :: hbar=6.350780668d0

c     conversion factor for pressure from internal units to katm

      real(8), parameter :: prsunt=0.163882576d0

c     main input channel

      integer, parameter :: ncontrol=5

c     main output channel

      integer, parameter :: nrite=6

c     force field input channel

      integer, parameter :: nfield=9

c     configuration file input channel

      integer,parameter :: nconfig=10

c     history file input channel
 
      integer,parameter :: nhist=11

c     statistics file input channel

      integer,parameter :: nstats=12

      integer,parameter :: nrev=555


      save xxx,yyy,zzz
      save framwkxxx,framwkyyy,framwkzzz,frambuff
      save origframwkxxx,origframwkyyy,origframwkzzz
      save atmcharge,atmchg,guestx,guesty,guestz
      save atmweight, atmwght
      save atomname,atmname
      save lentry,list
      save lfreezesite,lfzsite
      save nummols,numatoms,molnam
      save cell,volum,nsite,chainstats,dbuff,gasmol,gaseng
      save gstlentry,gstlist,locguest,locfram,statbuff
      save avgwindow, varwindow, sumwindowav
      save gstpress,angdist,node_avg,node_std,nodeweight
      save mcinsf, mcdelf, mcdisf, mcjmpf, mcflxf, mcswpf
      save mctraf, mcrotf, mcmjpf, mcmvnorm
      save xdf,ydf,zdf,rsqdf
      save ins,del,dis,jmp,flx,swp,ewald3en
      save lexatm,nexatm,noxatm,ilist,jlist,unqatm,ltpsit
      save mxlist,ltype,mxgrid,mxegrd,lexsit,nexsit
      save ckc,cks,clm,slm,elc,emc,enc,els,ems,ens,erc,fer
      save numtyp,numfrz,dens,ind,newx,newy,newz,despre,packf
      save ckcsum,ckssum,ckcsnew,ckssnew,ckcsorig,ckssorig
      save nprob,nprobsites,lprobsites,grid,selsum,seltot,selwin2
      save grid_norm
      save accept_disp, disp_count, accept_tran, tran_count,adspres
      save delrdisp, delr, disp_ratio, tran_ratio, tran_delr,adsguen
      save delE,energy,total_pressure,gstfuga,gstmolfract,selvar
      save Acc_factor,P_crit,T_crit,K_fug,selectivity,selwin

      contains

      subroutine initscan
     &(idnode,imcon,volm,keyfce,rcut,eps,alpha,kmax1,kmax2,kmax3,lprob,
     &initdelr,rvdw,ntpguest,ntprob,ntpsite,ntpvdw,maxmls,mxatm,mxatyp,
     &griddim,gridfactor,maxguest,maxatm,despre)
c**********************************************************************
c
c     scans input files for relevant maximums to allocate to arrays
c     communicate to all nodes??
c**********************************************************************
      
      implicit none
      integer, parameter :: mmk=1000

      logical loop,loop2,loop3,loop4,loop5,safe,lewald,lcut
      logical lrvdw,check,ldelr,kill,lprob,lmaxg,lmaxf
      character*8 keyword
      character*8 name, chr(mmk)
      integer imcon,keyfce,idnode,idum,ntpvdw,maxmls,mxatm
      integer n,nummls,numsit,kmax1,kmax2,kmax3,mcsteps,maxatm
      integer mxatyp,nrept,ifrz,nneu,ksite,isite,eqsteps,maxguest
      integer i,j,ntpguest,ntprob,ntpsite,temp,gsite,iprob,qprob
      real(8) alpha,initdelr,rvdw,ppp,width
      real(8) fac,tol,tol1,rcut,eps,volm,despre
      real(8), dimension(10) :: celprp
      real(8) griddim
      integer, dimension(3) :: gridfactor
      mxatm=0
      mxatyp=0
      keyfce=0
      ntprob=0
      ntpsite=0
      ntpguest=0
      maxguest=3000
      maxatm=15000
      data loop/.true./,loop2/.false./,loop4/.false./,lewald/.false./
      data lrvdw/.false./,ldelr/.false./,kill/.false./,loop5/.false./
      data lmaxg/.false./

      if(idnode.eq.0)open(nfield,file='FIELD',status='old')

      call getrec(safe,idnode,nfield)
      if(.not.safe)call abort_field_read(1,idnode,nfield)
      do while(loop)
        call getrec(safe,idnode,nfield)
        if(.not.safe)call abort_field_read(1,idnode,nfield)
        call lowcase(record, lenrec)
        call strip(record, lenrec)
        if(record(1).eq.'#'.or.record(1).eq.' ')then
        elseif (findstring('molecu',record,idum))then
           maxmls=intstr(record,lenrec,idum)

           do n=1,maxmls
             loop2=.true.
            
             do while(loop2)
               call getrec(safe,idnode,nfield)
               if(.not.safe)call abort_field_read(1,idnode,nfield)
               call lowcase(record,lenrec)
               call strip(record,lenrec)
               ksite=0

               if(record(1).eq.'#'.or.record(1).eq.' ')then
               elseif (findstring('nummol',record,idum))then
                 nummls=intstr(record,lenrec,idum)
               elseif (findstring('atoms',record,idum))then
                 numsit=intstr(record,lenrec,idum)
                 mxatm=mxatm+numsit*nummls
                 ksite=0

                 do isite=1,numsit
                   if (ksite.lt.numsit)then
                      call getrec(safe,idnode,nfield)
                      if(.not.safe)call abort_field_read
     &                    (1,idnode,nfield)
                      call getword(name,record,8,lenrec)
                      ppp=dblstr(record,lenrec,idum)
                      ppp=dblstr(record,lenrec,idum)
                      nrept=intstr(record,lenrec,idum)
                      ifrz=intstr(record,lenrec,idum)
                      nneu=intstr(record,lenrec,idum)
                      if(nrept.eq.0)nrept=1
                      ksite=ksite+nrept
                      if(mxatyp.eq.0)then
                         mxatyp=1
                         chr(1)=name
                      else
                         check=.true.
                         do j=1,mxatyp
                           if(name.eq.chr(j))check=.false.
                         enddo
                         if(check)then
                           mxatyp=mxatyp+1
                           if(mxatyp.lt.mmk)chr(mxatyp)=name
                         endif
                      endif
                   endif
                 enddo
               elseif (findstring('finish',record,idum))then
                 loop2=.false.
               endif
             enddo
           enddo
        elseif(findstring('vdw',record,idum))then
           ntpvdw=intstr(record,lenrec,idum)
        elseif(findstring('close',record,idum))then
           loop=.false.
        endif
      enddo
c     grab the cell vectors - need these to allocate
c     some of the arrays.
      if(idnode.eq.0)open(nconfig,file='CONFIG',status='old')
      call getrec(safe,idnode,nconfig)
      call getrec(safe,idnode,nconfig)
      imcon=intstr(record,lenrec,idum)
      imcon=intstr(record,lenrec,idum)
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
      call dcell(cell,celprp)
      if(idnode.eq.0)open(ncontrol,file='CONTROL',status='old')
      loop3=.true.
      ntpguest=0
      do while(loop3)
        call getrec(safe,idnode,ncontrol)
        if(.not.safe)call abort_control_read(1,idnode,ncontrol)
        call lowcase(record,lenrec)
        call strip(record,lenrec)
        if (record(1).eq.'#'.or.record(1).eq.' ')then

        elseif(findstring('cut',record,idum))then
          rcut=dblstr(record,lenrec,idum)
          lcut=.true.
        elseif (findstring('max guest atoms',record,idum))then
          maxguest=intstr(record,lenrec,idum)
          lmaxg=.true.
        elseif (findstring('max framework atoms',record,idum))then
          maxatm=intstr(record,lenrec,idum)
          lmaxf=.true.
        elseif (findstring('&guest',record,idum))then
           ntpguest=ntpguest+1
           gsite=0 
           loop4=.true.
           
           do while(loop4)
         
             call getrec(safe,idnode,ncontrol)
             call lowcase(record,lenrec)
             call strip(record,lenrec)
             if(record(1).eq.'#'.or.record(1).eq.' ')then
c            record is commented out
             elseif(findstring('probability',record,idum))then
               iprob=0
               lprob=.true.
               qprob=intstr(record,lenrec,idum)
               ntprob=ntprob+qprob
               loop5=.true.
               do while(loop5) 
                 call getrec(safe,idnode,ncontrol)
                 if(.not.safe)call abort_control_read(1,idnode,ncontrol)
                 call strip(record,lenrec)
                 if(record(1).eq.'#'.or.record(1).eq.' ')then
                 else
                   iprob=iprob+1
                   temp=intstr(record,lenrec,idum)
                   gsite=gsite+temp
                 endif
                 if(iprob.eq.qprob)loop5=.false.
               enddo
               
             elseif(findstring('&end',record,idum))then
               loop4=.false.
               ntpsite=max(gsite,ntpsite)
             endif    
           enddo 

        elseif(findstring('finishment',record,idum))then
          despre=dblstr(record,lenrec,idum)
        elseif(findstring('delr',record,idum))then
          initdelr=dblstr(record,lenrec,idum)
          ldelr=.true.
        elseif(findstring('finish',record,idum))then
          loop3=.false.
        elseif(findstring('ewald',record,idum))then
          lewald=.true.
          keyfce=2
          if (findstring('precision',record,idum))then
             eps=dblstr(record,lenrec,idum)
             if (.not.lcut)then
                call error(idnode,-433)
                kill=.true. 
             else
               if(lewald)then
      
                  eps=min(abs(eps),0.5d0)
                  tol=sqrt(abs(log(eps*rcut)))
                  alpha=sqrt(abs(log(eps*rcut*tol)))/rcut
                  tol1=sqrt(-log(eps*rcut*(2.d0*tol*alpha)**2))
                  fac=1.d0
                  if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7)then
                     fac=2.d0**(1.d0/3.d0)
                  endif
                  kmax1=nint(0.25d0+fac*celprp(1)*alpha*tol1/pi)
                  kmax2=nint(0.25d0+fac*celprp(2)*alpha*tol1/pi)
                  kmax3=nint(0.25d0+fac*celprp(3)*alpha*tol1/pi)
               endif
             endif
          else
              alpha=dblstr(record,lenrec,idum)
              kmax1=intstr(record,lenrec,idum)
              kmax2=intstr(record,lenrec,idum)
              kmax3=intstr(record,lenrec,idum)
          endif
        elseif(findstring('rvdw',record,idum))then
          rvdw=dblstr(record,lenrec,idum)
          lrvdw=.true.
c Need to find the grid parameters in first scan, before allocating
        elseif (findstring('grid',record,idum))then
          if (findstring('spa',record,idum))then
            griddim = dblstr(record,lenrec,idum)
          elseif (findstring('fac',record,idum))then
            gridfactor(1) = intstr(record,lenrec,idum)
            gridfactor(2) = intstr(record,lenrec,idum)
            gridfactor(3) = intstr(record,lenrec,idum)
          endif

        endif
      enddo
      if(.not.ldelr)then
        call error(idnode,-433)
        kill=.true.
      endif

      width=min(celprp(7),celprp(8),celprp(9))/2.d0
      if(imcon.eq.4)width=sqrt(3.d0)*cell(1)/4.d0
      if(imcon.eq.5)width=cell(1)/2.d0
      if(imcon.eq.6)width=min(celprp(7),celprp(8))/2.d0

c     halt program if cutoff exceeds cell width

      if(rcut.gt.width)call error(idnode,95)
      
      call dcell(cell,celprp)
      if(imcon.eq.0) then
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

      if(.not.lrvdw)then
        if(rcut.gt.0d0)then
           rvdw=rcut
        endif
      endif
      if(.not.lewald)then
        call error(idnode,2311)
      endif
      if(idnode.eq.0)then
        close(ncontrol)  
        close(nconfig) 
        close(nfield)
      endif
      return
      end subroutine initscan 

      subroutine alloc_prob_arrays
     &(idnode,ntpguest,ntpsite,ntprob,gridsize)
c*********************************************************************
c
c      subroutine to allocate probability arrays
c
c*********************************************************************
      implicit none
      integer, parameter :: np=5
      integer ntpguest,ntprob,ntpsite,gridsize
      integer i,j, idnode
      integer, dimension(np) :: fail 
      do i=1,np
        fail(i) = 0
      enddo

      allocate(nprob(ntpguest),stat=fail(1))
      allocate(nprobsites(ntprob),stat=fail(2))
      allocate(lprobsites(ntprob,ntpsite),stat=fail(3))
      allocate(grid(ntprob,gridsize),stat=fail(4))
      allocate(grid_norm(ntprob,gridsize),stat=fail(5))

      do i=1,np
        if(fail(i).gt.0)then
            if(idnode.eq.0)write(nrite,'(10i5)')fail
            call error(idnode, 1001)
        endif
      enddo
      do i=1,ntprob
        nprobsites(i)=0
        do j=1,gridsize
          grid(i,j)=0
          grid_norm(i,j)=0
        enddo
      enddo

      return
      end subroutine alloc_prob_arrays

      subroutine storeprob
     &(ntpguest,rcell,ngrida,ngridb,ngridc,lnorm_bin)
c*********************************************************************
c
c      subroutine to store guest positions in a grid
c
c*********************************************************************
      implicit none
      logical lnorm_bin
      integer mol,np,iprob,iguest,i,maxprob,itprob
      integer jj,itmols,nmols,natms,itatm,jatm,val
      integer ngrida,ngridb,ngridc,ntpguest
      real(8) comx,comy,comz
      real(8), dimension(9) :: rcell 

      itprob=0

c     iterate over all molecules (outer loop) and atoms
c     (inner loop), store each guest in a temp array,
c     get fractional coordinates.
      do iguest=1,ntpguest
        np=nprob(iguest)
        mol=locguest(iguest)
        nmols=nummols(mol)
        natms=numatoms(mol)
        jj=0
        do itmols=1,nmols
          do itatm=1,natms
            jj=jj+1
            newx(itatm)=framwkxxx(mol,jj)
            newy(itatm)=framwkyyy(mol,jj)
            newz(itatm)=framwkzzz(mol,jj)
          enddo

c         center of mass - calculated regardless of COM request for
c         now..

          call com(natms,mol,newx,newy,newz,comx,comy,comz)

c         convert centre of mass to fractional coordinates

          call frac2(comx,comy,comz,rcell)

c         convert guest atoms to fractional coordinates
          call frac(newx,newy,newz,natms,rcell)

c       iterate over number of probability plots

          do iprob=1,np
            if(.not.lnorm_bin)then
              do i=1,nprobsites(itprob+iprob)
                if(lprobsites(itprob+iprob,i).eq.0)then
                  call equitable_bin(comx,comy,comz,ngrida,ngridb,
     &ngridc,itprob+iprob)
                else
                  jatm=lprobsites(itprob+iprob,i)
                  call equitable_bin(newx(jatm),newy(jatm),newz(jatm),
     &ngrida,ngridb,ngridc,itprob+iprob)
                endif
              enddo
            else
              do i=1,nprobsites(itprob+iprob)
                if(lprobsites(itprob+iprob,i).eq.0)then
                  call equitable_bin(comx,comy,comz,ngrida,ngridb,
     &ngridc,itprob+iprob)
                  call normal_bin(comx,comy,comz,ngrida,ngridb,
     &ngridc,itprob+iprob)
                else
                  jatm=lprobsites(itprob+iprob,i)
                  call equitable_bin(newx(jatm),newy(jatm),newz(jatm),
     &ngrida,ngridb,ngridc,itprob+iprob)
                  call normal_bin(newx(jatm),newy(jatm),newz(jatm),
     &ngrida,ngridb,ngridc,itprob+iprob)
                endif
              enddo
            endif
          enddo
        enddo
        itprob=itprob+np
      enddo

      return
      end subroutine storeprob

      subroutine writeprob
     &(iguest,itprob,iprob,cell,ntpguest,ntpfram,gridsize,
     &ngrida,ngridb,ngridc,steps)
c*********************************************************************
c
c     write gaussian .cube files for visualization
c
c*********************************************************************
      implicit none
      character*25 filename
      real(8) vol
      integer itprob,ntpguest,iguest,i,j,ii
      integer igrid,k,l,ngrida,ngridb,ngridc,ntpfram
      integer nfram,framol,nmol,natms,iatm,m,n,o,p
      integer gridsize,steps,iprob,ip
      real(8), dimension(9) :: cell

      vol=(cell(2)*cell(6)-cell(3)*cell(5))/(ngrida*ngridb)*
     &cell(7)/ngridc+
     &(cell(3)*cell(4)-cell(1)*cell(6))/(ngrida*ngridb)*
     &cell(8)/ngridc+
     &(cell(1)*cell(5)-cell(2)*cell(4))/(ngrida*ngridb)*
     &cell(9)/ngridc
c     get number of framework atoms
      nfram=0
      do i=1,ntpfram
        framol=locfram(i)
        nfram=nfram+nummols(framol)*numatoms(framol)
      enddo
 
      write(filename,"('prob_guest',i2.2,'_prob_',i2.2,'.cube')")
     &iguest,iprob
      open(i,file=filename)

c     write header stuff for cube
      write(i,"('probability cube file',/,'outer loop a, middle loop
     & b, inner loop c')")
      write(i,'(i6, 3f12.6)')nfram,0.d0,0.d0,0.d0
      write(i,'(i6,3f12.6)')ngrida,cell(1)/ngrida,cell(2)/ngrida,
     &cell(3)/ngrida
      write(i,'(i6,3f12.6)')ngridb,cell(4)/ngridb,cell(5)/ngridb,
     &cell(6)/ngridb
      write(i,'(i6,3f12.6)')ngridc,cell(7)/ngridc,cell(8)/ngridc,
     &cell(9)/ngridc
c     write cartesians of framework
      do m=1,ntpfram
        framol=locfram(m)
        nmol=nummols(framol)
        natms=numatoms(framol)
        iatm=0
        do n=1,nmol
          do o=1,natms
            iatm=iatm+1
            write(i,'(i6,4f12.6)')atmnumber(atmwght(framol,o)),0.d0,
     &               framwkxxx(framol,iatm)*angs2bohr,
     &               framwkyyy(framol,iatm)*angs2bohr,
     &               framwkzzz(framol,iatm)*angs2bohr

          enddo
        enddo
      enddo
      igrid=0

c     inner loop is c vec, middle loop is b vec, outter loop is a vec.
c     this has been intrinsically stored in the array grid.
c     the cube file writes six entries per line *except* when it reaches
c     the end of the number of grid points in the c direction (fastest
c     loop).  In this case it will end the line before reaching the
c     sixth entry.  A nuisance.
      ip=0
      do j=1,gridsize
        ip=ip+1 
        write(i,'(e15.6)',advance='no')
     &(dble(grid(itprob,j))/vol/dble(steps))

c       conditions for writing to a new line 
        if((mod(j,ngridc).eq.0).or.(ip.eq.6))then
           write(i,'(x)')
           ip=0
        endif
      enddo
      close(i)
  

      return
      end subroutine writeprob

      subroutine writeprobnormalbin
     &(iguest,itprob,iprob,cell,ntpguest,ntpfram,gridsize,
     &ngrida,ngridb,ngridc,steps)
c*********************************************************************
c
c     write gaussian .cube files for visualization
c
c*********************************************************************
      implicit none
      character*34 filenamenormbin
      real(8) vol
      integer itprob,ntpguest,iguest,i,j,ii
      integer igrid,k,l,ngrida,ngridb,ngridc,ntpfram
      integer nfram,framol,nmol,natms,iatm,m,n,o,p
      integer gridsize,steps,iprob,ip
      real(8), dimension(9) :: cell

      vol=(cell(2)*cell(6)-cell(3)*cell(5))/(ngrida*ngridb)*
     &cell(7)/ngridc+
     &(cell(3)*cell(4)-cell(1)*cell(6))/(ngrida*ngridb)*
     &cell(8)/ngridc+
     &(cell(1)*cell(5)-cell(2)*cell(4))/(ngrida*ngridb)*
     &cell(9)/ngridc
c     get number of framework atoms
      nfram=0
      do i=1,ntpfram
        framol=locfram(i)
        nfram=nfram+nummols(framol)*numatoms(framol)
      enddo

      write(filenamenormbin,"('prob_guest',i2.2,'_prob_',i2.2,
     &'_norm_bin.cube')")
     &iguest,iprob
      open(i,file=filenamenormbin)

c     write header stuff for cube
      write(i,"('probability cube file',/,'outer loop a, middle loop
     & b, inner loop c')")
      write(i,'(i6, 3f12.6)')nfram,0.d0,0.d0,0.d0
      write(i,'(i6,3f12.6)')ngrida,cell(1)/ngrida,cell(2)/ngrida,
     &cell(3)/ngrida
      write(i,'(i6,3f12.6)')ngridb,cell(4)/ngridb,cell(5)/ngridb,
     &cell(6)/ngridb
      write(i,'(i6,3f12.6)')ngridc,cell(7)/ngridc,cell(8)/ngridc,
     &cell(9)/ngridc
c     write cartesians of framework
      do m=1,ntpfram
        framol=locfram(m)
        nmol=nummols(framol)
        natms=numatoms(framol)
        iatm=0
        do n=1,nmol
          do o=1,natms
            iatm=iatm+1
            write(i,'(i6,4f12.6)')atmnumber(atmwght(framol,o)),0.d0,
     &               framwkxxx(framol,iatm)*angs2bohr,
     &               framwkyyy(framol,iatm)*angs2bohr,
     &               framwkzzz(framol,iatm)*angs2bohr

          enddo
        enddo
      enddo
      igrid=0

c     inner loop is c vec, middle loop is b vec, outter loop is a vec.
c     this has been intrinsically stored in the array grid.
c     the cube file writes six entries per line *except* when it reaches
c     the end of the number of grid points in the c direction (fastest
c     loop).  In this case it will end the line before reaching the
c     sixth entry.  A nuisance.
      ip=0
      do j=1,gridsize
        ip=ip+1 
        write(i,'(e15.6)',advance='no')
     &(dble(grid_norm(itprob,j))/vol/dble(steps))

c       conditions for writing to a new line 
        if((mod(j,ngridc).eq.0).or.(ip.eq.6))then
          write(i,'(x)')
          ip=0
        endif
      enddo
      close(i)

      return
      end subroutine writeprobnormalbin

      subroutine frac(x,y,z,natms,rcell)
c*********************************************************************
c
c     convert cartesian coordinates to fractional coordinates 
c     haven't checked if it works with other cell types   
c     than parallelpiped (imcon.eq.3) and rectangular (imcon.eq.1 or 2)
c     PB
c 
c*********************************************************************
      implicit none
      integer i,natms
      real(8) ssx,ssy,ssz
      real(8), dimension(natms) :: x,y,z
      real(8), dimension(9) :: rcell


      do i=1,natms

        ssx=x(i)*rcell(1)+y(i)*rcell(4)+z(i)*rcell(7)
        ssy=x(i)*rcell(2)+y(i)*rcell(5)+z(i)*rcell(8)
        ssz=x(i)*rcell(3)+y(i)*rcell(6)+z(i)*rcell(9)
        x(i)=modulo(ssx, 1.0)
        y(i)=modulo(ssy, 1.0)
        z(i)=modulo(ssz, 1.0)
      enddo

      return
      end subroutine frac

      subroutine frac2(x,y,z,rcell)
c*********************************************************************
c
c     convert cartesian coordinates to fractional coordinates
c     for one atom only (this is a cheap fix for the COM calcluation
c     in the subroutine storeprob)
c     haven't checked if it works with other cell types   
c     than parallelpiped (imcon.eq.3) and rectangular (imcon.eq.1 or 2)
c     PB
c 
c*********************************************************************
      implicit none
      integer i,natms
      real(8) ssx,ssy,ssz
      real(8) x,y,z
      real(8), dimension(9) :: rcell

      ssx=x*rcell(1)+y*rcell(4)+z*rcell(7)
      ssy=x*rcell(2)+y*rcell(5)+z*rcell(8)
      ssz=x*rcell(3)+y*rcell(6)+z*rcell(9)
      x=modulo(ssx, 1.0)
      y=modulo(ssy, 1.0)
      z=modulo(ssz, 1.0)
 
      return
      end subroutine frac2

      subroutine alloc_config_arrays
     &(idnode,mxnode,maxmls,mxatm,mxatyp,volm,kmax2,kmax3,
     &mxebuf,mxewld,ntpguest,rcut,rvdw,tempdelr,maxguest)
c*********************************************************************
c
c     allocation of some arrays
c   
c*********************************************************************
      implicit none
      integer, parameter :: na = 114
      integer maxmls,mxatm,maxalloc,iatm,ntpguest,i,kk,j,k
      integer kmax2,kmax3,mxatyp,mxebuf,mxewld,idnode,mxnode,maxguest
      real(8) density,ratio,cut,volm,rcut,rvdw,tempdelr
      integer, dimension(na) :: fail

c     initialize fail check array
      do i=1,na
        fail(i) = 0
      enddo
      mxgrid=max(1000,int(rvdw/0.01d0+0.5d0)+4)
      mxegrd=mxgrid 

      maxalloc=mxatm+maxguest
c     the following line was added to test the size of maxalloc
c     which seemed to be giving insufficient virutal memory errors
c      if(idnode.eq.0)write(nrite,"('maxalloc: ', i9)")maxalloc
      density=dble(maxalloc)/volm
      cut=rcut+tempdelr
      ratio=1.5d0*density*(4.d0*pi/3.d0)*cut**3
      mxlist=min(nint(ratio),(maxalloc+1)/2)
      allocate(locguest(maxmls),stat=fail(1))
      allocate(guestx(ntpguest,mxguestsite),stat=fail(2))
      allocate(guesty(ntpguest,mxguestsite),stat=fail(3))
      allocate(guestz(ntpguest,mxguestsite),stat=fail(4))
      allocate(locfram(maxmls),stat=fail(5))
      allocate(list(maxalloc,mxlist),stat=fail(6))
      allocate(lentry(maxalloc),stat=fail(7))
      allocate(noxatm(maxalloc),stat=fail(8))
      allocate(nexatm(maxalloc),stat=fail(9))
      allocate(lexatm(maxalloc,mxexcl),stat=fail(10))
      allocate(nexsit(mxguestsite,mxexcl),stat=fail(11))
      allocate(lexsit(mxguestsite,mxguestsite,mxexcl),stat=fail(12))
      allocate(xxx(maxalloc),stat=fail(13))
      allocate(yyy(maxalloc),stat=fail(14))
      allocate(zzz(maxalloc),stat=fail(15))
      allocate(molnam(40,maxmls),stat=fail(16))
      allocate(nummols(maxmls),stat=fail(17))
      allocate(numatoms(maxmls),stat=fail(18))
      allocate(atmname(maxmls,maxalloc),stat=fail(19))
      allocate(atmchg(maxmls,maxalloc),stat=fail(20))
      allocate(atmwght(maxmls,maxalloc),stat=fail(21))
      allocate(lfzsite(maxmls,maxalloc),stat=fail(22))
      allocate(framwkxxx(maxmls,maxalloc),stat=fail(23))
      allocate(framwkyyy(maxmls,maxalloc),stat=fail(24))
      allocate(framwkzzz(maxmls,maxalloc),stat=fail(25))
      allocate(origframwkxxx(maxmls,maxalloc),stat=fail(26))
      allocate(origframwkyyy(maxmls,maxalloc),stat=fail(27))
      allocate(origframwkzzz(maxmls,maxalloc),stat=fail(28))
      allocate(frambuff(maxmls,maxalloc),stat=fail(29))
      allocate(lfreezesite(maxalloc),stat=fail(30))
      allocate(atmcharge(maxalloc),stat=fail(31))
      allocate(atomname(maxalloc),stat=fail(32))
      allocate(atmweight(maxalloc),stat=fail(33))
      allocate(gstlentry(mxguestsite),stat=fail(34))
      allocate(gstpress(ntpguest),stat=fail(35))
      allocate(mcinsf(ntpguest), stat=fail(36))
      allocate(mcdelf(ntpguest), stat=fail(37))
      allocate(mcdisf(ntpguest), stat=fail(38))
      allocate(mcjmpf(ntpguest), stat=fail(39))
      allocate(mcflxf(ntpguest), stat=fail(40))
      allocate(mcswpf(ntpguest), stat=fail(41))
      allocate(mctraf(ntpguest), stat=fail(42))
      allocate(mcrotf(ntpguest), stat=fail(43))
      allocate(mcmjpf(ntpguest), stat=fail(44))
      allocate(mcmvnorm(ntpguest), stat=fail(45))
      allocate(accept_disp(ntpguest), stat=fail(46))
      allocate(disp_count(ntpguest), stat=fail(47))
      allocate(accept_tran(ntpguest), stat=fail(48))
      allocate(tran_count(ntpguest), stat=fail(49))
      allocate(delrdisp(ntpguest), stat=fail(50))
      allocate(delr(ntpguest), stat=fail(51))
      allocate(disp_ratio(ntpguest), stat=fail(52))
      allocate(tran_ratio(ntpguest), stat=fail(53))
      allocate(tran_delr(ntpguest), stat=fail(54))
      allocate(statbuff(1+ntpguest*14),stat=fail(55))
      allocate(chainstats(1+ntpguest*14),stat=fail(56))
      allocate(gstlist(mxguestsite,maxalloc),stat=fail(57))
      allocate(ckc(maxalloc),stat=fail(58))
      allocate(cks(maxalloc),stat=fail(59))
      allocate(ckcsum(mxebuf),stat=fail(60))
      allocate(ckssum(mxebuf),stat=fail(61))
      allocate(ckcsorig(mxebuf),stat=fail(62))
      allocate(ckssorig(mxebuf),stat=fail(63))
      allocate(ckcsnew(mxebuf),stat=fail(64))
      allocate(ckssnew(mxebuf),stat=fail(65))
      allocate(clm(maxalloc),stat=fail(66))
      allocate(slm(maxalloc),stat=fail(67))
      allocate(elc(maxalloc,0:1),stat=fail(68))
      allocate(els(maxalloc,0:1),stat=fail(69))
      allocate(emc(maxalloc,0:kmax2),stat=fail(70))
      allocate(ems(maxalloc,0:kmax2),stat=fail(71))
      allocate(enc(maxalloc,0:kmax3),stat=fail(72))
      allocate(ens(maxalloc,0:kmax3),stat=fail(73))
      allocate(erc(mxegrd),stat=fail(74))
      allocate(fer(mxegrd),stat=fail(75))
      allocate(ilist(maxalloc),stat=fail(76))
      allocate(unqatm(maxalloc),stat=fail(77))
      allocate(jlist(maxalloc),stat=fail(78))
      allocate(ltpsit(maxmls,maxalloc),stat=fail(79))
      allocate(ltype(maxalloc),stat=fail(80))
      allocate(xdf(maxalloc),stat=fail(81))
      allocate(ydf(maxalloc),stat=fail(82))
      allocate(zdf(maxalloc),stat=fail(83))
      allocate(rsqdf(maxalloc),stat=fail(84))
      allocate(numtyp(mxatyp),stat=fail(85))
      allocate(numfrz(mxatyp),stat=fail(86))
      allocate(dens(mxatyp),stat=fail(87))
      allocate(newx(mxguestsite),stat=fail(88))
      allocate(newy(mxguestsite),stat=fail(89))
      allocate(newz(mxguestsite),stat=fail(90))
      allocate(ind(mxguestsite),stat=fail(91))
      allocate(energy(ntpguest),stat=fail(92))
      allocate(delE(ntpguest),stat=fail(93))
      allocate(avgwindow(ntpguest*8),stat=fail(94))
      allocate(sumwindowav(ntpguest*8),stat=fail(95))
      allocate(varwindow(ntpguest*8),stat=fail(96))
      allocate(nodeweight(mxnode),stat=fail(97))
      allocate(node_avg(mxnode,ntpguest*8),stat=fail(98))
      allocate(node_std(mxnode,ntpguest*8),stat=fail(99))
      allocate(ins(ntpguest),stat=fail(100))
      allocate(del(ntpguest),stat=fail(101))
      allocate(dis(ntpguest),stat=fail(102))
      allocate(jmp(ntpguest),stat=fail(103))
      allocate(flx(ntpguest),stat=fail(104))
      allocate(swp(ntpguest),stat=fail(105))
      allocate(ewald3en(ntpguest),stat=fail(106))
      allocate(gstfuga(ntpguest),stat=fail(107))
      allocate(gstmolfract(ntpguest),stat=fail(108))
      allocate(Acc_factor(ntpguest),stat=fail(109))
      allocate(T_crit(ntpguest),stat=fail(110))
      allocate(P_crit(ntpguest),stat=fail(111))
      allocate(K_fug(ntpguest, ntpguest),stat=fail(112))
      allocate(gasmol(ntpguest*500000),stat=fail(113))
      allocate(gaseng(ntpguest*500000),stat=fail(114))
      allocate(selectivity(ntpguest,ntpguest),stat=fail(115))
      allocate(selwin(ntpguest,ntpguest),stat=fail(116))
      allocate(selvar(ntpguest,ntpguest),stat=fail(117))
      allocate(selsum(ntpguest,ntpguest),stat=fail(118))
      allocate(seltot(ntpguest,ntpguest),stat=fail(119))
      allocate(selwin2(ntpguest,ntpguest),stat=fail(120))
      allocate(adspres(ntpguest),stat=fail(121))
      allocate(adsguen(ntpguest),stat=fail(122))
      do i=1,na
        if(fail(i).gt.0)then
            if(idnode.eq.0)write(nrite,'(10i5)')fail(i)
            call error(idnode, 1000)
        endif
      enddo

c     initialize statistic arrays
c     probably put all this initialization stuff in a separate
c     module
      do i=1,mxnode
        nodeweight(i)=0.d0
        do j=1,ntpguest
          node_avg(i,1:8) = 0.d0
          node_std(i,1:8) = 0.d0
        enddo
      enddo
      do i = 1, ntpguest
        do j = 1, ntpguest
            K_fug(i,j) = 0.d0
            selectivity(i,j) = 0.d0
            selwin(i,j) = 0.d0
            selvar(i,j) = 0.d0
            selsum(i,j) = 0.d0
            seltot(i,j) = 0.d0
            selwin2(j,i) = 0.d0
        enddo
      enddo
      do i=1,ntpguest*8
        avgwindow(i)=0.d0
        varwindow(i)=0.d0
        sumwindowav(i)=0.d0
      enddo
c     initialize the stat arrays
      chainstats(1:1+ntpguest*14)=0.d0
      ins(1:ntpguest) = 0.d0
      del(1:ntpguest) = 0.d0
      dis(1:ntpguest) = 0.d0
      jmp(1:ntpguest) = 0.d0
      flx(1:ntpguest) = 0.d0
      swp(1:ntpguest) = 0.d0
      energy(1:ntpguest) = 0.d0
      delE(1:ntpguest) = 0.d0
      ewald3en(1:ntpguest) = 0.d0
      gasmol(1:ntpguest*500000)=0.d0
      gaseng(1:ntpguest*500000)=0.d0
      return
      end subroutine alloc_config_arrays


      subroutine condensefram(totfram,ntpfram)
c**********************************************************************
c
c     subroutine to allocate all x,y,z coordinates and relevant
c     atomic data to 1-d arrays.  This is necessary for easier
c     calculations of neighbour lists ewald sums etc.
c*********************************************************************

      implicit none
      integer iatm,ntpfram,mol
      integer i,j,k,imol,nmols,totfram
      iatm=0
      imol=0
      totfram=0

      do i=1,ntpfram
        mol=locfram(i)
        nmols=nummols(mol)
        do j=1,nmols
          do k=1,numatoms(mol)
             totfram=totfram+1
             iatm=iatm+1
             imol=imol+1
             xxx(iatm)=framwkxxx(mol,imol)
             yyy(iatm)=framwkyyy(mol,imol)
             zzz(iatm)=framwkzzz(mol,imol)
             atmcharge(iatm)=atmchg(mol,k)
             lfreezesite(iatm)=lfzsite(mol,k)
             atomname(iatm)=atmname(mol,k)
             atmweight(iatm)=atmwght(mol,k)
             ltype(iatm)=ltpsit(mol,k)
          enddo
        enddo

        imol=0
      enddo

      return

      end subroutine condensefram

      subroutine floorimages
     x  (imcon,natms,cell,xxx,yyy,zzz)
      
c***********************************************************************
c     this will adjust the coordinates to fit within the original
c     box.
c
c     for
c     imcon=0 no boundary conditions apply
c     imcon=1 standard cubic boundaries apply
c     imcon=2 orthorhombic boundaries apply
c     imcon=3 parallelepiped boundaries apply
c     imcon=4 truncated octahedron boundaries apply
c     imcon=5 rhombic dodecahedron boundaries apply
c     imcon=6 x-y parallelogram boundary conditions : no periodicity in z
c     imcon=7 hexagonal prism boundaries apply
c     
c     
c***********************************************************************
      
      implicit none

      integer imcon,natms,i
      real(8) aaa,bbb,ccc,det,rt2,rt3,ssx
      real(8) ssy,ssz,ddd,xss,yss,zss

      real(8), dimension(natms):: xxx,yyy,zzz,sxx,syy,szz
      real(8), dimension(9) :: cell,rcell

      data rt2/1.41421356623d0/,rt3/1.7320508075d0/
      if(imcon.eq.1)then

c     standard cubic boundary conditions
        
        
        aaa=1.d0/cell(1)

        do i=1,natms
          xxx(i)=xxx(i)-cell(1)*floor(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*floor(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*floor(aaa*zzz(i))
        enddo
        
      else if(imcon.eq.2)then

c     rectangular (slab) boundary conditions
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(5)
        ccc=1.d0/cell(9)
        
        do i=1,natms
          
          xxx(i)=xxx(i)-cell(1)*floor(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(5)*floor(bbb*yyy(i))
          zzz(i)=zzz(i)-cell(9)*floor(ccc*zzz(i))
          
        enddo
        
      else if(imcon.eq.3)then

c     parallelpiped boundary conditions

        call invert(cell,rcell,det)
        do i=1,natms
          ssx=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))
          ssy=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))
          ssz=(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i))
          
          xss=ssx-floor(ssx)
          yss=ssy-floor(ssy)
          zss=ssz-floor(ssz)
          
          xxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
          yyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
          zzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
          
        enddo
      else if(imcon.eq.4)then

c     truncated octahedral boundary conditions
        
c        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
c     x    abs(cell(5)-cell(9)).lt.1.d-6)) call error(idnode,130)
        
        aaa=1.d0/cell(1)
        
        do i=1,natms
          
          xxx(i)=xxx(i)-cell(1)*floor(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*floor(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*floor(aaa*zzz(i))
          
          if((abs(sxx(i))+abs(syy(i))+abs(szz(i))).ge.
     x      (0.75d0*cell(1)))then
            
            xxx(i)=sxx(i)-0.5d0*sign(cell(1),sxx(i))
            yyy(i)=syy(i)-0.5d0*sign(cell(1),syy(i))
            zzz(i)=szz(i)-0.5d0*sign(cell(1),szz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.5)then

c     rhombic dodecahedral boundary conditions
        
c        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
c     x    abs(cell(9)-cell(1)*rt2).lt.1.d-6)) 
c     x    call error(idnode,140)
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(9)
        
        do i=1,natms
          
          xxx(i)=xxx(i)-cell(1)*floor(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*floor(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(9)*floor(bbb*zzz(i))
          
          if((abs(sxx(i))+abs(syy(i))+abs(rt2*szz(i))).ge.
     x      cell(1))then
            
            xxx(i)=sxx(i)-0.5d0*sign(cell(1),sxx(i))
            yyy(i)=syy(i)-0.5d0*sign(cell(1),syy(i))
            zzz(i)=szz(i)-0.5d0*sign(cell(9),szz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.6) then

c     x-y boundary conditions 

        det = cell(1)*cell(5) - cell(2)*cell(4)

c        if(abs(det).lt.1.d-6)call error(idnode,120)
        
        det = 1.d0/det

        rcell(1) =  det*cell(5)
        rcell(2) = -det*cell(2)
        rcell(4) = -det*cell(4)
        rcell(5) =  det*cell(1)
        
        do i=1,natms

          ssx = rcell(1)*xxx(i) + rcell(4)*yyy(i)
          ssy = rcell(2)*xxx(i) + rcell(5)*yyy(i)

          xss = ssx - floor(ssx)
          yss = ssy - floor(ssy)

          xxx(i)=cell(1)*xss + cell(4)*yss
          yyy(i)=cell(2)*xss + cell(5)*yss

        enddo

      else if(imcon.eq.7) then

c     hexagonal prism boundary conditions
        
c        if(abs(cell(1)-rt3*cell(5)).ge.1.d-6)
c     x    call error(idnode,135)
        
        aaa=cell(1)/(rt3*2.d0)
        bbb=cell(1)/rt3
        ccc=rt3/cell(1)
        ddd=1.d0/cell(9)
        
        do i=1,natms
          
          yyy(i)=yyy(i)-bbb*floor(ccc*yyy(i))
          zzz(i)=zzz(i)-cell(9)*floor(ddd*zzz(i))
          
          if((abs(syy(i))+abs(rt3*sxx(i))).ge.bbb)then
            
            xxx(i)=sxx(i)-rt3*sign(aaa,sxx(i))
            yyy(i)=syy(i)-sign(aaa,syy(i))
            
          endif
          
        enddo
        
      endif
      return
      end subroutine floorimages
      subroutine images
     x  (imcon,natms,cell,sxxx,syyy,szzz)
      
c***********************************************************************
c     
c     subroutine for calculating the minimum image
c     of atom pairs within a specified MD cell
c
c     note the following internal units apply everywhere
c     
c     unit of time      (to)    =          1 x 10**(-12) seconds
c     unit of length    (lo)    =          1 x 10**(-10) metres
c     unit of mass      (mo)    = 1.6605402  x 10**(-27) kilograms
c     unit of charge    (qo)    = 1.60217733 x 10**(-19) coulombs
c     unit of energy    (eo)    = 1.6605402  x 10**(-23) joules
c     unit of pressure  (po)    = 1.6605402  x 10**(  7) pascals
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     T3D optimised version. t.forester july 1994
c     
c     for
c     imcon=0 no boundary conditions apply
c     imcon=1 standard cubic boundaries apply
c     imcon=2 orthorhombic boundaries apply
c     imcon=3 parallelepiped boundaries apply
c     imcon=4 truncated octahedron boundaries apply
c     imcon=5 rhombic dodecahedron boundaries apply
c     imcon=6 x-y parallelogram boundary conditions : no periodicity in z
c     imcon=7 hexagonal prism boundaries apply
c     
c     note: in all cases the centre of the cell is at (0,0,0)
c     warning - replicated data version: does not re-merge 
c     coordinate arrays
c     
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c     
c***********************************************************************
      
      implicit none

      integer imcon,natms,i
      real(8) cell,sxxx,syyy,szzz,aaa,bbb,ccc,det,rt2,rt3,ssx
      real(8) ssy,ssz,ddd,xss,yss,zss,rcell

      dimension sxxx(*),syyy(*),szzz(*)
      dimension cell(9),rcell(9)

      data rt2/1.41421356623d0/,rt3/1.7320508075d0/
      if(imcon.eq.1)then

c     standard cubic boundary conditions
        
        
        aaa=1.d0/cell(1)

        do i=1,natms
          sxxx(i)=sxxx(i)-cell(1)*nint(aaa*sxxx(i))
          syyy(i)=syyy(i)-cell(1)*nint(aaa*syyy(i))
          szzz(i)=szzz(i)-cell(1)*nint(aaa*szzz(i))
        enddo
        
      else if(imcon.eq.2)then

c     rectangular (slab) boundary conditions
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(5)
        ccc=1.d0/cell(9)
        
        do i=1,natms
          
          sxxx(i)=sxxx(i)-cell(1)*nint(aaa*sxxx(i))
          syyy(i)=syyy(i)-cell(5)*nint(bbb*syyy(i))
          szzz(i)=szzz(i)-cell(9)*nint(ccc*szzz(i))
          
        enddo
        
      else if(imcon.eq.3)then

c     parallelpiped boundary conditions
        call invert(cell,rcell,det)
        do i=1,natms
          ssx=(rcell(1)*sxxx(i)+rcell(4)*syyy(i)+rcell(7)*szzz(i))
          ssy=(rcell(2)*sxxx(i)+rcell(5)*syyy(i)+rcell(8)*szzz(i))
          ssz=(rcell(3)*sxxx(i)+rcell(6)*syyy(i)+rcell(9)*szzz(i))
          
          xss=ssx-nint(ssx)
          yss=ssy-nint(ssy)
          zss=ssz-nint(ssz)
          
          sxxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
          syyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
          szzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
          
        enddo
      else if(imcon.eq.4)then

c     truncated octahedral boundary conditions
        
c        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
c     x    abs(cell(5)-cell(9)).lt.1.d-6)) call error(idnode,130)
        
        aaa=1.d0/cell(1)
        
        do i=1,natms
          
          sxxx(i)=sxxx(i)-cell(1)*nint(aaa*sxxx(i))
          syyy(i)=syyy(i)-cell(1)*nint(aaa*syyy(i))
          szzz(i)=szzz(i)-cell(1)*nint(aaa*szzz(i))
          
          if((abs(sxxx(i))+abs(syyy(i))+abs(szzz(i))).ge.
     x      (0.75d0*cell(1)))then
            
            sxxx(i)=sxxx(i)-0.5d0*sign(cell(1),sxxx(i))
            syyy(i)=syyy(i)-0.5d0*sign(cell(1),syyy(i))
            szzz(i)=szzz(i)-0.5d0*sign(cell(1),szzz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.5)then

c     rhombic dodecahedral boundary conditions
        
c        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
c     x    abs(cell(9)-cell(1)*rt2).lt.1.d-6)) 
c     x    call error(idnode,140)
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(9)
        
        do i=1,natms
          
          sxxx(i)=sxxx(i)-cell(1)*nint(aaa*sxxx(i))
          syyy(i)=syyy(i)-cell(1)*nint(aaa*syyy(i))
          szzz(i)=szzz(i)-cell(9)*nint(bbb*szzz(i))
          
          if((abs(sxxx(i))+abs(syyy(i))+abs(rt2*szzz(i))).ge.
     x      cell(1))then
            
            sxxx(i)=sxxx(i)-0.5d0*sign(cell(1),sxxx(i))
            syyy(i)=syyy(i)-0.5d0*sign(cell(1),syyy(i))
            szzz(i)=szzz(i)-0.5d0*sign(cell(9),szzz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.6) then

c     x-y boundary conditions 

        det = cell(1)*cell(5) - cell(2)*cell(4)

c        if(abs(det).lt.1.d-6)call error(idnode,120)
        
        det = 1.d0/det

        rcell(1) =  det*cell(5)
        rcell(2) = -det*cell(2)
        rcell(4) = -det*cell(4)
        rcell(5) =  det*cell(1)
        
        do i=1,natms

          ssx = rcell(1)*sxxx(i) + rcell(4)*syyy(i)
          ssy = rcell(2)*sxxx(i) + rcell(5)*syyy(i)

          xss = ssx - nint(ssx)
          yss = ssy - nint(ssy)

          sxxx(i)=cell(1)*xss + cell(4)*yss
          syyy(i)=cell(2)*xss + cell(5)*yss

        enddo

      else if(imcon.eq.7) then

c     hexagonal prism boundary conditions
        
c        if(abs(cell(1)-rt3*cell(5)).ge.1.d-6)
c     x    call error(idnode,135)
        
        aaa=cell(1)/(rt3*2.d0)
        bbb=cell(1)/rt3
        ccc=rt3/cell(1)
        ddd=1.d0/cell(9)
        
        do i=1,natms
          
          syyy(i)=syyy(i)-bbb*nint(ccc*syyy(i))
          szzz(i)=szzz(i)-cell(9)*nint(ddd*szzz(i))
          
          if((abs(syyy(i))+abs(rt3*sxxx(i))).ge.bbb)then
            
            sxxx(i)=sxxx(i)-rt3*sign(aaa,sxxx(i))
            syyy(i)=syyy(i)-sign(aaa,syyy(i))
            
          endif
          
        enddo
        
      endif
      return
      end subroutine images

      subroutine dcell(aaa,bbb)

c***********************************************************************
c     
c     dl_poly subroutine to calculate the dimensional properties of
c     a simulation cell specified by the input matrix aaa.
c     the results are returned in the array bbb, with :
c     
c     bbb(1 to 3) - lengths of cell vectors
c     bbb(4 to 6) - cosines of cell angles
c     bbb(7 to 9) - perpendicular cell widths
c     bbb(10)     - cell volume
c     
c     copyright daresbury laboratory 1992
c     author - w. smith         july 1992
c     
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c     
c***********************************************************************

      implicit none

      real(8) aaa,bbb,axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3

      dimension aaa(9),bbb(10)

c     calculate lengths of cell vectors

      bbb(1)=sqrt(aaa(1)*aaa(1)+aaa(2)*aaa(2)+aaa(3)*aaa(3))
      bbb(2)=sqrt(aaa(4)*aaa(4)+aaa(5)*aaa(5)+aaa(6)*aaa(6))
      bbb(3)=sqrt(aaa(7)*aaa(7)+aaa(8)*aaa(8)+aaa(9)*aaa(9))

c     calculate cosines of cell angles

      bbb(4)=(aaa(1)*aaa(4)+aaa(2)*aaa(5)+aaa(3)*aaa(6))/(bbb(1)*bbb(2))
      bbb(5)=(aaa(1)*aaa(7)+aaa(2)*aaa(8)+aaa(3)*aaa(9))/(bbb(1)*bbb(3))
      bbb(6)=(aaa(4)*aaa(7)+aaa(5)*aaa(8)+aaa(6)*aaa(9))/(bbb(2)*bbb(3))

c     calculate vector products of cell vectors

      axb1=aaa(2)*aaa(6)-aaa(3)*aaa(5)
      axb2=aaa(3)*aaa(4)-aaa(1)*aaa(6)
      axb3=aaa(1)*aaa(5)-aaa(2)*aaa(4)
      bxc1=aaa(5)*aaa(9)-aaa(6)*aaa(8)
      bxc2=aaa(6)*aaa(7)-aaa(4)*aaa(9)
      bxc3=aaa(4)*aaa(8)-aaa(5)*aaa(7)
      cxa1=aaa(8)*aaa(3)-aaa(2)*aaa(9)
      cxa2=aaa(1)*aaa(9)-aaa(3)*aaa(7)
      cxa3=aaa(2)*aaa(7)-aaa(1)*aaa(8)

c     calculate volume of cell

      bbb(10)=abs(aaa(1)*bxc1+aaa(2)*bxc2+aaa(3)*bxc3)

c     calculate cell perpendicular widths

      bbb(7)=bbb(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
      bbb(8)=bbb(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
      bbb(9)=bbb(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)

      return
      end subroutine dcell

      subroutine invert(a,b,d)

c***********************************************************************
c     
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       april 1992
c     
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c     
c***********************************************************************

      implicit none

      real(8) a,b,d,r

      dimension a(9),b(9)

c     calculate adjoint matrix
      b(1)=a(5)*a(9)-a(6)*a(8)
      b(2)=a(3)*a(8)-a(2)*a(9)
      b(3)=a(2)*a(6)-a(3)*a(5)
      b(4)=a(6)*a(7)-a(4)*a(9)
      b(5)=a(1)*a(9)-a(3)*a(7)
      b(6)=a(3)*a(4)-a(1)*a(6)
      b(7)=a(4)*a(8)-a(5)*a(7)
      b(8)=a(2)*a(7)-a(1)*a(8)
      b(9)=a(1)*a(5)-a(2)*a(4)

c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d

c     complete inverse matrix
      b(1)=r*b(1)
      b(2)=r*b(2)
      b(3)=r*b(3)
      b(4)=r*b(4)
      b(5)=r*b(5)
      b(6)=r*b(6)
      b(7)=r*b(7)
      b(8)=r*b(8)
      b(9)=r*b(9)

      return
      end subroutine invert

      subroutine jacobi(a,v,n)

c***********************************************************************
c     
c     diagonalisation of real symmetric matices by jacobi method
c     
c     input parameters:
c     
c     a(n,n) is the matrix to be diagonalised
c     v(n,n) is the eigenvector matrix
c     n   is the dimension of the matrices
c     
c     jacobi processes lower triangle only (upper triangle unchanged)
c     
c     variable rho sets absolute tolerance on convergence
c     variable tes is a moving tolerance that diminishes
c     on each pass until at true convergence tes<rho
c     
c     author w.smith 1993
c     
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c     
c***********************************************************************

      implicit none

      logical pass
      integer n,i,j,k
      real(8) a,v,rho,tes,scl,v1,v2,v3,u,omg,s,c,tem

      dimension a(n,n),v(n,n)

      rho=1.0d-16
      tes=0.0d0
      scl=0.0d0

c     initialize eigenvectors

      do i=1,n
        do j=1,n
          v(i,j)=0.0d0
        enddo
        v(i,i)=1.0d0
      enddo

c     rescale matrix for optimal accuracy

      do i=1,n
        if(abs(a(i,i)).gt.scl)scl=abs(a(i,i))
      enddo
      do i=1,n
        do j=1,i
          a(i,j)=a(i,j)/scl
        enddo
      enddo

c     set initial value of moving tolerance

      do i=2,n
        do j=1,i-1
          tes=tes+2.0d0*a(i,j)*a(i,j)
        enddo
      enddo
      tes=sqrt(tes)

c     recycle until absolute tolerance satisfied

      do while(tes.gt.rho)

        tes=tes/dble(n)
        if(tes.lt.rho)tes=rho
        
c     jacobi diagonalisation
        
        pass=.true.
        
c     recycle until moving tolerance satisfied
        
        do while(pass)
          
          pass=.false.
          
          do i=2,n
            
            do j=1,i-1
              
              if(abs(a(i,j)).ge.tes)then
                pass=.true.
                v1=a(j,j)
                v2=a(i,j)
                v3=a(i,i)
                u=0.5d0*(v1-v3)
                if(abs(u).lt.rho)then
                  omg=-1.0d0
                else
                  omg=-v2/sqrt(v2*v2+u*u)
                  if(u.lt.0.0d0)omg=-omg
                endif
                s=omg/sqrt(2.0d0*(1.0d0+sqrt(1.0d0-omg*omg)))
                c=sqrt(1.0d0-s*s)
                do k=1,n
                  if(k.ge.i)then
                    tem=a(k,j)*c-a(k,i)*s
                    a(k,i)=a(k,j)*s+a(k,i)*c
                    a(k,j)=tem
                  else if(k.lt.j)then
                    tem=a(j,k)*c-a(i,k)*s
                    a(i,k)=a(j,k)*s+a(i,k)*c
                    a(j,k)=tem
                  else
                    tem=a(k,j)*c-a(i,k)*s
                    a(i,k)=a(k,j)*s+a(i,k)*c
                    a(k,j)=tem
                  endif
                  tem=v(k,j)*c-v(k,i)*s
                  v(k,i)=v(k,j)*s+v(k,i)*c
                  v(k,j)=tem
                enddo
                a(j,j)=v1*c*c+v3*s*s-2.0d0*v2*s*c
                a(i,i)=v1*s*s+v3*c*c+2.0d0*v2*s*c
                a(i,j)=(v1-v3)*s*c+v2*(c*c-s*s)
              endif
              
            enddo
            
          enddo
          
        enddo

      enddo

c     rescale matrix

      do i=1,n
        do j=1,i
          a(i,j)=scl*a(i,j)
        enddo
      enddo

      return
      end subroutine jacobi

      subroutine set_block(nnn,ccc,aaa)

c**********************************************************************
c
c     dl_poly subroutine to initialise an array to a single value
c
c     copyright daresbury laboratory 1998
c     author w.smith july 1998
c
c**********************************************************************

      implicit none

      integer i,nnn
      real(8) ccc,aaa(nnn)

      do i=1,nnn,2

        aaa(i)=ccc
        aaa(i+1)=ccc

      enddo
      
      return
      end subroutine set_block

      subroutine matmul(aaa,bbb,ccc)

c***********************************************************************
c     
c     dlpoly utility to multiply 3x3 matrices
c
c     copyright daresbury laboratory
c     author      w.smith  oct  2005
c     
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c     
c**********************************************************************

      implicit none

      integer i
      real(8) aaa(9),bbb(9),ccc(9),tmp(9)

      tmp(1)=aaa(1)*bbb(1)+aaa(4)*bbb(2)+aaa(7)*bbb(3)
      tmp(2)=aaa(2)*bbb(1)+aaa(5)*bbb(2)+aaa(8)*bbb(3)
      tmp(3)=aaa(3)*bbb(1)+aaa(6)*bbb(2)+aaa(9)*bbb(3)

      tmp(4)=aaa(1)*bbb(4)+aaa(4)*bbb(5)+aaa(7)*bbb(6)
      tmp(5)=aaa(2)*bbb(4)+aaa(5)*bbb(5)+aaa(8)*bbb(6)
      tmp(6)=aaa(3)*bbb(4)+aaa(6)*bbb(5)+aaa(9)*bbb(6)

      tmp(7)=aaa(1)*bbb(7)+aaa(4)*bbb(8)+aaa(7)*bbb(9)
      tmp(8)=aaa(2)*bbb(7)+aaa(5)*bbb(8)+aaa(8)*bbb(9)
      tmp(9)=aaa(3)*bbb(7)+aaa(6)*bbb(8)+aaa(9)*bbb(9)
      
      do i=1,9
        ccc(i)=tmp(i)
      enddo
      
      return
      end subroutine matmul 

      subroutine debugging(idnode,natms,levcfg,imcon,cfgname,
     &eng,outdir,ins,del,dis,pass,delE,ewld1eng,ewld2sum,ewld3sum,
     &vdwsum,delrc,elrc,engunit,iguest)
c*********************************************************************
c     subroutine writes files necessary to get DL_poly energies
c     from insertions,deletions and displacements
c*********************************************************************
      implicit none
      logical ins,del,dis
      character*1 cfgname(80)
      character*8 outdir
      real(8) engunit,ewld1eng,ewld2sum,ewld3sum,vdwsum
      real(8) delrc,elrc
      real(8) eng,delE(*)
      integer pass,iguest,levcfg,imcon,idnode,natms

      if(pass.eq.1)then
        call revive_debug
     &(idnode,natms,levcfg,imcon,cfgname,eng,outdir,ins,del,dis,pass)
      else

        call revive_debug
     &(idnode,natms,levcfg,imcon,cfgname,eng,outdir,ins,del,dis,pass)
        if(ins)then          
          open(39,file=outdir//'/debug_ins')
          write(39,'(3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a13,f20.10,/,
     &      3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a13,f20.10)') 
     &      "Total eng: ",delE(iguest),"ewald 1 sum:",ewld1eng,
     &      "ewald 2 sum: ",ewld2sum,"ewald 3 sum: ",ewld3sum,
     &      "sum of vdw: ",vdwsum,
     &      "lrng correct:",delrc/engunit
          close(39)
        elseif(del)then
          open(39,file=outdir//'/debug_del')
           
          write(39,'(3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a13,f20.10,/,
     &      3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a13,f20.10)') 
     &      "Total eng: ",delE(iguest),"ewald 1 sum:",ewld1eng,
     &      "ewald 2 sum: ",ewld2sum,"ewald 3 sum: ",ewld3sum,
     &      "sum of vdw: ",vdwsum,
     &      "lrng correct:",delrc/engunit
          close(39) 
        elseif(dis)then
          open(39,file=outdir//'/debug_dis')
           
          write(39,'(3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a13,f20.10,/,
     &      3x,a13,f20.10,/,3x,a13,f20.10,/,3x,a19,f20.10)') 
     &      "Total eng: ",delE(iguest),"ewald 1 sum:",ewld1eng,
     &      "ewald 2 sum: ",ewld2sum,"ewald 3 sum: ",ewld3sum,
     &      "sum of vdw: ",vdwsum,
     &      "lrng correct total:",elrc/engunit
          close(39) 
        endif
      endif
      end subroutine debugging

      subroutine revive_debug
     &(idnode,natms,levcfg,imcon,cfgname,eng,outdir,ins,del,dis,pass)
c**********************************************************************
c 
c     subroutine to write the atom information to a file
c
c**********************************************************************
      implicit none
      integer, parameter :: nconfig=23
      logical ins,del,dis
      character*1 cfgname(80)
      character*8 outdir
      character*25 outfile
      integer i,natms,idnode,levcfg,imcon,pass
      real(8) eng


      if(ins)
     &write(outfile,'(a8,a11,i1)')outdir,'/REVCON_ins',pass
      if(del)
     &write(outfile,'(a8,a11,i1)')outdir,'/REVCON_del',pass
      if(dis)
     &write(outfile,'(a8,a11,i1)')outdir,'/REVCON_dis',pass
    
      open(nconfig,file=outfile,form='formatted')
      write(nconfig,'(80a1)')cfgname
      write(nconfig,'(2i10,e20.10)')levcfg,imcon,eng
      if(imcon.gt.0)write(nconfig,'(3f20.12)')cell

      do i=1,natms
        write(nconfig,'(a8,i10)')atomname(i),i
        write(nconfig,'(3g20.10)')xxx(i),yyy(i),zzz(i)
      enddo
      close(nconfig)
      return
      end subroutine revive_debug
      subroutine revive
     &(idnode,natms,levcfg,production,ntpguest,ntpmls,
     &imcon,cfgname,eng,outdir)
c**********************************************************************
c 
c     subroutine to write the atom information to a file
c
c**********************************************************************
      implicit none
      logical production
      character*1 cfgname(80)
      character*8 outdir
      character*25 outfile
      integer i,j,natms,idnode,levcfg,imcon,ntpguest,ntpmls,nsites
      real(8) eng

      write(outfile,'(a8,a7)')outdir,'/REVCON'
      open(nconfig,file=outfile,form='formatted')
      write(nconfig,'(80a1)')cfgname
      write(nconfig,'(2i10,e20.10)')levcfg,imcon,eng
      if(imcon.gt.0)write(nconfig,'(3f20.12)')cell

      do i=1,natms
        write(nconfig,'(a8,i10)')atomname(i),i
        write(nconfig,'(3g20.10)')xxx(i),yyy(i),zzz(i)
      enddo
      close(nconfig)

      write(outfile,'(a8,a7)')outdir,'/REVIVE'
      open(nrev,file=outfile,form='unformatted')
      write(nrev)production
      write(nrev)(nummols(i),i=1,ntpmls)
      write(nrev)(chainstats(i),i=1,ntpguest*14+1)
      write(nrev)natms
      do i=1,ntpmls
        nsites=nummols(i)*numatoms(i)
        write(nrev)(framwkxxx(i,j),j=1,nsites)
        write(nrev)(framwkyyy(i,j),j=1,nsites)
        write(nrev)(framwkzzz(i,j),j=1,nsites)
      enddo
      
   
      write(nrev)(energy(i),i=1,ntpguest)
      close(nrev)
      return
      end subroutine revive

      subroutine revscan(idnode,prevnodes)
c*********************************************************************
c
c     determine how many nodes were done in a previous run 
c
c*********************************************************************
      implicit none
      logical safe
      integer idnode,prevnodes,idum

      safe=.true.
      if(idnode.eq.0)then
         call system("ls | grep -c branch > 0110111")
         open(666,file="0110111")
      endif
      call getrec(safe,idnode,666)
      if(idnode.eq.0)close(666)
      prevnodes=intstr(record,lenrec,idum)
      if(idnode.eq.0)call system("rm 0110111")

      end subroutine revscan
      subroutine revread(localdir,production,ntpmls,totatm,ntpguest)
c*********************************************************************
c
c     read revive file and revcon file to get all the related data
c     to restart a gcmc simulation
c
c*********************************************************************
      implicit none
      logical production
      character*8 localdir
      character*15 outfile
      integer i,j,nsites,totatm,ntpguest,ios,ntpmls


      write(outfile,'(a8,a7)')localdir,'/REVIVE'
      open(nrev,file=outfile,form='unformatted',status='old',iostat=ios)

      if(ios.eq.29)then
c       value 29 means status=old is invalid ie. there is no existing
c       file called REVIVE.  Thus the restart will just start with the 
c       conditions in the CONFIG and FIELD file.

      else
        read(nrev)production
        read(nrev)(nummols(i),i=1,ntpmls)
        read(nrev)(chainstats(i),i=1,ntpguest*14+1)
        read(nrev)totatm
        do i=1,ntpmls
          nsites=nummols(i)*numatoms(i)
          read(nrev)(framwkxxx(i,j),j=1,nsites)
          read(nrev)(framwkyyy(i,j),j=1,nsites)
          read(nrev)(framwkzzz(i,j),j=1,nsites)
        enddo
   
        read(nrev)(energy(i),i=1,ntpguest)
        close(nrev)
      endif

      end subroutine revread
      
      subroutine normal_bin(fposa,fposb,fposc,ngrida,ngridb,ngridc,
     &subgrid)
c***********************************************************************
c
c     Update the global grid(subgrid, ...) with the configuration
c     Assigns the entire configuration to a single grid cell
c     
c***********************************************************************
      implicit none

      integer ngrida,ngridb,ngridc,subgrid
      real(8) fposa,fposb,fposc
c     local variables
      integer bin,aidx,bidx,cidx

c     Calculate the grid indices directly from the fractional positions
      aidx=int(fposa*ngrida)+1
      bidx=int(fposb*ngridb)+1
      cidx=int(fposc*ngridc)+1

c     Ensure indices stay within the grid bounds (periodic wrapping)
      if(aidx>ngrida)aidx=aidx-ngrida
      if(bidx>ngridb)bidx=bidx-ngridb
      if(cidx>ngridc)cidx=cidx-ngridc
      if(aidx<1)aidx=aidx+ngrida
      if(bidx<1)bidx=bidx+ngridb
      if(cidx<1)cidx=cidx+ngridc

c     Convert 3D indices to the 1D bin index
      bin=(aidx-1)*ngridb*ngridc+(bidx-1)*ngridc+cidx

c     Update the probability grid
      grid_norm(subgrid,bin)=grid(subgrid,bin)+1.0

      end subroutine normal_bin

      subroutine equitable_bin(fposa,fposb,fposc,ngrida,ngridb,ngridc,
     &subgrid)
c***********************************************************************
c
c     update the global grid(subgrid,...) with the equitable binned 
c     atom position
c     
c**********************************************************************
      implicit none
c arguments
      integer ngrida,ngridb,ngridc,subgrid
      real(8) fposa,fposb,fposc
c internal variables
      integer bin,aidx,bidx,cidx
      integer alidx,aridx,blidx,bridx,clidx,cridx
      real(8) gpta,gptb,gptc,apart,bpart,cpart
      real(8) al,ar,bl,br,cl,cr

c NOTE: since fortran arrays start at 1 and probabilities are 1D
c a?idx and b?idx are 0 based index and c?idx is 1 based
      gpta = ngrida*(fposa)
      aidx = ceiling(gpta)
      apart = modulo(gpta, 1.d0)
      if(apart.gt.0.5)then
        if(aidx.ge.ngrida)then
          alidx = (ngrida-1)*ngridb*ngridc
          aridx = 0
        else
          alidx = (aidx-1)*ngridb*ngridc
          aridx = (aidx)*ngridb*ngridc
        endif
        al = 1.5 - apart
        ar = apart - 0.5
      else
        if(aidx.le.1)then
          alidx = (ngrida-1)*ngridb*ngridc
          aridx = 0
        else
          alidx = (aidx-2)*ngridb*ngridc
          aridx = (aidx-1)*ngridb*ngridc
        endif
        al = 0.5 - apart
        ar = 0.5 + apart
      endif
c equitable spread for b vector
      gptb = ngridb*(fposb)
      bidx = ceiling(gptb)
      bpart = modulo(gptb, 1.d0)
      if(bpart.gt.0.5)then
        if(bidx.ge.ngridb)then
          blidx = (ngridb-1)*ngridc
          bridx = 0
        else
          blidx = (bidx-1)*ngridc
          bridx = (bidx)*ngridc
        endif
        bl = 1.5 - bpart
        br = bpart - 0.5
      else
        if(bidx.le.1)then
          blidx = (ngridb-1)*ngridc
          bridx = 0
        else
          blidx = (bidx-2)*ngridc
          bridx = (bidx-1)*ngridc
        endif
        bl = 0.5 - bpart
        br = 0.5 + bpart
      endif
c equitable spread for c vector
      gptc = ngridc*(fposc)
      cidx = ceiling(gptc)
      cpart = modulo(gptc, 1.d0)
      if(cpart.gt.0.5)then
        if(cidx.ge.ngridc)then
          clidx = ngridc
          cridx = 1
        else
          clidx = cidx
          cridx = cidx + 1
        endif
        cl = 1.5 - cpart
        cr = cpart - 0.5
      else
        if(cidx.le.1)then
          clidx = ngridc
          cridx = 1
        else
          clidx = cidx-1
          cridx = cidx
        endif
        cl = 0.5 - cpart
        cr = 0.5 + cpart
      endif
      bin = alidx+blidx+clidx
      grid(subgrid,bin) = grid(subgrid,bin)+(al*bl*cl)
      bin = aridx+blidx+clidx
      grid(subgrid,bin) = grid(subgrid,bin)+(ar*bl*cl)
      bin = alidx+bridx+clidx
      grid(subgrid,bin) = grid(subgrid,bin)+(al*br*cl)
      bin = alidx+blidx+cridx
      grid(subgrid,bin) = grid(subgrid,bin)+(al*bl*cr)
      bin = aridx+bridx+clidx
      grid(subgrid,bin) = grid(subgrid,bin)+(ar*br*cl)
      bin = aridx+blidx+cridx
      grid(subgrid,bin) = grid(subgrid,bin)+(ar*bl*cr)
      bin = alidx+bridx+cridx
      grid(subgrid,bin) = grid(subgrid,bin)+(al*br*cr)
      bin = aridx+bridx+cridx
      grid(subgrid,bin) = grid(subgrid,bin)+(ar*br*cr)
      end subroutine equitable_bin

      integer function atmnumber(i) 
c*******************************************************************
c     generate an atomic number from an atomic mass 
c     EDIT (pb 09/01/13): this function reads the mass reported
c     on a standard periodic table and assigns an atomic number.
c     You will run into problems if you are using atomic masses
c     of isotopes in the FIELD file.
c*******************************************************************
      implicit none
      real(8) i
      if ((i.ge.0.0).and.(i.le.1.5))then
        atmnumber=1
      elseif((i.ge.3.9).and.(i.le.4.5))then
        atmnumber=2
      elseif((i.ge.6.5).and.(i.le.7.1))then
        atmnumber=3
      elseif((i.ge.8.9).and.(i.le.9.5))then
        atmnumber=4
      elseif((i.ge.10.5).and.(i.le.11.1))then
        atmnumber=5
      elseif((i.ge.11.9).and.(i.le.12.5))then
        atmnumber=6
      elseif((i.ge.13.9).and.(i.le.14.5))then
        atmnumber=7
      elseif((i.ge.15.5).and.(i.le.16.1))then
        atmnumber=8
      elseif((i.ge.18.5).and.(i.le.19.1))then
        atmnumber=9
      elseif((i.ge.19.9).and.(i.le.20.5))then
        atmnumber=10
      elseif((i.ge.22.5).and.(i.le.23.1))then
        atmnumber=11
      elseif((i.ge.23.9).and.(i.le.24.5))then
        atmnumber=12
      elseif((i.ge.26.5).and.(i.le.27.1))then
        atmnumber=13
      elseif((i.ge.27.9).and.(i.le.28.5))then
        atmnumber=14
      elseif((i.ge.30.5).and.(i.le.31.1))then
        atmnumber=15
      elseif((i.ge.31.9).and.(i.le.32.5))then
        atmnumber=16
      elseif((i.ge.34.9).and.(i.le.36.1))then
        atmnumber=17
c     Ar (18) has mass range that overlaps with Ca (20). Be careful 
c     with mass rounding here!
      elseif((i.ge.39.5).and.(i.le.39.9999))then
        atmnumber=18
      elseif((i.ge.38.9).and.(i.le.39.4))then
        atmnumber=19
      elseif((i.ge.40.0).and.(i.le.40.5))then
        atmnumber=20
      elseif((i.ge.44.5).and.(i.le.45.1))then
        atmnumber=21
      elseif((i.ge.47.5).and.(i.le.48.1))then
        atmnumber=22
      elseif((i.ge.50.5).and.(i.le.51.1))then
        atmnumber=23
      elseif((i.ge.51.5).and.(i.le.52.1))then
        atmnumber=24
      elseif((i.ge.54.5).and.(i.le.55.1))then
        atmnumber=25
      elseif((i.ge.55.5).and.(i.le.56.1))then
        atmnumber=26
c     Co (27) and Ni (28) have very close mass ranges
      elseif((i.ge.58.76).and.(i.le.59.1))then
        atmnumber=27
      elseif((i.ge.58.5).and.(i.le.59.75))then
        atmnumber=28
      elseif((i.ge.62.9).and.(i.le.64.1))then
        atmnumber=29
      elseif((i.ge.64.9).and.(i.le.66.1))then
        atmnumber=30
      elseif((i.ge.69.5).and.(i.le.70.1))then
        atmnumber=31
      elseif((i.ge.72.5).and.(i.le.73.1))then
        atmnumber=32
      elseif((i.ge.74.5).and.(i.le.75.1))then
        atmnumber=33
      elseif((i.ge.78.5).and.(i.le.79.1))then
        atmnumber=34
      elseif((i.ge.79.5).and.(i.le.80.1))then
        atmnumber=35
      elseif((i.ge.83.5).and.(i.le.84.1))then
        atmnumber=36
      elseif((i.ge.84.9).and.(i.le.86.1))then
        atmnumber=37
      elseif((i.ge.87.5).and.(i.le.88.1))then
        atmnumber=38
      elseif((i.ge.88.5).and.(i.le.89.1))then
        atmnumber=39
      elseif((i.ge.90.9).and.(i.le.91.5))then
        atmnumber=40
      elseif((i.ge.92.5).and.(i.le.93.1))then
        atmnumber=41
      elseif((i.ge.95.5).and.(i.le.96.1))then
        atmnumber=42
      elseif((i.ge.97.9).and.(i.le.98.1))then
        atmnumber=43
      elseif((i.ge.109.9).and.(i.le.101.5))then
        atmnumber=44
      elseif((i.ge.102.5).and.(i.le.103.1))then
        atmnumber=45
      elseif((i.ge.105.9).and.(i.le.106.5))then
        atmnumber=46
      elseif((i.ge.107.5).and.(i.le.108.1))then
        atmnumber=47
      elseif((i.ge.111.9).and.(i.le.112.5))then
        atmnumber=48
      elseif((i.ge.114.5).and.(i.le.115.1))then
        atmnumber=49
      elseif((i.ge.118.5).and.(i.le.119.1))then
        atmnumber=50
      elseif((i.ge.121.5).and.(i.le.122.1))then
        atmnumber=51
      elseif((i.ge.127.5).and.(i.le.128.1))then
        atmnumber=52
      elseif((i.ge.126.5).and.(i.le.127.1))then
        atmnumber=53
      elseif((i.ge.130.9).and.(i.le.131.5))then
        atmnumber=54
      elseif((i.ge.132.5).and.(i.le.133.1))then
        atmnumber=55
      elseif((i.ge.136.9).and.(i.le.137.5))then
        atmnumber=56
      elseif((i.ge.138.5).and.(i.le.139.1))then
        atmnumber=57
      elseif((i.ge.139.9).and.(i.le.140.5))then
        atmnumber=58
      elseif((i.ge.140.6).and.(i.le.141.1))then
        atmnumber=59
      elseif((i.ge.144.0).and.(i.le.144.5))then
        atmnumber=60
      elseif((i.ge.144.9).and.(i.le.145.1))then
        atmnumber=61
      elseif((i.ge.150.0).and.(i.le.150.6))then
        atmnumber=62
      elseif((i.ge.151.5).and.(i.le.152.1))then
        atmnumber=63
      elseif((i.ge.156.9).and.(i.le.157.5))then
        atmnumber=64
      elseif((i.ge.158.5).and.(i.le.159.1))then
        atmnumber=65
      elseif((i.ge.162.0).and.(i.le.163.1))then
        atmnumber=66
      elseif((i.ge.164.5).and.(i.le.165.1))then
        atmnumber=67
      elseif((i.ge.166.5).and.(i.le.167.9))then
        atmnumber=68
      elseif((i.ge.168.0).and.(i.le.169.1))then
        atmnumber=69
      elseif((i.ge.172.9).and.(i.le.173.5))then
        atmnumber=70
      elseif((i.ge.174.0).and.(i.le.175.1))then
        atmnumber=71
      elseif((i.ge.178.0).and.(i.le.179.1))then
        atmnumber=72
      elseif((i.ge.180.0).and.(i.le.181.1))then
        atmnumber=73
      elseif((i.ge.183.0).and.(i.le.184.1))then
        atmnumber=74
      elseif((i.ge.185.9).and.(i.le.186.5))then
        atmnumber=75
      elseif((i.ge.189.9).and.(i.le.190.5))then
        atmnumber=76
      elseif((i.ge.191.9).and.(i.le.192.5))then
        atmnumber=77
      elseif((i.ge.194.9).and.(i.le.195.5))then
        atmnumber=78
      elseif((i.ge.196.5).and.(i.le.197.1))then
        atmnumber=79
      elseif((i.ge.200.0).and.(i.le.201.1))then
        atmnumber=80
      elseif((i.ge.203.9).and.(i.le.204.6))then
        atmnumber=81
      elseif((i.ge.206.9).and.(i.le.207.6))then
        atmnumber=82
      elseif((i.ge.208.5).and.(i.le.209.1))then
        atmnumber=83
c     Po atomic number 84 has the same mass range as Bi (83)
      elseif((i.ge.209.9).and.(i.le.210.1))then
        atmnumber=85
      elseif((i.ge.221.9).and.(i.le.222.1))then
        atmnumber=86
      elseif((i.ge.222.9).and.(i.le.223.1))then
        atmnumber=87
      elseif((i.ge.225.9).and.(i.le.226.1))then
        atmnumber=88
      elseif((i.ge.226.9).and.(i.le.227.1))then
        atmnumber=89
      elseif((i.ge.231.9).and.(i.le.232.1))then
        atmnumber=90
      elseif((i.ge.230.9).and.(i.le.231.1))then
        atmnumber=91
      elseif((i.ge.237.9).and.(i.le.238.1))then
        atmnumber=92
c     Np atomic number 93 has the same mass range as U (92)
      elseif((i.ge.243.9).and.(i.le.244.1))then
        atmnumber=94
      elseif((i.ge.242.9).and.(i.le.243.1))then
        atmnumber=95
      elseif((i.ge.246.9).and.(i.le.247.1))then
        atmnumber=96
c     Bk atomic number 97 has the same mass range as Cm (96)
      elseif((i.ge.250.9).and.(i.le.251.1))then
        atmnumber=98
      elseif((i.ge.251.9).and.(i.le.252.1))then
        atmnumber=99
      elseif((i.ge.256.9).and.(i.le.257.1))then
        atmnumber=100
      elseif((i.ge.257.9).and.(i.le.258.1))then
        atmnumber=101
      elseif((i.ge.258.9).and.(i.le.259.1))then
        atmnumber=102
      elseif((i.ge.261.9).and.(i.le.262.1))then
        atmnumber=103
      elseif((i.ge.266.9).and.(i.le.267.1))then
        atmnumber=104
      elseif((i.ge.267.9).and.(i.le.268.1))then
        atmnumber=105
      elseif((i.ge.270.9).and.(i.le.271.1))then
        atmnumber=106
      elseif((i.ge.269.9).and.(i.le.270.1))then
        atmnumber=107
      elseif((i.ge.268.9).and.(i.le.269.1))then
        atmnumber=108
      elseif((i.ge.277.9).and.(i.le.278.1))then
        atmnumber=109
      elseif((i.ge.280.9).and.(i.le.281.1))then
        atmnumber=110
c     Rg atomic number 112 has the same mass range as Ds (110)
      elseif((i.ge.284.9).and.(i.le.285.1))then
        atmnumber=112
      elseif((i.ge.285.9).and.(i.le.286.1))then
        atmnumber=113
      elseif((i.ge.288.9).and.(i.le.289.1))then
        atmnumber=114
c     Uup atomic number 115 has the same mass range as Fl (114)
      elseif((i.ge.292.9).and.(i.le.293.1))then
        atmnumber=116
      elseif((i.ge.293.9).and.(i.le.294.1))then
        atmnumber=117
c     Uuo atomic number 118 has the same mass range as Uus (117)
      endif
      return
      end function atmnumber

      function errorq(w,delw,x,delx,y,dely,z,delz)
c***********************************************************************
c
c     calculates the standard error based on the theory that
c     Q=f(w,x,y,z) has an error of
c
c    dQ^2 = (df/dw)^2 * delw^2 + (df/dx)^2 * delx^2 + (df/dy)^2 * dely^2
c            + (df/dz)^2 * delz^2
c
c
c***********************************************************************
      implicit none
      real(8) w,delw,x,delx,y,dely,z,delz
      real(8) numerator, denominator,denomsq
      real(8) dfdw,dfdx,dfdy,dfdz,errorq

      numerator=w-x*y
      denominator=z-y*y

c      denomsq=denominator*denominator

c      dfdw=1.d0/denominator
      
c      dfdx=-y/denominator

c      dfdy=2.d0*y*numerator/(denomsq)-
c     & x/denominator

c      dfdz=-1.d0*numerator/denomsq

c      errorq=sqrt((dfdw*delw/w)**2+(dfdx*delx/x)**2+
c     & (dfdy*dely/y)**2+(dfdz*delz/z)**2)

      errorq=abs(numerator/denominator)*
     &sqrt((delw/w)**2+(delx/x)**2+(dely/y)**2+(delz/z)**2)
      

      return

      end function errorq

      function calc_Qst(E2, E, N, N2, EN, temp)
c***********************************************************************
c
c     computes the isosteric heat of adsorption 
c
c***********************************************************************
      implicit none
      real(8) calc_Qst, E2, E, N, N2, EN, temp

      calc_Qst = -1.d0*(EN - E*N)/(N2-N*N) + Rgas*temp
      end function calc_Qst

      function old_cv(E2, E, N, N2, EN, temp)
c***********************************************************************
c
c     computes the heat capacity from gcmc 75 paper 
c
c***********************************************************************
      implicit none
      real(8) old_cv, E2, E, N, N2, EN, temp

      old_cv =(E2 - (E*E) - (EN - E*N)**2
     &   / (N2 - N*N)) / (N * kboltz * temp*temp) 
      end function old_cv

      function calc_Cv(E2, E, N, N2, EN, temp)
c***********************************************************************
c
c     computes the heat capacity at constant volume from Tildesley
c
c***********************************************************************
      implicit none
      real(8) calc_Cv, E2, E, N, N2, EN, temp

      calc_Cv = (3.d0/2.d0 * Rgas) + (E2 - (E*E) - (EN - E*N)**2
     &   / (N2 - N*N)) / (N * kboltz * temp*temp) 
      end function calc_Cv

      function errorcv(w, delw, x, delx, y, dely, z, delz, q, delq, 
     & temp)
c***********************************************************************
c
c     calculates the standard error based on the theory that
c     Q=f(w,x) has an error of
c
c    dQ^2 = (df/dw)^2 * delw^2 + (df/dx)^2 * delx^2 
c
c    w = avgE2, x = avgE, y = avgN, z = avgN2, q = avgEN
c***********************************************************************
      implicit none
      real(8) w, delw, x, delx, temp, numerator, denominator
      real(8) y, dely, z, delz, q, delq
      real(8) errorcv
      numerator=w-x*x - (q-x*y)**2/(z - y*y)
      denominator=kboltz*temp*temp*y
      errorcv=abs(numerator/denominator) *
     &sqrt((delw/w)**2+(delx/x)**2 + (dely/y)**2 + (delz/z)**2 +
     & (delq/q)**2) 
      return
      end function errorcv

      function sdot0(n,aaa,bbb)

c***********************************************************************
c     
c     dlpoly utility to calculate scalar product of two arrays
c
c     copyright daresbury laboratory
c     author      w.smith  july 2005
c     
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c     
c**********************************************************************

      implicit none

      integer n,i
      real(8) sdot0,aaa,bbb

      dimension aaa(*),bbb(*)

      sdot0=0.d0

      do i=1,n
        sdot0=sdot0+aaa(i)*bbb(i)
      enddo

      return
      end function sdot0

      integer function loc8(i,j)

c*********************************************************************
c
c     calculates array reference for a triangular matrix
c
c     copyright daresbury laboratory
c     author w.smith november 2005
c
c*********************************************************************
      
      integer i,j
      
      loc8=(max(i,j)*(max(i,j)-1))/2+min(i,j)
      
      return
      end function loc8

      character*3 function intstr3(nnn)

c*********************************************************************
c
c     converts a 3 digit integer to a string "001" etc.
c
c     copyright daresbury laboratory
c     author w.smith november 2005
c
c*********************************************************************

      implicit none

      integer nnn

      write(intstr3,'(i3.3)')nnn

      return
      end function intstr3

c      function duni()

c*********************************************************************
c     
c     dl_poly random number generator removed, new one in place 
c     to allow for time related seeds
c     
c*********************************************************************

c      implicit none

c      real(8) duni

c      call random_number(duni)

c      end function duni

      subroutine abort_field_read(kode,idnode,nfield)

c***********************************************************************
c     
c     dl_poly subroutine for aborting FIELD file read
c     
c     copyright - daresbury laboratory 
c     author    - w. smith    aug 2003
c     
c     wl
c     2006/12/07 13:05:20
c     1.11
c     Exp
c     
c***********************************************************************
      
      implicit none

      integer kode,idnode,nfield

      if(idnode.eq.0) close (nfield)

      if(kode.eq.1)then

c     end of field file error exit

        call error(idnode,52)

      else if(kode.eq.2)then

c     unrecognised directive in field file

        call error(idnode,4)

      endif

      return
      end subroutine abort_field_read

      subroutine abort_control_read(kode,idnode,nread)

c***********************************************************************
c     
c     dl_poly subroutine for aborting CONTROL file read
c     
c     copyright - daresbury laboratory 
c     author    - w. smith    aug 2003
c     
c     wl
c     2006/12/07 13:05:20
c     1.11
c     Exp
c     
c***********************************************************************
      
      implicit none

      integer kode,idnode,nread

      if(idnode.eq.0) close (nread)

      if(kode.eq.1)then

c     end of control file error exit

        call error(idnode,53)

      else if(kode.eq.2)then

c     general error exit from field file processing

        call error(idnode,0)

      endif

      return
      end subroutine abort_control_read

      subroutine abort_config_read(kode,idnode,nconf)

c***********************************************************************
c     
c     dl_poly subroutine for aborting CONTROL file read
c     
c     copyright - daresbury laboratory 
c     author    - w. smith    aug 2003
c     
c     wl
c     2006/12/07 13:05:20
c     1.11
c     Exp
c     
c***********************************************************************
      
      implicit none

      integer kode,idnode,nconf

      if(idnode.eq.0) close (nconf)

      if(kode.eq.1)then

c     general error exit from field file processing

        call error(idnode,54)

      else if(kode.eq.2)then

c     end of config file error exit

        call error(idnode,55)

      endif

      return
      end subroutine abort_config_read
      
      subroutine random_ins(imcon,natms,totatm,iguest,rcut,delr,
     &randa,randb,randc,rand1,rand2,rand3,rand4)
c*****************************************************************************
c
c     routine to randomly insert a guest into the cell 
c     ignore insertions placed too close to framework atoms?
c     the radial scan might be a bit costly
c
c*****************************************************************************
      implicit none
      logical good,insert,done
      integer i,ii,iguest,jatm,imcon,maxguest
      integer mol,natms,j,jj,k,totatm,counter,idnode
      real(8) randa,randb,randc,rmin,rcut,rclim,delr
      real(8) xc,yc,zc,comx,comy,comz,theta
      real(8) u1,u2,u3,q1,q2,q3,q4,u1sqrt,u1m1sqrt
      real(8) sintheta,beta1,beta2,beta3,norm,rand1,rand2,rand3,rand4

      good=.false.
      insert=.true.
      done=.false.
      xc=0.d0
      yc=0.d0
      zc=0.d0
      rmin=0.1d0
      rclim=(rcut+delr)**2

c     convert rand fractions to cartesian coordinates in the cell 
      call cartesian(randa,randb,randc,xc,yc,zc)
c     xc,yc,zc are the coordinates of the com guestx,guesty,guestz 
c     are the positions of atom i relative to the com
c      if(totatm+natms.ge.maxguest) call error(idnode,75)
      do i=1, natms
        newx(i)=guestx(iguest,i)
        newy(i)=guesty(iguest,i)
        newz(i)=guestz(iguest,i)
      enddo
c     the guestx,guesty, and guestz positions are centered at the
c     origin, this is done when the values are read from the 
c     FIELD file
      comx=0.d0
      comy=0.d0
      comz=0.d0 
c     rotate

c    setup quaternion which will uniformly sample the sphere
      
c      u1=duni()
c      u2=duni()
c      u3=duni()

c      u1sqrt=sqrt(u1)
c      u1m1sqrt=sqrt(1.d0-u1)
c     not sure if this representation is correct.
c     ie. the Q matrix is assuming the scalar term comes first
c     q(w,x,y,z) -- taken from wikipedia entry for "Rotation matrix"
c     the values for w,x,y,z were taken from 
c     "planning.cs.uiuc.edu/node198.html".. no mention where the 
c     scalar "w" lies in this representation... assuming it's the first
c     entry.
c      q1=u1m1sqrt*sin(2.d0*pi*u2)
c      q2=u1m1sqrt*cos(2.d0*pi*u2)
c      q3=u1sqrt*sin(2.d0*pi*u3)
c      q4=u1sqrt*cos(2.d0*pi*u3)

      theta=rand1*2.d0*pi
      beta1=2.d0*rand2-1.d0
      beta2=2.d0*rand3-1.d0
      beta3=2.d0*rand4-1.d0

c      beta1=duni()
c      beta2=duni()
c      beta3=duni()
      sintheta=sin(theta/2.d0)

      norm=sqrt(beta1*beta1+beta2*beta2+beta3*beta3)
      u1=beta1/norm
      u2=beta2/norm
      u3=beta3/norm

c     quaternion below
      q1=cos(theta/2.d0)
      q2=sintheta*u1
      q3=sintheta*u2
      q4=sintheta*u3
      call rotation(newx,newy,newz,natms,q1,q2,q3,q4) 
c      call rotationeuler(newx,newy,newz,natms,2.d0*pi) 
c     add the com to each atom in the guest

      do i=1, natms
        newx(i)=newx(i)+xc
        newy(i)=newy(i)+yc
        newz(i)=newz(i)+zc
      enddo


c     test to see if radial overlap with all other atoms. 
c      call images(imcon,natms,cell,newx,newy,newz)
c      really messed right now, commented out for later
c      investigation.  Uncomment the do while loop at the
c      beginning once this radial evaluation is working
c      call radial_eval(imcon,totatm,insert,natms,newx,newy,newz,good)
      
c      enddo
      return
      end subroutine random_ins

      subroutine random_del(nmols)
c*****************************************************************************
c
c     routine to randomly delete a guest from the cell 
c
c*****************************************************************************
      implicit none
      logical flag
      integer nmols,randguest,nloc,i
      real(8) rand

      data randguest/0/,flag/.false./
c     randdel should work like this: 
c     first do an ewald calc for the atoms being deleted
c     add to ewald1, which should've been calculated previously
c     do vdw calc etc..
c     move the xyz's to the end of
c     framwkxxx,framwkyyy,framwkzzz and shift everything back by
c     the number of atoms
c     change nummols to nummols-1
c     the deltaE will be -ve ewald2,vdw and the difference
c     between ewald1 before and ewald1 after.

      return
      end subroutine random_del


      subroutine random_disp
     &(imcon,natms,mol,newx,newy,newz,delr,cell,rcell,rotangle,
     &randa,randb,randc,rand1,rand2,rand3,rand4)
c*****************************************************************************
c
c     routine to randomly displace newx,newy,newz atoms
c
c*****************************************************************************
      implicit none
      logical done
      integer imcon,natms,mol,it,i
      real(8), dimension(natms) :: newx,newy,newz
c     these are the unit cell vectors, random moves
c     will be done in the dimensions of the cell.
c     only need to be calculated once at the beginning
c     of the simulation.
      real(8), dimension(9) :: cell,rcell 
      real(8) comx,comy,comz,delr,tol,movx,movy,movz,det
      real(8) randa,randb,randc,halfdelr,rotangle,x,y,z
      real(8) ssx,ssy,ssz,xss,yss,zss,shiftx,shifty,shiftz
      real(8) rand1,rand2,rand3,rand4
      done=.false.

      halfdelr=delr*0.5d0
c      tol=-3.d0/2.d0*delr+1.d-5
c     generate random numbers in the a,b, and c directions
c      do while(.not.done)
      randa=(2.d0*randa-1.d0)*halfdelr
      randb=(2.d0*randb-1.d0)*halfdelr
      randc=(2.d0*randc-1.d0)*halfdelr
c     some kind of tolerance here?
c        if(randa+randb+randc.gt.tol)done=.true.

c      enddo
c      movx=randa*iabc(1)+randb*iabc(4)+randc*iabc(7)
c      movy=randa*iabc(2)+randb*iabc(5)+randc*iabc(8)
c      movz=randa*iabc(3)+randb*iabc(6)+randc*iabc(9)

c     update the new coordinates with the shift vector
c     used to be shifts along the cell vectors..
      do i=1,natms
        newx(i)=randa+newx(i)
        newy(i)=randb+newy(i)
        newz(i)=randc+newz(i)
c        newx(i)=movx+newx(i)
c        newy(i)=movy+newy(i)
c        newz(i)=movz+newz(i)
      enddo
c     check pbc for com (keep the molecule intact)
      call com(natms,mol,newx,newy,newz,comx,comy,comz)

c     all this to make sure the COM lies within the boundary of 
c     the simulation... sheeit.
      x=comx
      y=comy
      z=comz
 
      ssx=(rcell(1)*x+rcell(4)*y+rcell(7)*z)
      ssy=(rcell(2)*x+rcell(5)*y+rcell(8)*z)
      ssz=(rcell(3)*x+rcell(6)*y+rcell(9)*z)
          
      xss=ssx-floor(ssx)
      yss=ssy-floor(ssy)
      zss=ssz-floor(ssz)
       
      x=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
      y=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
      z=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
 
      shiftx=x-comx
      shifty=y-comy
      shiftz=z-comz
 
      comx=x
      comy=y
      comz=z
      newx=newx+shiftx
      newy=newy+shifty
      newz=newz+shiftz

      call random_rot(rotangle,newx,newy,newz,comx,comy,comz,natms,mol,
     &rand1,rand2,rand3,rand4)

      return
      end subroutine random_disp


      subroutine random_jump
     &(natms,mol,newx,newy,newz,randa,randb,randc,
     &rand1,rand2,rand3,rand4)
c*****************************************************************************
c
c     place molecule at a completely random position in the cell
c     should be similar to a simultaneous insertion and deletion
c
c*****************************************************************************
      implicit none
      logical done
      integer natms, mol, i
      real(8), dimension(natms) :: newx,newy,newz
      real(8) comx,comy,comz
      real(8) randa,randb,randc,jumpangle
      real(8) rand1, rand2, rand3, rand4
      real(8) xdisp, ydisp, zdisp

c rotate randomly up to 2 pies 
c      jumpangle = 8.D0*datan(1.D0)*duni()
      jumpangle = 8.D0*datan(1.D0)
c random new position anywhere in fractional cell
c      xduni=duni()
c      yduni=duni()
c      zduni=duni()
      call cartesian(randa,randb,randc,xdisp,ydisp,zdisp)

      call com(natms,mol,newx,newy,newz,comx,comy,comz)

c make relative to com
      xdisp = xdisp - comx
      ydisp = ydisp - comy
      zdisp = zdisp - comz

c move molecule to new position
      do i=1,natms
        newx(i)=xdisp+newx(i)
        newy(i)=ydisp+newy(i)
        newz(i)=zdisp+newz(i)
      enddo

c rotate at the new position
      call com(natms,mol,newx,newy,newz,comx,comy,comz)
      call random_rot(jumpangle,newx,newy,newz,comx,comy,comz,natms,mol,
     &rand1,rand2,rand3,rand4)

      return
      end subroutine random_jump

      subroutine rotationeuler(xxr,yyr,zzr,natms,rotangle)
c*********************************************************************      
c     subroutine to rotate by the old euler angles.  singularities
c     are a problem with this, implemented to test robustness of 
c     quaternion rotation
c*********************************************************************
      implicit none
      integer i,natms
      real(8) rotangle,alpha,beta,kappa,randmult,xtmp,ytmp,ztmp
      real(8), dimension(9) :: R1,R2,R3,mult1,mult2
      real(8), dimension(natms) :: xxr,yyr,zzr
      real(8) rand1,rand2,rand3,rand4
c     alpha beta kappa generated from rotangle
      alpha=(rand1*2.d0-1.d0)*rotangle
      beta=(rand2*2.d0-1.d0)*rotangle
      kappa=(rand3*2.d0-1.d0)*rotangle

c     R1, rotation about the x axis by angle alpha

      R1(1)=1.d0
      R1(2)=0.d0
      R1(3)=0.d0
      R1(4)=0.d0
      R1(5)=cos(alpha)
      R1(6)=-1*sin(alpha)
      R1(7)=0.d0
      R1(8)=sin(alpha)
      R1(9)=cos(alpha)

c     R2, rotation about the y axis by angle beta

      R2(1)=cos(beta)
      R2(2)=0.d0
      R2(3)=sin(beta)
      R2(4)=0.d0
      R2(5)=1.d0
      R2(6)=0.d0
      R2(7)=-1*sin(beta)
      R2(8)=0.d0
      R2(9)=cos(beta)

c     R3, rotation about the z axis by angle kappa

      R3(1)=cos(kappa)
      R3(2)=-1*sin(kappa)
      R3(3)=0.d0
      R3(4)=sin(kappa)
      R3(5)=cos(kappa)
      R3(6)=0.d0
      R3(7)=0.d0
      R3(8)=0.d0
      R3(9)=1.d0

c     do some matrix multiplication.  This may have a significant
c     impact on the speed of each gcmc step

      if(rand4.lt.0.5)then
        call matmul(R2,R3,mult2)
      else
        call matmul(R3,R2,mult2)
      endif

      do i=1,natms
        xtmp=xxr(i)*mult2(1)+yyr(i)*mult2(4)+zzr(i)*mult2(7)
        ytmp=xxr(i)*mult2(2)+yyr(i)*mult2(5)+zzr(i)*mult2(8)
        ztmp=xxr(i)*mult2(3)+yyr(i)*mult2(6)+zzr(i)*mult2(9)
        xxr(i)=xtmp
        yyr(i)=ytmp
        zzr(i)=ztmp
      enddo

      end subroutine rotationeuler
      subroutine rotation(xxr,yyr,zzr,natms,q1,q2,q3,q4)
c*****************************************************************************
c
c     quaternion rotation
c     the quaternion is chosen randomly 
c
c*****************************************************************************
      implicit none
      integer i,natms
      real(8) q1,q2,q3,q4
      real(8), dimension(natms) :: xxr,yyr,zzr,tmpx,tmpy,tmpz
      real(8), dimension(9) :: Q

      do i=1,9
        Q(i)=0.d0
      enddo
c     taken from wikipedia entry "Rotation matrix"
      Q(1)=2.d0*q2*q2+2.d0*q1*q1-1.d0
      Q(2)=2.d0*q2*q3-2.d0*q1*q4
      Q(3)=2.d0*q2*q4+2.d0*q1*q3
      Q(4)=2.d0*q2*q3+2.d0*q1*q4
      Q(5)=2.d0*q1*q1+2.d0*q3*q3-1.d0
      Q(6)=2.d0*q3*q4-2.d0*q1*q2
      Q(7)=2.d0*q2*q4-2.d0*q1*q3
      Q(8)=2.d0*q3*q4+2.d0*q1*q2
      Q(9)=2.d0*q4*q4+2.d0*q1*q1-1.d0

      do i=1,natms
        tmpx(i)=Q(1)*xxr(i)+Q(2)*yyr(i)+Q(3)*zzr(i)
        tmpy(i)=Q(4)*xxr(i)+Q(5)*yyr(i)+Q(6)*zzr(i)
        tmpz(i)=Q(7)*xxr(i)+Q(8)*yyr(i)+Q(9)*zzr(i)
      enddo
      do i=1,natms
        xxr(i)=tmpx(i)
        yyr(i)=tmpy(i)
        zzr(i)=tmpz(i)
      enddo

      return
      end subroutine rotation


      subroutine random_rot
     &(rotangle,newx,newy,newz,comx,comy,comz,natms,mol,
     &rand1,rand2,rand3,rand4)
c*****************************************************************************
c
c     routine to randomly rotate natms at the index ik
c
c*****************************************************************************
      implicit none
      logical done
      integer i,natms,ik,mol
      real(8), dimension(4) :: quatern
      real(8) comx,comy,comz,rotangle,theta,beta1,beta2,beta3,norm
      real(8) sintheta,q1,q2,q3,q4,lenq,u1,u2,u3,u1sqrt,u1m1sqrt
      real(8), dimension(natms) :: newx,newy,newz,xtemp,ytemp,ztemp
      real(8) rand1,rand2,rand3,rand4
      done=.false.
 
c     com is subtracted to rotate about it
      do i=1,natms     
        newx(i)=newx(i)-comx
        newy(i)=newy(i)-comy
        newz(i)=newz(i)-comz
      enddo 
c    setup quaternion which will uniformly sample the sphere
        
c      u1=duni()
c      u2=duni()
c      u3=duni()

c      u1sqrt=sqrt(u1)
c      u1m1sqrt=sqrt(1.d0-u1)
c      q1=u1m1sqrt*sin(rotangle/2.d0*u2)
c      q2=u1m1sqrt*cos(rotangle/2.d0*u2)
c      q3=u1sqrt*sin(rotangle/2.d0*u3)
c      q4=u1sqrt*cos(rotangle/2.d0*u3)

c      call rotation(newx,newy,newz,natms,q1,q2,q3,q4) 
c      theta=(2.d0*duni()-1.d0)*rotangle*0.5d0
      theta=rand1*rotangle
      beta1=2.d0*rand2-1.d0
      beta2=2.d0*rand3-1.d0
      beta3=2.d0*rand4-1.d0

c      beta1=duni()
c      beta2=duni()
c      beta3=duni()
      sintheta=sin(theta/2.d0)

      norm=sqrt(beta1*beta1+beta2*beta2+beta3*beta3)
      u1=beta1/norm
      u2=beta2/norm
      u3=beta3/norm

c     quaternion below
      q1=cos(theta/2.d0)
      q2=sintheta*u1
      q3=sintheta*u2
      q4=sintheta*u3

      call rotation(newx,newy,newz,natms,q1,q2,q3,q4)
c      call rotationeuler(newx,newy,newz,natms,rotangle)
c     com re-added to the xyz terms once the rotation is done

      do i=1,natms
        newx(i)=newx(i)+comx
        newy(i)=newy(i)+comy
        newz(i)=newz(i)+comz 
      enddo

      return
      end subroutine random_rot


      subroutine com(natms,mol,newx,newy,newz,comx,comy,comz)
c*****************************************************************************
c
c     routine to calculate the centre of mass 
c     given an index and number of atoms to iterate over
c
c*****************************************************************************
      implicit none
      real(8) weight,comx,comy,comz,newx,newy,newz,mass
      integer natms,j,i,atomindex,mol
      dimension newx(*),newy(*),newz(*)


c     initialized values
      comx=0.d0
      comy=0.d0
      comz=0.d0 
      weight=0.d0
c     calculate the centre of mass of the molecule
      do i=1,natms
        mass=atmwght(mol,i)
        weight=weight+mass
        comx=comx+newx(i)*mass
        comy=comy+newy(i)*mass
        comz=comz+newz(i)*mass
      enddo
      
      comx=comx/weight
      comy=comy/weight
      comz=comz/weight
      return
      end subroutine com


      subroutine round(val,dum,base)
c*********************************************************************
c     rounds a quad precision number to the nearest base value
c*********************************************************************
      implicit none
      real(8) val,dum,mult
      integer base

      mult=10.d0**base
      dum=dble(anint(val*mult))/mult

      end subroutine round
      subroutine energy_eval
     &(eng,rande,volm,press,ngsts,temp,beta,
     &displace,insert,delete,accepted)
c*****************************************************************************
c
c     subroutine evaluates whether to accept or reject 
c     a gcmc move
c
c*****************************************************************************
      implicit none
      logical accepted,displace,insert,delete
      integer ngsts
      real(8) rande,test,edummy
      real(8) volm,press,temp,beta
      real(8) eng
      accepted=.false.
      test=0.d0

c     TODO(pboyd): include biased MC sampling here.
c     Define a probability of acceptance associated with the inserted COM
c     ENERGY BASED?

c     round the energy to the nearest second decimal place.
c     *** WARNING - the higher the value of eng, the worse this
c        rounding is ***
c      call round(duni(),rande,2)
      edummy=eng

c     first if statement prevents exp overflow
      if (-1*beta*edummy.gt.700.d0)then
        test=dexp(700.d0)
      elseif(displace)then
        test=dexp(-1.d0*beta*edummy)
      elseif(insert)then
        test=volm*press/((ngsts+1)*boltz*temp)*dexp(-1.d0*beta*edummy)
      elseif(delete)then
        if(abs(press).lt.1.d-7)then
            test = 1.d0
        else
            test=ngsts*boltz*temp/(volm*press)*dexp(-1.d0*beta*edummy)
        endif
      endif
      if(rande.lt.test)then
         accepted=.true.
      endif
      return
      end subroutine energy_eval
      subroutine cartesian(fraca,fracb,fracc,cartx,carty,cartz)
c*****************************************************************************
c
c     subroutine takes 3 fractional coordinates and returns 
c     3 cartesian coordinates based on the cell vectors
c
c*****************************************************************************
      implicit none
      real(8) fraca,fracb,fracc,cartx,carty,cartz

      cartx=fraca*cell(1)+fracb*cell(4)+fracc*cell(7)
      carty=fraca*cell(2)+fracb*cell(5)+fracc*cell(8)
      cartz=fraca*cell(3)+fracb*cell(6)+fracc*cell(9)

      return
      end subroutine cartesian
      
      subroutine guestlistgen
     &(imcon,iguest,totatm,rcut,delr,
     &natms,newx,newy,newz)
c*****************************************************************************
c    
c     builds temporary neighbour list.. area for 
c     improvement.. especially finding the index
c     of the xnew,ynew,and znew components.    
c 
c*****************************************************************************
      implicit none
      logical chk
      integer imcon,iguest,natms,i,ii,j,totatm
      integer itatms,atmadd,n,mol,at
      real(8) rsq,rmin,rclim,delr,rcut
      real(8),dimension(natms) :: newx,newy,newz
   
      rclim=(rcut+delr)**2
c     initialize the counter array
      do i=1,natms
        gstlentry(i)=0
      enddo 

      do i=1,natms
        ii=0
       
        do j=1,totatm
          chk=.true.
          do itatms=1,natms
            if (j.eq.ind(itatms))chk=.false.
          enddo
          
          if(chk)then
            ii=ii+1
            xdf(ii)=newx(i)-xxx(j)
            ydf(ii)=newy(i)-yyy(j)
            zdf(ii)=newz(i)-zzz(j)
          endif
        enddo
        call images(imcon,ii,cell,xdf,ydf,zdf)
        ii=0
        do j=1,totatm
          chk=.true.
          do itatms=1,natms
            if (j.eq.ind(itatms))chk=.false.
          enddo
         
          if(chk)then
            ii=ii+1
            if(imcon.eq.6)then
               rsq=xdf(ii)*xdf(ii)+ydf(ii)*ydf(ii)
            else 
               rsq=xdf(ii)*xdf(ii)+ydf(ii)*ydf(ii)+zdf(ii)*zdf(ii)
            endif
            if(rsq.le.rclim)then
              gstlentry(i)=gstlentry(i)+1
              gstlist(i,gstlentry(i))=j
            endif
          endif
        enddo
    
      enddo

      return
      end subroutine guestlistgen

      subroutine guest_exclude
     &(ntpguest)
c*****************************************************************************
c
c     subroutine populates the exclude lists for a guest
c     (excludes the other atoms of the guest.)
c     
c*****************************************************************************
      implicit none
      logical lchk
      integer nmols,natms,indx,last,mpm2,npm2,m,i,j,k,ii
      integer ntpguest,mol,jj,jz
 
c     make a general list of exclusion sites for all guest atom types
c     then we can build the neighbour lists based on these general lists

c     initialize arrays

      do i=1,ntpguest
        mol=locguest(i)
        natms=numatoms(mol)
        do ii=1,natms
          nexsit(i,ii)=0
          do jj=1,natms
            lexsit(i,ii,jj)=0
          enddo
        enddo
      enddo
      do i=1,ntpguest
        mol=locguest(i)
        natms=numatoms(mol)

        do m=1,natms-1
          ii=0

          do k=m+1,natms
            lchk=.true.
            do jz=1,min(nexsit(i,m),mxexcl)
              if(lexsit(i,m,jz).eq.k)lchk=.false.
            enddo
            if(lchk)then
              nexsit(i,m)=nexsit(i,m)+1
              nexsit(i,k)=nexsit(i,k)+1
              lexsit(i,m,nexsit(i,m))=k
              lexsit(i,k,nexsit(i,k))=m
            endif
            
          enddo
        enddo
      enddo
c     once the sites are built, then build the atom list of excluded
c     atoms, this can be found in subroutine exclude_atom
c     in exclude_terms.f

c     this is the old code which populates a list for immediate
c     calculation

c      last=natms
c      mpm2=natms/2
c      npm2=(natms-1)/2
c      ii=0
c 
c      do i=1,natms
c        ii=ii+1
c        nexatm(ii)=0
c      enddo
c 
c 
c      do m=1,mpm2
c        if(m.gt.npm2)last=mpm2
c 
c        ii=0
c 
c        do i=1,last
c          j=i+m
c          if(j.gt.natms)j=j-natms
c 
c          ii=ii+1
c          nexatm(ii)=nexatm(ii)+1
c          lexatm(ii,nexatm(ii))=j+indx
c 
c        enddo
c      enddo
      return
      end subroutine guest_exclude

      subroutine radial_eval
     &(imcon,totatm,insert,ngstatm,newx,newy,newz,good)
c*****************************************************************************
c
c     subroutine evaluates the radial distance between  
c     the cartesian coordinates in xnew,ynew,and znew
c     if they are above the minimum tolerance then returns
c     good=.true.
c     also builds temporary neighbour list.. area for 
c     improvement..
c
c*****************************************************************************
      implicit none
      logical good,done,insert
      integer imcon,natms,i,j,ngstatm,jatm,totatm
      real(8) rsq,rmin
      real(8),dimension(ngstatm) :: newx,newy,newz
      good=.true.
      done=.false.
      rmin=2.25d0
c     the do while loop is so that if an atom overlap is found
c     the loop termninates prematurely
c     so the gcmc move can be re-started

      do while(.not.done)
        do i=1,ngstatm
          if(.not.insert)then
            natms=gstlentry(i)
          else
            natms=totatm
          endif
          do j=1,natms
            if(.not.insert)then
              jatm=gstlist(i,j)
            else
              jatm=j
            endif
            xdf(j)=newx(i)-xxx(jatm)
            ydf(j)=newy(i)-yyy(jatm)
            zdf(j)=newz(i)-zzz(jatm)
          enddo
          call images(imcon,natms,cell,xdf,ydf,zdf)
          do j=1,natms
            if(imcon.eq.6)then
               rsq=xdf(j)*xdf(j)+ydf(j)*ydf(j)
            else 
               rsq=xdf(j)*xdf(j)+ydf(j)*ydf(j)+zdf(j)*zdf(j)
            endif
            if(rsq.lt.rmin)then
              good=.false.
              done=.true.
            endif
            if(i.eq.ngstatm.and.j.eq.natms)done=.true.
          enddo
        enddo
      enddo
      return

      end subroutine radial_eval

      subroutine condense(imcon,totatm,ntpmls,ntpfram,ntpguest)
c*****************************************************************************
c     
c     subroutine takes the current list of atoms and condenses it
c     to a 1-d array.  This makes energy calculations and neighbour list
c     generations easier and quicker to access. 
c
c*****************************************************************************
      implicit none
      integer i,n,itmls,ntpmls,ntpguest,itgst,itmol,imcon,ntpfram
      integer mol,indnam,nmol,molind,ii,jatm,isite,nsite,newatm
      integer j,k,kk,lsite,jj,indexsite,natms,totatm,latom
c     guests are placed first
      jj=0
      indexsite=0
      lsite=0 
      do itgst=1,ntpguest
        mol=locguest(itgst)
        natms=numatoms(mol)
        nmol=nummols(mol)
        ii=0
        
        do itmol=1,nmol 
          do isite=1,natms
            ii=ii+1
            jj=jj+1
            xxx(jj)=framwkxxx(mol,ii)
            yyy(jj)=framwkyyy(mol,ii)
            zzz(jj)=framwkzzz(mol,ii)
            atmcharge(jj)=atmchg(mol,isite)
            lfreezesite(jj)=lfzsite(mol,isite)
            atomname(jj)=atmname(mol,isite)
            atmweight(jj)=atmwght(mol,isite)
            ltype(jj)=ltpsit(mol,isite)
            kk=0
c      set up exclusion sites for the atoms
c      this is based on the assumption that 
c      guest atoms are RIGID and pairwise
c      interactions WILL NOT be calculated
c      between them.
            
            
            do k=1,nexsit(itgst,isite)
              newatm=lexsit(itgst,isite,k)+lsite

              if(((newatm.gt.jj).and.
     &          (newatm-jj.le.totatm/2)).or.
     &          ((newatm.lt.jj)
     &          .and.(newatm+totatm-jj.le.(totatm-1)/2)
     &          ))then
                kk=kk+1
                

                lexatm(jj,kk)=newatm

                if(kk.gt.1)then
                  do j=kk,2,-1
                    if(lexatm(jj,j).lt.lexatm(jj,j-1))
     &                then
                      latom=lexatm(jj,j)
                      lexatm(jj,j)=lexatm(jj,j-1)
                      lexatm(jj,j-1)=latom
                    endif
                  enddo
                endif
              endif
            enddo

            nexatm(jj)=kk

          enddo
          lsite=lsite+natms
        enddo
      enddo    
c     framework atoms are added after
      indnam=0
      do itmls=1,ntpfram
        ii=0
        mol=locfram(itmls)
        natms=numatoms(mol)
        nmol=nummols(mol)
        do i=1,nmol
          do k=1,natms
            jj=jj+1
            ii=ii+1
            xxx(jj)=framwkxxx(mol,ii)
            yyy(jj)=framwkyyy(mol,ii)
            zzz(jj)=framwkzzz(mol,ii)
            atmcharge(jj)=atmchg(mol,k)
            lfreezesite(jj)=lfzsite(mol,k)
            atomname(jj)=atmname(mol,k)
            atmweight(jj)=atmwght(mol,k)
            ltype(jj)=ltpsit(mol,k)
      
          enddo
        enddo
      enddo 
      return
      end subroutine condense

      subroutine angular_dist(mol,newx,newy,newz,natms,maxanglegrid)
c*****************************************************************************
c     
c     calculate the angular distribution of all guests at a particular 
c     step
c
c*****************************************************************************
      implicit none
      integer i,j,nmol,mol
      integer natms,jatm,iangle,maxanglegrid
      real(8) angle,comx,comy,comz
      real(8), dimension(3) :: zaxis,vector1
      real(8), dimension(natms) :: newx,newy,newz

      zaxis=(/0.d0,0.d0,1.d0/)
      angle=0.d0 
  
      comx=0.d0
      comy=0.d0
      comz=0.d0

c      call com(natms,mol,newx,newy,newz,comx,comy,comz)

c     the first vector with which the angle is calculated 
c     is between the second atom of the molecule and the com

      vector1=(/newx(2)-newx(1),newy(2)-newy(1),newz(2)-newz(1)/)

      call calc_angle(angle,vector1,zaxis)

      iangle=nint(angle/(pi/dble(maxanglegrid)))
   
      if(iangle.eq.0)iangle=1
      angdist(iangle)=angdist(iangle)+1

      end subroutine angular_dist
      
      subroutine jobcheck(idnode,jobsafe,ljob,production)
c*****************************************************************************
c     
c     checks a file called jobcontrol.in for job information 
c
c*****************************************************************************
      implicit none
      logical jobsafe, ljob, production, check, safe
      integer idum,idnode,i

      check=.true.
      safe=.true.

      if(idnode.eq.0)open(205,file='jobcontrol.in',status='old')

      call getrec(safe,idnode,205)
c     below is just a non mpi version of 'getrec' from parse_module.f
c     ...why??
c      read(205,'(a150)',end=100)line

c      do i=1,lenrec
c        lrec(i)=line(i:i)
c      enddo 

      call lowcase(record,lenrec)
      if(findstring('terminate',record,idum))then
        jobsafe=.false.
      elseif(findstring('start averaging',record,idum))then
        if((ljob).and.(.not.production))then 
           production=.true.
           if(idnode.eq.0)write(nrite,'(/,3x,a,/)')
     &"'START AVERAGING' found in jobcontrol.in, starting production
     & averaging"
        endif
      endif 
c     this is just nothing
c100   check=.false.
      if(idnode.eq.0)close(205)
      return
      end subroutine jobcheck


      subroutine calc_angle(angle,vect1,vect2) 
c*****************************************************************************
c     
c     calculate the angle between two vectors of dimension 3
c
c*****************************************************************************
      implicit none
      real(8) numer,denom,angle,dotprod
      real(8), dimension(3) :: vect1,vect2,crossprod

c     this is the numerator of the angle calc (dot product of vectors)
      angle=acos(dot_product(vect1,vect2)/
     &sqrt(dot_product(vect1,vect1))/sqrt(dot_product(vect2,vect2)))

      end subroutine calc_angle

      end module utility_pack

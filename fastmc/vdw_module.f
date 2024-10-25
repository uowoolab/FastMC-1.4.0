      module vdw_module

c***********************************************************************
c     
c     dl_poly module for defining van der waald potential arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c     wl
c     2006/11/28 16:32:48
c     1.4
c     Exp
c     
c***********************************************************************
      use utility_pack

      implicit none

      integer, allocatable :: ltpvdw(:),lstvdw(:)
      real(8), allocatable :: vvv(:,:),prmvdw(:,:)
      real(8), allocatable :: steadd(:)
c     maximum number of vdw parameters

      integer, parameter :: mxpvdw=5

      save ltpvdw,lstvdw,prmvdw,vvv

      contains
      
      subroutine alloc_vdw_arrays(idnode,maxvdw)
      implicit none
      integer, parameter :: nv=5
      integer maxvdw,i,idnode
      integer, dimension(nv) :: fail

      do i=1,nv
        fail(i) = 0
      enddo
      allocate (ltpvdw(maxvdw),stat=fail(1))
      allocate (lstvdw(maxvdw),stat=fail(2))
      allocate (steadd(maxvdw),stat=fail(3))
      allocate (prmvdw(maxvdw,mxpvdw),stat=fail(4))
      allocate (vvv(mxgrid,maxvdw),stat=fail(5))

      do i=1,nv
        if(fail(i).gt.0)then
            if(idnode.eq.0)write(nrite,'(10i5)')fail(i)
            call error(idnode, 1002)
        endif
      enddo
      end subroutine alloc_vdw_arrays

      end module vdw_module

!     -*- f90 -*-
!     This file is autogenerated with f2py (version:2)
!     It contains Fortran 90 wrappers to fortran functions.

      subroutine f2py_readr_getdims_real_table(r,s,f2pysetdata,flag)
      use readr, only: d => real_table

      integer flag
      external f2pysetdata
      logical ns
      integer r,i
      integer(8) s(*)
      ns = .FALSE.
      if (allocated(d)) then
         do i=1,r
            if ((size(d,i).ne.s(i)).and.(s(i).ge.0)) then
               ns = .TRUE.
            end if
         end do
         if (ns) then
            deallocate(d)
         end if
      end if
      if ((.not.allocated(d)).and.(s(1).ge.1)) then
       allocate(d(s(1),s(2)))
      end if
      if (allocated(d)) then
         do i=1,r
            s(i) = size(d,i)
         end do
      end if
      flag = 1
      call f2pysetdata(d,allocated(d))
      end subroutine f2py_readr_getdims_real_table
      subroutine f2py_readr_getdims_integer_table(r,s,f2pysetdata,flag)
      use readr, only: d => integer_table

      integer flag
      external f2pysetdata
      logical ns
      integer r,i
      integer(8) s(*)
      ns = .FALSE.
      if (allocated(d)) then
         do i=1,r
            if ((size(d,i).ne.s(i)).and.(s(i).ge.0)) then
               ns = .TRUE.
            end if
         end do
         if (ns) then
            deallocate(d)
         end if
      end if
      if ((.not.allocated(d)).and.(s(1).ge.1)) then
       allocate(d(s(1),s(2)))
      end if
      if (allocated(d)) then
         do i=1,r
            s(i) = size(d,i)
         end do
      end if
      flag = 1
      call f2pysetdata(d,allocated(d))
      end subroutine f2py_readr_getdims_integer_table
      subroutine f2py_readr_getdims_long_table(r,s,f2pysetdata,flag)
      use readr, only: d => long_table

      integer flag
      external f2pysetdata
      logical ns
      integer r,i
      integer(8) s(*)
      ns = .FALSE.
      if (allocated(d)) then
         do i=1,r
            if ((size(d,i).ne.s(i)).and.(s(i).ge.0)) then
               ns = .TRUE.
            end if
         end do
         if (ns) then
            deallocate(d)
         end if
      end if
      if ((.not.allocated(d)).and.(s(1).ge.1)) then
       allocate(d(s(1),s(2)))
      end if
      if (allocated(d)) then
         do i=1,r
            s(i) = size(d,i)
         end do
      end if
      flag = 1
      call f2pysetdata(d,allocated(d))
      end subroutine f2py_readr_getdims_long_table
      subroutine f2py_readr_getdims_byte_table(r,s,f2pysetdata,flag)
      use readr, only: d => byte_table

      integer flag
      external f2pysetdata
      logical ns
      integer r,i
      integer(8) s(*)
      ns = .FALSE.
      if (allocated(d)) then
         do i=1,r
            if ((size(d,i).ne.s(i)).and.(s(i).ge.0)) then
               ns = .TRUE.
            end if
         end do
         if (ns) then
            deallocate(d)
         end if
      end if
      if ((.not.allocated(d)).and.(s(1).ge.1)) then
       allocate(d(s(1),s(2)))
      end if
      if (allocated(d)) then
         do i=1,r
            s(i) = size(d,i)
         end do
      end if
      flag = 1
      call f2pysetdata(d,allocated(d))
      end subroutine f2py_readr_getdims_byte_table
      subroutine f2pywrap_readr_read_part (repo, iout, cpu_list, mode, x&
     &lims, verbose, longint, f2py_cpu_list_d0, f2py_xlims_d0)
      use readr, only : read_part
      integer iout
      logical verbose
      logical longint
      integer f2py_cpu_list_d0
      integer f2py_xlims_d0
      character(len=128) repo
      integer cpu_list(f2py_cpu_list_d0)
      character(len=10) mode
      real(kind=8) xlims(f2py_xlims_d0)
      call read_part(repo, iout, cpu_list, mode, xlims, verbose, longint&
     &)
      end subroutine f2pywrap_readr_read_part
      subroutine f2pywrap_readr_charind (charindf2pywrap, iout)
      use readr, only : charind
      integer iout
      character(len=5) charindf2pywrap
      charindf2pywrap = charind(iout)
      end subroutine f2pywrap_readr_charind
      subroutine f2pywrap_readr_part_filename (part_filenamef2pywrap, re&
     &po, iout, icpu, mode)
      use readr, only : part_filename
      integer iout
      integer icpu
      character(len=128) repo
      character(len=10) mode
      character(len=128) part_filenamef2pywrap
      part_filenamef2pywrap = part_filename(repo, iout, icpu, mode)
      end subroutine f2pywrap_readr_part_filename
      subroutine f2pywrap_readr_sinkprop_filename (sinkprop_filenamef2py&
     &wrap, repo, iout)
      use readr, only : sinkprop_filename
      integer iout
      character(len=128) repo
      character(len=128) sinkprop_filenamef2pywrap
      sinkprop_filenamef2pywrap = sinkprop_filename(repo, iout)
      end subroutine f2pywrap_readr_sinkprop_filename
      subroutine f2pywrap_readr_sink_filename (sink_filenamef2pywrap, re&
     &po, iout, icpu)
      use readr, only : sink_filename
      integer iout
      integer icpu
      character(len=128) repo
      character(len=128) sink_filenamef2pywrap
      sink_filenamef2pywrap = sink_filename(repo, iout, icpu)
      end subroutine f2pywrap_readr_sink_filename
      
      subroutine f2pyinitreadr(f2pysetupfunc)
      use readr, only : real_table
      use readr, only : integer_table
      use readr, only : long_table
      use readr, only : byte_table
      use readr, only : nhvar
      use readr, only : aexp
      use readr, only : read_sinkprop
      use readr, only : read_sink
      use readr, only : close
      use readr, only : skip_read
      use readr, only : progress_bar
      interface 
      subroutine f2pywrap_readr_read_part (repo, iout, cpu_list, mode, x&
     &lims, verbose, longint, f2py_cpu_list_d0, f2py_xlims_d0)
      integer iout
      logical verbose
      logical longint
      integer f2py_cpu_list_d0
      integer f2py_xlims_d0
      character(len=128) repo
      integer cpu_list(f2py_cpu_list_d0)
      character(len=10) mode
      real(kind=8) xlims(f2py_xlims_d0)
      end subroutine f2pywrap_readr_read_part 
      subroutine f2pywrap_readr_charind (charindf2pywrap, charind, iout)
      integer iout
      character(len=5) charind
      character(len=5) charindf2pywrap
      end subroutine f2pywrap_readr_charind 
      subroutine f2pywrap_readr_part_filename (part_filenamef2pywrap, pa&
     &rt_filename, repo, iout, icpu, mode)
      integer iout
      integer icpu
      character(len=128) repo
      character(len=10) mode
      character(len=128) part_filename
      character(len=128) part_filenamef2pywrap
      end subroutine f2pywrap_readr_part_filename 
      subroutine f2pywrap_readr_sinkprop_filename (sinkprop_filenamef2py&
     &wrap, sinkprop_filename, repo, iout)
      integer iout
      character(len=128) repo
      character(len=128) sinkprop_filename
      character(len=128) sinkprop_filenamef2pywrap
      end subroutine f2pywrap_readr_sinkprop_filename 
      subroutine f2pywrap_readr_sink_filename (sink_filenamef2pywrap, si&
     &nk_filename, repo, iout, icpu)
      integer iout
      integer icpu
      character(len=128) repo
      character(len=128) sink_filename
      character(len=128) sink_filenamef2pywrap
      end subroutine f2pywrap_readr_sink_filename
      end interface
      external f2pysetupfunc
      external f2py_readr_getdims_real_table
      external f2py_readr_getdims_integer_table
      external f2py_readr_getdims_long_table
      external f2py_readr_getdims_byte_table
      call f2pysetupfunc(f2py_readr_getdims_real_table,f2py_readr_getdim&
     &s_integer_table,f2py_readr_getdims_long_table,f2py_readr_getdims_b&
     &yte_table,nhvar,aexp,f2pywrap_readr_read_part,read_sinkprop,read_s&
     &ink,close,skip_read,f2pywrap_readr_charind,f2pywrap_readr_part_fil&
     &ename,f2pywrap_readr_sinkprop_filename,f2pywrap_readr_sink_filenam&
     &e,progress_bar)
      end subroutine f2pyinitreadr


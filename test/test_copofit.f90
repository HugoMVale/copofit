module test_copofit
!! Test for module 'break1d' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use copofit, only: dp, fit, dtype_mayo, dtype_fp1x
   implicit none
   private

   public :: collect_tests_fit

   logical, parameter :: verbose = .false.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_fit(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fit F1(f1) data", test_fit_mayo_data), &
                  new_unittest("fit F1(x, f10) data", test_fit_drift_data)]

   end subroutine

   subroutine test_fit_mayo_data(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n=7
      real(dp), dimension(n) :: f1, fp1, x, f10, sx, sf1, sfp1
      real(dp) :: rcopo(2)
      integer, dimension(n) :: dsid, dtype

      x = 0
      f10 = 0
      sx = 0
      dtype = dtype_mayo
      sf1 = 1e-10_dp
      sfp1 = 0.05_dp
      rcopo = [1.0_dp, 1.0_dp]
      
      f1 = [0.100_dp, 0.300_dp, 0.400_dp, 0.500_dp, 0.600_dp, 0.700_dp, 0.800_dp]
      fp1 = [0.059_dp, 0.243_dp, 0.364_dp, 0.486_dp, 0.583_dp, 0.721_dp, 0.824_dp]

      call fit(x, f1, fp1, f10, sx, sf1, sfp1, dtype, rcopo)
      call check(error, rcopo, [1.43_dp, 1.67_dp], rel=.true., thr=1e-2_dp)

   end subroutine test_fit_mayo_data

   subroutine test_fit_drift_data(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n=7
      real(dp), dimension(n) :: f1, fp1, x, f10, sx, sf1, sfp1
      real(dp) :: rcopo(2)
      integer, dimension(n) :: dsid, dtype

      f1 = 0
      sf1 = 0
      fp1 = [0.059_dp, 0.243_dp, 0.364_dp, 0.486_dp, 0.583_dp, 0.721_dp, 0.824_dp]
      sfp1 = 0.05_dp
      x = 0.01_dp
      sx = 1e-10_dp
      f10 = [0.100_dp, 0.300_dp, 0.400_dp, 0.500_dp, 0.600_dp, 0.700_dp, 0.800_dp]     
      dtype = dtype_fp1x
      
      rcopo = [1.0_dp, 1.0_dp]    

      call fit(x, f1, fp1, f10, sx, sf1, sfp1, dtype, rcopo)
      call check(error, rcopo, [1.43_dp, 1.67_dp], rel=.true., thr=1e-2_dp)

   end subroutine test_fit_drift_data

end module test_copofit

module copofit
   use iso_fortran_env, only: dp => real64
   use odrpack, only: odr
   implicit none
   private

   public :: dp
   public :: dtype_mayo, dtype_f1x, dtype_fp1x
   public :: fit

   real(dp), parameter :: zero = 0.0_dp, &
                          one = 1.0_dp, &
                          eps = epsilon(one)

   integer, parameter :: npar = 2, &
                         dtype_mayo = 1, &
                         dtype_f1x = 2, &
                         dtype_fp1x = 3

   abstract interface
      pure real(dp) function integrand(x, y)
         import :: dp
         real(dp), intent(in) :: x
         real(dp), intent(in) :: y
      end function
   end interface

contains

   subroutine fit(x, f1, fp1, f10, sx, sf1, sfp1, dtype, rcopo)

      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: f1(:)
      real(dp), intent(in) :: fp1(:)
      real(dp), intent(in), target :: f10(:)
      real(dp), intent(in) :: sx(:)
      real(dp), intent(in) :: sf1(:)
      real(dp), intent(in) :: sfp1(:)
      integer, intent(in), target :: dtype(:)
      real(dp), intent(inout) :: rcopo(npar)

      common/odr_common/f10_, dtype_

      integer, parameter :: np = npar, m = 1, nq = 1
      real(dp), dimension(size(x), m) :: xdata
      real(dp), dimension(size(x), nq) :: ydata
      real(dp), dimension(size(x), nq, nq) :: we
      real(dp), dimension(size(x), m, m) :: wd
      real(dp), pointer :: f10_(:)
      integer, pointer :: dtype_(:)
      integer :: i, n, info

      f10_ => f10
      dtype_ => dtype
      n = size(x)

      do concurrent(i=1:n)
         select case (dtype(i))
         case (dtype_mayo)
            xdata(i, 1) = f1(i)
            ydata(i, 1) = fp1(i)
            wd(i, 1, 1) = sf1(i)
            we(i, 1, 1) = sfp1(i)
         case (dtype_f1x)
            xdata(i, 1) = x(i)
            ydata(i, 1) = f1(i)
            wd(i, 1, 1) = sx(i)
            we(i, 1, 1) = sf1(i)
         case (dtype_fp1x)
            xdata(i, 1) = x(i)
            ydata(i, 1) = fp1(i)
            wd(i, 1, 1) = sx(i)
            we(i, 1, 1) = sfp1(i)
         end select
      end do

      we = 1/we**2
      wd = 1/wd**2

      call odr(fcn, n, m, np, nq, rcopo, ydata, xdata, &
               we=we, wd=wd, &
               lower=[zero, zero], upper=[1e2_dp, 1e2_dp], &
               info=info)

   end subroutine

   pure subroutine fcn(n, m, np, nq, ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, &
                       ldifx, ideval, f, fjacb, fjacd, istop)
   !! User-supplied subroutine for evaluating the model.
      integer, intent(in) :: ideval, ldifx, ldm, ldn, ldnp, m, n, np, nq
      integer, intent(in) :: ifixb(np), ifixx(ldifx, m)
      real(dp), intent(in) :: beta(np), xplusd(ldn, m)
      real(dp), intent(out) :: f(ldn, nq), fjacb(ldn, ldnp, nq), fjacd(ldn, ldm, nq)
      integer, intent(out) :: istop

      real(dp), pointer :: f10_(:)
      integer, pointer :: dtype_(:)
      integer :: i

      common/odr_common/f10_, dtype_

      istop = 0
      fjacb(:, :, :) = zero
      fjacd(:, :, :) = zero

      if (mod(ideval, 10) >= 1) then
         do concurrent(i=1:n)
            if (dtype_(i) == dtype_mayo) then
               f(i, 1) = mayo_equation(xplusd(i, 1), beta(1), beta(2))
            else
               f(i, 1) = f1x_equation(xplusd(i, 1), f10_(i), beta(1), beta(2))
               if (dtype_(i) == dtype_fp1x) then
                  ! Convert f1 to F1
                  f(i, 1) = f(i, 1) + (f10_(i) - f(i, 1))/(xplusd(i, 1) + eps)
               end if
            end if
         end do
      end if

   end subroutine

   pure elemental function mayo_equation(f1, r1, r2) result(fp1)
   !! Calculate the instantaneous copolymer composition using the Mayo-Lewis equation.
      real(dp), intent(in) :: f1
         !! Molar fraction of M1.
      real(dp), intent(in) :: r1
         !! Reactivity ratio of M1.
      real(dp), intent(in) :: r2
         !! Reactivity ratio of M2.
      real(dp) :: fp1
         !! Instantaneous copolymer composition, F1.
      real(dp) :: f2
      f2 = one - f1
      fp1 = (r1*f1**2 + f1*f2)/(r1*f1**2 + 2*f1*f2 + r2*f2**2)
   end function

   pure function f1x_equation(x, f10, r1, r2) result(f1x)
      real(dp), intent(in) :: x
         !! Total monomer conversion.
      real(dp), intent(in) :: f10
         !! Initial molar fraction of M1.
      real(dp), intent(in) :: r1
         !! Reactivity ratio of M1.
      real(dp), intent(in) :: r2
         !! Reactivity ratio of M2.
      real(dp) :: f1x
         !! Monomer fraction of M1 at a conversion `x`.

      real(dp), parameter :: dx = 1e-2_dp ! to be adjusted
      real(dp) :: xf

      xf = x
      call rkint(zero, xf, f10, f1x, dx, df1dx)

   contains

      pure real(dp) function df1dx(x_, f1_)
      !! Skeist equation for a binary system.
         real(dp), intent(in) :: x_, f1_
         df1dx = (f1_ - mayo_equation(f1_, r1, r2))/(one - x_ + eps)
      end function

   end function

   pure subroutine rkint(x0, xf, y0, yf, dx, f)
   !! Explicit, constant-step, Runga-Kutta integration.
      real(dp), intent(in) :: x0
         !! Initial value of `x`.
      real(dp), intent(inout) :: xf
         !! Final value of `x`.
      real(dp), intent(in) :: y0
         !! Initial value of `y`, i.e. `y(x0)`.
      real(dp), intent(out) :: yf
         !! Value of `y(xf)`.
      real(dp), intent(in) :: dx
         !! Step size.
      procedure(integrand) :: f
         !! Function to be integrated.

      integer :: i, nsteps
      real(dp) :: dxf

      nsteps = floor((xf - x0)/dx)
      dxf = (xf - x0) - nsteps*dx

      xf = x0
      yf = y0
      do i = 1, nsteps
         call rk_step(xf, yf, dx, f)
      end do

      call rk_step(xf, yf, dxf, f)

   end subroutine

   pure subroutine rk_step(x, y, dx, f)
   !! 2nd-order explicit mid-point step.
      real(dp), intent(inout) :: x
         !! Independent variable. `x(n)` on input. `x(n+1)` on output.
      real(dp), intent(inout) :: y
         !! Dependent variable. `y(n)` on input. `y(n+1)` on output.
      real(dp), intent(in) :: dx
         !! Step size.
      procedure(integrand) :: f
         !! Function to be integrated.
      y = y + f(x + dx/2, y + f(x, y)*dx/2)*dx
      x = x + dx
   end subroutine rk_step

end module copofit

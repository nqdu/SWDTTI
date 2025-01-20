!===========================================================================
!============================= AUTO CODE FROM SYMPY ========================
!===========================================================================

subroutine get_comp1 (NGL,A,C,F,L,N,theta0,dphi,jaco, &
                        weight,hp,hpT, K0U,K0V,K0W,K1U,K1V,K1W,&
                        K2U,K2V,K2W) bind(c,name='get_comp1_')
  use iso_c_binding,only: c_int,dp => c_double ,dcp => c_double_complex
  implicit none

  integer(c_int),value,intent(in)         :: NGL
  real(dp),dimension(NGL),intent(in)  :: A,C,F,L,N,theta0,dphi,weight
  real(dp),dimension(NGL,NGL),intent(in)  :: hp,hpT
  real(dp),value :: jaco 
  complex(dcp),dimension(NGL,NGL),intent(out) :: K1U,K1V,K1W,K2U,K2V,K2W,K0U,K0V,K0W 
  
  !local vars
  integer(c_int)        :: i,j 
  complex(dcp),dimension(NGL)  :: temp
  real(dp),dimension(NGL) :: costh0, sinth0,cosphi,sinphi
  complex(dcp),parameter :: imag_i = cmplx(0.0,1.0,kind=dcp)

  !init matrices
  K0U(:,:) = (0.0,0.0); K0V(:,:) = (0.0,0.0); K0W(:,:) = (0.0,0.0);
  K1U(:,:) = (0.0,0.0); K1V(:,:) = (0.0,0.0); K1W(:,:) = (0.0,0.0);
  K2U(:,:) = (0.0,0.0); K2V(:,:) = (0.0,0.0); K2W(:,:) = (0.0,0.0);
  costh0 = cos(theta0); sinth0 = sin(theta0);
  cosphi = cos(dphi); sinphi = sin(dphi);

  ! k^0, Udot psidot
  temp(:) = A*sinth0**2*cosphi**2*costh0**2 + C*sinth0**2*cosphi**2*costh0**2 &
      - 2*F*sinth0**2*cosphi**2*costh0**2 &
      - 4*L*sinth0**2*cosphi**2*costh0**2 + L*sinth0**2*cosphi**2 &
      + L*costh0**2 - N*sinth0**2*cosphi**2 - N*costh0**2 &
      + N
  do j=1,NGL; do i=1,NGL; 
    K0U(i,j) = K0U(i,j) + sum(weight * temp / jaco * hp(:,j) * hp(:,i))
  enddo; enddo; 

  ! k^0, Wdot psidot
  temp(:) = A*sinphi*sinth0**2*cosphi*costh0**2 &
      + C*sinphi*sinth0**2*cosphi*costh0**2 - 2*F*sinphi*sinth0**2*cosphi*costh0**2 &
      - 4*L*sinphi*sinth0**2*cosphi*costh0**2 &
      + L*sinphi*sinth0**2*cosphi &
      - N*sinphi*sinth0**2*cosphi
  do j=1,NGL; do i=1,NGL; 
    K0W(i,j) = K0W(i,j) + sum(weight * temp / jaco * hp(:,j) * hp(:,i))
  enddo; enddo; 

  ! k^0, Vdot psidot
  temp(:) = A*sinth0*cosphi*costh0**3 - A*sinth0*cosphi*costh0 &
      + C*sinth0*cosphi*costh0**3 - 2*F*sinth0*cosphi*costh0**3 &
      + F*sinth0*cosphi*costh0 &
      - 4*L*sinth0*cosphi*costh0**3 + 2*L*sinth0*cosphi*costh0
  do j=1,NGL; do i=1,NGL; 
    K0V(i,j) = K0V(i,j) + sum(weight * temp / jaco * hp(:,j) * hp(:,i))
  enddo; enddo; 

  ! k^1, U psidot
  temp(:) = - imag_i*A*sinth0**3*cosphi**3*costh0 &
      + imag_i*A*sinth0*cosphi*costh0 - imag_i*C*sinth0**3*cosphi**3*costh0 &
      + 2*imag_i*F*sinth0**3*cosphi**3*costh0 &
      - imag_i*F*sinth0*cosphi*costh0 &
      + 4*imag_i*L*sinth0**3*cosphi**3*costh0 &
      - 2*imag_i*L*sinth0*cosphi*costh0
  do j=1,NGL; do i=1,NGL; 
    K1U(i,j) = K1U(i,j) + temp(j) * weight(j) * hpT(i,j)
  enddo; enddo; 

  ! k^1, Udot psi
  temp(:) = imag_i*A*sinth0**3*cosphi**3*costh0 &
      - imag_i*A*sinth0*cosphi*costh0 + imag_i*C*sinth0**3*cosphi**3*costh0 &
      - 2*imag_i*F*sinth0**3*cosphi**3*costh0 &
      + imag_i*F*sinth0*cosphi*costh0 &
      - 4*imag_i*L*sinth0**3*cosphi**3*costh0 &
      + 2*imag_i*L*sinth0*cosphi*costh0
  do j=1,NGL; do i=1,NGL; 
    K1U(i,j) = K1U(i,j) + temp(i) * weight(i) * hp(i,j)
  enddo; enddo; 

  ! k^1, W psidot
  temp(:) = - imag_i*A*sinphi*sinth0**3*cosphi**2*costh0 &
      - imag_i*C*sinphi*sinth0**3*cosphi**2*costh0 &
      + 2*imag_i*F*sinphi*sinth0**3*cosphi**2*costh0 &
      + 4*imag_i*L*sinphi*sinth0**3*cosphi**2*costh0 &
      - imag_i*L*sinphi*sinth0*costh0 &
      + imag_i*N*sinphi*sinth0*costh0
  do j=1,NGL; do i=1,NGL; 
    K1W(i,j) = K1W(i,j) + temp(j) * weight(j) * hpT(i,j)
  enddo; enddo; 

  ! k^1, Wdot psi
  temp(:) = imag_i*A*sinphi*sinth0**3*cosphi**2*costh0 &
      - imag_i*A*sinphi*sinth0*costh0 &
      + imag_i*C*sinphi*sinth0**3*cosphi**2*costh0 - 2*imag_i*F*sinphi*sinth0**3*cosphi**2*costh0 &
      + imag_i*F*sinphi*sinth0*costh0 &
      - 4*imag_i*L*sinphi*sinth0**3*cosphi**2*costh0 &
      + 2*imag_i*N*sinphi*sinth0*costh0
  do j=1,NGL; do i=1,NGL; 
    K1W(i,j) = K1W(i,j) + temp(i) * weight(i) * hp(i,j)
  enddo; enddo; 

  ! k^1, V psidot
  temp(:) = - imag_i*A*sinth0**2*cosphi**2*costh0**2 &
      - imag_i*C*sinth0**2*cosphi**2*costh0**2 &
      + 2*imag_i*F*sinth0**2*cosphi**2*costh0**2 + 4*imag_i*L*sinth0**2*cosphi**2*costh0**2 &
      - imag_i*L*sinth0**2*cosphi**2 &
      - imag_i*L*costh0**2 + imag_i*N*sinth0**2*cosphi**2 &
      + imag_i*N*costh0**2 - imag_i*N
  do j=1,NGL; do i=1,NGL; 
    K1V(i,j) = K1V(i,j) + temp(j) * weight(j) * hpT(i,j)
  enddo; enddo; 

  ! k^1, Vdot psi
  temp(:) = imag_i*A*sinth0**2*cosphi**2*costh0**2 &
      - imag_i*A*sinth0**2*cosphi**2 - imag_i*A*costh0**2 + imag_i*A &
      + imag_i*C*sinth0**2*cosphi**2*costh0**2 &
      - 2*imag_i*F*sinth0**2*cosphi**2*costh0**2 + imag_i*F*sinth0**2*cosphi**2 &
      + imag_i*F*costh0**2 &
      - 4*imag_i*L*sinth0**2*cosphi**2*costh0**2 + 2*imag_i*N*sinth0**2*cosphi**2 &
      + 2*imag_i*N*costh0**2 - 2*imag_i*N
  do j=1,NGL; do i=1,NGL; 
    K1V(i,j) = K1V(i,j) + temp(i) * weight(i) * hp(i,j)
  enddo; enddo; 

  ! k^2, U psi
  temp(:) = A*sinth0**4*cosphi**4 - 2*A*sinth0**2*cosphi**2 + A &
      + C*sinth0**4*cosphi**4 - 2*F*sinth0**4*cosphi**4 &
      + 2*F*sinth0**2*cosphi**2 - 4*L*sinth0**4*cosphi**4 &
      + 4*L*sinth0**2*cosphi**2
  do i=1,NGL; K2U(i,i) = K2U(i,i) + temp(i) * weight(i) * jaco; enddo

  ! k^2, W psi
  temp(:) = A*sinphi*sinth0**4*cosphi**3 - A*sinphi*sinth0**2*cosphi &
      + C*sinphi*sinth0**4*cosphi**3 &
      - 2*F*sinphi*sinth0**4*cosphi**3 + F*sinphi*sinth0**2*cosphi &
      - 4*L*sinphi*sinth0**4*cosphi**3 &
      + 2*L*sinphi*sinth0**2*cosphi
  do i=1,NGL; K2W(i,i) = K2W(i,i) + temp(i) * weight(i) * jaco; enddo

  ! k^2, V psi
  temp(:) = A*sinth0**3*cosphi**3*costh0 - A*sinth0*cosphi*costh0 &
      + C*sinth0**3*cosphi**3*costh0 &
      - 2*F*sinth0**3*cosphi**3*costh0 + F*sinth0*cosphi*costh0 &
      - 4*L*sinth0**3*cosphi**3*costh0 &
      + 2*L*sinth0*cosphi*costh0
  do i=1,NGL; K2V(i,i) = K2V(i,i) + temp(i) * weight(i) * jaco; enddo

  !transpose
  K0U = transpose(K0U); K0V = transpose(K0V); K0W = transpose(K0W);
  K1U = transpose(K1U); K1V = transpose(K1V); K1W = transpose(K1W);
  K2U = transpose(K2U); K2V = transpose(K2V); K2W = transpose(K2W);

end subroutine get_comp1

subroutine get_comp2 (NGL,A,C,F,L,N,theta0,dphi,jaco, &
                        weight,hp,hpT, K0U,K0V,K0W,K1U,K1V,K1W,&
                        K2U,K2V,K2W) bind(c,name='get_comp2_')
  use iso_c_binding,only: c_int,dp => c_double ,dcp => c_double_complex
  implicit none

  integer(c_int),value,intent(in)         :: NGL
  real(dp),dimension(NGL),intent(in)  :: A,C,F,L,N,theta0,dphi,weight
  real(dp),dimension(NGL,NGL),intent(in)  :: hp,hpT
  real(dp),value :: jaco 
  complex(dcp),dimension(NGL,NGL),intent(out) :: K1U,K1V,K1W,K2U,K2V,K2W,K0U,K0V,K0W 
  
  !local vars
  integer(c_int)        :: i,j 
  complex(dcp),dimension(NGL)  :: temp
  real(dp),dimension(NGL) :: costh0, sinth0,cosphi,sinphi
  complex(dcp),parameter :: imag_i = cmplx(0.0,1.0,kind=dcp)

  !init matrices
  K0U(:,:) = (0.0,0.0); K0V(:,:) = (0.0,0.0); K0W(:,:) = (0.0,0.0);
  K1U(:,:) = (0.0,0.0); K1V(:,:) = (0.0,0.0); K1W(:,:) = (0.0,0.0);
  K2U(:,:) = (0.0,0.0); K2V(:,:) = (0.0,0.0); K2W(:,:) = (0.0,0.0);
  costh0 = cos(theta0); sinth0 = sin(theta0);
  cosphi = cos(dphi); sinphi = sin(dphi);

  ! k^0, Udot psidot
  temp(:) = A*sinphi*sinth0**2*cosphi*costh0**2 &
      + C*sinphi*sinth0**2*cosphi*costh0**2 - 2*F*sinphi*sinth0**2*cosphi*costh0**2 &
      - 4*L*sinphi*sinth0**2*cosphi*costh0**2 &
      + L*sinphi*sinth0**2*cosphi &
      - N*sinphi*sinth0**2*cosphi
  do j=1,NGL; do i=1,NGL; 
    K0U(i,j) = K0U(i,j) + sum(weight * temp / jaco * hp(:,j) * hp(:,i))
  enddo; enddo; 

  ! k^0, Wdot psidot
  temp(:) = A*sinphi**2*sinth0**2*costh0**2 + C*sinphi**2*sinth0**2*costh0**2 &
      - 2*F*sinphi**2*sinth0**2*costh0**2 &
      - 4*L*sinphi**2*sinth0**2*costh0**2 + L*sinphi**2*sinth0**2 &
      + L*costh0**2 - N*sinphi**2*sinth0**2 - N*costh0**2 &
      + N
  do j=1,NGL; do i=1,NGL; 
    K0W(i,j) = K0W(i,j) + sum(weight * temp / jaco * hp(:,j) * hp(:,i))
  enddo; enddo; 

  ! k^0, Vdot psidot
  temp(:) = A*sinphi*sinth0*costh0**3 - A*sinphi*sinth0*costh0 &
      + C*sinphi*sinth0*costh0**3 - 2*F*sinphi*sinth0*costh0**3 &
      + F*sinphi*sinth0*costh0 &
      - 4*L*sinphi*sinth0*costh0**3 + 2*L*sinphi*sinth0*costh0
  do j=1,NGL; do i=1,NGL; 
    K0V(i,j) = K0V(i,j) + sum(weight * temp / jaco * hp(:,j) * hp(:,i))
  enddo; enddo; 

  ! k^1, U psidot
  temp(:) = - imag_i*A*sinphi*sinth0**3*cosphi**2*costh0 &
      + imag_i*A*sinphi*sinth0*costh0 &
      - imag_i*C*sinphi*sinth0**3*cosphi**2*costh0 + 2*imag_i*F*sinphi*sinth0**3*cosphi**2*costh0 &
      - imag_i*F*sinphi*sinth0*costh0 &
      + 4*imag_i*L*sinphi*sinth0**3*cosphi**2*costh0 &
      - 2*imag_i*N*sinphi*sinth0*costh0
  do j=1,NGL; do i=1,NGL; 
    K1U(i,j) = K1U(i,j) + temp(j) * weight(j) * hpT(i,j)
  enddo; enddo; 

  ! k^1, Udot psi
  temp(:) = imag_i*A*sinphi*sinth0**3*cosphi**2*costh0 &
      + imag_i*C*sinphi*sinth0**3*cosphi**2*costh0 &
      - 2*imag_i*F*sinphi*sinth0**3*cosphi**2*costh0 &
      - 4*imag_i*L*sinphi*sinth0**3*cosphi**2*costh0 &
      + imag_i*L*sinphi*sinth0*costh0 &
      - imag_i*N*sinphi*sinth0*costh0
  do j=1,NGL; do i=1,NGL; 
    K1U(i,j) = K1U(i,j) + temp(i) * weight(i) * hp(i,j)
  enddo; enddo; 

  ! k^1, W psidot
  temp(:) = - imag_i*A*sinphi**2*sinth0**3*cosphi*costh0 &
      - imag_i*C*sinphi**2*sinth0**3*cosphi*costh0 &
      + 2*imag_i*F*sinphi**2*sinth0**3*cosphi*costh0 &
      + 4*imag_i*L*sinphi**2*sinth0**3*cosphi*costh0 &
      - imag_i*L*sinth0*cosphi*costh0 &
      + imag_i*N*sinth0*cosphi*costh0
  do j=1,NGL; do i=1,NGL; 
    K1W(i,j) = K1W(i,j) + temp(j) * weight(j) * hpT(i,j)
  enddo; enddo; 

  ! k^1, Wdot psi
  temp(:) = imag_i*A*sinphi**2*sinth0**3*cosphi*costh0 &
      + imag_i*C*sinphi**2*sinth0**3*cosphi*costh0 &
      - 2*imag_i*F*sinphi**2*sinth0**3*cosphi*costh0 &
      - 4*imag_i*L*sinphi**2*sinth0**3*cosphi*costh0 &
      + imag_i*L*sinth0*cosphi*costh0 &
      - imag_i*N*sinth0*cosphi*costh0
  do j=1,NGL; do i=1,NGL; 
    K1W(i,j) = K1W(i,j) + temp(i) * weight(i) * hp(i,j)
  enddo; enddo; 

  ! k^1, V psidot
  temp(:) = - imag_i*A*sinphi*sinth0**2*cosphi*costh0**2 &
      - imag_i*C*sinphi*sinth0**2*cosphi*costh0**2 &
      + 2*imag_i*F*sinphi*sinth0**2*cosphi*costh0**2 &
      + 4*imag_i*L*sinphi*sinth0**2*cosphi*costh0**2 &
      - imag_i*L*sinphi*sinth0**2*cosphi &
      + imag_i*N*sinphi*sinth0**2*cosphi
  do j=1,NGL; do i=1,NGL; 
    K1V(i,j) = K1V(i,j) + temp(j) * weight(j) * hpT(i,j)
  enddo; enddo; 

  ! k^1, Vdot psi
  temp(:) = imag_i*A*sinphi*sinth0**2*cosphi*costh0**2 &
      - imag_i*A*sinphi*sinth0**2*cosphi &
      + imag_i*C*sinphi*sinth0**2*cosphi*costh0**2 - 2*imag_i*F*sinphi*sinth0**2*cosphi*costh0**2 &
      + imag_i*F*sinphi*sinth0**2*cosphi &
      - 4*imag_i*L*sinphi*sinth0**2*cosphi*costh0**2 &
      + 2*imag_i*N*sinphi*sinth0**2*cosphi
  do j=1,NGL; do i=1,NGL; 
    K1V(i,j) = K1V(i,j) + temp(i) * weight(i) * hp(i,j)
  enddo; enddo; 

  ! k^2, U psi
  temp(:) = A*sinphi*sinth0**4*cosphi**3 - A*sinphi*sinth0**2*cosphi &
      + C*sinphi*sinth0**4*cosphi**3 &
      - 2*F*sinphi*sinth0**4*cosphi**3 + F*sinphi*sinth0**2*cosphi &
      - 4*L*sinphi*sinth0**4*cosphi**3 &
      + 2*L*sinphi*sinth0**2*cosphi
  do i=1,NGL; K2U(i,i) = K2U(i,i) + temp(i) * weight(i) * jaco; enddo

  ! k^2, W psi
  temp(:) = A*sinphi**2*sinth0**4*cosphi**2 + C*sinphi**2*sinth0**4*cosphi**2 &
      - 2*F*sinphi**2*sinth0**4*cosphi**2 &
      - 4*L*sinphi**2*sinth0**4*cosphi**2 + L*sinphi**2*sinth0**2 &
      + L*sinth0**2*cosphi**2 - N*sinphi**2*sinth0**2 &
      - N*sinth0**2*cosphi**2 + N
  do i=1,NGL; K2W(i,i) = K2W(i,i) + temp(i) * weight(i) * jaco; enddo

  ! k^2, V psi
  temp(:) = A*sinphi*sinth0**3*cosphi**2*costh0 &
      + C*sinphi*sinth0**3*cosphi**2*costh0 - 2*F*sinphi*sinth0**3*cosphi**2*costh0 &
      - 4*L*sinphi*sinth0**3*cosphi**2*costh0 &
      + L*sinphi*sinth0*costh0 &
      - N*sinphi*sinth0*costh0
  do i=1,NGL; K2V(i,i) = K2V(i,i) + temp(i) * weight(i) * jaco; enddo

  !transpose
  K0U = transpose(K0U); K0V = transpose(K0V); K0W = transpose(K0W);
  K1U = transpose(K1U); K1V = transpose(K1V); K1W = transpose(K1W);
  K2U = transpose(K2U); K2V = transpose(K2V); K2W = transpose(K2W);

end subroutine get_comp2

subroutine get_comp3 (NGL,A,C,F,L,N,theta0,dphi,jaco, &
                        weight,hp,hpT, K0U,K0V,K0W,K1U,K1V,K1W,&
                        K2U,K2V,K2W) bind(c,name='get_comp3_')
  use iso_c_binding,only: c_int,dp => c_double ,dcp => c_double_complex
  implicit none

  integer(c_int),value,intent(in)         :: NGL
  real(dp),dimension(NGL),intent(in)  :: A,C,F,L,N,theta0,dphi,weight
  real(dp),dimension(NGL,NGL),intent(in)  :: hp,hpT
  real(dp),value :: jaco 
  complex(dcp),dimension(NGL,NGL),intent(out) :: K1U,K1V,K1W,K2U,K2V,K2W,K0U,K0V,K0W 
  
  !local vars
  integer(c_int)        :: i,j 
  complex(dcp),dimension(NGL)  :: temp
  real(dp),dimension(NGL) :: costh0, sinth0,cosphi,sinphi
  complex(dcp),parameter :: imag_i = cmplx(0.0,1.0,kind=dcp)

  !init matrices
  K0U(:,:) = (0.0,0.0); K0V(:,:) = (0.0,0.0); K0W(:,:) = (0.0,0.0);
  K1U(:,:) = (0.0,0.0); K1V(:,:) = (0.0,0.0); K1W(:,:) = (0.0,0.0);
  K2U(:,:) = (0.0,0.0); K2V(:,:) = (0.0,0.0); K2W(:,:) = (0.0,0.0);
  costh0 = cos(theta0); sinth0 = sin(theta0);
  cosphi = cos(dphi); sinphi = sin(dphi);

  ! k^0, Udot psidot
  temp(:) = A*sinth0*cosphi*costh0**3 - A*sinth0*cosphi*costh0 &
      + C*sinth0*cosphi*costh0**3 - 2*F*sinth0*cosphi*costh0**3 &
      + F*sinth0*cosphi*costh0 &
      - 4*L*sinth0*cosphi*costh0**3 + 2*L*sinth0*cosphi*costh0
  do j=1,NGL; do i=1,NGL; 
    K0U(i,j) = K0U(i,j) + sum(weight * temp / jaco * hp(:,j) * hp(:,i))
  enddo; enddo; 

  ! k^0, Wdot psidot
  temp(:) = A*sinphi*sinth0*costh0**3 - A*sinphi*sinth0*costh0 &
      + C*sinphi*sinth0*costh0**3 - 2*F*sinphi*sinth0*costh0**3 &
      + F*sinphi*sinth0*costh0 &
      - 4*L*sinphi*sinth0*costh0**3 + 2*L*sinphi*sinth0*costh0
  do j=1,NGL; do i=1,NGL; 
    K0W(i,j) = K0W(i,j) + sum(weight * temp / jaco * hp(:,j) * hp(:,i))
  enddo; enddo; 

  ! k^0, Vdot psidot
  temp(:) = A*costh0**4 - 2*A*costh0**2 + A + C*costh0**4 &
      - 2*F*costh0**4 + 2*F*costh0**2 - 4*L*costh0**4 &
      + 4*L*costh0**2
  do j=1,NGL; do i=1,NGL; 
    K0V(i,j) = K0V(i,j) + sum(weight * temp / jaco * hp(:,j) * hp(:,i))
  enddo; enddo; 

  ! k^1, U psidot
  temp(:) = - imag_i*A*sinth0**2*cosphi**2*costh0**2 &
      + imag_i*A*sinth0**2*cosphi**2 + imag_i*A*costh0**2 - imag_i*A &
      - imag_i*C*sinth0**2*cosphi**2*costh0**2 &
      + 2*imag_i*F*sinth0**2*cosphi**2*costh0**2 &
      - imag_i*F*sinth0**2*cosphi**2 - imag_i*F*costh0**2 &
      + 4*imag_i*L*sinth0**2*cosphi**2*costh0**2 - 2*imag_i*N*sinth0**2*cosphi**2 &
      - 2*imag_i*N*costh0**2 + 2*imag_i*N
  do j=1,NGL; do i=1,NGL; 
    K1U(i,j) = K1U(i,j) + temp(j) * weight(j) * hpT(i,j)
  enddo; enddo; 

  ! k^1, Udot psi
  temp(:) = imag_i*A*sinth0**2*cosphi**2*costh0**2 &
      + imag_i*C*sinth0**2*cosphi**2*costh0**2 &
      - 2*imag_i*F*sinth0**2*cosphi**2*costh0**2 - 4*imag_i*L*sinth0**2*cosphi**2*costh0**2 &
      + imag_i*L*sinth0**2*cosphi**2 &
      + imag_i*L*costh0**2 - imag_i*N*sinth0**2*cosphi**2 &
      - imag_i*N*costh0**2 + imag_i*N
  do j=1,NGL; do i=1,NGL; 
    K1U(i,j) = K1U(i,j) + temp(i) * weight(i) * hp(i,j)
  enddo; enddo; 

  ! k^1, W psidot
  temp(:) = - imag_i*A*sinphi*sinth0**2*cosphi*costh0**2 &
      + imag_i*A*sinphi*sinth0**2*cosphi &
      - imag_i*C*sinphi*sinth0**2*cosphi*costh0**2 + 2*imag_i*F*sinphi*sinth0**2*cosphi*costh0**2 &
      - imag_i*F*sinphi*sinth0**2*cosphi &
      + 4*imag_i*L*sinphi*sinth0**2*cosphi*costh0**2 &
      - 2*imag_i*N*sinphi*sinth0**2*cosphi
  do j=1,NGL; do i=1,NGL; 
    K1W(i,j) = K1W(i,j) + temp(j) * weight(j) * hpT(i,j)
  enddo; enddo; 

  ! k^1, Wdot psi
  temp(:) = imag_i*A*sinphi*sinth0**2*cosphi*costh0**2 &
      + imag_i*C*sinphi*sinth0**2*cosphi*costh0**2 &
      - 2*imag_i*F*sinphi*sinth0**2*cosphi*costh0**2 &
      - 4*imag_i*L*sinphi*sinth0**2*cosphi*costh0**2 &
      + imag_i*L*sinphi*sinth0**2*cosphi &
      - imag_i*N*sinphi*sinth0**2*cosphi
  do j=1,NGL; do i=1,NGL; 
    K1W(i,j) = K1W(i,j) + temp(i) * weight(i) * hp(i,j)
  enddo; enddo; 

  ! k^1, V psidot
  temp(:) = - imag_i*A*sinth0*cosphi*costh0**3 &
      + imag_i*A*sinth0*cosphi*costh0 - imag_i*C*sinth0*cosphi*costh0**3 &
      + 2*imag_i*F*sinth0*cosphi*costh0**3 &
      - imag_i*F*sinth0*cosphi*costh0 + 4*imag_i*L*sinth0*cosphi*costh0**3 &
      - 2*imag_i*L*sinth0*cosphi*costh0
  do j=1,NGL; do i=1,NGL; 
    K1V(i,j) = K1V(i,j) + temp(j) * weight(j) * hpT(i,j)
  enddo; enddo; 

  ! k^1, Vdot psi
  temp(:) = imag_i*A*sinth0*cosphi*costh0**3 &
      - imag_i*A*sinth0*cosphi*costh0 + imag_i*C*sinth0*cosphi*costh0**3 &
      - 2*imag_i*F*sinth0*cosphi*costh0**3 &
      + imag_i*F*sinth0*cosphi*costh0 - 4*imag_i*L*sinth0*cosphi*costh0**3 &
      + 2*imag_i*L*sinth0*cosphi*costh0
  do j=1,NGL; do i=1,NGL; 
    K1V(i,j) = K1V(i,j) + temp(i) * weight(i) * hp(i,j)
  enddo; enddo; 

  ! k^2, U psi
  temp(:) = A*sinth0**3*cosphi**3*costh0 - A*sinth0*cosphi*costh0 &
      + C*sinth0**3*cosphi**3*costh0 &
      - 2*F*sinth0**3*cosphi**3*costh0 + F*sinth0*cosphi*costh0 &
      - 4*L*sinth0**3*cosphi**3*costh0 &
      + 2*L*sinth0*cosphi*costh0
  do i=1,NGL; K2U(i,i) = K2U(i,i) + temp(i) * weight(i) * jaco; enddo

  ! k^2, W psi
  temp(:) = A*sinphi*sinth0**3*cosphi**2*costh0 &
      + C*sinphi*sinth0**3*cosphi**2*costh0 - 2*F*sinphi*sinth0**3*cosphi**2*costh0 &
      - 4*L*sinphi*sinth0**3*cosphi**2*costh0 &
      + L*sinphi*sinth0*costh0 &
      - N*sinphi*sinth0*costh0
  do i=1,NGL; K2W(i,i) = K2W(i,i) + temp(i) * weight(i) * jaco; enddo

  ! k^2, V psi
  temp(:) = A*sinth0**2*cosphi**2*costh0**2 + C*sinth0**2*cosphi**2*costh0**2 &
      - 2*F*sinth0**2*cosphi**2*costh0**2 &
      - 4*L*sinth0**2*cosphi**2*costh0**2 + L*sinth0**2*cosphi**2 &
      + L*costh0**2 - N*sinth0**2*cosphi**2 - N*costh0**2 &
      + N
  do i=1,NGL; K2V(i,i) = K2V(i,i) + temp(i) * weight(i) * jaco; enddo

  !transpose
  K0U = transpose(K0U); K0V = transpose(K0V); K0W = transpose(K0W);
  K1U = transpose(K1U); K1V = transpose(K1V); K1W = transpose(K1W);
  K2U = transpose(K2U); K2V = transpose(K2V); K2W = transpose(K2W);

end subroutine get_comp3

subroutine get_kernels (NGL,k,A,C,F,L,N,theta0,dphi, &
                        U,V,W,Udot,Vdot,Wdot,Kwvnm,&
                        KA,KC,KF,KL,KN,Ktheta0,Kdphi,dL_dkv) bind(c,name='get_kernels_')
  use iso_c_binding,only: c_int,dp => c_double ,dcp => c_double_complex
  implicit none

  integer(c_int),value,intent(in)         :: NGL
  real(dp),value,intent(in)               :: k
  real(dp),dimension(NGL),intent(in)      :: A,C,F,L,N,theta0,dphi
  complex(dcp),dimension(NGL),intent(in)  :: U,V,W,Udot,Vdot,Wdot
  real(dcp),dimension(NGL),intent(out)    :: Kwvnm,KA,KC,KF,KL,KN,Ktheta0,Kdphi
  real(dcp),dimension(NGL,2),intent(out)  :: dL_dkv
  
  !local vars
  complex(dcp),dimension(NGL)  :: temp
  real(dp),dimension(NGL) :: costh0, sinth0,cosphi,sinphi
  complex(dcp),parameter :: imag_i = cmplx(0.0,1.0,kind=dcp)

  !init matrices
  costh0 = cos(theta0); sinth0 = sin(theta0);
  cosphi = cos(dphi); sinphi = sin(dphi);

  temp(:) = - U*k**2*(sinphi*sinth0**4*cosphi**3 &
      - sinphi*sinth0**2*cosphi)*conjg(W) - U*k**2*(sinth0**3*cosphi**3*costh0 &
      - sinth0*cosphi*costh0)*conjg(V) - U*k**2*(sinth0**4*cosphi**4 &
      - 2*sinth0**2*cosphi**2 + 1)*conjg(U) &
      + imag_i*U*k*(sinth0**3*cosphi**3*costh0 &
      - sinth0*cosphi*costh0)*conjg(Udot) + imag_i*U*k*(sinphi*sinth0**3*cosphi**2*costh0 &
      - sinphi*sinth0*costh0)*conjg(Wdot) &
      + imag_i*U*k*(sinth0**2*cosphi**2*costh0**2 &
      - sinth0**2*cosphi**2 - costh0**2 + 1)*conjg(Vdot) &
      - imag_i*Udot*k*(sinth0**3*cosphi**3*costh0 - sinth0*cosphi*costh0)*conjg(U) &
      - imag_i*Udot*k*sinphi*sinth0**3*cosphi**2*costh0*conjg(W) &
      - imag_i*Udot*k*sinth0**2*cosphi**2*costh0**2*conjg(V) &
      - Udot*(sinth0*cosphi*costh0**3 &
      - sinth0*cosphi*costh0)*conjg(Vdot) &
      - Udot*sinphi*sinth0**2*cosphi*costh0**2*conjg(Wdot) &
      - Udot*sinth0**2*cosphi**2*costh0**2*conjg(Udot) - V*k**2*(sinth0**3*cosphi**3*costh0 &
      - sinth0*cosphi*costh0)*conjg(U) &
      - V*k**2*sinphi*sinth0**3*cosphi**2*costh0*conjg(W) &
      - V*k**2*sinth0**2*cosphi**2*costh0**2*conjg(V) &
      + imag_i*V*k*(sinth0*cosphi*costh0**3 - sinth0*cosphi*costh0)*conjg(Vdot) &
      + imag_i*V*k*sinphi*sinth0**2*cosphi*costh0**2*conjg(Wdot) &
      + imag_i*V*k*sinth0**2*cosphi**2*costh0**2*conjg(Udot) &
      - imag_i*Vdot*k*(sinth0*cosphi*costh0**3 &
      - sinth0*cosphi*costh0)*conjg(V) &
      - imag_i*Vdot*k*(sinphi*sinth0**2*cosphi*costh0**2 &
      - sinphi*sinth0**2*cosphi)*conjg(W) - imag_i*Vdot*k*(sinth0**2*cosphi**2*costh0**2 &
      - sinth0**2*cosphi**2 - costh0**2 &
      + 1)*conjg(U) - Vdot*(sinphi*sinth0*costh0**3 &
      - sinphi*sinth0*costh0)*conjg(Wdot) - Vdot*(sinth0*cosphi*costh0**3 &
      - sinth0*cosphi*costh0)*conjg(Udot) &
      - Vdot*(costh0**4 - 2*costh0**2 + 1)*conjg(Vdot) &
      - W*k**2*(sinphi*sinth0**4*cosphi**3 - sinphi*sinth0**2*cosphi)*conjg(U) &
      - W*k**2*sinphi**2*sinth0**4*cosphi**2*conjg(W) &
      - W*k**2*sinphi*sinth0**3*cosphi**2*costh0*conjg(V) &
      + imag_i*W*k*(sinphi*sinth0**2*cosphi*costh0**2 &
      - sinphi*sinth0**2*cosphi)*conjg(Vdot) &
      + imag_i*W*k*sinphi**2*sinth0**3*cosphi*costh0*conjg(Wdot) &
      + imag_i*W*k*sinphi*sinth0**3*cosphi**2*costh0*conjg(Udot) &
      - imag_i*Wdot*k*(sinphi*sinth0**3*cosphi**2*costh0 &
      - sinphi*sinth0*costh0)*conjg(U) - imag_i*Wdot*k*sinphi**2*sinth0**3*cosphi*costh0*conjg(W) &
      - imag_i*Wdot*k*sinphi*sinth0**2*cosphi*costh0**2*conjg(V) - Wdot*(sinphi*sinth0*costh0**3 &
      - sinphi*sinth0*costh0)*conjg(Vdot) &
      - Wdot*sinphi**2*sinth0**2*costh0**2*conjg(Wdot) &
      - Wdot*sinphi*sinth0**2*cosphi*costh0**2*conjg(Udot)
  KA(:) = real(temp(:),kind=dp)

  temp(:) = - U*k**2*sinphi*sinth0**4*cosphi**3*conjg(W) &
      - U*k**2*sinth0**4*cosphi**4*conjg(U) - U*k**2*sinth0**3*cosphi**3*costh0*conjg(V) &
      + imag_i*U*k*sinphi*sinth0**3*cosphi**2*costh0*conjg(Wdot) &
      + imag_i*U*k*sinth0**3*cosphi**3*costh0*conjg(Udot) &
      + imag_i*U*k*sinth0**2*cosphi**2*costh0**2*conjg(Vdot) &
      - imag_i*Udot*k*sinphi*sinth0**3*cosphi**2*costh0*conjg(W) &
      - imag_i*Udot*k*sinth0**3*cosphi**3*costh0*conjg(U) &
      - imag_i*Udot*k*sinth0**2*cosphi**2*costh0**2*conjg(V) &
      - Udot*sinphi*sinth0**2*cosphi*costh0**2*conjg(Wdot) &
      - Udot*sinth0**2*cosphi**2*costh0**2*conjg(Udot) &
      - Udot*sinth0*cosphi*costh0**3*conjg(Vdot) &
      - V*k**2*sinphi*sinth0**3*cosphi**2*costh0*conjg(W) &
      - V*k**2*sinth0**3*cosphi**3*costh0*conjg(U) &
      - V*k**2*sinth0**2*cosphi**2*costh0**2*conjg(V) &
      + imag_i*V*k*sinphi*sinth0**2*cosphi*costh0**2*conjg(Wdot) &
      + imag_i*V*k*sinth0**2*cosphi**2*costh0**2*conjg(Udot) &
      + imag_i*V*k*sinth0*cosphi*costh0**3*conjg(Vdot) &
      - imag_i*Vdot*k*sinphi*sinth0**2*cosphi*costh0**2*conjg(W) &
      - imag_i*Vdot*k*sinth0**2*cosphi**2*costh0**2*conjg(U) &
      - imag_i*Vdot*k*sinth0*cosphi*costh0**3*conjg(V) &
      - Vdot*sinphi*sinth0*costh0**3*conjg(Wdot) - Vdot*sinth0*cosphi*costh0**3*conjg(Udot) &
      - Vdot*costh0**4*conjg(Vdot) &
      - W*k**2*sinphi**2*sinth0**4*cosphi**2*conjg(W) &
      - W*k**2*sinphi*sinth0**4*cosphi**3*conjg(U) - W*k**2*sinphi*sinth0**3*cosphi**2*costh0*conjg(V) &
      + imag_i*W*k*sinphi**2*sinth0**3*cosphi*costh0*conjg(Wdot) &
      + imag_i*W*k*sinphi*sinth0**3*cosphi**2*costh0*conjg(Udot) &
      + imag_i*W*k*sinphi*sinth0**2*cosphi*costh0**2*conjg(Vdot) &
      - imag_i*Wdot*k*sinphi**2*sinth0**3*cosphi*costh0*conjg(W) &
      - imag_i*Wdot*k*sinphi*sinth0**3*cosphi**2*costh0*conjg(U) &
      - imag_i*Wdot*k*sinphi*sinth0**2*cosphi*costh0**2*conjg(V) &
      - Wdot*sinphi**2*sinth0**2*costh0**2*conjg(Wdot) &
      - Wdot*sinphi*sinth0**2*cosphi*costh0**2*conjg(Udot) &
      - Wdot*sinphi*sinth0*costh0**3*conjg(Vdot)
  KC(:) = real(temp(:),kind=dp)

  temp(:) = - U*k**2*( - 2*sinth0**4*cosphi**4 + 2*sinth0**2*cosphi**2)*conjg(U) &
      - U*k**2*( - 2*sinphi*sinth0**4*cosphi**3 &
      + sinphi*sinth0**2*cosphi)*conjg(W) - U*k**2*( &
      - 2*sinth0**3*cosphi**3*costh0 + sinth0*cosphi*costh0)*conjg(V) &
      + imag_i*U*k*( - 2*sinth0**3*cosphi**3*costh0 &
      + sinth0*cosphi*costh0)*conjg(Udot) + imag_i*U*k*( &
      - 2*sinphi*sinth0**3*cosphi**2*costh0 + sinphi*sinth0*costh0)*conjg(Wdot) &
      + imag_i*U*k*( - 2*sinth0**2*cosphi**2*costh0**2 &
      + sinth0**2*cosphi**2 + costh0**2)*conjg(Vdot) &
      - imag_i*Udot*k*( - 2*sinth0**3*cosphi**3*costh0 &
      + sinth0*cosphi*costh0)*conjg(U) &
      + 2*imag_i*Udot*k*sinphi*sinth0**3*cosphi**2*costh0*conjg(W) &
      + 2*imag_i*Udot*k*sinth0**2*cosphi**2*costh0**2*conjg(V) - Udot*( &
      - 2*sinth0*cosphi*costh0**3 + sinth0*cosphi*costh0)*conjg(Vdot) &
      + 2*Udot*sinphi*sinth0**2*cosphi*costh0**2*conjg(Wdot) &
      + 2*Udot*sinth0**2*cosphi**2*costh0**2*conjg(Udot) &
      - V*k**2*( - 2*sinth0**3*cosphi**3*costh0 &
      + sinth0*cosphi*costh0)*conjg(U) + 2*V*k**2*sinphi*sinth0**3*cosphi**2*costh0*conjg(W) &
      + 2*V*k**2*sinth0**2*cosphi**2*costh0**2*conjg(V) &
      + imag_i*V*k*( - 2*sinth0*cosphi*costh0**3 &
      + sinth0*cosphi*costh0)*conjg(Vdot) &
      - 2*imag_i*V*k*sinphi*sinth0**2*cosphi*costh0**2*conjg(Wdot) &
      - 2*imag_i*V*k*sinth0**2*cosphi**2*costh0**2*conjg(Udot) &
      - imag_i*Vdot*k*( - 2*sinth0*cosphi*costh0**3 &
      + sinth0*cosphi*costh0)*conjg(V) - imag_i*Vdot*k*( &
      - 2*sinphi*sinth0**2*cosphi*costh0**2 + sinphi*sinth0**2*cosphi)*conjg(W) &
      - imag_i*Vdot*k*( - 2*sinth0**2*cosphi**2*costh0**2 &
      + sinth0**2*cosphi**2 + costh0**2)*conjg(U) &
      - Vdot*( - 2*sinphi*sinth0*costh0**3 &
      + sinphi*sinth0*costh0)*conjg(Wdot) - Vdot*( - 2*sinth0*cosphi*costh0**3 &
      + sinth0*cosphi*costh0)*conjg(Udot) - Vdot*( &
      - 2*costh0**4 + 2*costh0**2)*conjg(Vdot) - W*k**2*( &
      - 2*sinphi*sinth0**4*cosphi**3 + sinphi*sinth0**2*cosphi)*conjg(U) &
      + 2*W*k**2*sinphi**2*sinth0**4*cosphi**2*conjg(W) &
      + 2*W*k**2*sinphi*sinth0**3*cosphi**2*costh0*conjg(V) + imag_i*W*k*( &
      - 2*sinphi*sinth0**2*cosphi*costh0**2 &
      + sinphi*sinth0**2*cosphi)*conjg(Vdot) &
      - 2*imag_i*W*k*sinphi**2*sinth0**3*cosphi*costh0*conjg(Wdot) &
      - 2*imag_i*W*k*sinphi*sinth0**3*cosphi**2*costh0*conjg(Udot) - imag_i*Wdot*k*( &
      - 2*sinphi*sinth0**3*cosphi**2*costh0 &
      + sinphi*sinth0*costh0)*conjg(U) &
      + 2*imag_i*Wdot*k*sinphi**2*sinth0**3*cosphi*costh0*conjg(W) &
      + 2*imag_i*Wdot*k*sinphi*sinth0**2*cosphi*costh0**2*conjg(V) - Wdot*( &
      - 2*sinphi*sinth0*costh0**3 + sinphi*sinth0*costh0)*conjg(Vdot) &
      + 2*Wdot*sinphi**2*sinth0**2*costh0**2*conjg(Wdot) &
      + 2*Wdot*sinphi*sinth0**2*cosphi*costh0**2*conjg(Udot)
  KF(:) = real(temp(:),kind=dp)

  temp(:) = - U*k**2*( - 4*sinth0**4*cosphi**4 + 4*sinth0**2*cosphi**2)*conjg(U) &
      - U*k**2*( - 4*sinphi*sinth0**4*cosphi**3 &
      + 2*sinphi*sinth0**2*cosphi)*conjg(W) - U*k**2*( &
      - 4*sinth0**3*cosphi**3*costh0 + 2*sinth0*cosphi*costh0)*conjg(V) &
      + imag_i*U*k*( - 4*sinth0**3*cosphi**3*costh0 &
      + 2*sinth0*cosphi*costh0)*conjg(Udot) &
      - 4*imag_i*U*k*sinphi*sinth0**3*cosphi**2*costh0*conjg(Wdot) &
      - 4*imag_i*U*k*sinth0**2*cosphi**2*costh0**2*conjg(Vdot) - imag_i*Udot*k*( &
      - 4*sinth0**3*cosphi**3*costh0 &
      + 2*sinth0*cosphi*costh0)*conjg(U) - imag_i*Udot*k*( &
      - 4*sinphi*sinth0**3*cosphi**2*costh0 + sinphi*sinth0*costh0)*conjg(W) &
      - imag_i*Udot*k*( - 4*sinth0**2*cosphi**2*costh0**2 &
      + sinth0**2*cosphi**2 + costh0**2)*conjg(V) - Udot*( &
      - 4*sinth0*cosphi*costh0**3 + 2*sinth0*cosphi*costh0)*conjg(Vdot) &
      - Udot*( - 4*sinphi*sinth0**2*cosphi*costh0**2 &
      + sinphi*sinth0**2*cosphi)*conjg(Wdot) - Udot*( &
      - 4*sinth0**2*cosphi**2*costh0**2 + sinth0**2*cosphi**2 &
      + costh0**2)*conjg(Udot) - V*k**2*( &
      - 4*sinth0**3*cosphi**3*costh0 + 2*sinth0*cosphi*costh0)*conjg(U) &
      - V*k**2*( - 4*sinphi*sinth0**3*cosphi**2*costh0 &
      + sinphi*sinth0*costh0)*conjg(W) - V*k**2*( &
      - 4*sinth0**2*cosphi**2*costh0**2 + sinth0**2*cosphi**2 + costh0**2)*conjg(V) &
      + imag_i*V*k*( - 4*sinth0*cosphi*costh0**3 &
      + 2*sinth0*cosphi*costh0)*conjg(Vdot) + imag_i*V*k*( &
      - 4*sinphi*sinth0**2*cosphi*costh0**2 &
      + sinphi*sinth0**2*cosphi)*conjg(Wdot) + imag_i*V*k*( - 4*sinth0**2*cosphi**2*costh0**2 &
      + sinth0**2*cosphi**2 + costh0**2)*conjg(Udot) &
      - imag_i*Vdot*k*( - 4*sinth0*cosphi*costh0**3 &
      + 2*sinth0*cosphi*costh0)*conjg(V) &
      + 4*imag_i*Vdot*k*sinphi*sinth0**2*cosphi*costh0**2*conjg(W) &
      + 4*imag_i*Vdot*k*sinth0**2*cosphi**2*costh0**2*conjg(U) - Vdot*( &
      - 4*sinphi*sinth0*costh0**3 + 2*sinphi*sinth0*costh0)*conjg(Wdot) &
      - Vdot*( - 4*sinth0*cosphi*costh0**3 &
      + 2*sinth0*cosphi*costh0)*conjg(Udot) - Vdot*( - 4*costh0**4 &
      + 4*costh0**2)*conjg(Vdot) - W*k**2*( - 4*sinphi*sinth0**4*cosphi**3 &
      + 2*sinphi*sinth0**2*cosphi)*conjg(U) &
      - W*k**2*( - 4*sinphi*sinth0**3*cosphi**2*costh0 &
      + sinphi*sinth0*costh0)*conjg(V) - W*k**2*( - 4*sinphi**2*sinth0**4*cosphi**2 &
      + sinphi**2*sinth0**2 + sinth0**2*cosphi**2)*conjg(W) &
      + imag_i*W*k*( - 4*sinphi*sinth0**3*cosphi**2*costh0 &
      + sinphi*sinth0*costh0)*conjg(Udot) &
      + imag_i*W*k*( - 4*sinphi**2*sinth0**3*cosphi*costh0 &
      + sinth0*cosphi*costh0)*conjg(Wdot) &
      - 4*imag_i*W*k*sinphi*sinth0**2*cosphi*costh0**2*conjg(Vdot) - imag_i*Wdot*k*( &
      - 4*sinphi*sinth0**2*cosphi*costh0**2 &
      + sinphi*sinth0**2*cosphi)*conjg(V) - imag_i*Wdot*k*( &
      - 4*sinphi**2*sinth0**3*cosphi*costh0 + sinth0*cosphi*costh0)*conjg(W) &
      + 4*imag_i*Wdot*k*sinphi*sinth0**3*cosphi**2*costh0*conjg(U) &
      - Wdot*( - 4*sinphi*sinth0*costh0**3 &
      + 2*sinphi*sinth0*costh0)*conjg(Vdot) - Wdot*( &
      - 4*sinphi*sinth0**2*cosphi*costh0**2 &
      + sinphi*sinth0**2*cosphi)*conjg(Udot) - Wdot*( - 4*sinphi**2*sinth0**2*costh0**2 &
      + sinphi**2*sinth0**2 + costh0**2)*conjg(Wdot)
  KL(:) = real(temp(:),kind=dp)

  temp(:) = imag_i*U*k*(2*sinth0**2*cosphi**2 + 2*costh0**2 &
      - 2)*conjg(Vdot) + 2*imag_i*U*k*sinphi*sinth0*costh0*conjg(Wdot) &
      - imag_i*Udot*k*( - sinth0**2*cosphi**2 - costh0**2 &
      + 1)*conjg(V) + imag_i*Udot*k*sinphi*sinth0*costh0*conjg(W) &
      - Udot*( - sinth0**2*cosphi**2 - costh0**2 &
      + 1)*conjg(Udot) + Udot*sinphi*sinth0**2*cosphi*conjg(Wdot) &
      - V*k**2*( - sinth0**2*cosphi**2 - costh0**2 + 1)*conjg(V) &
      + V*k**2*sinphi*sinth0*costh0*conjg(W) + imag_i*V*k*( &
      - sinth0**2*cosphi**2 - costh0**2 + 1)*conjg(Udot) &
      - imag_i*V*k*sinphi*sinth0**2*cosphi*conjg(Wdot) &
      - imag_i*Vdot*k*(2*sinth0**2*cosphi**2 + 2*costh0**2 - 2)*conjg(U) &
      - 2*imag_i*Vdot*k*sinphi*sinth0**2*cosphi*conjg(W) - W*k**2*( &
      - sinphi**2*sinth0**2 - sinth0**2*cosphi**2 + 1)*conjg(W) &
      + W*k**2*sinphi*sinth0*costh0*conjg(V) &
      + 2*imag_i*W*k*sinphi*sinth0**2*cosphi*conjg(Vdot) &
      - imag_i*W*k*sinphi*sinth0*costh0*conjg(Udot) &
      - imag_i*W*k*sinth0*cosphi*costh0*conjg(Wdot) + imag_i*Wdot*k*sinphi*sinth0**2*cosphi*conjg(V) &
      - 2*imag_i*Wdot*k*sinphi*sinth0*costh0*conjg(U) &
      + imag_i*Wdot*k*sinth0*cosphi*costh0*conjg(W) &
      - Wdot*( - sinphi**2*sinth0**2 - costh0**2 &
      + 1)*conjg(Wdot) + Wdot*sinphi*sinth0**2*cosphi*conjg(Udot)
  KN(:) = real(temp(:),kind=dp)

  temp(:) = - U*k**2*(8*(L - N)*sinth0*cosphi**2*costh0 + 4*( - A + F &
      + 2*N)*sinth0*cosphi**2*costh0 + 4*(A + C - 2*F &
      - 4*L)*sinth0**3*cosphi**4*costh0)*conjg(U) - U*k**2*(4*(L &
      - N)*sinphi*sinth0*cosphi*costh0 + 2*( - A + F &
      + 2*N)*sinphi*sinth0*cosphi*costh0 + 4*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**3*costh0)*conjg(W) - U*k**2*( - 2*(L &
      - N)*sinth0**2*cosphi + 2*(L - N)*cosphi*costh0**2 - ( &
      - A + F + 2*N)*sinth0**2*cosphi + ( - A + F + 2*N)*cosphi*costh0**2 &
      - (A + C - 2*F - 4*L)*sinth0**4*cosphi**3 + 3*(A &
      + C - 2*F - 4*L)*sinth0**2*cosphi**3*costh0**2)*conjg(V) &
      + imag_i*U*k*((2*sinth0*cosphi**2*costh0 &
      - 2*sinth0*costh0)*( - A + F + 2*N) - 2*(A + C - 2*F &
      - 4*L)*sinth0**3*cosphi**2*costh0 + 2*(A + C - 2*F &
      - 4*L)*sinth0*cosphi**2*costh0**3)*conjg(Vdot) + imag_i*U*k*( - ( - A + F &
      + 2*N)*sinphi*sinth0**2 + ( - A + F + 2*N)*sinphi*costh0**2 - (A &
      + C - 2*F - 4*L)*sinphi*sinth0**4*cosphi**2 + 3*(A + C &
      - 2*F - 4*L)*sinphi*sinth0**2*cosphi**2*costh0**2)*conjg(Wdot) &
      + imag_i*U*k*( - 2*(L - N)*sinth0**2*cosphi + 2*(L &
      - N)*cosphi*costh0**2 - ( - A + F + 2*N)*sinth0**2*cosphi &
      + ( - A + F + 2*N)*cosphi*costh0**2 - (A + C - 2*F &
      - 4*L)*sinth0**4*cosphi**3 + 3*(A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**3*costh0**2)*conjg(Udot) - imag_i*Udot*k*((L &
      - N)*(2*sinth0*cosphi**2*costh0 - 2*sinth0*costh0) - 2*(A &
      + C - 2*F - 4*L)*sinth0**3*cosphi**2*costh0 + 2*(A + C &
      - 2*F - 4*L)*sinth0*cosphi**2*costh0**3)*conjg(V) &
      - imag_i*Udot*k*( - (L - N)*sinphi*sinth0**2 + (L &
      - N)*sinphi*costh0**2 - (A + C - 2*F - 4*L)*sinphi*sinth0**4*cosphi**2 &
      + 3*(A + C - 2*F - 4*L)*sinphi*sinth0**2*cosphi**2*costh0**2)*conjg(W) &
      - imag_i*Udot*k*( - 2*(L - N)*sinth0**2*cosphi &
      + 2*(L - N)*cosphi*costh0**2 - ( - A + F &
      + 2*N)*sinth0**2*cosphi + ( - A + F + 2*N)*cosphi*costh0**2 - (A &
      + C - 2*F - 4*L)*sinth0**4*cosphi**3 + 3*(A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**3*costh0**2)*conjg(U) - Udot*((L &
      - N)*(2*sinth0*cosphi**2*costh0 - 2*sinth0*costh0) &
      - 2*(A + C - 2*F - 4*L)*sinth0**3*cosphi**2*costh0 + 2*(A &
      + C - 2*F - 4*L)*sinth0*cosphi**2*costh0**3)*conjg(Udot) &
      - Udot*(2*(L - N)*sinphi*sinth0*cosphi*costh0 - 2*(A &
      + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi*costh0 + 2*(A &
      + C - 2*F - 4*L)*sinphi*sinth0*cosphi*costh0**3)*conjg(Wdot) &
      - Udot*( - 2*(L - N)*sinth0**2*cosphi + 2*(L &
      - N)*cosphi*costh0**2 - ( - A + F + 2*N)*sinth0**2*cosphi + ( &
      - A + F + 2*N)*cosphi*costh0**2 - 3*(A + C - 2*F &
      - 4*L)*sinth0**2*cosphi*costh0**2 + (A + C - 2*F &
      - 4*L)*cosphi*costh0**4)*conjg(Vdot) - V*k**2*((L - N)*(2*sinth0*cosphi**2*costh0 &
      - 2*sinth0*costh0) - 2*(A + C - 2*F &
      - 4*L)*sinth0**3*cosphi**2*costh0 + 2*(A + C - 2*F &
      - 4*L)*sinth0*cosphi**2*costh0**3)*conjg(V) - V*k**2*( - (L &
      - N)*sinphi*sinth0**2 + (L - N)*sinphi*costh0**2 - (A + C - 2*F &
      - 4*L)*sinphi*sinth0**4*cosphi**2 + 3*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi**2*costh0**2)*conjg(W) - V*k**2*( &
      - 2*(L - N)*sinth0**2*cosphi + 2*(L - N)*cosphi*costh0**2 &
      - ( - A + F + 2*N)*sinth0**2*cosphi + ( - A + F &
      + 2*N)*cosphi*costh0**2 - (A + C - 2*F - 4*L)*sinth0**4*cosphi**3 &
      + 3*(A + C - 2*F - 4*L)*sinth0**2*cosphi**3*costh0**2)*conjg(U) &
      + imag_i*V*k*((L - N)*(2*sinth0*cosphi**2*costh0 &
      - 2*sinth0*costh0) - 2*(A + C - 2*F &
      - 4*L)*sinth0**3*cosphi**2*costh0 + 2*(A + C - 2*F &
      - 4*L)*sinth0*cosphi**2*costh0**3)*conjg(Udot) + imag_i*V*k*(2*(L &
      - N)*sinphi*sinth0*cosphi*costh0 - 2*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi*costh0 + 2*(A + C - 2*F &
      - 4*L)*sinphi*sinth0*cosphi*costh0**3)*conjg(Wdot) + imag_i*V*k*( &
      - 2*(L - N)*sinth0**2*cosphi + 2*(L - N)*cosphi*costh0**2 &
      - ( - A + F + 2*N)*sinth0**2*cosphi + ( - A + F &
      + 2*N)*cosphi*costh0**2 - 3*(A + C - 2*F - 4*L)*sinth0**2*cosphi*costh0**2 &
      + (A + C - 2*F - 4*L)*cosphi*costh0**4)*conjg(Vdot) &
      - imag_i*Vdot*k*((2*sinth0*cosphi**2*costh0 &
      - 2*sinth0*costh0)*( - A + F + 2*N) - 2*(A + C - 2*F &
      - 4*L)*sinth0**3*cosphi**2*costh0 + 2*(A + C - 2*F &
      - 4*L)*sinth0*cosphi**2*costh0**3)*conjg(U) - imag_i*Vdot*k*(2*( - A &
      + F + 2*N)*sinphi*sinth0*cosphi*costh0 - 2*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi*costh0 + 2*(A + C - 2*F &
      - 4*L)*sinphi*sinth0*cosphi*costh0**3)*conjg(W) &
      - imag_i*Vdot*k*( - 2*(L - N)*sinth0**2*cosphi + 2*(L &
      - N)*cosphi*costh0**2 - ( - A + F + 2*N)*sinth0**2*cosphi + ( &
      - A + F + 2*N)*cosphi*costh0**2 - 3*(A + C - 2*F &
      - 4*L)*sinth0**2*cosphi*costh0**2 + (A + C - 2*F - 4*L)*cosphi*costh0**4)*conjg(V) &
      - Vdot*( - 8*(L - N)*sinth0*costh0 - 4*( &
      - A + F + 2*N)*sinth0*costh0 - 4*(A + C - 2*F &
      - 4*L)*sinth0*costh0**3)*conjg(Vdot) - Vdot*( - 2*(L - N)*sinphi*sinth0**2 &
      + 2*(L - N)*sinphi*costh0**2 - ( - A + F &
      + 2*N)*sinphi*sinth0**2 + ( - A + F + 2*N)*sinphi*costh0**2 - 3*(A &
      + C - 2*F - 4*L)*sinphi*sinth0**2*costh0**2 + (A + C &
      - 2*F - 4*L)*sinphi*costh0**4)*conjg(Wdot) - Vdot*( - 2*(L &
      - N)*sinth0**2*cosphi + 2*(L - N)*cosphi*costh0**2 - ( &
      - A + F + 2*N)*sinth0**2*cosphi + ( - A + F + 2*N)*cosphi*costh0**2 &
      - 3*(A + C - 2*F - 4*L)*sinth0**2*cosphi*costh0**2 &
      + (A + C - 2*F - 4*L)*cosphi*costh0**4)*conjg(Udot) &
      - W*k**2*((L - N)*(2*sinphi**2*sinth0*costh0 &
      + 2*sinth0*cosphi**2*costh0) + 4*(A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**3*cosphi**2*costh0)*conjg(W) - W*k**2*(4*(L &
      - N)*sinphi*sinth0*cosphi*costh0 + 2*( - A + F &
      + 2*N)*sinphi*sinth0*cosphi*costh0 + 4*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**3*costh0)*conjg(U) - W*k**2*( - (L &
      - N)*sinphi*sinth0**2 + (L - N)*sinphi*costh0**2 - (A + C - 2*F &
      - 4*L)*sinphi*sinth0**4*cosphi**2 + 3*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi**2*costh0**2)*conjg(V) &
      + imag_i*W*k*(2*( - A + F + 2*N)*sinphi*sinth0*cosphi*costh0 &
      - 2*(A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi*costh0 &
      + 2*(A + C - 2*F - 4*L)*sinphi*sinth0*cosphi*costh0**3)*conjg(Vdot) &
      + imag_i*W*k*( - (L - N)*sinphi*sinth0**2 + (L &
      - N)*sinphi*costh0**2 - (A + C - 2*F &
      - 4*L)*sinphi*sinth0**4*cosphi**2 + 3*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi**2*costh0**2)*conjg(Udot) + imag_i*W*k*( - (L &
      - N)*sinth0**2*cosphi + (L - N)*cosphi*costh0**2 - (A + C &
      - 2*F - 4*L)*sinphi**2*sinth0**4*cosphi + 3*(A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**2*cosphi*costh0**2)*conjg(Wdot) &
      - imag_i*Wdot*k*(2*(L - N)*sinphi*sinth0*cosphi*costh0 &
      - 2*(A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi*costh0 &
      + 2*(A + C - 2*F - 4*L)*sinphi*sinth0*cosphi*costh0**3)*conjg(V) &
      - imag_i*Wdot*k*( - (L - N)*sinth0**2*cosphi &
      + (L - N)*cosphi*costh0**2 - (A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**4*cosphi + 3*(A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**2*cosphi*costh0**2)*conjg(W) - imag_i*Wdot*k*( &
      - ( - A + F + 2*N)*sinphi*sinth0**2 + ( - A + F &
      + 2*N)*sinphi*costh0**2 - (A + C - 2*F - 4*L)*sinphi*sinth0**4*cosphi**2 &
      + 3*(A + C - 2*F - 4*L)*sinphi*sinth0**2*cosphi**2*costh0**2)*conjg(U) &
      - Wdot*((L - N)*(2*sinphi**2*sinth0*costh0 &
      - 2*sinth0*costh0) - 2*(A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**3*costh0 + 2*(A + C - 2*F &
      - 4*L)*sinphi**2*sinth0*costh0**3)*conjg(Wdot) - Wdot*(2*(L &
      - N)*sinphi*sinth0*cosphi*costh0 - 2*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi*costh0 + 2*(A + C - 2*F &
      - 4*L)*sinphi*sinth0*cosphi*costh0**3)*conjg(Udot) - Wdot*( - 2*(L &
      - N)*sinphi*sinth0**2 + 2*(L - N)*sinphi*costh0**2 - ( &
      - A + F + 2*N)*sinphi*sinth0**2 + ( - A + F + 2*N)*sinphi*costh0**2 &
      - 3*(A + C - 2*F - 4*L)*sinphi*sinth0**2*costh0**2 &
      + (A + C - 2*F - 4*L)*sinphi*costh0**4)*conjg(Vdot)
  Ktheta0(:) = real(temp(:),kind=dp)

  temp(:) = - U*k**2*( - 2*(L - N)*sinphi*sinth0*costh0 - ( - A + F &
      + 2*N)*sinphi*sinth0*costh0 - 3*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(V) - U*k**2*( - 8*(L &
      - N)*sinphi*sinth0**2*cosphi - 4*( - A + F &
      + 2*N)*sinphi*sinth0**2*cosphi - 4*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**4*cosphi**3)*conjg(U) - U*k**2*( - 2*(L - N)*sinphi**2*sinth0**2 &
      + 2*(L - N)*sinth0**2*cosphi**2 - ( - A + F &
      + 2*N)*sinphi**2*sinth0**2 + ( - A + F + 2*N)*sinth0**2*cosphi**2 &
      - 3*(A + C - 2*F - 4*L)*sinphi**2*sinth0**4*cosphi**2 &
      + (A + C - 2*F - 4*L)*sinth0**4*cosphi**4)*conjg(W) &
      + imag_i*U*k*( - 2*( - A + F + 2*N)*sinphi*sinth0**2*cosphi - 2*(A &
      + C - 2*F - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(Vdot) &
      + imag_i*U*k*( - 2*(L - N)*sinphi*sinth0*costh0 &
      - ( - A + F + 2*N)*sinphi*sinth0*costh0 - 3*(A + C &
      - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(Udot) &
      + imag_i*U*k*(( - A + F + 2*N)*sinth0*cosphi*costh0 &
      - 2*(A + C - 2*F - 4*L)*sinphi**2*sinth0**3*cosphi*costh0 &
      + (A + C - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(Wdot) &
      - imag_i*Udot*k*( - 2*(L - N)*sinphi*sinth0**2*cosphi &
      - 2*(A + C - 2*F - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(V) &
      - imag_i*Udot*k*( - 2*(L &
      - N)*sinphi*sinth0*costh0 - ( - A + F + 2*N)*sinphi*sinth0*costh0 &
      - 3*(A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(U) &
      - imag_i*Udot*k*((L - N)*sinth0*cosphi*costh0 &
      - 2*(A + C - 2*F - 4*L)*sinphi**2*sinth0**3*cosphi*costh0 &
      + (A + C - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(W) &
      - Udot*( - 2*(L - N)*sinphi*sinth0**2*cosphi &
      - 2*(A + C - 2*F - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(Udot) &
      - Udot*( - 2*(L - N)*sinphi*sinth0*costh0 &
      - ( - A + F + 2*N)*sinphi*sinth0*costh0 - (A &
      + C - 2*F - 4*L)*sinphi*sinth0*costh0**3)*conjg(Vdot) &
      - Udot*( - (L - N)*sinphi**2*sinth0**2 + (L - N)*sinth0**2*cosphi**2 &
      - (A + C - 2*F - 4*L)*sinphi**2*sinth0**2*costh0**2 &
      + (A + C - 2*F - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(Wdot) &
      - V*k**2*( - 2*(L - N)*sinphi*sinth0**2*cosphi &
      - 2*(A + C - 2*F - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(V) &
      - V*k**2*( - 2*(L - N)*sinphi*sinth0*costh0 &
      - ( - A + F + 2*N)*sinphi*sinth0*costh0 - 3*(A + C &
      - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(U) &
      - V*k**2*((L - N)*sinth0*cosphi*costh0 - 2*(A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**3*cosphi*costh0 + (A + C &
      - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(W) &
      + imag_i*V*k*( - 2*(L - N)*sinphi*sinth0**2*cosphi - 2*(A + C &
      - 2*F - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(Udot) &
      + imag_i*V*k*( - 2*(L - N)*sinphi*sinth0*costh0 &
      - ( - A + F + 2*N)*sinphi*sinth0*costh0 - (A + C - 2*F &
      - 4*L)*sinphi*sinth0*costh0**3)*conjg(Vdot) + imag_i*V*k*( &
      - (L - N)*sinphi**2*sinth0**2 + (L - N)*sinth0**2*cosphi**2 &
      - (A + C - 2*F - 4*L)*sinphi**2*sinth0**2*costh0**2 &
      + (A + C - 2*F - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(Wdot) &
      - imag_i*Vdot*k*( - 2*( - A + F + 2*N)*sinphi*sinth0**2*cosphi &
      - 2*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(U) - imag_i*Vdot*k*( - 2*(L &
      - N)*sinphi*sinth0*costh0 - ( - A + F + 2*N)*sinphi*sinth0*costh0 &
      - (A + C - 2*F - 4*L)*sinphi*sinth0*costh0**3)*conjg(V) &
      - imag_i*Vdot*k*( - ( - A + F + 2*N)*sinphi**2*sinth0**2 &
      + ( - A + F + 2*N)*sinth0**2*cosphi**2 - (A + C &
      - 2*F - 4*L)*sinphi**2*sinth0**2*costh0**2 + (A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(W) - Vdot*( &
      - 2*(L - N)*sinphi*sinth0*costh0 - ( - A + F &
      + 2*N)*sinphi*sinth0*costh0 - (A + C - 2*F &
      - 4*L)*sinphi*sinth0*costh0**3)*conjg(Udot) - Vdot*(2*(L - N)*sinth0*cosphi*costh0 &
      + ( - A + F + 2*N)*sinth0*cosphi*costh0 + (A &
      + C - 2*F - 4*L)*sinth0*cosphi*costh0**3)*conjg(Wdot) &
      - W*k**2*( - 2*(A + C - 2*F - 4*L)*sinphi**3*sinth0**4*cosphi &
      + 2*(A + C - 2*F - 4*L)*sinphi*sinth0**4*cosphi**3)*conjg(W) &
      - W*k**2*((L - N)*sinth0*cosphi*costh0 - 2*(A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**3*cosphi*costh0 + (A + C &
      - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(V) - W*k**2*( &
      - 2*(L - N)*sinphi**2*sinth0**2 + 2*(L - N)*sinth0**2*cosphi**2 &
      - ( - A + F + 2*N)*sinphi**2*sinth0**2 + ( - A + F &
      + 2*N)*sinth0**2*cosphi**2 - 3*(A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**4*cosphi**2 + (A + C - 2*F - 4*L)*sinth0**4*cosphi**4)*conjg(U) &
      + imag_i*W*k*( - (L - N)*sinphi*sinth0*costh0 &
      - (A + C - 2*F - 4*L)*sinphi**3*sinth0**3*costh0 &
      + 2*(A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(Wdot) &
      + imag_i*W*k*((L - N)*sinth0*cosphi*costh0 &
      - 2*(A + C - 2*F - 4*L)*sinphi**2*sinth0**3*cosphi*costh0 &
      + (A + C - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(Udot) &
      + imag_i*W*k*( - ( - A + F + 2*N)*sinphi**2*sinth0**2 &
      + ( - A + F + 2*N)*sinth0**2*cosphi**2 - (A + C &
      - 2*F - 4*L)*sinphi**2*sinth0**2*costh0**2 + (A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(Vdot) &
      - imag_i*Wdot*k*( - (L - N)*sinphi*sinth0*costh0 - (A + C &
      - 2*F - 4*L)*sinphi**3*sinth0**3*costh0 + 2*(A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(W) &
      - imag_i*Wdot*k*(( - A + F + 2*N)*sinth0*cosphi*costh0 &
      - 2*(A + C - 2*F - 4*L)*sinphi**2*sinth0**3*cosphi*costh0 &
      + (A + C - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(U) &
      - imag_i*Wdot*k*( - (L - N)*sinphi**2*sinth0**2 + (L &
      - N)*sinth0**2*cosphi**2 - (A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**2*costh0**2 + (A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(V) - Wdot*(2*(L - N)*sinphi*sinth0**2*cosphi &
      + 2*(A + C - 2*F - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(Wdot) &
      - Wdot*(2*(L - N)*sinth0*cosphi*costh0 &
      + ( - A + F + 2*N)*sinth0*cosphi*costh0 + (A &
      + C - 2*F - 4*L)*sinth0*cosphi*costh0**3)*conjg(Vdot) &
      - Wdot*( - (L - N)*sinphi**2*sinth0**2 + (L - N)*sinth0**2*cosphi**2 &
      - (A + C - 2*F - 4*L)*sinphi**2*sinth0**2*costh0**2 &
      + (A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(Udot)
  Kdphi(:) = real(temp(:),kind=dp)

  temp(:) = - 2*U*k*(2*(L - N)*sinphi*sinth0**2*cosphi + ( - A + F &
      + 2*N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**4*cosphi**3)*conjg(W) - 2*U*k*(2*(L &
      - N)*sinth0*cosphi*costh0 + ( - A + F + 2*N)*sinth0*cosphi*costh0 &
      + (A + C - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(V) &
      - 2*U*k*(A + 4*(L - N)*sinth0**2*cosphi**2 + 2*( - A + F &
      + 2*N)*sinth0**2*cosphi**2 + (A + C - 2*F &
      - 4*L)*sinth0**4*cosphi**4)*conjg(U) + imag_i*U*(( - A + F + 2*N)*sinphi*sinth0*costh0 &
      + (A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(Wdot) &
      + imag_i*U*(2*(L &
      - N)*sinth0*cosphi*costh0 + ( - A + F + 2*N)*sinth0*cosphi*costh0 &
      + (A + C - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(Udot) &
      + imag_i*U*(A - 2*N + (sinth0**2*cosphi**2 &
      + costh0**2)*( - A + F + 2*N) + (A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(Vdot) - imag_i*Udot*((L &
      - N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(W) - imag_i*Udot*(N + (L &
      - N)*(sinth0**2*cosphi**2 + costh0**2) + (A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(V) - imag_i*Udot*(2*(L &
      - N)*sinth0*cosphi*costh0 + ( - A + F &
      + 2*N)*sinth0*cosphi*costh0 + (A + C - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(U) &
      - 2*V*k*((L - N)*sinphi*sinth0*costh0 &
      + (A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(W) &
      - 2*V*k*(N + (L - N)*(sinth0**2*cosphi**2 &
      + costh0**2) + (A + C - 2*F - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(V) &
      - 2*V*k*(2*(L - N)*sinth0*cosphi*costh0 &
      + ( - A + F + 2*N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(U) + imag_i*V*((L &
      - N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(Wdot) &
      + imag_i*V*(N + (L - N)*(sinth0**2*cosphi**2 + costh0**2) + (A &
      + C - 2*F - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(Udot) &
      + imag_i*V*(2*(L - N)*sinth0*cosphi*costh0 + ( - A &
      + F + 2*N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinth0*cosphi*costh0**3)*conjg(Vdot) - imag_i*Vdot*(( &
      - A + F + 2*N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(W) &
      - imag_i*Vdot*(2*(L - N)*sinth0*cosphi*costh0 + ( - A + F &
      + 2*N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinth0*cosphi*costh0**3)*conjg(V) - imag_i*Vdot*(A - 2*N &
      + (sinth0**2*cosphi**2 + costh0**2)*( - A + F + 2*N) + (A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(U) - 2*W*k*((L &
      - N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(V) - 2*W*k*(N + (L &
      - N)*(sinphi**2*sinth0**2 + sinth0**2*cosphi**2) + (A + C &
      - 2*F - 4*L)*sinphi**2*sinth0**4*cosphi**2)*conjg(W) &
      - 2*W*k*(2*(L - N)*sinphi*sinth0**2*cosphi + ( - A + F &
      + 2*N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**4*cosphi**3)*conjg(U) + imag_i*W*((L &
      - N)*sinphi*sinth0*costh0 + (A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(Udot) &
      + imag_i*W*((L &
      - N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(Wdot) + imag_i*W*(( - A + F &
      + 2*N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(Vdot) - imag_i*Wdot*((L &
      - N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(V) - imag_i*Wdot*((L &
      - N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(W) - imag_i*Wdot*(( &
      - A + F + 2*N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(U)
  Kwvnm(:) = real(temp(:),kind=dp)

  temp(:) = - 2*U*k*(2*(L - N)*sinphi*sinth0**2*cosphi + ( - A + F &
      + 2*N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**4*cosphi**3)*conjg(W) - 2*U*k*(2*(L &
      - N)*sinth0*cosphi*costh0 + ( - A + F + 2*N)*sinth0*cosphi*costh0 &
      + (A + C - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(V) &
      - 2*U*k*(A + 4*(L - N)*sinth0**2*cosphi**2 + 2*( - A + F &
      + 2*N)*sinth0**2*cosphi**2 + (A + C - 2*F &
      - 4*L)*sinth0**4*cosphi**4)*conjg(U) + imag_i*U*(( - A + F + 2*N)*sinphi*sinth0*costh0 &
      + (A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(Wdot) &
      + imag_i*U*(2*(L &
      - N)*sinth0*cosphi*costh0 + ( - A + F + 2*N)*sinth0*cosphi*costh0 &
      + (A + C - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(Udot) &
      + imag_i*U*(A - 2*N + (sinth0**2*cosphi**2 &
      + costh0**2)*( - A + F + 2*N) + (A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(Vdot) - imag_i*Udot*((L &
      - N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(W) - imag_i*Udot*(N + (L &
      - N)*(sinth0**2*cosphi**2 + costh0**2) + (A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(V) - imag_i*Udot*(2*(L &
      - N)*sinth0*cosphi*costh0 + ( - A + F &
      + 2*N)*sinth0*cosphi*costh0 + (A + C - 2*F - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(U) &
      - 2*V*k*((L - N)*sinphi*sinth0*costh0 &
      + (A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(W) &
      - 2*V*k*(N + (L - N)*(sinth0**2*cosphi**2 &
      + costh0**2) + (A + C - 2*F - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(V) &
      - 2*V*k*(2*(L - N)*sinth0*cosphi*costh0 &
      + ( - A + F + 2*N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinth0**3*cosphi**3*costh0)*conjg(U) + imag_i*V*((L &
      - N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(Wdot) &
      + imag_i*V*(N + (L - N)*(sinth0**2*cosphi**2 + costh0**2) + (A &
      + C - 2*F - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(Udot) &
      + imag_i*V*(2*(L - N)*sinth0*cosphi*costh0 + ( - A &
      + F + 2*N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinth0*cosphi*costh0**3)*conjg(Vdot) - imag_i*Vdot*(( &
      - A + F + 2*N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(W) &
      - imag_i*Vdot*(2*(L - N)*sinth0*cosphi*costh0 + ( - A + F &
      + 2*N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinth0*cosphi*costh0**3)*conjg(V) - imag_i*Vdot*(A - 2*N &
      + (sinth0**2*cosphi**2 + costh0**2)*( - A + F + 2*N) + (A + C - 2*F &
      - 4*L)*sinth0**2*cosphi**2*costh0**2)*conjg(U) - 2*W*k*((L &
      - N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(V) - 2*W*k*(N + (L &
      - N)*(sinphi**2*sinth0**2 + sinth0**2*cosphi**2) + (A + C &
      - 2*F - 4*L)*sinphi**2*sinth0**4*cosphi**2)*conjg(W) &
      - 2*W*k*(2*(L - N)*sinphi*sinth0**2*cosphi + ( - A + F &
      + 2*N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**4*cosphi**3)*conjg(U) + imag_i*W*((L &
      - N)*sinphi*sinth0*costh0 + (A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(Udot) &
      + imag_i*W*((L &
      - N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(Wdot) + imag_i*W*(( - A + F &
      + 2*N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(Vdot) - imag_i*Wdot*((L &
      - N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(V) - imag_i*Wdot*((L &
      - N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(W) - imag_i*Wdot*(( &
      - A + F + 2*N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(U)
  dL_dkv(:,1) = real(temp(:),kind=dp)
  temp(:) = - U*k*((L - N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(V) - U*k*(( &
      - A + F + 2*N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(V) - U*k*(N &
      + (L - N)*(sinphi**2*sinth0**2 + sinth0**2*cosphi**2) &
      + (A + C - 2*F - 4*L)*sinphi**2*sinth0**4*cosphi**2)*conjg(W) &
      - 2*U*k*(2*(L - N)*sinphi*sinth0**2*cosphi + ( - A + F &
      + 2*N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**4*cosphi**3)*conjg(U) - U*k*(A - 2*N &
      + (sinphi**2*sinth0**2 + sinth0**2*cosphi**2)*( - A + F + 2*N) + (A + C &
      - 2*F - 4*L)*sinphi**2*sinth0**4*cosphi**2)*conjg(W) &
      + imag_i*U*((L - N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(Udot) &
      + imag_i*U*((L - N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(Wdot) &
      + imag_i*U*(( - A + F + 2*N)*sinphi*sinth0**2*cosphi + (A &
      + C - 2*F - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(Vdot) &
      - imag_i*Udot*((L - N)*sinphi*sinth0*costh0 &
      + (A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(U) &
      - imag_i*Udot*((L - N)*sinphi*sinth0**2*cosphi &
      + (A + C - 2*F - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(V) &
      - imag_i*Udot*(( - A + F + 2*N)*sinth0*cosphi*costh0 &
      + (A + C - 2*F - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(W) &
      - V*k*((L - N)*sinphi*sinth0*costh0 &
      + (A + C - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(U) &
      - 2*V*k*((L - N)*sinphi*sinth0**2*cosphi &
      + (A + C - 2*F - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(V) &
      - V*k*((L - N)*sinth0*cosphi*costh0 + (A + C &
      - 2*F - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(W) &
      - V*k*(( - A + F + 2*N)*sinphi*sinth0*costh0 + (A + C &
      - 2*F - 4*L)*sinphi*sinth0**3*cosphi**2*costh0)*conjg(U) &
      - V*k*(( - A + F + 2*N)*sinth0*cosphi*costh0 + (A + C &
      - 2*F - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(W) &
      + imag_i*V*((L - N)*sinphi*sinth0**2*cosphi + (A + C &
      - 2*F - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(Udot) &
      + imag_i*V*(N + (L - N)*(sinphi**2*sinth0**2 + costh0**2) &
      + (A + C - 2*F - 4*L)*sinphi**2*sinth0**2*costh0**2)*conjg(Wdot) &
      + imag_i*V*(2*(L - N)*sinphi*sinth0*costh0 &
      + ( - A + F + 2*N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi*sinth0*costh0**3)*conjg(Vdot) - imag_i*Vdot*(( &
      - A + F + 2*N)*sinphi*sinth0**2*cosphi + (A + C &
      - 2*F - 4*L)*sinphi*sinth0**2*cosphi*costh0**2)*conjg(U) &
      - imag_i*Vdot*(2*(L - N)*sinphi*sinth0*costh0 + ( - A &
      + F + 2*N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi*sinth0*costh0**3)*conjg(V) - imag_i*Vdot*(A - 2*N &
      + (sinphi**2*sinth0**2 + costh0**2)*( - A + F + 2*N) + (A &
      + C - 2*F - 4*L)*sinphi**2*sinth0**2*costh0**2)*conjg(W) &
      - W*k*((L - N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(V) - W*k*(( &
      - A + F + 2*N)*sinth0*cosphi*costh0 + (A + C - 2*F &
      - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(V) - W*k*(N &
      + (L - N)*(sinphi**2*sinth0**2 + sinth0**2*cosphi**2) &
      + (A + C - 2*F - 4*L)*sinphi**2*sinth0**4*cosphi**2)*conjg(U) &
      - 2*W*k*(2*(L - N)*sinphi*sinth0**2*cosphi + ( - A + F &
      + 2*N)*sinphi*sinth0**2*cosphi + (A + C - 2*F &
      - 4*L)*sinphi**3*sinth0**4*cosphi)*conjg(W) - W*k*(A - 2*N &
      + (sinphi**2*sinth0**2 + sinth0**2*cosphi**2)*( - A + F + 2*N) + (A + C &
      - 2*F - 4*L)*sinphi**2*sinth0**4*cosphi**2)*conjg(U) &
      + imag_i*W*(( - A + F + 2*N)*sinth0*cosphi*costh0 + (A + C &
      - 2*F - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(Udot) &
      + imag_i*W*(2*(L - N)*sinphi*sinth0*costh0 + ( &
      - A + F + 2*N)*sinphi*sinth0*costh0 + (A + C - 2*F &
      - 4*L)*sinphi**3*sinth0**3*costh0)*conjg(Wdot) + imag_i*W*(A &
      - 2*N + (sinphi**2*sinth0**2 + costh0**2)*( - A + F + 2*N) &
      + (A + C - 2*F - 4*L)*sinphi**2*sinth0**2*costh0**2)*conjg(Vdot) &
      - imag_i*Wdot*((L - N)*sinth0*cosphi*costh0 &
      + (A + C - 2*F - 4*L)*sinphi**2*sinth0**3*cosphi*costh0)*conjg(U) &
      - imag_i*Wdot*(N + (L - N)*(sinphi**2*sinth0**2 &
      + costh0**2) + (A + C - 2*F - 4*L)*sinphi**2*sinth0**2*costh0**2)*conjg(V) &
      - imag_i*Wdot*(2*(L - N)*sinphi*sinth0*costh0 &
      + ( - A + F + 2*N)*sinphi*sinth0*costh0 + (A &
      + C - 2*F - 4*L)*sinphi**3*sinth0**3*costh0)*conjg(W)
  dL_dkv(:,2) = real(temp(:),kind=dp)

end subroutine get_kernels


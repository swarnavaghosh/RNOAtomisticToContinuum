!DEC$ FREEFORM
!DIR$ FREEFORM
! Abaqus 6.14-compatible thread-safe RNO25 VUMAT.
! VEXTERNALDB performs DAT I/O before parallel VUMAT calls.
! Abaqus mutex 1 protects one-time initialization.
! Auxiliary weight-file I/O uses unit 201 (>100) to avoid Abaqus units.
! The authoritative external weights are read from rno_weight.dat in the Abaqus job working directory.

module rno25_checkpoint_module
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: nx = 7, nh = 25, nw = 100, ny = 6
  real(dp), save :: x_min(nx), x_max(nx), y_min(ny), y_max(ny)
  real(dp), save :: fw1(nw,nx+nh), fb1(nw)
  real(dp), save :: fw2(nw,nw), fb2(nw)
  real(dp), save :: fw3(nh,nw), fb3(nh)
  real(dp), save :: gw1(nw,2*nx+nh), gb1(nw)
  real(dp), save :: gw2(nw,nw), gb2(nw)
  real(dp), save :: gw3(ny,nw), gb3(ny)
  real(dp), save :: checkpoint_dt = 0.0_dp
  real(dp), save :: hidden_rate_limit = 0.0_dp
  logical, save :: weights_loaded = .false.
contains
  elemental real(dp) function selu(value)
    real(dp), intent(in) :: value
    real(dp), parameter :: lambda_selu=1.0507009873554805_dp
    real(dp), parameter :: alpha_selu=1.6732632423543772_dp
    if (value > 0.0_dp) then
      selu=lambda_selu*value
    else
      selu=lambda_selu*alpha_selu*(exp(value)-1.0_dp)
    end if
  end function selu

  subroutine weight_file_failure(code)
    integer, intent(in) :: code
    write(*,*) "RNO25 WEIGHT FILE ERROR. CODE=", code
    write(*,*) "Expected file: rno_weight.dat"
    call xplb_exit
  end subroutine weight_file_failure

  subroutine load_rno25_weights_unlocked()
    integer, parameter :: u = 201
    integer :: ios,hf,wf,xf,yf,i,j
    character(len=128), parameter :: weight_file = &
         "rno_weight.dat"

    if (weights_loaded) return

    open(unit=u,file=weight_file,status="old",form="formatted", &
         action="read",iostat=ios)
    if (ios /= 0) call weight_file_failure(9101)

    read(u,*,iostat=ios) hf,wf,xf,yf
    if (ios /= 0) call weight_file_failure(9102)
    if (hf /= nh .or. wf /= nw .or. xf /= nx .or. yf /= ny) then
      call weight_file_failure(9103)
    end if

    read(u,*,iostat=ios) checkpoint_dt,hidden_rate_limit
    if (ios /= 0) call weight_file_failure(9104)

    do i=1,nx
      read(u,*,iostat=ios) x_min(i)
      if (ios /= 0) call weight_file_failure(9110)
    end do
    do i=1,nx
      read(u,*,iostat=ios) x_max(i)
      if (ios /= 0) call weight_file_failure(9111)
    end do
    do i=1,ny
      read(u,*,iostat=ios) y_min(i)
      if (ios /= 0) call weight_file_failure(9112)
    end do
    do i=1,ny
      read(u,*,iostat=ios) y_max(i)
      if (ios /= 0) call weight_file_failure(9113)
    end do

    do i=1,nw
      do j=1,nx+nh
        read(u,*,iostat=ios) fw1(i,j)
        if (ios /= 0) call weight_file_failure(9120)
      end do
    end do
    do i=1,nw
      read(u,*,iostat=ios) fb1(i)
      if (ios /= 0) call weight_file_failure(9121)
    end do
    do i=1,nw
      do j=1,nw
        read(u,*,iostat=ios) fw2(i,j)
        if (ios /= 0) call weight_file_failure(9122)
      end do
    end do
    do i=1,nw
      read(u,*,iostat=ios) fb2(i)
      if (ios /= 0) call weight_file_failure(9123)
    end do
    do i=1,nh
      do j=1,nw
        read(u,*,iostat=ios) fw3(i,j)
        if (ios /= 0) call weight_file_failure(9124)
      end do
    end do
    do i=1,nh
      read(u,*,iostat=ios) fb3(i)
      if (ios /= 0) call weight_file_failure(9125)
    end do

    do i=1,nw
      do j=1,2*nx+nh
        read(u,*,iostat=ios) gw1(i,j)
        if (ios /= 0) call weight_file_failure(9130)
      end do
    end do
    do i=1,nw
      read(u,*,iostat=ios) gb1(i)
      if (ios /= 0) call weight_file_failure(9131)
    end do
    do i=1,nw
      do j=1,nw
        read(u,*,iostat=ios) gw2(i,j)
        if (ios /= 0) call weight_file_failure(9132)
      end do
    end do
    do i=1,nw
      read(u,*,iostat=ios) gb2(i)
      if (ios /= 0) call weight_file_failure(9133)
    end do
    do i=1,ny
      do j=1,nw
        read(u,*,iostat=ios) gw3(i,j)
        if (ios /= 0) call weight_file_failure(9134)
      end do
    end do
    do i=1,ny
      read(u,*,iostat=ios) gb3(i)
      if (ios /= 0) call weight_file_failure(9135)
    end do

    close(unit=u,iostat=ios)
    if (ios /= 0) call weight_file_failure(9140)

    weights_loaded=.true.
  end subroutine load_rno25_weights_unlocked

  subroutine normalize_x(xp,xn)
    real(dp),intent(in)::xp(nx); real(dp),intent(out)::xn(nx)
    xn=(xp-x_min)/max(x_max-x_min,1.0e-12_dp)*2.0_dp-1.0_dp
  end subroutine normalize_x

  subroutine decode_y(yn,yp)
    real(dp),intent(in)::yn(ny); real(dp),intent(out)::yp(ny)
    yp=(yn+1.0_dp)*0.5_dp*max(y_max-y_min,1.0e-12_dp)+y_min
  end subroutine decode_y

  subroutine evaluate_f(x_previous,h_previous,h_rate)
    real(dp),intent(in)::x_previous(nx),h_previous(nh)
    real(dp),intent(out)::h_rate(nh)
    real(dp)::network_input(nx+nh),h1(nw),h2(nw),raw(nh)
    network_input(1:nx)=x_previous
    network_input(nx+1:nx+nh)=h_previous
    h1=selu(matmul(fw1,network_input)+fb1)
    h2=selu(matmul(fw2,h1)+fb2)
    raw=matmul(fw3,h2)+fb3
    if(hidden_rate_limit>0.0_dp) then
      h_rate=hidden_rate_limit*tanh(raw/hidden_rate_limit)
    else
      h_rate=raw
    end if
  end subroutine evaluate_f

  subroutine evaluate_g(x_current,x_rate,h_current,stress_n)
    real(dp),intent(in)::x_current(nx),x_rate(nx),h_current(nh)
    real(dp),intent(out)::stress_n(ny)
    real(dp)::network_input(2*nx+nh),h1(nw),h2(nw)
    network_input(1:nx)=x_current
    network_input(nx+1:2*nx)=x_rate
    network_input(2*nx+1:2*nx+nh)=h_current
    h1=selu(matmul(gw1,network_input)+gb1)
    h2=selu(matmul(gw2,h1)+gb2)
    stress_n=matmul(gw3,h2)+gb3
  end subroutine evaluate_g

  subroutine rno25_step(x_previous_physical,x_current_physical,h_previous,delta_t_normalized,h_current,stress_physical)
    real(dp),intent(in)::x_previous_physical(nx),x_current_physical(nx)
    real(dp),intent(in)::h_previous(nh),delta_t_normalized
    real(dp),intent(out)::h_current(nh),stress_physical(ny)
    real(dp)::xp(nx),xc(nx),xr(nx),hr(nh),sn(ny),effective_dt
    call load_rno25_weights_unlocked()
    effective_dt=max(delta_t_normalized,1.0e-12_dp)
    call normalize_x(x_previous_physical,xp)
    call normalize_x(x_current_physical,xc)
    call evaluate_f(xp,h_previous,hr)
    h_current=h_previous+hr*effective_dt
    xr=(xc-xp)/effective_dt
    call evaluate_g(xc,xr,h_current,sn)
    call decode_y(sn,stress_physical)
  end subroutine rno25_step

  subroutine ensure_rno25_weights()
    if (weights_loaded) return
    call MutexLock(1)
    if (.not. weights_loaded) then
      call load_rno25_weights_unlocked()
    end if
    call MutexUnlock(1)
  end subroutine ensure_rno25_weights

end module rno25_checkpoint_module


subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)
  use rno25_checkpoint_module
  implicit none
  integer, intent(in) :: lOp, niArray, nrArray
  integer, intent(inout) :: i_Array(niArray)
  real(dp), intent(inout) :: r_Array(nrArray)
  integer, parameter :: j_int_StartAnalysis = 0

  if (lOp == j_int_StartAnalysis) then
    call MutexInit(1)
    call MutexLock(1)
    if (.not. weights_loaded) then
      call load_rno25_weights_unlocked()
    end if
    call MutexUnlock(1)
  end if

  return
end subroutine vexternaldb

subroutine vumat( &
     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
     steptime, totaltime, dt, cmname, coordmp, charlength, &
     props, density, straininc, relspininc, &
     tempold, stretchold, defgradold, fieldold, &
     stressold, stateold, enerinternold, enerinelasold, &
     tempnew, stretchnew, defgradnew, fieldnew, &
     stressnew, statenew, enerinternnew, enerinelasnew)

  use rno25_checkpoint_module
  implicit none
  integer,intent(in)::nblock,ndir,nshr,nstatev,nfieldv,nprops,lanneal
  character(len=80),intent(in)::cmname
  real(dp),intent(in)::steptime,totaltime,dt
  real(dp),intent(in)::coordmp(nblock,*),charlength(nblock),props(nprops)
  real(dp),intent(in)::density(nblock),straininc(nblock,ndir+nshr)
  real(dp),intent(in)::relspininc(nblock,nshr),tempold(nblock)
  real(dp),intent(in)::stretchold(nblock,ndir+nshr)
  real(dp),intent(in)::defgradold(nblock,ndir+2*nshr)
  real(dp),intent(in)::fieldold(nblock,nfieldv),stressold(nblock,ndir+nshr)
  real(dp),intent(in)::stateold(nblock,nstatev),enerinternold(nblock)
  real(dp),intent(in)::enerinelasold(nblock)
  real(dp),intent(inout)::tempnew(nblock)
  real(dp),intent(in)::stretchnew(nblock,ndir+nshr)
  real(dp),intent(in)::defgradnew(nblock,ndir+2*nshr)
  real(dp),intent(in)::fieldnew(nblock,nfieldv)
  real(dp),intent(out)::stressnew(nblock,ndir+nshr)
  real(dp),intent(out)::statenew(nblock,nstatev)
  real(dp),intent(out)::enerinternnew(nblock),enerinelasnew(nblock)

  integer::point,component,substep,number_substeps
  real(dp)::previous_x(nx),current_x(nx),interpolated_x(nx)
  real(dp)::substep_previous_x(nx),previous_hidden(nh)
  real(dp)::current_hidden(nh),next_hidden(nh),predicted_stress(ny)
  real(dp)::reference_duration,normalized_dt,substep_dt
  real(dp)::interpolation_fraction,shear_factor,stress_power
  real(dp)::bulk_modulus,two_shear_modulus,elastic_modulus
  real(dp)::poisson_ratio,trace_increment
  if(nstatev<33) then
    write(*,*) "RNO25 VUMAT ERROR: NSTATEV must be at least 33."
    call xplb_exit
  end if
  elastic_modulus=2.82e9_dp
  poisson_ratio=0.286_dp
  if(nprops>=1 .and. props(1)>0.0_dp) elastic_modulus=props(1)
  if(nprops>=2) poisson_ratio=props(2)
  if(nprops<9) then
    write(*,*) "RNO25 VUMAT ERROR: provide 9 constants."
    write(*,*) "PROPS(8)=physical duration; PROPS(9)=shear factor."
    call xplb_exit
  end if
  reference_duration=props(8)
  shear_factor=props(9)
  if(reference_duration<=0.0_dp) call xplb_exit
  if(abs(shear_factor-1.0_dp)>1.0e-12_dp .and. &
     abs(shear_factor-2.0_dp)>1.0e-12_dp) call xplb_exit
  normalized_dt=dt/reference_duration

  if(totaltime==0.0_dp .and. steptime==0.0_dp) then
    two_shear_modulus=elastic_modulus/(1.0_dp+poisson_ratio)
    bulk_modulus=elastic_modulus/(3.0_dp*(1.0_dp-2.0_dp*poisson_ratio))
    do point=1,nblock
      trace_increment=straininc(point,1)+straininc(point,2)+straininc(point,3)
      stressnew(point,1)=stressold(point,1) + &
           (bulk_modulus-two_shear_modulus/3.0_dp)*trace_increment + &
           two_shear_modulus*straininc(point,1)
      stressnew(point,2)=stressold(point,2) + &
           (bulk_modulus-two_shear_modulus/3.0_dp)*trace_increment + &
           two_shear_modulus*straininc(point,2)
      stressnew(point,3)=stressold(point,3) + &
           (bulk_modulus-two_shear_modulus/3.0_dp)*trace_increment + &
           two_shear_modulus*straininc(point,3)
      do component=4,ndir+nshr
        stressnew(point,component)=stressold(point,component) + &
             two_shear_modulus*straininc(point,component)
      end do
      statenew(point,:)=stateold(point,:)
      enerinternnew(point)=enerinternold(point)
      enerinelasnew(point)=enerinelasold(point)
      tempnew(point)=tempold(point)
    end do
    return
  end if

  call ensure_rno25_weights()

  do point=1,nblock
    statenew(point,:)=stateold(point,:)
    previous_hidden=stateold(point,1:25)
    previous_x(1:6)=stateold(point,26:31)
    if(stateold(point,33)<0.5_dp) then
      previous_x(1:6)=0.0_dp
      previous_x(7)=tempold(point)
    else
      previous_x(7)=stateold(point,32)
    end if
    current_x=previous_x
    current_x(1:3)=previous_x(1:3)+straininc(point,1:3)
    if(nshr>=1) current_x(4)=previous_x(4)+shear_factor*straininc(point,4)
    if(nshr>=2) current_x(5)=previous_x(5)+shear_factor*straininc(point,5)
    if(nshr>=3) current_x(6)=previous_x(6)+shear_factor*straininc(point,6)
    current_x(7)=tempold(point)
    number_substeps=max(1,ceiling(normalized_dt/max(checkpoint_dt,1.0e-12_dp)))
    substep_dt=normalized_dt/real(number_substeps,dp)
    substep_previous_x=previous_x
    current_hidden=previous_hidden
    do substep=1,number_substeps
      interpolation_fraction=real(substep,dp)/real(number_substeps,dp)
      interpolated_x=previous_x+interpolation_fraction*(current_x-previous_x)
      call rno25_step(substep_previous_x, interpolated_x, &
           current_hidden, substep_dt, next_hidden, predicted_stress)
      current_hidden=next_hidden
      substep_previous_x=interpolated_x
    end do
    stressnew(point,1:6)=predicted_stress
    statenew(point,1:25)=current_hidden
    statenew(point,26:31)=current_x(1:6)
    statenew(point,32)=current_x(7)
    statenew(point,33)=1.0_dp
    stress_power=0.0_dp
    do component=1,min(3,ndir)
      stress_power=stress_power + 0.5_dp * &
           (stressold(point,component)+stressnew(point,component)) * &
           straininc(point,component)
    end do
    do component=4,ndir+nshr
      stress_power=stress_power + &
           (stressold(point,component)+stressnew(point,component)) * &
           straininc(point,component)
    end do
    enerinternnew(point)=enerinternold(point) + &
         stress_power/max(density(point),1.0e-30_dp)
    enerinelasnew(point)=enerinelasold(point)
    tempnew(point)=tempold(point)
  end do
end subroutine vumat

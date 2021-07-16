! *** Finite Element Analysis
!     Alternating-Direction Explicit Method
!     Simulating Heat Diffusion of a 9x9 cm^2 2D Square
!     dT/dt = D(Del)^2T = D(d2T/dx2 + d2T/dy2)

      Program FEA

      Real, Dimension(9,9,5000000) :: T=0.0
      Open(1,FILE='FEA.out',STATUS='UNKNOWN')


      D = 0.002 !cm^2/s - Thermal Diffusion Coefficient


      dt = 1.0 !s
      dx = .1 !cm
      dy = .1 !cm - lattice def by 1x1 cm squares
      
      CFL = D*((1./dx**2)+(1./dy**2))*dt
      !If(CFL.GT.0.5)Then
      !STOP 'CFL'
      !End If

      Write(1,*)"  Total_Time ","       Gradient ","    dt=",dt

      !Write(*,*)"Enter Minutes of Temperature Ramping"
      !Read(*,*)MINUTES

      HOURS = 10
      !Do 5 HOURS = 10,11,1
      Write(1,*)"  999999","   ",HOURS

      T_SLOPE = 15.0/(HOURS*3600.) !deg C/s... our variable

      !Initial Conditions - all at zero deg celsius
      STEP = 1
      T(:,:,STEP) = 0.0    !degrees celsius

      TOTAL_STEPS = INT(HOURS*3600./dt)

      Do 10 STEP = 2,TOTAL_STEPS
      T(1,:,STEP) = (STEP-1)*dt*T_SLOPE
      T(:,1,STEP) = (STEP-1)*dt*T_SLOPE
      T(9,:,STEP) = (STEP-1)*dt*T_SLOPE
      T(:,9,STEP) = (STEP-1)*dt*T_SLOPE

      Do 20 i=2,8
      Do 30 j=2,8
      
      GROUP = (1/dx**2)*(T(i-1,j,STEP-1) + T(i+1,j,STEP-1)
     & -2.*T(i,j,STEP-1))
     & + (1/dy**2)*(T(i,j-1,STEP-1) + T(i,j+1,STEP-1)
     & -2.*T(i,j,STEP-1))

      T(i,j,STEP) = D*dt*GROUP + T(i,j,STEP-1)

   30 Continue
   20 Continue

      GRADIENT = MAXVAL(T(:,:,STEP))-MINVAL(T(:,:,STEP))
      Write(1,*)(STEP-1)*dt, GRADIENT

   10 Continue
      !Write(*,*)"Grad:",GRADIENT
      !Write(*,*)"Gradient Time:",STEP_GRADIENT*dt,"seconds"
    !5 Continue
      STOP
      End Program FEA
